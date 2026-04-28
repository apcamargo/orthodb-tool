from __future__ import annotations

import json
import random
import threading
import time
import zlib
from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime, timezone
from email.utils import parsedate_to_datetime
from typing import Any, Callable, TypeVar

import httpx
from loguru import logger

API_BASE_URL = "https://data.orthodb.org/v12"
OG_PAIRS_URL = (
    "https://data.orthodb.org/current/download/odb_data_dump/odb12v2_OG_pairs.tab.gz"
)
RETRIABLE_HTTP_CODES = {429, 500, 502, 503, 504}
SEARCH_PAGE_SIZE = 1_000_000
MAX_RETRY_DELAY_SEC = 60.0

_T = TypeVar("_T")


class ApiRequestError(Exception):
    """Raised for network or HTTP transport failures."""

    def __init__(self, message: str, *, status_code: int | None = None) -> None:
        super().__init__(message)
        self.status_code = status_code


class ApiResponseError(Exception):
    """Raised for malformed or error payloads."""


@dataclass(frozen=True)
class HttpConfig:
    timeout_sec: float = 30.0
    retries: int = 4
    retry_delay_base_sec: float = 1.0
    request_interval_sec: float = 0.0


@dataclass(frozen=True)
class _RetryPolicy:
    max_attempts: int
    base_delay_sec: float
    max_delay_sec: float


@dataclass
class _RetryState:
    attempt_number: int
    start_monotonic_sec: float
    next_sleep_sec: float = 0.0

    @property
    def seconds_since_start(self) -> float:
        return max(0.0, time.monotonic() - self.start_monotonic_sec)


@dataclass(frozen=True)
class SearchFilters:
    level_taxid: str
    universal: float | None = None
    single_copy: float | None = None
    min_genes: int | None = None


class RequestPacer:
    """Adaptive request pacing with a user-provided minimum interval."""

    def __init__(self, min_interval_sec: float) -> None:
        self.min_interval_sec = max(0.0, float(min_interval_sec))
        self.current_interval_sec = self.min_interval_sec
        self.max_interval_sec = max(self.min_interval_sec, 60.0)
        self.last_request_start: float | None = None
        self._lock = threading.Lock()

    def before_request(self) -> None:
        with self._lock:
            if self.last_request_start is not None:
                elapsed = time.monotonic() - self.last_request_start
                sleep_for = self.current_interval_sec - elapsed
                if sleep_for > 0:
                    time.sleep(sleep_for)
            self.last_request_start = time.monotonic()

    def on_success(self) -> None:
        with self._lock:
            if self.current_interval_sec <= self.min_interval_sec:
                self.current_interval_sec = self.min_interval_sec
                return
            self.current_interval_sec = max(
                self.min_interval_sec, self.current_interval_sec * 0.9
            )

    def on_retryable_error(self) -> None:
        with self._lock:
            base = self.current_interval_sec
            if base <= 0:
                base = self.min_interval_sec if self.min_interval_sec > 0 else 0.1
            self.current_interval_sec = min(self.max_interval_sec, base * 1.8)


def to_int(value: Any, default: int = 0) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def format_filter_cutoff(value: float) -> str:
    if abs(value - 1.0) < 1e-9:
        return "1"
    return f"{value:.1f}"


def build_search_params(filters: SearchFilters) -> dict[str, Any]:
    params: dict[str, Any] = {"level": int(filters.level_taxid)}
    if filters.universal is not None:
        params["universal"] = format_filter_cutoff(filters.universal)
    if filters.single_copy is not None:
        params["singlecopy"] = format_filter_cutoff(filters.single_copy)
    return params


def _extract_orthodb_error(payload: Any) -> str | None:
    if payload is None:
        return "API returned null payload."
    if not isinstance(payload, dict):
        return None

    status = payload.get("status")
    message = payload.get("message")
    error = payload.get("error")

    if message in (None, 0, "0"):
        message = None
    if status == "error":
        return str(message or error or "status=error")
    if message is not None:
        return str(message)
    if error not in (None, ""):
        return str(error)
    return None


def _iter_gzip_lines(chunks: Iterator[bytes]) -> Iterator[bytes]:
    decompressor = zlib.decompressobj(16 + zlib.MAX_WBITS)
    pending = b""
    for chunk in chunks:
        data = decompressor.decompress(chunk)
        if not data:
            continue
        pending += data
        parts = pending.splitlines(keepends=True)
        if parts and not parts[-1].endswith((b"\n", b"\r")):
            pending = parts.pop()
        else:
            pending = b""
        for line in parts:
            yield line

    pending += decompressor.flush()
    if pending:
        for line in pending.splitlines(keepends=True):
            yield line


class OrthodbClient:
    """HTTP client for OrthoDB and related supporting endpoints."""

    def __init__(self, config: HttpConfig):
        self.config = config
        self.pacer = RequestPacer(config.request_interval_sec)
        self._retry_policy = _RetryPolicy(
            max_attempts=config.retries + 1,
            base_delay_sec=max(0.0, config.retry_delay_base_sec),
            max_delay_sec=MAX_RETRY_DELAY_SEC,
        )
        self._client = httpx.Client(timeout=config.timeout_sec, follow_redirects=True)

    def close(self) -> None:
        self._client.close()

    def __enter__(self) -> OrthodbClient:
        return self

    def __exit__(self, exc_type: Any, exc: Any, tb: Any) -> None:
        self.close()

    def get_tree(self) -> dict[str, Any]:
        payload = self._request_json("GET", f"{API_BASE_URL}/tree", check_orthodb=True)
        if not isinstance(payload, dict):
            raise ApiResponseError("Malformed /tree payload.")
        return payload

    def search(
        self,
        filters: SearchFilters,
        *,
        skip: int = 0,
        take: int = SEARCH_PAGE_SIZE,
        counts_only: bool = False,
    ) -> dict[str, Any]:
        params = build_search_params(filters)
        if counts_only:
            params["counts_only"] = 1
        else:
            params["skip"] = skip
            params["take"] = take

        payload = self._request_json(
            "GET", f"{API_BASE_URL}/search", params=params, check_orthodb=True
        )
        if not isinstance(payload, dict):
            raise ApiResponseError(
                f"Malformed /search payload for taxid={filters.level_taxid}."
            )
        return payload

    def iter_search_records(
        self,
        filters: SearchFilters,
    ) -> Iterator[list[dict[str, Any]]]:
        skip = 0
        while True:
            payload = self.search(filters, skip=skip, take=SEARCH_PAGE_SIZE)
            records = payload.get("bigdata")
            if records is None:
                records = payload.get("data")

            if records is None:
                total_count = to_int(payload.get("count"), default=-1)
                if total_count >= 0 and skip >= total_count:
                    break
                raise ApiResponseError(
                    "Malformed /search payload: missing bigdata/data array."
                )

            if not isinstance(records, list):
                raise ApiResponseError(
                    "Malformed /search payload: bigdata/data is not a list."
                )
            if not records:
                total_count = to_int(payload.get("count"), default=-1)
                if total_count >= 0 and skip < total_count:
                    raise ApiResponseError(
                        f"Unexpected empty search page at skip={skip} "
                        f"before reported count={total_count}."
                    )
                break

            normalized = [record for record in records if isinstance(record, dict)]
            if normalized:
                yield normalized
            skip += len(records)

            total_count = to_int(payload.get("count"), default=-1)
            if total_count >= 0 and skip >= total_count:
                break

    def count_search_results(self, filters: SearchFilters) -> int:
        payload = self.search(filters, counts_only=True)
        return to_int(payload.get("count"), default=0)

    def get_external_json(self, url: str) -> Any:
        return self._request_json("GET", url, check_orthodb=False)

    @contextmanager
    def stream_og_pairs(self) -> Iterator[Iterator[bytes]]:
        with self._stream_request(OG_PAIRS_URL) as response:
            yield _iter_gzip_lines(response.iter_bytes())

    @contextmanager
    def stream_og_description(self, taxid: str) -> Iterator[Iterator[bytes]]:
        with self._stream_request(
            f"{API_BASE_URL}/og_description", params={"clade": int(taxid)}
        ) as response:
            yield response.iter_bytes()

    @contextmanager
    def stream_fasta(
        self,
        og_id: str,
        *,
        cds: bool = False,
    ) -> Iterator[Iterator[bytes]]:
        params = {
            "id": og_id,
            "seqtype": "cds" if cds else "protein",
        }
        with self._stream_request(f"{API_BASE_URL}/fasta", params=params) as response:
            yield response.iter_bytes()

    def _request_json(
        self,
        method: str,
        url: str,
        *,
        params: dict[str, Any] | None = None,
        check_orthodb: bool,
    ) -> Any:
        def run_request() -> Any:
            response = self._client.request(method, url, params=params)
            response.raise_for_status()
            payload = response.json()
            if check_orthodb:
                message = _extract_orthodb_error(payload)
                if message is not None:
                    raise ApiResponseError(message)
            return payload

        return self._run_with_retry(
            run_request,
            method=method,
            url=url,
            failure_label=f"calling {method} {url}",
        )

    @contextmanager
    def _stream_request(
        self,
        url: str,
        *,
        params: dict[str, Any] | None = None,
    ) -> Iterator[httpx.Response]:
        def run_request() -> httpx.Response:
            request = self._client.build_request("GET", url, params=params)
            response = self._client.send(request, stream=True)
            try:
                response.raise_for_status()
            except Exception:
                response.close()
                raise
            return response

        response = self._run_with_retry(
            run_request,
            method="GET",
            url=url,
            failure_label=f"opening stream {url}",
        )
        try:
            yield response
        finally:
            response.close()

    def _run_with_retry(
        self,
        run_request: Callable[[], _T],
        *,
        method: str,
        url: str,
        failure_label: str,
    ) -> _T:
        state = _RetryState(attempt_number=0, start_monotonic_sec=time.monotonic())

        while True:
            state.attempt_number += 1
            try:
                self.pacer.before_request()
                result = run_request()
                self.pacer.on_success()
                return result
            except ApiResponseError:
                raise
            except (httpx.HTTPError, json.JSONDecodeError) as exc:
                if not self._is_retryable_exception(exc) or self._should_stop_retry(
                    state
                ):
                    raise self._build_request_error(
                        exc,
                        failure_label=failure_label,
                    ) from exc

                self.pacer.on_retryable_error()
                retry_after_sec = self._parse_retry_after_sec(exc)
                state.next_sleep_sec = self._compute_wait_sec(
                    state,
                    retry_after_sec=retry_after_sec,
                )
                self._log_retry(
                    state,
                    exc,
                    method=method,
                    url=url,
                    retry_after_sec=retry_after_sec,
                )
                if state.next_sleep_sec > 0:
                    time.sleep(state.next_sleep_sec)

    def _is_retryable_exception(self, exc: BaseException) -> bool:
        if isinstance(exc, httpx.HTTPStatusError):
            return exc.response.status_code in RETRIABLE_HTTP_CODES
        return isinstance(exc, (httpx.HTTPError, json.JSONDecodeError))

    def _should_stop_retry(self, state: _RetryState) -> bool:
        return state.attempt_number >= self._retry_policy.max_attempts

    def _compute_wait_sec(
        self,
        state: _RetryState,
        *,
        retry_after_sec: float | None,
    ) -> float:
        if retry_after_sec is not None:
            return retry_after_sec

        if self._retry_policy.base_delay_sec <= 0:
            return 0.0

        raw_delay_sec = min(
            self._retry_policy.base_delay_sec * (2 ** (state.attempt_number - 1)),
            self._retry_policy.max_delay_sec,
        )
        return random.uniform(0.0, raw_delay_sec)

    def _parse_retry_after_sec(self, exc: BaseException) -> float | None:
        if not isinstance(exc, httpx.HTTPStatusError):
            return None

        header_value = exc.response.headers.get("Retry-After")
        if header_value is None:
            return None

        raw_value = header_value.strip()
        if not raw_value:
            return None

        try:
            return max(0.0, float(raw_value))
        except ValueError:
            pass

        try:
            retry_at = parsedate_to_datetime(raw_value)
        except (TypeError, ValueError, IndexError, OverflowError):
            return None

        if retry_at.tzinfo is None:
            retry_at = retry_at.replace(tzinfo=timezone.utc)
        return max(0.0, (retry_at - datetime.now(timezone.utc)).total_seconds())

    def _log_retry(
        self,
        state: _RetryState,
        exc: BaseException,
        *,
        method: str,
        url: str,
        retry_after_sec: float | None,
    ) -> None:
        status_code: int | None = None
        if isinstance(exc, httpx.HTTPStatusError):
            status_code = exc.response.status_code

        reason = (
            f"HTTP {status_code}" if status_code is not None else type(exc).__name__
        )
        retry_after_detail = (
            f", Retry-After={retry_after_sec:.1f}s"
            if retry_after_sec is not None
            else ""
        )
        logger.debug(
            "Request failed with {}, retrying {} {} in {:.1f}s "
            "(attempt {}/{}{}, elapsed {:.1f}s).",
            reason,
            method,
            url,
            state.next_sleep_sec,
            state.attempt_number + 1,
            self._retry_policy.max_attempts,
            retry_after_detail,
            state.seconds_since_start,
        )

    def _build_request_error(
        self,
        exc: BaseException,
        *,
        failure_label: str,
    ) -> ApiRequestError:
        if isinstance(exc, httpx.HTTPStatusError):
            detail = exc.response.text.strip() or str(exc)
            return ApiRequestError(
                f"HTTP {exc.response.status_code} while {failure_label}: {detail}",
                status_code=exc.response.status_code,
            )
        return ApiRequestError(f"Failed {failure_label}: {exc}")
