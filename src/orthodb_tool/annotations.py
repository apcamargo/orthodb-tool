from __future__ import annotations

import csv
import sqlite3
import tempfile
import time
from collections.abc import Sequence
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Iterable

from loguru import logger

from orthodb_tool.api import OrthodbClient

ANNOTATION_VALUE_KEYS = [
    "parent_og_id",
    "cog_category_ids",
    "go_molecular_function_ids",
    "go_biological_process_ids",
    "ec_ids",
    "kegg_ids",
    "interpro_ids",
]
SQLITE_IN_CLAUSE_CHUNK_SIZE = 900


class AnnotationStore:
    """Store annotation and parent OG values for selected OGs."""

    def __init__(self, db_path: Path):
        self.conn = sqlite3.connect(str(db_path))
        self.conn.execute("PRAGMA synchronous=NORMAL")
        self.conn.execute(
            "CREATE TABLE IF NOT EXISTS selected_ogs ("
            "og_id TEXT PRIMARY KEY,"
            "parent_og_id TEXT,"
            "cog_category_ids TEXT,"
            "go_molecular_function_ids TEXT,"
            "go_biological_process_ids TEXT,"
            "ec_ids TEXT,"
            "kegg_ids TEXT,"
            "interpro_ids TEXT"
            ") WITHOUT ROWID"
        )
        self.conn.commit()
        self._selected_og_ids: set[str] = set()

    def add_needed_og_ids(self, og_ids: list[str]) -> None:
        if not og_ids:
            return
        self._selected_og_ids.update(og_id for og_id in og_ids if og_id)
        self.conn.executemany(
            "INSERT OR IGNORE INTO selected_ogs(og_id) VALUES (?)",
            [(og_id,) for og_id in og_ids],
        )

    def commit(self) -> None:
        self.conn.commit()

    def apply_pairs_stream(self, *, stream: Iterable[bytes]) -> tuple[int, int]:
        processed_lines = 0
        updated_rows = 0
        batch: list[tuple[str, str]] = []
        commit_every = 250_000

        self.commit()
        self.conn.execute("BEGIN")
        try:
            for raw_line in stream:
                processed_lines += 1
                line = _decode_utf8_line(raw_line, context="OG pairs stream")
                if not line:
                    continue
                parts = line.split("\t")
                if len(parts) < 2:
                    continue
                if parts[0] not in self._selected_og_ids:
                    continue
                batch.append((parts[1], parts[0]))
                if len(batch) >= 5000:
                    updated_rows += self._update_parent_batch(batch)
                    batch.clear()
                if processed_lines % commit_every == 0:
                    if batch:
                        updated_rows += self._update_parent_batch(batch)
                        batch.clear()
                    self.conn.commit()
                    self.conn.execute("BEGIN")

            if batch:
                updated_rows += self._update_parent_batch(batch)
            self.conn.commit()
            return processed_lines, updated_rows
        except Exception:
            self.conn.rollback()
            raise

    def apply_og_description_stream(
        self, *, stream: Iterable[bytes]
    ) -> tuple[int, int]:
        processed_lines = 0
        updated_rows = 0
        batch: list[tuple[str, str, str, str, str, str, str]] = []
        header_index: dict[str, int] | None = None
        commit_every = 250_000

        self.commit()
        self.conn.execute("BEGIN")
        try:
            for raw_line in stream:
                line = _decode_utf8_line(raw_line, context="OG description stream")
                if not line:
                    continue
                processed_lines += 1
                parts = line.split("\t")
                if header_index is None:
                    header_index = _parse_og_description_header(parts)
                    continue

                og_id = _safe_field(parts, header_index["cluster_id"])
                if not og_id:
                    continue
                if og_id not in self._selected_og_ids:
                    continue

                batch.append(
                    (
                        _normalize_og_description_value(
                            _safe_field(parts, header_index["cog"])
                        ),
                        _normalize_og_description_value(
                            _safe_field(parts, header_index["molfunction_go"])
                        ),
                        _normalize_og_description_value(
                            _safe_field(parts, header_index["bioprocess_go"])
                        ),
                        _normalize_og_description_value(
                            _safe_field(parts, header_index["ec"])
                        ),
                        _normalize_og_description_value(
                            _safe_field(parts, header_index["kegg"])
                        ),
                        _normalize_og_description_value(
                            _safe_field(parts, header_index["interpro"])
                        ),
                        og_id,
                    )
                )

                if len(batch) >= 5000:
                    updated_rows += self._update_annotation_batch(batch)
                    batch.clear()
                if processed_lines % commit_every == 0:
                    if batch:
                        updated_rows += self._update_annotation_batch(batch)
                        batch.clear()
                    self.conn.commit()
                    self.conn.execute("BEGIN")

            if batch:
                updated_rows += self._update_annotation_batch(batch)
            self.conn.commit()
            return processed_lines, updated_rows
        except Exception:
            self.conn.rollback()
            raise

    def get_annotation_values_batch(
        self,
        og_ids: Sequence[str],
    ) -> dict[str, dict[str, str]]:
        unique_og_ids = list(dict.fromkeys(og_id for og_id in og_ids if og_id))
        if not unique_og_ids:
            return {}

        results: dict[str, dict[str, str]] = {}
        select_prefix = (
            "SELECT og_id, parent_og_id, cog_category_ids, go_molecular_function_ids, "
            "go_biological_process_ids, ec_ids, kegg_ids, interpro_ids "
            "FROM selected_ogs WHERE og_id IN ({placeholders})"
        )
        for offset in range(0, len(unique_og_ids), SQLITE_IN_CLAUSE_CHUNK_SIZE):
            chunk = unique_og_ids[offset : offset + SQLITE_IN_CLAUSE_CHUNK_SIZE]
            placeholders = ", ".join("?" for _ in chunk)
            rows = self.conn.execute(
                select_prefix.format(placeholders=placeholders),
                chunk,
            ).fetchall()
            for row in rows:
                og_id = str(row[0])
                results[og_id] = {
                    key: (str(value) if value is not None else "")
                    for key, value in zip(ANNOTATION_VALUE_KEYS, row[1:])
                }
        return results

    def close(self) -> None:
        self.conn.close()

    def _update_parent_batch(self, rows: list[tuple[str, str]]) -> int:
        before = self.conn.total_changes
        self.conn.executemany(
            "UPDATE selected_ogs SET parent_og_id = ? WHERE og_id = ?",
            rows,
        )
        return self.conn.total_changes - before

    def _update_annotation_batch(
        self,
        rows: list[tuple[str, str, str, str, str, str, str]],
    ) -> int:
        before = self.conn.total_changes
        self.conn.executemany(
            "UPDATE selected_ogs "
            "SET cog_category_ids = ?, go_molecular_function_ids = ?, "
            "go_biological_process_ids = ?, ec_ids = ?, kegg_ids = ?, interpro_ids = ? "
            "WHERE og_id = ?",
            rows,
        )
        return self.conn.total_changes - before


def create_annotation_store_if_needed(
    *,
    work_dir: Path,
    name_hint: str,
    with_parent_og: bool,
    with_annotation_ids: bool,
) -> AnnotationStore | None:
    if not with_parent_og and not with_annotation_ids:
        return None
    handle = tempfile.NamedTemporaryFile(
        prefix=f".{name_hint}.",
        suffix=".parents.sqlite3",
        dir=str(work_dir),
        delete=False,
    )
    handle.close()
    db_path = Path(handle.name)
    return AnnotationStore(db_path)


def flush_annotation_ids(
    *,
    annotation_store: AnnotationStore | None,
    pending_ids: list[str],
) -> None:
    if annotation_store is None or not pending_ids:
        return
    annotation_store.add_needed_og_ids(pending_ids)
    pending_ids.clear()


def render_with_annotations(
    *,
    base_tsv_path: Path,
    final_tsv_path: Path,
    annotation_store: AnnotationStore,
    output_columns: list[str],
) -> None:
    with base_tsv_path.open("r", encoding="utf-8", newline="") as in_file:
        reader = csv.DictReader(in_file, delimiter="\t")
        with final_tsv_path.open("w", encoding="utf-8", newline="") as out_file:
            writer = csv.DictWriter(
                out_file,
                fieldnames=output_columns,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            pending_rows: list[dict[str, str]] = []
            chunk_size = 5000
            for row in reader:
                pending_rows.append(row)
                if len(pending_rows) >= chunk_size:
                    _write_rows_with_annotations(
                        writer=writer,
                        rows=pending_rows,
                        annotation_store=annotation_store,
                        output_columns=output_columns,
                    )
                    pending_rows.clear()
            if pending_rows:
                _write_rows_with_annotations(
                    writer=writer,
                    rows=pending_rows,
                    annotation_store=annotation_store,
                    output_columns=output_columns,
                )


def apply_annotation_mapping(
    *,
    client: OrthodbClient,
    annotation_store: AnnotationStore,
    taxids_with_rows: list[str],
    work_dir: Path,
    name_hint: str,
    workers: int,
) -> tuple[int, int]:
    if not taxids_with_rows:
        return 0, 0

    chunk_prefix = f".{name_hint}.og_description.{time.monotonic_ns()}"
    annotation_stream_rows = 0
    annotation_matches = 0

    try:
        if workers <= 1 or len(taxids_with_rows) <= 1:
            for taxid in taxids_with_rows:
                chunk_path = _download_og_description_chunk(
                    client=client,
                    taxid=taxid,
                    temp_dir=work_dir,
                    chunk_prefix=chunk_prefix,
                )
                rows, matched = _apply_annotation_chunk(
                    annotation_store=annotation_store,
                    chunk_path=chunk_path,
                )
                annotation_stream_rows += rows
                annotation_matches += matched
            return annotation_stream_rows, annotation_matches

        with ThreadPoolExecutor(
            max_workers=min(workers, len(taxids_with_rows))
        ) as executor:
            futures = {
                executor.submit(
                    _download_og_description_chunk,
                    client=client,
                    taxid=taxid,
                    temp_dir=work_dir,
                    chunk_prefix=chunk_prefix,
                ): taxid
                for taxid in taxids_with_rows
            }
            for future in as_completed(futures):
                chunk_path = future.result()
                rows, matched = _apply_annotation_chunk(
                    annotation_store=annotation_store,
                    chunk_path=chunk_path,
                )
                annotation_stream_rows += rows
                annotation_matches += matched
        return annotation_stream_rows, annotation_matches
    finally:
        for orphan_path in work_dir.glob(f"{chunk_prefix}*.tsv"):
            orphan_path.unlink(missing_ok=True)


def _download_og_description_chunk(
    *,
    client: OrthodbClient,
    taxid: str,
    temp_dir: Path,
    chunk_prefix: str,
) -> Path:
    handle = tempfile.NamedTemporaryFile(
        prefix=f"{chunk_prefix}.{taxid}.",
        suffix=".tsv",
        dir=str(temp_dir),
        delete=False,
    )
    handle.close()
    chunk_path = Path(handle.name)
    try:
        with client.stream_og_description(taxid) as stream:
            with chunk_path.open("wb") as out_file:
                for chunk in stream:
                    out_file.write(chunk)
        logger.debug("Downloaded annotation data for taxid {}.", taxid)
        return chunk_path
    except Exception:
        chunk_path.unlink(missing_ok=True)
        raise


def _apply_annotation_chunk(
    *,
    annotation_store: AnnotationStore,
    chunk_path: Path,
) -> tuple[int, int]:
    try:
        with chunk_path.open("rb") as chunk_file:
            return annotation_store.apply_og_description_stream(stream=chunk_file)
    finally:
        chunk_path.unlink(missing_ok=True)


def _normalize_og_description_value(value: str) -> str:
    cleaned = value.strip()
    return "" if cleaned == "-" else cleaned


def _decode_utf8_line(raw_line: bytes, *, context: str) -> str:
    try:
        return raw_line.decode("utf-8").rstrip("\r\n")
    except UnicodeDecodeError as exc:
        raise ValueError(f"{context} contains invalid UTF-8 data.") from exc


def _safe_field(parts: list[str], index: int) -> str:
    if index < 0 or index >= len(parts):
        return ""
    return parts[index]


def _parse_og_description_header(parts: list[str]) -> dict[str, int]:
    header_index = {key: idx for idx, key in enumerate(parts)}
    required = [
        "cluster_id",
        "cog",
        "molfunction_go",
        "bioprocess_go",
        "ec",
        "kegg",
        "interpro",
    ]
    missing = [key for key in required if key not in header_index]
    if missing:
        raise ValueError(
            "Malformed og_description header; missing columns: " + ", ".join(missing)
        )
    return header_index


def _write_rows_with_annotations(
    *,
    writer: csv.DictWriter,
    rows: list[dict[str, str]],
    annotation_store: AnnotationStore,
    output_columns: list[str],
) -> None:
    annotations_by_og_id = annotation_store.get_annotation_values_batch(
        [row.get("og_id", "") for row in rows]
    )
    for row in rows:
        annotations = annotations_by_og_id.get(row.get("og_id", ""), {})
        for key in output_columns:
            if key in annotations:
                row[key] = annotations[key]
        writer.writerow({key: row.get(key, "") for key in output_columns})
