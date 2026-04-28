from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from urllib.parse import urlencode

from loguru import logger

from orthodb_tool.api import ApiRequestError, ApiResponseError, OrthodbClient

NCBI_MAX_URL_LENGTH = 1800
NCBI_TAXONOMY_ESUMMARY_URL = (
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
)


@dataclass(frozen=True)
class TaxonInfo:
    name: str


@dataclass(frozen=True)
class TaxonomyIndex:
    taxa_by_id: dict[str, TaxonInfo]
    children_by_id: dict[str, list[str]]
    organism_count_by_taxid: dict[str, int]


def normalize_rank(value: str) -> str:
    return value.strip().lower().replace("_", " ")


def taxid_sort_key(taxid: str) -> int:
    return int(taxid)


def iter_taxid_batches(
    taxids: list[str],
    *,
    max_url_length: int = NCBI_MAX_URL_LENGTH,
) -> list[tuple[int, list[str]]]:
    batches: list[tuple[int, list[str]]] = []
    current_batch: list[str] = []
    current_offset = 0

    for offset, taxid in enumerate(taxids):
        if not current_batch:
            current_batch = [taxid]
            current_offset = offset
            continue

        candidate_batch = [*current_batch, taxid]
        if len(_build_taxonomy_esummary_url(candidate_batch)) <= max_url_length:
            current_batch = candidate_batch
            continue

        batches.append((current_offset, current_batch))
        current_batch = [taxid]
        current_offset = offset

    if current_batch:
        batches.append((current_offset, current_batch))

    return batches


def load_taxonomy_index(client: OrthodbClient) -> TaxonomyIndex:
    payload = client.get_tree()
    roots = payload.get("data")
    if not isinstance(roots, list):
        raise ApiResponseError("Malformed /tree payload: missing data list.")

    taxa_by_id: dict[str, TaxonInfo] = {}
    children_by_id: dict[str, list[str]] = {}
    all_children_by_id: dict[str, list[str]] = {}
    stack: list[dict[str, Any]] = [node for node in roots if isinstance(node, dict)]

    while stack:
        node = stack.pop()
        key = str(node.get("key", "")).strip()
        if key:
            parent_taxid = str(node.get("parent", "")).strip()
            if parent_taxid:
                all_children_by_id.setdefault(parent_taxid, []).append(key)
            if "_" not in key and key.isdigit():
                taxa_by_id[key] = TaxonInfo(name=str(node.get("name", "")).strip())
                if parent_taxid and parent_taxid.isdigit():
                    children_by_id.setdefault(parent_taxid, []).append(key)
        children = node.get("children")
        if isinstance(children, list):
            stack.extend(child for child in children if isinstance(child, dict))

    if not taxa_by_id:
        raise ApiResponseError("No taxonomic-level nodes found in /tree.")

    return TaxonomyIndex(
        taxa_by_id=taxa_by_id,
        children_by_id=children_by_id,
        organism_count_by_taxid=_build_organism_count_index(
            taxa_by_id=taxa_by_id,
            all_children_by_id=all_children_by_id,
        ),
    )


def _build_organism_count_index(
    *,
    taxa_by_id: dict[str, TaxonInfo],
    all_children_by_id: dict[str, list[str]],
) -> dict[str, int]:
    organism_count_by_taxid: dict[str, int] = {}

    def count_organisms(node_id: str) -> int:
        cached = organism_count_by_taxid.get(node_id)
        if cached is not None:
            return cached
        if node_id not in taxa_by_id:
            return 1

        total = 0
        for child_id in all_children_by_id.get(node_id, []):
            total += count_organisms(child_id)
        organism_count_by_taxid[node_id] = total
        return total

    for taxid in taxa_by_id:
        count_organisms(taxid)

    return organism_count_by_taxid


def collect_descendant_taxids(parent_taxid: str, index: TaxonomyIndex) -> list[str]:
    descendants: set[str] = set()
    stack: list[str] = list(index.children_by_id.get(parent_taxid, []))
    while stack:
        taxid = stack.pop()
        if taxid in descendants:
            continue
        descendants.add(taxid)
        stack.extend(index.children_by_id.get(taxid, []))
    return sorted(descendants, key=taxid_sort_key)


def expand_descendants(seed_taxids: set[str], index: TaxonomyIndex) -> set[str]:
    expanded = set(seed_taxids)
    stack = list(seed_taxids)
    while stack:
        current = stack.pop()
        for child in index.children_by_id.get(current, []):
            if child not in expanded:
                expanded.add(child)
                stack.append(child)
    return expanded


def parse_taxids_csv(raw: str | None) -> list[str]:
    if not raw:
        return []
    parsed: list[str] = []
    for token in raw.split(","):
        value = token.strip()
        if not value:
            continue
        if not value.isdigit():
            raise ValueError(f"Invalid taxid: {value}")
        parsed.append(value)
    return parsed


def parse_taxids_file(path: Path | None) -> list[str]:
    if path is None:
        return []
    parsed: list[str] = []
    for line_no, raw_line in enumerate(
        path.read_text(encoding="utf-8").splitlines(), start=1
    ):
        line = raw_line.split("#", 1)[0].strip()
        if not line:
            continue
        if not line.isdigit():
            raise ValueError(
                f"Invalid taxid in {path} at line {line_no}: {raw_line.strip()}"
            )
        parsed.append(line)
    return parsed


def resolve_selected_taxids(
    index: TaxonomyIndex,
    *,
    taxids_raw: str | None,
    taxids_file: Path | None,
    include_descendant_taxids: bool,
) -> list[str]:
    requested = set(parse_taxids_csv(taxids_raw))
    requested.update(parse_taxids_file(taxids_file))

    if requested:
        missing = sorted(taxid for taxid in requested if taxid not in index.taxa_by_id)
        if missing:
            raise ValueError(
                "Requested taxids are not present in OrthoDB tree: "
                + ", ".join(missing)
            )
        selected = set(requested)
    else:
        selected = set(index.taxa_by_id)

    if include_descendant_taxids:
        selected = expand_descendants(selected, index)

    selected.discard("1")
    if not selected:
        raise ValueError("No taxids selected after filtering.")

    return sorted(selected, key=taxid_sort_key)


def resolve_target_taxids_by_rank(
    index: TaxonomyIndex,
    *,
    parent_taxid: str,
    target_rank: str,
    rank_by_taxid: dict[str, str],
) -> list[str]:
    return filter_taxids_by_rank(
        collect_descendant_taxids(parent_taxid, index),
        target_rank=target_rank,
        rank_by_taxid=rank_by_taxid,
    )


def filter_taxids_by_rank(
    taxids: list[str],
    *,
    target_rank: str,
    rank_by_taxid: dict[str, str],
) -> list[str]:
    target_rank_normalized = normalize_rank(target_rank)
    return [
        taxid
        for taxid in taxids
        if normalize_rank(rank_by_taxid.get(taxid, "")) == target_rank_normalized
    ]


def fetch_taxid_rank_map(
    client: OrthodbClient,
    taxids: list[str],
    *,
    workers: int,
    strict: bool,
) -> dict[str, str]:
    if not taxids:
        return {}

    unique_taxids = sorted({str(taxid) for taxid in taxids}, key=taxid_sort_key)
    rank_by_taxid: dict[str, str] = {taxid: "" for taxid in unique_taxids}
    batches = iter_taxid_batches(unique_taxids)

    if workers <= 1 or len(batches) <= 1:
        for offset, batch in batches:
            _apply_rank_batch(
                rank_by_taxid=rank_by_taxid,
                batch=batch,
                offset=offset,
                client=client,
                strict=strict,
            )
        return rank_by_taxid

    with ThreadPoolExecutor(max_workers=min(workers, len(batches))) as executor:
        futures = {
            executor.submit(_fetch_taxid_rank_batch, client, batch): (offset, batch)
            for offset, batch in batches
        }
        for future in as_completed(futures):
            offset, batch = futures[future]
            try:
                rank_by_taxid.update(future.result())
            except (ApiRequestError, ApiResponseError) as exc:
                if strict:
                    raise
                logger.warning(
                    "Could not look up taxonomy ranks for {} taxids, those rank "
                    "values will be left blank.",
                    len(batch),
                )
                logger.debug(
                    "Taxonomy rank lookup failed for batch starting at {}: {}",
                    offset,
                    exc,
                )
                for taxid in batch:
                    rank_by_taxid[taxid] = ""

    return rank_by_taxid


def _apply_rank_batch(
    *,
    rank_by_taxid: dict[str, str],
    batch: list[str],
    offset: int,
    client: OrthodbClient,
    strict: bool,
) -> None:
    try:
        rank_by_taxid.update(_fetch_taxid_rank_batch(client, batch))
    except (ApiRequestError, ApiResponseError) as exc:
        if strict:
            raise
        logger.warning(
            "Could not look up taxonomy ranks for {} taxids, those rank values will "
            "be left blank.",
            len(batch),
        )
        logger.debug(
            "Taxonomy rank lookup failed for batch starting at {}: {}",
            offset,
            exc,
        )
        for taxid in batch:
            rank_by_taxid[taxid] = ""


def _fetch_taxid_rank_batch(
    client: OrthodbClient,
    batch: list[str],
) -> dict[str, str]:
    try:
        payload = client.get_external_json(_build_taxonomy_esummary_url(batch))
    except ApiRequestError as exc:
        if exc.status_code == 414 and len(batch) > 1:
            midpoint = len(batch) // 2
            split_updates = _fetch_taxid_rank_batch(client, batch[:midpoint])
            split_updates.update(_fetch_taxid_rank_batch(client, batch[midpoint:]))
            return split_updates
        raise

    result_raw = payload.get("result") if isinstance(payload, dict) else None
    if not isinstance(result_raw, dict):
        raise ApiResponseError("Malformed NCBI taxonomy response.")

    result = {str(key): value for key, value in result_raw.items()}
    updates: dict[str, str] = {}
    for taxid in batch:
        entry = result.get(taxid)
        if not isinstance(entry, dict):
            updates[taxid] = ""
            continue
        rank = entry.get("rank")
        updates[taxid] = rank.strip() if isinstance(rank, str) else ""
    return updates


def _build_taxonomy_esummary_url(taxids: list[str]) -> str:
    query = urlencode({"db": "taxonomy", "id": ",".join(taxids), "retmode": "json"})
    return f"{NCBI_TAXONOMY_ESUMMARY_URL}?{query}"
