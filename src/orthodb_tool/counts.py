from __future__ import annotations

import csv
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path

from loguru import logger

from orthodb_tool.api import OrthodbClient, SearchFilters, to_int
from orthodb_tool.taxonomy import (
    TaxonInfo,
    TaxonomyIndex,
    collect_descendant_taxids,
    fetch_taxid_rank_map,
    filter_taxids_by_rank,
    load_taxonomy_index,
)


@dataclass(frozen=True)
class CountOptions:
    parent_taxid: int
    target_rank: str
    universal: float | None
    single_copy: float | None
    min_genes: int | None
    workers: int


@dataclass(frozen=True)
class CountRow:
    taxid: str
    taxon_name: str
    n_genomes: int
    n_og: int


def run_counts(client: OrthodbClient, options: CountOptions) -> list[CountRow]:
    logger.info("Loading OrthoDB taxonomy…")
    index = load_taxonomy_index(client)
    parent_taxid = str(options.parent_taxid)
    if parent_taxid not in index.taxa_by_id:
        raise ValueError(
            f"Parent taxid {options.parent_taxid} was not found in OrthoDB /tree."
        )

    parent_name = index.taxa_by_id[parent_taxid].name
    logger.info("Found parent taxid {} ({}).", parent_taxid, parent_name)
    descendant_taxids = collect_descendant_taxids(parent_taxid, index)
    logger.info(
        "Looking up taxonomy ranks for {} descendant taxids…",
        len(descendant_taxids),
    )
    rank_by_taxid = fetch_taxid_rank_map(
        client,
        descendant_taxids,
        workers=options.workers,
        strict=True,
    )
    descendants = filter_taxids_by_rank(
        descendant_taxids,
        target_rank=options.target_rank,
        rank_by_taxid=rank_by_taxid,
    )
    if descendants:
        logger.info(
            "Found {} descendant taxids at rank {}.",
            len(descendants),
            options.target_rank,
        )
    else:
        logger.info(
            "No descendant taxids matched rank {}, writing an empty result.",
            options.target_rank,
        )

    logger.info("Counting OGs for {} taxids…", len(descendants))
    rows = _count_target_taxa(client, index, descendants, options)
    rows.sort(key=lambda row: (-row.n_og, int(row.taxid)))
    return rows


def write_counts_tsv(rows: list[CountRow], output_path: Path | None) -> str:
    if output_path is None:
        writer = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
        writer.writerow(["taxid", "taxon_name", "n_genomes", "n_og"])
        for row in rows:
            writer.writerow([row.taxid, row.taxon_name, row.n_genomes, row.n_og])
        return "stdout"

    with output_path.open("w", encoding="utf-8", newline="") as file_obj:
        writer = csv.writer(file_obj, delimiter="\t", lineterminator="\n")
        writer.writerow(["taxid", "taxon_name", "n_genomes", "n_og"])
        for row in rows:
            writer.writerow([row.taxid, row.taxon_name, row.n_genomes, row.n_og])
    return str(output_path)


def _count_target_taxa(
    client: OrthodbClient,
    index: TaxonomyIndex,
    target_taxids: list[str],
    options: CountOptions,
) -> list[CountRow]:
    if not target_taxids:
        return []
    if options.workers <= 1:
        return [
            _count_single_taxon(client, index, taxid, options)
            for taxid in target_taxids
        ]

    rows: list[CountRow] = []
    with ThreadPoolExecutor(max_workers=options.workers) as executor:
        futures = {
            executor.submit(_count_single_taxon, client, index, taxid, options): taxid
            for taxid in target_taxids
        }
        for future in as_completed(futures):
            rows.append(future.result())
    return rows


def _count_single_taxon(
    client: OrthodbClient,
    index: TaxonomyIndex,
    taxid: str,
    options: CountOptions,
) -> CountRow:
    filters = SearchFilters(
        level_taxid=taxid,
        universal=options.universal,
        single_copy=options.single_copy,
        min_genes=options.min_genes,
    )
    count = _fetch_og_count(client, filters)
    taxon_name = index.taxa_by_id.get(taxid, TaxonInfo(name="")).name
    logger.debug("Counted {} OGs for taxid {} ({}).", count, taxid, taxon_name)
    return CountRow(
        taxid=taxid,
        taxon_name=taxon_name,
        n_genomes=index.organism_count_by_taxid.get(taxid, 0),
        n_og=count,
    )


def _fetch_og_count(client: OrthodbClient, filters: SearchFilters) -> int:
    if filters.min_genes is None or filters.min_genes <= 1:
        return client.count_search_results(filters)

    matched_count = 0
    for records in client.iter_search_records(filters):
        for record in records:
            if to_int(record.get("gene_count"), default=0) >= filters.min_genes:
                matched_count += 1
    return matched_count
