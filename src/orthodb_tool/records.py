from __future__ import annotations

import csv
import sys
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable

from loguru import logger

from orthodb_tool.annotations import (
    AnnotationStore,
    apply_annotation_mapping,
    create_annotation_store_if_needed,
    flush_annotation_ids,
    render_with_annotations,
)
from orthodb_tool.api import OrthodbClient, SearchFilters, to_int
from orthodb_tool.taxonomy import (
    TaxonInfo,
    fetch_taxid_rank_map,
    load_taxonomy_index,
    resolve_selected_taxids,
    taxid_sort_key,
)

BASE_OUTPUT_COLUMNS = [
    "og_id",
    "taxid",
    "taxon_name",
    "protein_count",
    "description",
]
ANNOTATION_OUTPUT_COLUMNS = [
    "cog_category_ids",
    "go_molecular_function_ids",
    "go_biological_process_ids",
    "ec_ids",
    "kegg_ids",
    "interpro_ids",
]


@dataclass
class ExportStats:
    exported_rows: int = 0
    filtered_rows: int = 0
    scanned_rows: int = 0
    search_pages: int = 0
    exported_taxids: set[str] = field(default_factory=set)


@dataclass
class TaxidChunkResult:
    chunk_path: Path
    stats: ExportStats


@dataclass(frozen=True)
class RecordsOptions:
    output_path: Path | None
    taxids: str | None
    taxids_file: Path | None
    include_descendant_taxids: bool
    min_genes: int | None
    universal: float | None
    single_copy: float | None
    workers: int
    with_parent_og: bool
    with_annotation_ids: bool
    with_taxid_rank: bool


@dataclass(frozen=True)
class RecordExportResult:
    output_label: str
    stats: ExportStats
    parent_mapping_applied: bool
    annotation_mapping_applied: bool


@dataclass(frozen=True)
class _RowExportContext:
    client: OrthodbClient
    options: RecordsOptions
    selected_taxids: list[str]
    taxa_by_id: dict[str, TaxonInfo]
    rank_by_taxid: dict[str, str]
    output_columns: list[str]


def run_records(client: OrthodbClient, options: RecordsOptions) -> RecordExportResult:
    logger.info("Loading OrthoDB taxonomy…")
    index = load_taxonomy_index(client)
    selected_taxids = resolve_selected_taxids(
        index,
        taxids_raw=options.taxids,
        taxids_file=options.taxids_file,
        include_descendant_taxids=options.include_descendant_taxids,
    )
    logger.info("Selected {} taxids for export.", len(selected_taxids))
    rank_by_taxid: dict[str, str] = {}
    if options.with_taxid_rank:
        logger.info("Looking up taxonomy ranks for {} taxids…", len(selected_taxids))
        rank_by_taxid = fetch_taxid_rank_map(
            client, selected_taxids, workers=options.workers, strict=False
        )
        rank_count = sum(1 for rank in rank_by_taxid.values() if rank)
        logger.info(
            "Added taxonomy ranks for {} taxids, {} ranks were unavailable.",
            rank_count,
            len(rank_by_taxid) - rank_count,
        )

    output_columns = build_output_columns(
        with_taxid_rank=options.with_taxid_rank,
        with_annotation_ids=options.with_annotation_ids,
        with_parent_og=options.with_parent_og,
    )

    temp_dir_handle: tempfile.TemporaryDirectory[str]
    if options.output_path is None:
        temp_dir_handle = tempfile.TemporaryDirectory(prefix="orthodb_tool.")
        name_hint = "stdout"
    else:
        temp_dir_handle = tempfile.TemporaryDirectory(
            prefix=f".{options.output_path.name}.",
            dir=str(options.output_path.parent),
        )
        name_hint = options.output_path.name

    work_dir = Path(temp_dir_handle.name)
    base_tmp_path = work_dir / "base.tsv"
    final_tmp_path = work_dir / "final.tsv"
    annotation_store: AnnotationStore | None = None
    stats = ExportStats()
    parent_mapping_applied = False
    annotation_mapping_applied = False

    try:
        annotation_store = create_annotation_store_if_needed(
            work_dir=work_dir,
            name_hint=name_hint,
            with_parent_og=options.with_parent_og,
            with_annotation_ids=options.with_annotation_ids,
        )
        row_context = _RowExportContext(
            client=client,
            options=options,
            selected_taxids=selected_taxids,
            taxa_by_id=index.taxa_by_id,
            rank_by_taxid=rank_by_taxid,
            output_columns=output_columns,
        )

        worker_detail = (
            f" using {options.workers} workers" if options.workers > 1 else ""
        )
        logger.info(
            "Searching OrthoDB and writing matching OG records{}…",
            worker_detail,
        )
        if options.workers <= 1:
            with base_tmp_path.open("w", encoding="utf-8", newline="") as out_file:
                writer = csv.writer(out_file, delimiter="\t", lineterminator="\n")
                writer.writerow(output_columns)
                stats = _write_base_rows(
                    context=row_context,
                    writer=writer,
                    annotation_store=annotation_store,
                )
        else:
            stats = _write_base_rows_parallel(
                context=row_context,
                base_tsv_path=base_tmp_path,
                annotation_store=annotation_store,
            )

        logger.info(
            "Found {} OG records after scanning {} records across {} search pages.",
            stats.exported_rows,
            stats.scanned_rows,
            stats.search_pages,
        )
        if stats.filtered_rows > 0:
            logger.info(
                "Skipped {} records below the gene-count filter.",
                stats.filtered_rows,
            )

        if (
            options.with_parent_og
            and annotation_store is not None
            and stats.exported_rows > 0
        ):
            logger.info("Adding parent OG IDs…")
            parent_mapping_applied = True
            annotation_store.commit()
            with client.stream_og_pairs() as stream:
                parent_pairs_scanned, parent_ogs_added = (
                    annotation_store.apply_pairs_stream(stream=stream)
                )
            logger.info(
                "Added parent OG IDs for {} of {} OG records.",
                parent_ogs_added,
                stats.exported_rows,
            )
            logger.debug("Scanned {} parent OG pair rows.", parent_pairs_scanned)
        elif options.with_parent_og:
            logger.info(
                "Skipping parent OG IDs because no OG records matched the filters."
            )

        if (
            options.with_annotation_ids
            and annotation_store is not None
            and stats.exported_rows > 0
        ):
            taxids_with_rows = sorted(stats.exported_taxids, key=taxid_sort_key)
            logger.info(
                "Adding annotation IDs for {} taxids with exported records…",
                len(taxids_with_rows),
            )
            annotation_mapping_applied = True
            annotation_store.commit()
            annotation_rows_scanned, annotation_ogs_added = apply_annotation_mapping(
                client=client,
                annotation_store=annotation_store,
                taxids_with_rows=taxids_with_rows,
                work_dir=work_dir,
                name_hint=name_hint,
                workers=options.workers,
            )
            logger.info(
                "Added annotation IDs for {} of {} OG records.",
                annotation_ogs_added,
                stats.exported_rows,
            )
            logger.debug("Scanned {} annotation rows.", annotation_rows_scanned)
        elif options.with_annotation_ids:
            logger.info(
                "Skipping annotation IDs because no OG records matched the filters."
            )

        if annotation_store is not None and (
            parent_mapping_applied or annotation_mapping_applied
        ):
            logger.info("Writing enriched records…")
            render_with_annotations(
                base_tsv_path=base_tmp_path,
                final_tsv_path=final_tmp_path,
                annotation_store=annotation_store,
                output_columns=output_columns,
            )

        result_tsv_path = (
            final_tmp_path
            if parent_mapping_applied or annotation_mapping_applied
            else base_tmp_path
        )
        output_target = (
            str(options.output_path) if options.output_path is not None else "stdout"
        )
        logger.info("Writing record table to {}…", output_target)
        output_label = _emit_output_tsv(result_tsv_path, options.output_path)
        return RecordExportResult(
            output_label=output_label,
            stats=stats,
            parent_mapping_applied=parent_mapping_applied,
            annotation_mapping_applied=annotation_mapping_applied,
        )
    finally:
        if annotation_store is not None:
            annotation_store.close()
        temp_dir_handle.cleanup()


def build_output_columns(
    *,
    with_taxid_rank: bool,
    with_annotation_ids: bool,
    with_parent_og: bool,
) -> list[str]:
    columns = ["og_id", "taxid", "taxon_name"]
    if with_taxid_rank:
        columns.append("taxid_rank")
    columns.extend(BASE_OUTPUT_COLUMNS[3:])
    if with_annotation_ids:
        columns.extend(ANNOTATION_OUTPUT_COLUMNS)
    if with_parent_og:
        columns.append("parent_og_id")
    return columns


def build_output_row_values(
    *,
    og_id: str,
    og_taxid: str,
    taxon_name: str,
    protein_count: str,
    description: str,
    rank_by_taxid: dict[str, str],
    output_columns: list[str],
) -> list[str]:
    values = {
        "og_id": og_id,
        "taxid": og_taxid,
        "taxon_name": taxon_name,
        "taxid_rank": rank_by_taxid.get(og_taxid, ""),
        "protein_count": protein_count,
        "description": description,
        "cog_category_ids": "",
        "go_molecular_function_ids": "",
        "go_biological_process_ids": "",
        "ec_ids": "",
        "kegg_ids": "",
        "interpro_ids": "",
        "parent_og_id": "",
    }
    return [values.get(column, "") for column in output_columns]


def extract_taxid_from_og(og_id: str) -> str:
    if "at" not in og_id:
        return ""
    suffix = og_id.rsplit("at", 1)[1]
    return suffix if suffix.isdigit() else ""


def _write_base_rows(
    *,
    context: _RowExportContext,
    writer: Any,
    annotation_store: AnnotationStore | None,
) -> ExportStats:
    total_stats = ExportStats()
    pending_parent_ids: list[str] = []

    for level_taxid in context.selected_taxids:

        def write_row(row: list[str]) -> None:
            writer.writerow(row)
            if annotation_store is not None:
                pending_parent_ids.append(row[0])
                if len(pending_parent_ids) >= 5000:
                    flush_annotation_ids(
                        annotation_store=annotation_store,
                        pending_ids=pending_parent_ids,
                    )

        _accumulate_stats(
            total_stats,
            _export_taxid_rows(
                context=context,
                level_taxid=level_taxid,
                write_row=write_row,
            ),
        )

    flush_annotation_ids(
        annotation_store=annotation_store,
        pending_ids=pending_parent_ids,
    )
    return total_stats


def _write_base_rows_parallel(
    *,
    context: _RowExportContext,
    base_tsv_path: Path,
    annotation_store: AnnotationStore | None,
) -> ExportStats:
    chunk_results: dict[str, TaxidChunkResult] = {}
    chunk_prefix = f".chunk.{time.monotonic_ns()}"
    try:
        with ThreadPoolExecutor(max_workers=context.options.workers) as executor:
            futures = {
                executor.submit(
                    _export_taxid_to_chunk,
                    context=context,
                    level_taxid=taxid,
                    chunk_prefix=chunk_prefix,
                    temp_dir=base_tsv_path.parent,
                ): taxid
                for taxid in context.selected_taxids
            }
            for future in as_completed(futures):
                taxid = futures[future]
                chunk_results[taxid] = future.result()

        return _merge_taxid_chunks(
            context=context,
            chunk_results=chunk_results,
            base_tsv_path=base_tsv_path,
            annotation_store=annotation_store,
        )
    finally:
        for orphan_path in base_tsv_path.parent.glob(f"{chunk_prefix}*.tsv"):
            orphan_path.unlink(missing_ok=True)


def _export_taxid_rows(
    *,
    context: _RowExportContext,
    level_taxid: str,
    write_row: Callable[[list[str]], None],
) -> ExportStats:
    stats = ExportStats()
    filters = SearchFilters(
        level_taxid=level_taxid,
        universal=context.options.universal,
        single_copy=context.options.single_copy,
        min_genes=context.options.min_genes,
    )
    for records in context.client.iter_search_records(filters):
        stats.search_pages += 1
        stats.scanned_rows += len(records)
        for record in records:
            og_id = str(record.get("id", "")).strip()
            if not og_id:
                continue
            protein_count_int = to_int(record.get("gene_count"), default=0)
            if (
                context.options.min_genes is not None
                and protein_count_int < context.options.min_genes
            ):
                stats.filtered_rows += 1
                continue

            og_taxid = extract_taxid_from_og(og_id)
            taxon_info = (
                context.taxa_by_id.get(og_taxid)
                or context.taxa_by_id.get(level_taxid)
                or TaxonInfo(name="")
            )
            row = build_output_row_values(
                og_id=og_id,
                og_taxid=og_taxid,
                taxon_name=str(record.get("level_name") or taxon_info.name),
                protein_count=str(record.get("gene_count", "")),
                description=str(record.get("name") or ""),
                rank_by_taxid=context.rank_by_taxid,
                output_columns=context.output_columns,
            )
            write_row(row)
            stats.exported_rows += 1
            stats.exported_taxids.add(level_taxid)
    taxon_name = context.taxa_by_id.get(level_taxid, TaxonInfo(name="")).name
    logger.debug(
        "Exported {} OG records for taxid {} ({}), scanned {}, skipped {}.",
        stats.exported_rows,
        level_taxid,
        taxon_name,
        stats.scanned_rows,
        stats.filtered_rows,
    )
    return stats


def _export_taxid_to_chunk(
    *,
    context: _RowExportContext,
    level_taxid: str,
    chunk_prefix: str,
    temp_dir: Path,
) -> TaxidChunkResult:
    handle = tempfile.NamedTemporaryFile(
        prefix=f"{chunk_prefix}.{level_taxid}.",
        suffix=".tsv",
        dir=str(temp_dir),
        delete=False,
    )
    handle.close()
    chunk_path = Path(handle.name)
    try:
        with chunk_path.open("w", encoding="utf-8", newline="") as chunk_file:
            writer = csv.writer(chunk_file, delimiter="\t", lineterminator="\n")
            stats = _export_taxid_rows(
                context=context,
                level_taxid=level_taxid,
                write_row=writer.writerow,
            )
        return TaxidChunkResult(chunk_path=chunk_path, stats=stats)
    except Exception:
        chunk_path.unlink(missing_ok=True)
        raise


def _merge_taxid_chunks(
    *,
    context: _RowExportContext,
    chunk_results: dict[str, TaxidChunkResult],
    base_tsv_path: Path,
    annotation_store: AnnotationStore | None,
) -> ExportStats:
    merged_stats = ExportStats()
    pending_parent_ids: list[str] = []

    with base_tsv_path.open("w", encoding="utf-8", newline="") as out_file:
        writer = csv.writer(out_file, delimiter="\t", lineterminator="\n")
        writer.writerow(context.output_columns)
        for taxid in context.selected_taxids:
            chunk_result = chunk_results[taxid]
            _accumulate_stats(merged_stats, chunk_result.stats)
            with chunk_result.chunk_path.open(
                "r", encoding="utf-8", newline=""
            ) as chunk_file:
                reader = csv.reader(chunk_file, delimiter="\t")
                for row in reader:
                    if not row:
                        continue
                    writer.writerow(row)
                    if annotation_store is not None:
                        pending_parent_ids.append(row[0])
                        if len(pending_parent_ids) >= 5000:
                            flush_annotation_ids(
                                annotation_store=annotation_store,
                                pending_ids=pending_parent_ids,
                            )
            chunk_result.chunk_path.unlink(missing_ok=True)

    flush_annotation_ids(
        annotation_store=annotation_store, pending_ids=pending_parent_ids
    )
    return merged_stats


def _accumulate_stats(target: ExportStats, source: ExportStats) -> None:
    target.exported_rows += source.exported_rows
    target.filtered_rows += source.filtered_rows
    target.scanned_rows += source.scanned_rows
    target.search_pages += source.search_pages
    target.exported_taxids.update(source.exported_taxids)


def _emit_output_tsv(
    source_path: Path,
    output_path: Path | None,
) -> str:
    if output_path is None:
        with source_path.open("r", encoding="utf-8") as in_file:
            for line in in_file:
                sys.stdout.write(line)
        sys.stdout.flush()
        return "stdout"

    with source_path.open("r", encoding="utf-8") as in_file:
        with output_path.open("w", encoding="utf-8") as out_file:
            for line in in_file:
                out_file.write(line)
    return str(output_path)
