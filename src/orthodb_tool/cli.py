from __future__ import annotations

import os
import sys
from functools import wraps
from pathlib import Path
from typing import Any, Callable

import rich_click as click
from loguru import logger

from orthodb_tool.api import (
    ApiRequestError,
    ApiResponseError,
    HttpConfig,
    OrthodbClient,
)
from orthodb_tool.counts import CountOptions, run_counts, write_counts_tsv
from orthodb_tool.fasta import FastaDownloadOptions, download_fasta
from orthodb_tool.records import RecordExportResult, RecordsOptions, run_records

ALLOWED_FILTER_CUTOFFS = (1.0, 0.9, 0.8)


class CutoffParamType(click.ParamType):
    name = "cutoff"

    def convert(
        self,
        value: Any,
        param: click.Parameter | None,
        ctx: click.Context | None,
    ) -> float:
        try:
            numeric = float(value)
        except (TypeError, ValueError):
            self.fail(f"{value!r} is not a valid float", param, ctx)

        for cutoff in ALLOWED_FILTER_CUTOFFS:
            if abs(numeric - cutoff) < 1e-9:
                return cutoff

        allowed_values = ", ".join(f"{cutoff:.1f}" for cutoff in ALLOWED_FILTER_CUTOFFS)
        self.fail(f"must be one of: {allowed_values}", param, ctx)


CUTOFF_TYPE = CutoffParamType()
LOG_FORMAT = (
    "<dim>{time:YYYY-MM-DD HH:mm:ss.SSS}</dim> | <level>{level: <8}</level> | {message}"
)


def _configure_logging(*, quiet: bool = False, verbose: bool = False) -> None:
    level = "ERROR" if quiet else "DEBUG" if verbose else "INFO"
    logger.remove()
    logger.enable("orthodb_tool")
    logger.add(
        sys.stderr,
        level=level,
        format=LOG_FORMAT,
    )


def _format_cutoff(value: float) -> str:
    return f"{value:.1f}"


def _format_filters(
    *,
    universal: float | None,
    single_copy: float | None,
    min_genes: int | None,
) -> str:
    filters: list[str] = []
    if universal is not None:
        filters.append(f"universal ≥ {_format_cutoff(universal)}")
    if single_copy is not None:
        filters.append(f"single-copy ≥ {_format_cutoff(single_copy)}")
    if min_genes is not None:
        filters.append(f"gene_count ≥ {min_genes}")
    return ", ".join(filters) if filters else "no filters"


def _format_output_target(output_path: Path | None) -> str:
    return str(output_path) if output_path is not None else "stdout"


def _format_worker_detail(workers: int) -> str:
    return f", using {workers} workers" if workers > 1 else ""


def _count_taxid_tokens(raw: str | None) -> int:
    if not raw:
        return 0
    return len([token for token in (part.strip() for part in raw.split(",")) if token])


def _format_records_selection(
    *,
    taxids: str | None,
    taxids_file: Path | None,
    include_descendant_taxids: bool,
) -> str:
    parts: list[str] = []
    taxid_count = _count_taxid_tokens(taxids)
    if taxid_count:
        parts.append(f"{taxid_count} requested taxids")
    if taxids_file is not None:
        parts.append(f"taxids from {taxids_file}")
    if not parts:
        parts.append("all taxids")
    if include_descendant_taxids:
        parts.append("including descendants")
    return ", ".join(parts)


def _format_requested_columns(
    *,
    with_taxid_rank: bool,
    with_parent_og: bool,
    with_annotation_ids: bool,
) -> str:
    columns: list[str] = []
    if with_taxid_rank:
        columns.append("taxonomy ranks")
    if with_parent_og:
        columns.append("parent OG IDs")
    if with_annotation_ids:
        columns.append("annotation IDs")
    return f", including {', '.join(columns)}" if columns else ""


def _format_enrichment_status(
    *,
    requested: bool,
    rows_exported: int,
    applied: bool,
) -> str:
    if not requested:
        return "not requested"
    if rows_exported == 0:
        return "skipped because no records matched"
    return "added" if applied else "not added"


def _validate_output_path(output_path: Path | None) -> None:
    if output_path is None:
        return
    if output_path.parent and not output_path.parent.exists():
        raise click.ClickException(
            f"Output directory does not exist: {output_path.parent}"
        )


def _http_options(func: Callable[..., Any]) -> Callable[..., Any]:
    for option in reversed(
        [
            click.option(
                "--timeout-sec",
                default=30.0,
                show_default=True,
                type=click.FloatRange(min=0.1),
                help="HTTP timeout in seconds.",
            ),
            click.option(
                "--retries",
                default=4,
                show_default=True,
                type=click.IntRange(min=0),
                help="Retry count for transient HTTP/network errors.",
            ),
            click.option(
                "--retry-delay-base-sec",
                default=1.0,
                show_default=True,
                type=click.FloatRange(min=0.0),
                help="Base delay in seconds between retry attempts.",
            ),
            click.option(
                "--request-interval-sec",
                default=0.0,
                show_default=True,
                type=click.FloatRange(min=0.0),
                help="Minimum seconds between API requests.",
            ),
        ]
    ):
        func = option(func)
    return func


def _common_worker_options(func: Callable[..., Any]) -> Callable[..., Any]:
    for option in reversed(
        [
            click.option(
                "--workers",
                default=6,
                show_default=True,
                type=click.IntRange(min=1),
                help="Number of worker threads.",
            ),
            click.option(
                "--verbose", is_flag=True, default=False, help="Print detailed logs."
            ),
            click.option(
                "--quiet", is_flag=True, default=False, help="Hide progress logs."
            ),
        ]
    ):
        func = option(func)
    return func


def _with_client(func: Callable[..., Any]) -> Callable[..., Any]:
    @wraps(func)
    def wrapper(*args: Any, **kwargs: Any) -> Any:
        verbose = bool(kwargs.pop("verbose", False))
        quiet = bool(kwargs.pop("quiet", False))
        _configure_logging(quiet=quiet, verbose=verbose)
        try:
            config = HttpConfig(
                timeout_sec=kwargs.pop("timeout_sec"),
                retries=kwargs.pop("retries"),
                retry_delay_base_sec=kwargs.pop("retry_delay_base_sec"),
                request_interval_sec=kwargs.pop("request_interval_sec"),
            )
            with OrthodbClient(config) as client:
                return func(*args, client=client, **kwargs)
        except BrokenPipeError:
            raise
        except (ApiRequestError, ApiResponseError, ValueError, OSError) as exc:
            raise click.ClickException(str(exc)) from exc

    return wrapper


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
def cli() -> None:
    """CLI toolkit for querying OrthoDB."""


@cli.command("counts")
@click.option(
    "--parent-taxid",
    required=True,
    type=click.IntRange(min=1),
    help="Parent NCBI taxid used to define the descendant subtree.",
)
@click.option(
    "--target-rank",
    required=True,
    type=str,
    help="Target taxonomic rank to evaluate.",
)
@click.option(
    "--universal",
    default=None,
    type=CUTOFF_TYPE,
    help="Universality cutoff; allowed values: 1.0, 0.9, 0.8.",
)
@click.option(
    "--single-copy",
    "single_copy",
    default=None,
    type=CUTOFF_TYPE,
    help="Single-copy cutoff; allowed values: 1.0, 0.9, 0.8.",
)
@click.option(
    "--min-genes",
    default=None,
    type=click.IntRange(min=1),
    help="Only count OGs with gene_count >= this value.",
)
@click.option(
    "--output",
    "output_path",
    default=None,
    type=click.Path(path_type=Path, dir_okay=False),
    help="Output TSV path; defaults to stdout.",
)
@_http_options
@_common_worker_options
@_with_client
def counts_command(
    *,
    client: OrthodbClient,
    parent_taxid: int,
    target_rank: str,
    universal: float | None,
    single_copy: float | None,
    min_genes: int | None,
    output_path: Path | None,
    workers: int,
) -> None:
    _validate_output_path(output_path)
    target = _format_output_target(output_path)
    logger.info(
        "Counting OGs below parent taxid {} at rank {} ({}, output: {}{}).",
        parent_taxid,
        target_rank,
        _format_filters(
            universal=universal,
            single_copy=single_copy,
            min_genes=min_genes,
        ),
        target,
        _format_worker_detail(workers),
    )
    rows = run_counts(
        client,
        CountOptions(
            parent_taxid=parent_taxid,
            target_rank=target_rank,
            universal=universal,
            single_copy=single_copy,
            min_genes=min_genes,
            workers=workers,
        ),
    )
    logger.info("Writing count table to {}…", target)
    output_label = write_counts_tsv(rows, output_path)
    logger.info("Wrote {} taxid rows to {}.", len(rows), output_label)


@cli.command("records")
@click.option(
    "--output",
    "output_path",
    default=None,
    type=click.Path(path_type=Path, dir_okay=False),
    help="Output TSV path; defaults to stdout.",
)
@click.option(
    "--taxids",
    default=None,
    help="Comma-separated list of numeric NCBI taxids to include.",
)
@click.option(
    "--taxids-file",
    default=None,
    type=click.Path(path_type=Path, dir_okay=False, exists=True, readable=True),
    help="File with one numeric NCBI taxid per line.",
)
@click.option(
    "--include-descendant-taxids",
    is_flag=True,
    default=False,
    help="Include all descendant taxonomic-level taxids for selected taxids.",
)
@click.option(
    "--min-genes",
    default=None,
    type=click.IntRange(min=1),
    help="Only export OGs with gene_count >= this value.",
)
@click.option(
    "--universal",
    default=None,
    type=CUTOFF_TYPE,
    help="Universality cutoff; allowed values: 1.0, 0.9, 0.8.",
)
@click.option(
    "--single-copy",
    "single_copy",
    default=None,
    type=CUTOFF_TYPE,
    help="Single-copy cutoff; allowed values: 1.0, 0.9, 0.8.",
)
@click.option(
    "--with-parent-og",
    is_flag=True,
    default=False,
    help="Fill parent_og_id from OG hierarchy pairs.",
)
@click.option(
    "--with-annotation-ids",
    is_flag=True,
    default=False,
    help="Retrieve and include OG annotation ID columns.",
)
@click.option(
    "--with-taxid-rank",
    is_flag=True,
    default=False,
    help="Enable NCBI taxonomy rank lookup to populate taxid_rank.",
)
@_http_options
@_common_worker_options
@_with_client
def records_command(
    *,
    client: OrthodbClient,
    output_path: Path | None,
    taxids: str | None,
    taxids_file: Path | None,
    include_descendant_taxids: bool,
    min_genes: int | None,
    universal: float | None,
    single_copy: float | None,
    workers: int,
    with_parent_og: bool,
    with_annotation_ids: bool,
    with_taxid_rank: bool,
) -> None:
    _validate_output_path(output_path)
    target = _format_output_target(output_path)
    logger.info(
        "Exporting OG records ({}, {}, output: {}{}{}).",
        _format_records_selection(
            taxids=taxids,
            taxids_file=taxids_file,
            include_descendant_taxids=include_descendant_taxids,
        ),
        _format_filters(
            universal=universal,
            single_copy=single_copy,
            min_genes=min_genes,
        ),
        target,
        _format_worker_detail(workers),
        _format_requested_columns(
            with_taxid_rank=with_taxid_rank,
            with_parent_og=with_parent_og,
            with_annotation_ids=with_annotation_ids,
        ),
    )
    result: RecordExportResult = run_records(
        client,
        RecordsOptions(
            output_path=output_path,
            taxids=taxids,
            taxids_file=taxids_file,
            include_descendant_taxids=include_descendant_taxids,
            min_genes=min_genes,
            universal=universal,
            single_copy=single_copy,
            workers=workers,
            with_parent_og=with_parent_og,
            with_annotation_ids=with_annotation_ids,
            with_taxid_rank=with_taxid_rank,
        ),
    )
    parent_status = _format_enrichment_status(
        requested=with_parent_og,
        rows_exported=result.stats.exported_rows,
        applied=result.parent_mapping_applied,
    )
    annotation_status = _format_enrichment_status(
        requested=with_annotation_ids,
        rows_exported=result.stats.exported_rows,
        applied=result.annotation_mapping_applied,
    )
    logger.info(
        "Wrote {} OG records to {}. Scanned {} records across {} search pages. "
        "Skipped {} below the gene-count filter. Parent OG IDs: {}. "
        "Annotation IDs: {}.",
        result.stats.exported_rows,
        result.output_label,
        result.stats.scanned_rows,
        result.stats.search_pages,
        result.stats.filtered_rows,
        parent_status,
        annotation_status,
    )


@cli.command("download-fasta")
@click.argument("og_id", type=str)
@click.option(
    "--cds",
    is_flag=True,
    default=False,
    help="Download CDS sequences instead of protein sequences.",
)
@click.option(
    "--output",
    "output_path",
    default=None,
    type=click.Path(path_type=Path, dir_okay=False),
    help="Output FASTA path; defaults to stdout.",
)
@_http_options
@click.option("--verbose", is_flag=True, default=False, help="Print detailed logs.")
@click.option("--quiet", is_flag=True, default=False, help="Hide progress logs.")
@_with_client
def download_fasta_command(
    *,
    client: OrthodbClient,
    og_id: str,
    cds: bool,
    output_path: Path | None,
) -> None:
    _validate_output_path(output_path)
    target = _format_output_target(output_path)
    logger.info(
        "Downloading {} FASTA for {} to {}…",
        "CDS" if cds else "protein",
        og_id,
        target,
    )
    download_fasta(
        client,
        FastaDownloadOptions(
            og_id=og_id,
            cds=cds,
            output_path=output_path,
        ),
    )


def _silence_broken_pipe() -> None:
    try:
        devnull = os.open(os.devnull, os.O_WRONLY)
        try:
            os.dup2(devnull, sys.stdout.fileno())
        finally:
            os.close(devnull)
    except OSError:
        return


def main(argv: list[str] | None = None) -> int:
    try:
        cli.main(args=argv, standalone_mode=False)
        return 0
    except BrokenPipeError:
        _silence_broken_pipe()
        return 0
    except click.ClickException as exc:
        _configure_logging()
        logger.error("Error: {}", exc.format_message())
        return exc.exit_code
    except click.Abort:
        _configure_logging()
        logger.error("Aborted!")
        return 1
