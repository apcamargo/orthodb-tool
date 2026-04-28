from __future__ import annotations

import sys
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path

from loguru import logger

from orthodb_tool.api import OrthodbClient


@dataclass(frozen=True)
class FastaDownloadOptions:
    og_id: str
    cds: bool
    output_path: Path | None


def normalize_fasta_target_id(target_id: str) -> str:
    normalized = target_id.strip()
    if not normalized:
        raise ValueError("A non-empty OrthoDB target ID is required.")
    return normalized


def download_fasta(
    client: OrthodbClient,
    options: FastaDownloadOptions,
) -> str:
    og_id = normalize_fasta_target_id(options.og_id)
    if og_id != options.og_id:
        logger.debug(
            "Normalized FASTA target ID from {!r} to {!r}.",
            options.og_id,
            og_id,
        )
    sequence_type = "CDS" if options.cds else "protein"
    output_label = (
        str(options.output_path) if options.output_path is not None else "stdout"
    )
    logger.debug("Requesting {} sequences from OrthoDB…", sequence_type)
    line_count = 0
    sequence_count = 0
    if options.output_path is None:
        with client.stream_fasta(og_id, cds=options.cds) as stream:
            for chunk in _iter_nonempty_lines(stream):
                line_count += 1
                if chunk.lstrip().startswith(b">"):
                    sequence_count += 1
                sys.stdout.buffer.write(chunk)
        sys.stdout.flush()
        logger.debug("Received {} non-empty FASTA lines for {}.", line_count, og_id)
        logger.info(
            "Wrote {} {} sequences for {} to {}.",
            sequence_count,
            sequence_type,
            og_id,
            output_label,
        )
        return output_label

    with client.stream_fasta(og_id, cds=options.cds) as stream:
        with options.output_path.open("wb") as out_file:
            for chunk in _iter_nonempty_lines(stream):
                line_count += 1
                if chunk.lstrip().startswith(b">"):
                    sequence_count += 1
                out_file.write(chunk)
    logger.debug("Received {} non-empty FASTA lines for {}.", line_count, og_id)
    logger.info(
        "Wrote {} {} sequences for {} to {}.",
        sequence_count,
        sequence_type,
        og_id,
        output_label,
    )
    return output_label


def _iter_nonempty_lines(chunks: Iterator[bytes]) -> Iterator[bytes]:
    pending = b""
    for chunk in chunks:
        if not chunk:
            continue
        pending += chunk
        parts = pending.splitlines(keepends=True)
        if parts and not parts[-1].endswith((b"\n", b"\r")):
            pending = parts.pop()
        else:
            pending = b""
        for line in parts:
            if line.strip():
                yield line

    if pending and pending.strip():
        yield pending
