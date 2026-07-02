import argparse
import json
import sys
from pathlib import Path

from phredline.aggregator import FastqAggregator
from phredline.parser import (
    CHUNK_SIZE,
    ScratchBuffers,
    compute_summary,
    parse_records,
    read_chunks,
    read_lines,
)


def passes_filters(
    seq: bytes,
    q_scores: list[int] | memoryview | tuple[int, ...],
    min_qual: float,
    min_length: int,
    max_n_fraction: float | None,
) -> bool:
    """Return whether a read satisfies the length, quality, and N-content filters.

    The helper is intentionally lightweight so the CLI can evaluate each record as it
    streams through the parser without allocating additional intermediate structures.
    It rejects reads that are too short, fall below the minimum average Phred score,
    or exceed the configured fraction of ambiguous `N` bases.
    """
    if len(seq) < min_length:
        return False

    if q_scores:
        avg_qual = sum(q_scores) / len(q_scores)
        if avg_qual < min_qual:
            return False

    if max_n_fraction is not None:
        n_count = sum(1 for base in seq if chr(base).upper() == "N")
        n_fraction = n_count / len(seq) if len(seq) > 0 else 1.0
        if n_fraction > max_n_fraction:
            return False

    return True


def main() -> None:
    """Parse command-line arguments, run FASTQ QC, and write MultiQC section files.

    The CLI streams input records through the parser and aggregator, optionally writes
    passing reads to a filtered FASTQ file, and finally emits a directory of
    `_mqc.json` files that MultiQC can discover directly from a directory scan.
    """
    parser = argparse.ArgumentParser(
        description="Memory-bounded FASTQ quality control and filtering."
    )

    parser.add_argument("input", help="Input FASTQ file (use '-' to read from stdin)")
    parser.add_argument(
        "output",
        help="Output directory for MultiQC *_mqc.json section files",
    )
    parser.add_argument(
        "--max-read-len",
        type=int,
        default=300,
        help="Maximum read length to track (default: 300)",
    )
    parser.add_argument("--verbose", action="store_true", help="Enable verbose output")

    parser.add_argument(
        "--min-qual", type=float, default=0, help="Minimum average Phred across read"
    )
    parser.add_argument(
        "--min-length",
        type=int,
        default=0,
        help="Minimum read length after trimming, or raw length if we are not trimming",
    )
    parser.add_argument(
        "--max-n-fraction", type=float, help="Maximum fraction of N bases allowed"
    )
    parser.add_argument(
        "--filtered-fastq",
        type=str,
        help=(
            "If set, write passing reads there; otherwise compute QC and drop "
            "failing reads from the stream."
        ),
    )

    args = parser.parse_args()

    filtered_handle = None

    aggregator = FastqAggregator(max_read_len=args.max_read_len)
    scratch = ScratchBuffers(initial_capacity=args.max_read_len)
    output_dir = Path(args.output)

    try:
        output_dir.mkdir(parents=True, exist_ok=True)

        if args.filtered_fastq:
            filtered_handle = open(args.filtered_fastq, "wb")

        aggregator = FastqAggregator(max_read_len=args.max_read_len)

        chunks = read_chunks(args.input, CHUNK_SIZE)
        lines = read_lines(chunks)
        records = parse_records(lines)

        for header, seq, plus, qual in records:
            q_scores, error_probs, _ = scratch.decode(qual)
            aggregator.add_record(seq, error_probs)

            passes = passes_filters(
                seq,
                q_scores,
                args.min_qual,
                args.min_length,
                args.max_n_fraction,
            )

            if filtered_handle is not None and passes:
                filtered_handle.write(header + b"\n")
                filtered_handle.write(seq + b"\n")
                filtered_handle.write(plus + b"\n")
                filtered_handle.write(qual + b"\n")

        summary = compute_summary(aggregator)

        report_sections = format_for_multiqc(summary, args.input)
        sample_name = Path(args.input).stem

        for section in report_sections:
            section_path = output_dir / f"{sample_name}_{section['id']}_mqc.json"
            with open(section_path, "w") as f:
                json.dump(section, f, indent=2)

        if args.verbose:
            print(f"Report written to {output_dir}")

    except FileNotFoundError as e:
        print(f"Phredline: Input file not found: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Phredline: Error {e}", file=sys.stderr)
        sys.exit(1)

    finally:
        if filtered_handle is not None:
            filtered_handle.close()


def format_for_multiqc(
    summary: dict[str, object], input_path: str
) -> list[dict[str, object]]:
    """Convert summary statistics into the MultiQC section layout.

    The output groups per-position quality values and aggregate sample metrics into the
    structure MultiQC expects, using the input filename stem as the sample identifier.
    """
    sample_name = Path(input_path).stem

    quality_data = {sample_name: {}}
    for stat in summary["per_position_stats"]:
        position = str(stat["position"] + 1)
        quality_data[sample_name][position] = round(stat["mean_quality"], 2)

    quality_section = {
        "id": "per_position_quality",
        "section_name": "Per-Position Quality Scores",
        "description": "Mean Phred quality score at each position in the read.",
        "plot_type": "linegraph",
        "pconfig": {
            "id": "quality_plot",
            "title": "Mean Quality per Position",
            "xlab": "Position (bp)",
            "ylab": "Mean Quality (Phred Score)",
        },
        "data": quality_data,
    }

    summary_section = {
        "id": "fastq_summary",
        "section_name": "Summary Statistics",
        "description": "Overall sequencing metrics.",
        "plot_type": "table",
        "data": {
            sample_name: {
                "total_reads": summary["total_reads"],
                "gc_ratio": round(summary["overall_gc_ratio"], 4),
            }
        },
    }

    return [quality_section, summary_section]


if __name__ == "__main__":
    main()
