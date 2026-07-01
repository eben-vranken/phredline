import argparse
import sys
import json
from phredline.parser import (
    read_chunks,
    read_lines,
    parse_records,
    decode_phred,
    compute_summary,
    CHUNK_SIZE,
)
from phredline.aggregator import FastqAggregator
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(
        description="Memory-bounded FASTQ quality control and filtering."
    )

    parser.add_argument("input", help="Input FASTQ file (use '-' to read from stdin)")
    parser.add_argument("output", help="Output JSON report file")
    parser.add_argument(
        "--max-read-len",
        type=int,
        default=300,
        help="Maximum read length to track (default: 300)",
    )
    parser.add_argument("--verbose", action="store_true", help="Enable verbose output")

    args = parser.parse_args()

    try:
        aggregator = FastqAggregator(max_read_len=args.max_read_len)

        chunks = read_chunks(args.input, CHUNK_SIZE)
        lines = read_lines(chunks)
        records = parse_records(lines)

        for header, seq, plus, qual in records:
            q_scores = decode_phred(qual)
            aggregator.add_record(seq, q_scores)

        summary = compute_summary(aggregator)

        report = format_for_multiqc(summary, args.input)

        with open(args.output, "w") as f:
            json.dump({"multiqc_data": report}, f, indent=2)

        if args.verbose:
            print(f"Report written to {args.output}")

    except FileNotFoundError as e:
        print(f"Phredline: Input file not found: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Phredline: Error {e}", file=sys.stderr)
        sys.exit(1)


def format_for_multiqc(summary, input_path):
    """Convert summary stats to MultiQC-compatible JSON format"""
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
