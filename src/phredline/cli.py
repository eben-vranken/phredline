import argparse
import sys

def main():
    parser = argparse.ArgumentParser(
        description="Memory-bounded FASTQ quality control and filtering."
    )

    parser.add_argument(
        "input",
        help="Input FASTQ file (use '-' to read from stdin)"
    )
    parser.add_argument(
        "output",
        help="Output JSON report file"
    )
    parser.add_argument(
        "--max-read-len",
        type=int,
        default=300,
        help="Maximum read length to track (default: 300)"
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose output"
    )

    args = parser.parse_args()

if __name__ == "__main__":
    main()