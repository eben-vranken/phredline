import gzip
import sys
import math

CHUNK_SIZE = 16 * 1024


class ParseError(Exception):
    """Error for when FATQ format is malformed"""

    pass


def read_chunks(file_path, chunk_size):
    if file_path is None or file_path == "-":
        stream = sys.stdin.buffer
        while True:
            bchunk = stream.read(chunk_size)
            if bchunk == b"":
                break
            yield bchunk
        return

    with open(file_path, "rb") as raw:
        magic = raw.read(2)
        raw.seek(0)

        if magic == b"\x1f\x8b":
            stream = gzip.open(raw, "rb")
        else:
            stream = raw

        while True:
            chunk = stream.read(chunk_size)
            if chunk == b"":
                break
            yield chunk


def read_lines(chunk_stream):
    tail = b""

    for chunk in chunk_stream:
        data = tail + chunk
        parts = data.split(b"\n")

        for split in parts[:-1]:
            yield split

        tail = parts[-1]

    if tail:
        yield tail


def parse_records(line_stream):
    record = []

    for line in line_stream:
        record.append(line)

        if len(record) == 4:
            header, seq, plus, qual = record

            if not header.startswith(b"@"):
                raise ParseError("FASTQ record does not start with @")
            if not plus.startswith(b"+"):
                raise ParseError("FASTQ record seperator does not start with +")
            if len(seq) != len(qual):
                raise ParseError("Sequence and quality lengths do not match")

            yield header, seq, plus, qual
            record = []

    if record:
        raise ParseError("Truncated FASTQ record at end of file")


def decode_phred(quality_string):
    return [byte - 33 for byte in quality_string]


def decode_error_probability(q_array):
    return [10 ** (-q / 10) for q in q_array]


def compute_summary(aggregator):
    """Computes final stats from the Aggregator's state"""
    summary = {
        "per_position_stats": [],
        "overall_gc_ratio": (
            aggregator.total_gc / aggregator.total_length
            if aggregator.total_length > 0
            else 0
        ),
        "total_reads": aggregator.total_reads,
    }

    for i in range(aggregator.max_read_len):
        count = aggregator.read_counts[i]
        if count == 0:
            break

        mean_error_prob = aggregator.error_prob_sums[i] / count
        mean_quality = (
            -10 * math.log10(mean_error_prob) if mean_error_prob > 0 else float("inf")
        )

        summary["per_position_stats"].append(
            {
                "position": i,
                "mean_quality": mean_quality,
                "frequencies": {
                    base: aggregator.base_counts[base][i] / count
                    for base in ["A", "T", "C", "G", "N"]
                },
            }
        )

    return summary
