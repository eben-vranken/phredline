import gzip
import sys
import math
import array
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

def decode_phred(quality_string: bytes) -> array.array:
    out = array.array("b", bytes(len(quality_string)))
    decode_phred_into(quality_string, out)
    return out

def decode_phred_into(q_buffer, count: int, out_buffer: array.array) -> int:
    for i in range(count):
        out_buffer[i] = 10 ** (-q_buffer[i] / 10)
    return count

def decode_error_probability(q_array) -> array.array:
    q_array = list(q_array)
    out = array.array("d", [0.0] * len(q_array))
    decode_error_probability_into(q_array, len(q_array), out)
    return out

def decode_error_probability_into(q_buffer, count: int, out_buffer: array.array) -> int:
    for i in range(count):
        out_buffer[i] = 10 ** (-q_buffer[i] / 10)
    return count

class ScratchBuffers:
    """Reusable, growable buffers for per-record Phred decoding.

    q_scores/error_probs returned by decode() are zero-copy memoryview
    slices into shared internal buffers. These are valid only until the NEXT call
    to decode(). Consume them before decoding the next record.
    """
    def __init__(self, initial_capacity: int = 300):
        self._q_scores = array.array("b", bytes(initial_capacity))
        self._error_probs = array.array("d", [0.0] * initial_capacity)

    def _ensure_capacity(self, n: int) -> None:
        if len(self._q_scores) < n:
            self._q_scores = array.array("b", bytes(n))
            self._error_probs = array.array("d", [0.0] * n)

    def decode(self, quality_string: bytes):
        n = len(quality_string)
        self._ensure_capacity(n)
        decode_phred_into(quality_string, self._q_scores)
        decode_error_probability_into(self._q_scores, n, self._error_probs)
        return memoryview(self._q_scores)[:n], memoryview(self._error_probs)[:n], n

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
