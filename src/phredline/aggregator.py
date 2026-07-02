import array


class FastqAggregator:
    def __init__(self, max_read_len: int = 300) -> None:
        """Initialize per-position counters and running totals for streamed FASTQ reads.

        The aggregator keeps fixed-width arrays for base counts, error probability sums,
        and read coverage so it can accumulate statistics in constant memory relative to
        the configured maximum read length.
        """
        self.max_read_len = max_read_len

        self.base_counts = {
            "A": array.array("I", [0] * max_read_len),
            "T": array.array("I", [0] * max_read_len),
            "C": array.array("I", [0] * max_read_len),
            "G": array.array("I", [0] * max_read_len),
            "N": array.array("I", [0] * max_read_len),
        }
        self.error_prob_sums = array.array("d", [0.0] * max_read_len)
        self.read_counts = array.array("I", [0] * max_read_len)
        self.total_gc = 0
        self.total_length = 0
        self.total_reads = 0

    def add_record(self, sequence: bytes, error_probabilities) -> None:
        """Update the running statistics in place for a single FASTQ record.

        Bases are counted per position, error probabilities are accumulated alongside
        coverage, and the overall GC and length totals are advanced using the sequence
        prefix that fits within the configured maximum read length.
        """
        self.total_reads += 1
        seq_len = min(len(sequence), self.max_read_len)

        for i in range(seq_len):
            base = chr(sequence[i]).upper()
            if base in self.base_counts:
                self.base_counts[base][i] += 1

            self.error_prob_sums[i] += error_probabilities[i]
            self.read_counts[i] += 1

            if base in ("G", "C"):
                self.total_gc += 1

        self.total_length += seq_len
