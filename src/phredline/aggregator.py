import array


class FastqAggregator:
    def __init__(self, max_read_len=300):
        self.max_read_len = max_read_len

        self.base_counts = {
            "A": array.array("I", [0] * max_read_len),
            "T": array.array("I", [0] * max_read_len),
            "C": array.array("I", [0] * max_read_len),
            "G": array.array("I", [0] * max_read_len),
            "N": array.array("I", [0] * max_read_len),
        }
        self.qual_sums = array.array("d", [0.0] * max_read_len)
        self.read_counts = array.array("I", [0] * max_read_len)
        self.total_gc = 0
        self.total_length = 0

    def add_record(self, sequence, qualities):
        """Updates internal arrays in-place using the record data."""
        seq_len = min(len(sequence), self.max_read_len)

        for i in range(seq_len):
            base = chr(sequence[i]).upper()
            if base in self.base_counts:
                self.base_counts[base][i] += 1

            self.qual_sums[i] += qualities[i]
            self.read_counts[i] += 1

            if base in ("G", "C"):
                self.total_gc += 1

        self.total_length += seq_len
