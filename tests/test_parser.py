import gzip
import math

import pytest
from hypothesis import given, settings, strategies as st

from phredline.aggregator import FastqAggregator
from phredline.parser import (
    ParseError,
    ScratchBuffers,
    compute_summary,
    decode_error_probability,
    decode_phred,
    read_chunks,
    read_lines,
    parse_records,
)


DNA_BASES = [b"A", b"C", b"G", b"T", b"N"]


def sequence_bytes_strategy(size: int):
    return st.lists(st.sampled_from(DNA_BASES), min_size=size, max_size=size).map(
        b"".join
    )


def quality_bytes_strategy(size: int):
    return st.lists(
        st.integers(min_value=33, max_value=73), min_size=size, max_size=size
    ).map(bytes)


def valid_fastq_records_strategy():
    read_length = st.integers(min_value=0, max_value=40)

    return st.lists(
        read_length.flatmap(
            lambda size: st.tuples(
                sequence_bytes_strategy(size),
                quality_bytes_strategy(size),
            )
        ),
        min_size=0,
        max_size=25,
    )


def malformed_fastq_bytes_strategy():
    good_size = st.integers(min_value=0, max_value=20)

    bad_header = good_size.flatmap(
        lambda size: st.tuples(
            sequence_bytes_strategy(size),
            quality_bytes_strategy(size),
        ).map(lambda pair: b"r1\n" + pair[0] + b"\n+\n" + pair[1] + b"\n")
    )

    bad_plus = good_size.flatmap(
        lambda size: st.tuples(
            sequence_bytes_strategy(size),
            quality_bytes_strategy(size),
        ).map(lambda pair: b"@r1\n" + pair[0] + b"\n-\n" + pair[1] + b"\n")
    )

    length_mismatch = good_size.flatmap(
        lambda size: st.tuples(
            sequence_bytes_strategy(size),
            quality_bytes_strategy(size + 1),
        ).map(lambda pair: b"@r1\n" + pair[0] + b"\n+\n" + pair[1] + b"\n")
    )

    truncated = st.just(b"@r1\nACGT\n+\n")

    return st.one_of(bad_header, bad_plus, length_mismatch, truncated)


def test_clean_record():
    data = b"@r1\nACGT\n+\n!!!!\n"
    records = list(parse_records(read_lines([data])))
    assert records == [(b"@r1", b"ACGT", b"+", b"!!!!")]


def test_split_exactly_on_chunk_boundary():
    chunks = [b"@r1\nAC", b"GT\n+\n!!", b"!!\n"]
    records = list(parse_records(read_lines(chunks)))
    assert records == [(b"@r1", b"ACGT", b"+", b"!!!!")]


def test_empty_file():
    records = list(parse_records(read_lines([])))
    assert records == []


def test_truncated_final_record_raises():
    data = b"@r1\nACGT\n+\n"
    with pytest.raises(ParseError):
        list(parse_records(read_lines([data])))


def test_malformed_seperator_raises():
    data = b"@r1\nACGT\n-\n!!!!\n"
    with pytest.raises(ParseError):
        list(parse_records(read_lines([data])))


def test_gzip_input(tmp_path):
    path = tmp_path / "sample.fastq.gz"
    with gzip.open(
        path,
        "wb",
    ) as f:
        f.write(b"@r1\nACGT\n+\n!!!!\n")

    chunks = read_chunks(str(path), chunk_size=4)
    records = list(parse_records(read_lines(chunks)))
    assert records == [(b"@r1", b"ACGT", b"+", b"!!!!")]


def test_plain_input_file(tmp_path):
    path = tmp_path / "sample.fastq"
    path.write_bytes(b"@r1\nACGT\n+\n!!!!\n")

    chunks = read_chunks(str(path), chunk_size=4)
    records = list(parse_records(read_lines(chunks)))
    assert records == [(b"@r1", b"ACGT", b"+", b"!!!!")]


def test_read_chunks_from_stdin(monkeypatch):
    class DummyBuffer:
        def __init__(self, data: bytes):
            self._data = data
            self._offset = 0

        def read(self, chunk_size: int) -> bytes:
            chunk = self._data[self._offset : self._offset + chunk_size]
            self._offset += len(chunk)
            return chunk

    class DummyStdin:
        def __init__(self, data: bytes):
            self.buffer = DummyBuffer(data)

    monkeypatch.setattr(
        "phredline.parser.sys.stdin", DummyStdin(b"@r1\nACGT\n+\n!!!!\n")
    )

    chunks = read_chunks("-", chunk_size=5)
    records = list(parse_records(read_lines(chunks)))
    assert records == [(b"@r1", b"ACGT", b"+", b"!!!!")]


def test_read_lines_preserves_final_unterminated_chunk():
    assert list(read_lines([b"@r1\nAC", b"GT\n+\n!!!!"])) == [
        b"@r1",
        b"ACGT",
        b"+",
        b"!!!!",
    ]


def test_decode_phred_and_scratch_buffers():
    assert list(decode_phred(b"!!!!")) == [0, 0, 0, 0]
    error_probs = decode_error_probability([40])
    assert error_probs[0] == pytest.approx(0.0001)

    assert list(decode_phred(b"I")) == [40]

    scratch = ScratchBuffers(initial_capacity=4)
    q_scores, error_probs, count = scratch.decode(b"IIII")

    assert count == 4
    assert list(q_scores) == [40, 40, 40, 40]
    assert list(error_probs) == [0.0001, 0.0001, 0.0001, 0.0001]


def test_malformed_header_line_raises_parse_error():
    data = b"r1\nACGT\n+\n!!!!\n"

    with pytest.raises(ParseError):
        list(parse_records(read_lines([data])))


def test_scratch_buffers_grow_for_longer_reads():
    scratch = ScratchBuffers(initial_capacity=2)

    short_q_scores, short_error_probs, short_count = scratch.decode(b"II")
    long_q_scores, long_error_probs, long_count = scratch.decode(b"IIIIII")

    assert short_count == 2
    assert list(short_q_scores) == [40, 40]
    assert list(short_error_probs) == [0.0001, 0.0001]
    assert long_count == 6
    assert list(long_q_scores) == [40, 40, 40, 40, 40, 40]
    assert list(long_error_probs) == [0.0001] * 6


def test_fastq_aggregator_add_record_matches_manual_counts():
    aggregator = FastqAggregator(max_read_len=4)

    aggregator.add_record(b"ACGT", [0.10, 0.20, 0.30, 0.40])
    aggregator.add_record(b"AGGT", [0.50, 0.60, 0.70, 0.80])

    summary = compute_summary(aggregator)

    assert aggregator.total_reads == 2
    assert aggregator.total_length == 8
    assert aggregator.total_gc == 4

    assert aggregator.read_counts.tolist() == [2, 2, 2, 2]
    assert aggregator.base_counts["A"].tolist() == [2, 0, 0, 0]
    assert aggregator.base_counts["C"].tolist() == [0, 1, 0, 0]
    assert aggregator.base_counts["G"].tolist() == [0, 1, 2, 0]
    assert aggregator.base_counts["T"].tolist() == [0, 0, 0, 2]

    first_position = summary["per_position_stats"][0]
    second_position = summary["per_position_stats"][1]

    assert first_position["frequencies"] == {
        "A": 1.0,
        "T": 0.0,
        "C": 0.0,
        "G": 0.0,
        "N": 0.0,
    }
    assert second_position["frequencies"] == {
        "A": 0.0,
        "T": 0.0,
        "C": 0.5,
        "G": 0.5,
        "N": 0.0,
    }

    expected_first_mean_quality = -10 * math.log10((0.10 + 0.50) / 2)
    expected_second_mean_quality = -10 * math.log10((0.20 + 0.60) / 2)

    assert first_position["mean_quality"] == pytest.approx(expected_first_mean_quality)
    assert second_position["mean_quality"] == pytest.approx(
        expected_second_mean_quality
    )


@given(valid_fastq_records_strategy())
@settings(max_examples=50)
def test_valid_fastq_parser_and_aggregator_preserve_counts(records):
    lines = []
    expected_records = []
    total_gc = 0
    total_length = 0

    for index, (sequence, quality) in enumerate(records):
        header = f"@read{index}".encode()
        expected_records.append((header, sequence, b"+", quality))
        lines.extend([header, sequence, b"+", quality])
        total_gc += sum(base in (ord("G"), ord("C")) for base in sequence)
        total_length += len(sequence)

    parsed_records = list(
        parse_records(read_lines([b"\n".join(lines) + (b"\n" if lines else b"")]))
    )
    assert parsed_records == expected_records

    max_read_len = max((len(sequence) for sequence, _ in records), default=0)
    aggregator = FastqAggregator(max_read_len=max_read_len or 1)
    for _, sequence, _, quality in parsed_records:
        q_scores = [10 ** (-(byte - 33) / 10) for byte in quality]
        aggregator.add_record(sequence, q_scores)

    assert aggregator.total_reads == len(records)
    assert aggregator.total_gc == total_gc
    assert aggregator.total_length == total_length

    for position in range(max_read_len):
        assert (
            sum(aggregator.base_counts[base][position] for base in "ATCGN")
            == aggregator.read_counts[position]
        )


@given(malformed_fastq_bytes_strategy())
@settings(max_examples=50)
def test_malformed_fastq_input_always_raises_parse_error(raw_fastq):
    with pytest.raises(ParseError):
        list(parse_records(read_lines([raw_fastq])))
