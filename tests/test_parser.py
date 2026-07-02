import io
import gzip
import pytest

from phredline.parser import (
    ParseError,
    ScratchBuffers,
    decode_phred,
    read_chunks,
    read_lines,
    parse_records,
)


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


def test_decode_phred_and_scratch_buffers():
    assert list(decode_phred(b"!!!!")) == [0, 0, 0, 0]

    scratch = ScratchBuffers(initial_capacity=4)
    q_scores, error_probs, count = scratch.decode(b"IIII")

    assert count == 4
    assert list(q_scores) == [40, 40, 40, 40]
    assert list(error_probs) == [0.0001, 0.0001, 0.0001, 0.0001]
