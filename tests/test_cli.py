import json
import sys

from phredline import cli


def test_passes_filters_helper():
    assert cli.passes_filters(b"ACGT", [40, 40, 40, 40], 30, 4, 0.5) is True
    assert cli.passes_filters(b"ACNT", [40, 40, 40, 40], 30, 4, 0.1) is False
    assert cli.passes_filters(b"ACGT", [10, 10, 10, 10], 20, 4, None) is False
    assert cli.passes_filters(b"ACGT", [40, 40, 40, 40], 30, 5, None) is False


def test_main_writes_only_passing_reads_to_filtered_fastq(tmp_path, monkeypatch):
    input_path = tmp_path / "input.fastq"
    output_path = tmp_path / "report.json"
    filtered_path = tmp_path / "filtered.fastq"

    records = [
        (b"@pass", b"ACGT", b"+", b"IIII"),
        (b"@fail", b"ACNT", b"+", b"!!!!"),
    ]

    monkeypatch.setattr(cli, "read_chunks", lambda *args, **kwargs: [b""])
    monkeypatch.setattr(cli, "read_lines", lambda chunks: [])
    monkeypatch.setattr(cli, "parse_records", lambda lines: iter(records))

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "phredline",
            str(input_path),
            str(output_path),
            "--min-qual",
            "30",
            "--max-n-fraction",
            "0.1",
            "--filtered-fastq",
            str(filtered_path),
        ],
    )

    cli.main()

    report = json.loads(output_path.read_text())
    summary_section = next(
        section
        for section in report["multiqc_data"]
        if section["id"] == "fastq_summary"
    )
    assert summary_section["data"][input_path.stem]["total_reads"] == 2

    filtered = filtered_path.read_text().splitlines()
    assert filtered == ["@pass", "ACGT", "+", "IIII"]
