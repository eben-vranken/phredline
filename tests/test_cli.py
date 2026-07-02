import json
import sys

import pytest

from phredline import cli


def test_passes_filters_helper():
    assert cli.passes_filters(b"ACGT", [40, 40, 40, 40], 30, 4, 0.5) is True
    assert cli.passes_filters(b"ACNT", [40, 40, 40, 40], 30, 4, 0.1) is False
    assert cli.passes_filters(b"ACGT", [10, 10, 10, 10], 20, 4, None) is False
    assert cli.passes_filters(b"ACGT", [40, 40, 40, 40], 30, 5, None) is False


def test_main_writes_multiqc_section_files_in_output_directory(tmp_path, monkeypatch):
    input_path = tmp_path / "input.fastq"
    output_dir = tmp_path / "report"
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
            str(output_dir),
            "--min-qual",
            "30",
            "--max-n-fraction",
            "0.1",
            "--verbose",
            "--filtered-fastq",
            str(filtered_path),
        ],
    )

    cli.main()

    assert output_dir.is_dir()

    report_files = sorted(output_dir.glob("*_mqc.json"))
    assert len(report_files) == 2

    parsed_sections = [json.loads(path.read_text()) for path in report_files]
    assert {section["id"] for section in parsed_sections} == {
        "per_position_quality",
        "fastq_summary",
    }

    summary_section = next(
        section for section in parsed_sections if section["id"] == "fastq_summary"
    )
    assert summary_section["data"][input_path.stem]["total_reads"] == 2

    filtered = filtered_path.read_text().splitlines()
    assert filtered == ["@pass", "ACGT", "+", "IIII"]


def test_main_reports_unexpected_error(tmp_path, monkeypatch, capsys):
    input_path = tmp_path / "input.fastq"
    output_dir = tmp_path / "report"

    monkeypatch.setattr(cli, "read_chunks", lambda *args, **kwargs: [b""])
    monkeypatch.setattr(cli, "read_lines", lambda chunks: [])
    monkeypatch.setattr(cli, "parse_records", lambda lines: iter([]))
    monkeypatch.setattr(
        cli,
        "compute_summary",
        lambda aggregator: (_ for _ in ()).throw(RuntimeError("boom")),
    )

    monkeypatch.setattr(
        sys,
        "argv",
        ["phredline", str(input_path), str(output_dir)],
    )

    with pytest.raises(SystemExit) as excinfo:
        cli.main()

    assert excinfo.value.code == 1
    assert "Error boom" in capsys.readouterr().err


def test_main_reports_missing_input_file(tmp_path, monkeypatch, capsys):
    input_path = tmp_path / "missing.fastq"
    output_dir = tmp_path / "report"

    monkeypatch.setattr(
        sys,
        "argv",
        ["phredline", str(input_path), str(output_dir)],
    )

    with pytest.raises(SystemExit) as excinfo:
        cli.main()

    assert excinfo.value.code == 1
    assert "Input file not found" in capsys.readouterr().err
