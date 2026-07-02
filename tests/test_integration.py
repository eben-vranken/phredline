import json
import subprocess
import sys


def build_synthetic_fastq(record_count: int = 3000, read_length: int = 80) -> bytes:
    chunks = []
    qualities = b"I" * read_length

    for index in range(record_count):
        sequence = (b"ACGT" * (read_length // 4 + 1))[:read_length]
        if index % 3 == 1:
            sequence = sequence.replace(b"A", b"N", 1)
        elif index % 3 == 2:
            sequence = sequence.replace(b"C", b"G", 1)

        chunks.extend([f"@read{index}".encode(), sequence, b"+", qualities])

    return b"\n".join(chunks) + b"\n"


def test_full_cli_pipeline_writes_valid_multiqc_section_files(tmp_path):
    input_path = tmp_path / "synthetic.fastq"
    output_dir = tmp_path / "multiqc_sections"

    input_path.write_bytes(build_synthetic_fastq())

    subprocess.run(
        [
            sys.executable,
            "-m",
            "phredline.cli",
            str(input_path),
            str(output_dir),
        ],
        check=True,
        capture_output=True,
        text=True,
    )

    assert output_dir.is_dir()

    section_files = sorted(output_dir.glob("*_mqc.json"))
    assert {path.name for path in section_files} == {
        f"{input_path.stem}_fastq_summary_mqc.json",
        f"{input_path.stem}_per_position_quality_mqc.json",
    }

    sections = [json.loads(path.read_text()) for path in section_files]
    assert {section["id"] for section in sections} == {
        "fastq_summary",
        "per_position_quality",
    }

    summary = next(section for section in sections if section["id"] == "fastq_summary")
    quality = next(
        section for section in sections if section["id"] == "per_position_quality"
    )

    assert summary["data"][input_path.stem]["total_reads"] == 3000
    assert "gc_ratio" in summary["data"][input_path.stem]
    assert quality["data"][input_path.stem]["1"] >= 0
