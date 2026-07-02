<h1 align="center">🧬 Phredline</h1>

<p align="center">
    A memory-bounded streaming FASTQ quality control and filtering engine.
</p>

## Output Contract

Phredline writes one output directory per run. The directory contains MultiQC custom-content
section files named `*_mqc.json`, one file per section, so `multiqc .` can discover the
results directly during its directory scan.