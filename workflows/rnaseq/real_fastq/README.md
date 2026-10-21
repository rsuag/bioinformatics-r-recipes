# Seagrass meta-transcriptome FASTQ excerpts

These four files are paired-end excerpts from two owner-provided seagrass meta-transcriptome samples. They contain **the first 250 read pairs** from each cleaned source pair, not a random sample. Source files on the owner's server:

- `SCS_M_0_1.R1.clean.fastq.gz` and `SCS_M_0_1.R2.clean.fastq.gz`
- `SSS_M_0_1.R1.clean.fastq.gz` and `SSS_M_0_1.R2.clean.fastq.gz`

The original instrument headers were replaced with `@<sample>_read_<number>/<mate>`. Nucleotide sequences and quality strings were retained. The repository owner explicitly authorized publication of these excerpts on 2026-10-21. The checksum list is in [SHA256SUMS](SHA256SUMS).

All four complete source `.fastq.gz` files passed `gzip -t` before publication of these excerpts. The four published excerpts also passed archive checks, record-count checks, and mate-pair checks.

Run basic paired-read QC from the repository root:

```bash
Rscript workflows/rnaseq/fastq_qc.R \
  workflows/rnaseq/real_fastq/SCS_M_0_1.R1.first250.fastq.gz \
  workflows/rnaseq/real_fastq/SCS_M_0_1.R2.first250.fastq.gz \
  fastq_qc_output 250
```

These short, first-read excerpts are intended for file-format and QC demonstrations only. They cannot be used as gene counts or support biological comparisons or DESeq2 inference. The sample IDs are not interpreted as treatment labels here. The repository's MIT license covers original code; no separate data reuse license is specified for these excerpts.
