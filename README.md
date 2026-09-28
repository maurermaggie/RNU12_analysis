# RNU12 analysis

RNU12-related RNA-seq processing and FRASER splicing outlier analysis (June 2025 cohort work).

## Oak path

Primary project directory on Sherlock/Oak:

`/oak/stanford/groups/smontgom/shared/UDN/Analysis/RNU12_analysis`

Previously:

`/oak/stanford/groups/smontgom/shared/UDN/Analysis/Transcriptome_Wide_Splicing_Analysis/Projects/RNU12/Jun_2025`

## Layout

| Directory | Purpose |
|-----------|---------|
| `config/` | YAML configs for alignment and FRASER |
| `cram_to_fastq_workflow/` | CRAM → FASTQ |
| `fastq_to_bam_workflow/`, `nygc_fastq_to_bam_workflow/` | Alignment pipelines |
| `realign_bams/` | BAM realignment |
| `FRASER/` | FRASER1/FRASER2 workflow and configs |
| `output/` | Pipeline and FRASER outputs |
| `symlinked_bams/`, `symlinked_fastqs/` | Sample symlinks |
| `metadata/` | Sample metadata |
| `French_approach/` | [intronaute](https://github.com/JMdeSteAgathe/intronaute) (git submodule) |

## GitHub

This tree is linked to the private repo [maurermaggie/RNU12_analysis](https://github.com/maurermaggie/RNU12_analysis). Large data (BAMs, FASTQs, FRASER count matrices, most logs) are excluded via `.gitignore` and remain on Oak only.

After cloning this repo, initialize submodules:

```bash
git clone --recurse-submodules git@github.com:maurermaggie/RNU12_analysis.git
# or, if already cloned:
git submodule update --init --recursive
```
