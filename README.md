# TTISS On- and Off-Target Insertion Analysis

Pipeline for determining R2 retrotransposon insertion site specificity using tagmentation-based tag integration site sequencing (TTISS). Adapted from [Edmonds *et al.*, 2025](https://doi.org/10.1038/s41467-025-61321-z) for use in plant genomes.

## Overview

TTISS identifies the genomic locations of payload insertions by enriching for read pairs where one end originates from the inserted sequence (detected via the R2Tg 3′ UTR) and the other end maps to the flanking genomic sequence. Reads are classified as:

- **On-target** — both reads of a pair map to a 25S rDNA site
- **Off-target** — both reads map to the genome but not to a 25S site
- **Plasmid-derived** — both reads map to the payload cassette (excluded from on/off-target analysis)

A payload-only negative control (no R2 protein) is run in parallel to estimate the false-positive rate from PCR chimeras and correct final counts.

## Dependencies and environment setup

| Tool | Version tested | Purpose |
|------|---------------|---------|
| [cutadapt](https://cutadapt.readthedocs.io/) | ≥4.0 | Adapter trimming and filtering |
| [bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/) | ≥2.4 | Short-read alignment |
| [samtools](http://www.htslib.org/) | ≥1.16 | BAM processing |
| [bedtools](https://bedtools.readthedocs.io/) | ≥2.30 | Genomic interval operations |
| [barrnap](https://github.com/tseemann/barrnap) | ≥0.9 | rRNA annotation |

### Apple Silicon (M-series) note

barrnap depends on HMMER 3.1b2, which runs on x86-64 only. On Apple Silicon, Rosetta 2 is required and barrnap must be installed in a separate x64 conda environment. Install Rosetta if not already present:

```bash
softwareupdate --install-rosetta
```

We provide two conda environment files: `r2_environment.yml` for all pipeline steps, and `barrnap_environment.yml` for the barrnap annotation step only.

```bash
conda env create -f r2_environment.yml                            # creates ttiss-workflow
CONDA_SUBDIR=osx-64 conda env create -f barrnap_environment.yml  # creates barrnap-workflow
```

Activate `ttiss-workflow` before running the pipeline (all steps except Step 7):

```bash
conda activate ttiss-workflow
```

For Step 7 only, switch to the barrnap environment:

```bash
conda activate barrnap-workflow
conda config --env --set subdir osx-64
```

## Required Input Files

```
data/
  <sample>/
    <sample>_R1.fastq.gz   # paired-end reads, R1 should contain R2Tg 3' UTR
    <sample>_R2.fastq.gz
genome/
  genome_plasmid.fasta     # N. benthamiana genome + payload cassette, concatenated
                           # payload contig must be named "plasmid"
index/
  benth_plasmid.*          # bowtie2 index built from genome_plasmid.fasta
```

Build the bowtie2 index once before running the pipeline:

```bash
bowtie2-build genome/genome_plasmid.fasta index/benth_plasmid
```

Set `SAMPLE` to the sample prefix (e.g., `rep1`, `neg_control`) at the start of each run. All commands below use this variable.

```bash
SAMPLE=rep1
mkdir -p final
```

---

## Pipeline

### Step 1 — Filter for R2Tg 3′ UTR and trim insert sequence

Retains only read pairs where R1 begins with the complete R2Tg 3′ UTR, then removes the UTR sequence from R1 to improve downstream mapping.

```bash
cutadapt -j 0 \
  -g 'GGAACATATATAATTTATGTGTGTTCGATAAATAGC;min_overlap=36;max_error_rate=0.05' \
  --no-indels \
  --discard-untrimmed \
  --pair-filter=first \
  --minimum-length 20 \
  -o final/${SAMPLE}_trim_R1.fastq.gz \
  -p final/${SAMPLE}_trim_R2.fastq.gz \
  data/${SAMPLE}/${SAMPLE}_R1.fastq.gz \
  data/${SAMPLE}/${SAMPLE}_R2.fastq.gz \
  > final/${SAMPLE}_cutadapt_1.log
```

Approximately 91% of R1 reads should contain the full UTR and pass this filter.

---

### Step 2 — Trim residual adaptor sequences

Trims the Tn5 ME motif from the 3′ end of R1 (`-a`) and its reverse complement from the 5′ end of R2 (`-G`), and removes the reverse complement of the R2Tg 3′ UTR from the 3′ end of R2 (`-A`). No reads are discarded in this step.

```bash
cutadapt -j 0 \
  -a 'CTGTCTCTTATACACATCT;min_overlap=12;max_error_rate=0.1' \
  -G 'AGATGTGTATAAGAGACAG;min_overlap=19;max_error_rate=0.1' \
  -A 'GCTATTTATCGAACACACATAAATTATATATGTTCC;min_overlap=15;max_error_rate=0.1' \
  --no-indels \
  --minimum-length 20 \
  -o final/${SAMPLE}_trim2_R1.fastq.gz \
  -p final/${SAMPLE}_trim2_R2.fastq.gz \
  final/${SAMPLE}_trim_R1.fastq.gz \
  final/${SAMPLE}_trim_R2.fastq.gz \
  > final/${SAMPLE}_cutadapt_2.log
```

---

### Step 3 — Map to genome + payload reference

Map trimmed read pairs to the combined genome/payload index. Expect this to take 15–45 minutes per sample.

```bash
bowtie2 -p 10 \
  -x index/benth_plasmid \
  -1 final/${SAMPLE}_trim2_R1.fastq.gz \
  -2 final/${SAMPLE}_trim2_R2.fastq.gz \
  -D 20 -R 3 -N 0 -L 30 -i S,1,0.50 \
  --no-unal \
  2> final/${SAMPLE}_bowtie.log \
  | samtools view -@ 4 -b -o final/${SAMPLE}_mapped.bam
```

---

### Step 4 — Sort and index BAM

```bash
samtools sort final/${SAMPLE}_mapped.bam -o final/${SAMPLE}_sorted.bam
samtools index -@ 10 final/${SAMPLE}_sorted.bam
```

---

### Step 5 — Count all concordantly mapped read pairs

Extracts read names where both members of a pair mapped successfully. This is the denominator for subsequent categorization steps.

```bash
samtools view -F 0x90C -f 0x1 final/${SAMPLE}_sorted.bam \
  | cut -f1 | sort -u > final/${SAMPLE}_all_pairs_mapped.txt

wc -l final/${SAMPLE}_all_pairs_mapped.txt
```

---

### Step 6 — Count reads mapping to the payload cassette

Both reads in a pair must map to the plasmid contig. The same concordance criterion is applied throughout to avoid biasing comparisons between categories (reads mapping to 25S sites may map to multiple loci, so strict concordance would undercount them).

```bash
samtools view -F 0x90C -f 0x1 final/${SAMPLE}_sorted.bam plasmid \
  | cut -f1 \
  | sort \
  | uniq -c \
  | awk '$1==2 {print $2}' > final/${SAMPLE}_plasmid.txt

wc -l final/${SAMPLE}_plasmid.txt
```

Approximately 98% of mapped read pairs are expected to come from the payload cassette and are excluded from on/off-target analysis.

---

### Step 7 — Annotate 25S rDNA sites and count on-target reads

First, generate a BED file of 25S (28S) rRNA features using barrnap (run once per genome). On Apple Silicon, activate the barrnap environment before running this step (see [Apple Silicon note](#apple-silicon-m-series-note) above). For the 2.88 Gb *N. benthamiana* genome this takes approximately 30 minutes with 10 threads on M4 Pro.

```bash
barrnap --kingdom euk --threads 10 genome/genome_plasmid.fasta > genome.gff

awk 'BEGIN{FS=OFS="\t"}
  $0 !~ /^#/ && $3 ~ /rRNA/ && tolower($9) ~ /28s/ {
    print $1, $4-1, $5, "25S_rDNA", ".", $7
  }' genome.gff > 25S_features.bed
```

Then count read pairs where both reads overlap a 25S site:

```bash
samtools view -b -F 0x90C -f 0x1 final/${SAMPLE}_sorted.bam \
  | bedtools intersect -abam - -b 25S_features.bed -u \
  | samtools view - \
  | awk '{
      q=$1; f=$2+0
      if (int(f/64)%2)  a[q]=1
      if (int(f/128)%2) b[q]=1
    }
    END{
      for (q in a) if (b[q]) print q
    }' \
  | sort -u > final/${SAMPLE}_25S.txt

wc -l final/${SAMPLE}_25S.txt
```

---

### Step 8 — Count off-target reads

Off-target insertions are all concordantly mapped read pairs that do not map to the payload or a 25S site.

```bash
# All reads touching the plasmid (either member of pair)
samtools view -F 0x90C -f 0x1 final/${SAMPLE}_sorted.bam plasmid \
  | cut -f1 | sort -u > final/${SAMPLE}_pairs_touching_plasmid.txt

# All reads touching a 25S site (either member of pair)
samtools view -b -F 0x90C -f 0x1 final/${SAMPLE}_sorted.bam \
  | bedtools intersect -abam - -b 25S_features.bed -u \
  | samtools view - \
  | cut -f1 | sort -u > final/${SAMPLE}_pairs_touching_25S.txt

# Union of plasmid- and 25S-touching reads
cat final/${SAMPLE}_pairs_touching_plasmid.txt \
    final/${SAMPLE}_pairs_touching_25S.txt \
  | sort -u > final/${SAMPLE}_pairs_touching_plasmid_or_25S.txt

# Off-target = all mapped pairs minus plasmid/25S-touching pairs
comm -23 \
  final/${SAMPLE}_all_pairs_mapped.txt \
  final/${SAMPLE}_pairs_touching_plasmid_or_25S.txt \
  > final/${SAMPLE}_pairs_other.txt

wc -l final/${SAMPLE}_pairs_other.txt
```

---

## Negative Control Correction

The payload-only negative control (no R2 protein) is processed through the full pipeline. Any on- or off-target reads recovered from this sample represent PCR chimeras, not true integration events. Use the negative control counts to calculate false-positive rates:

```
FP rate (on-target)  = neg_control_25S_reads  / neg_control_all_mapped_reads
FP rate (off-target) = neg_control_other_reads / neg_control_all_mapped_reads
```

Apply these rates to each experimental replicate:

```
corrected_on-target  = observed_on-target  - (FP_rate_on  × all_mapped_reads)
corrected_off-target = observed_off-target - (FP_rate_off × all_mapped_reads)
```

Final on-target rate = `corrected_on-target / (corrected_on-target + corrected_off-target)`.

---

## Output Files

Per-genome intermediates (generated once in Step 7, reusable across samples):

| File | Contents |
|------|----------|
| `genome.gff` | Full barrnap rRNA annotation |
| `25S_features.bed` | 25S rDNA intervals used for on/off-target classification |

Per-sample outputs:

| File | Contents |
|------|----------|
| `final/<sample>_all_pairs_mapped.txt` | Read names of all concordantly mapped pairs |
| `final/<sample>_plasmid.txt` | Read names mapping to payload cassette |
| `final/<sample>_25S.txt` | Read names mapping to 25S rDNA sites (on-target) |
| `final/<sample>_pairs_other.txt` | Read names of off-target insertions |
| `final/<sample>_cutadapt_1.log` | UTR filtering summary |
| `final/<sample>_cutadapt_2.log` | Adaptor trimming summary |
| `final/<sample>_bowtie.log` | Mapping summary |

---

## References

- Edmonds *et al.* (2025) *Nature Communications*. TTISS methodology. [doi:10.1038/s41467-025-61321-z](https://doi.org/10.1038/s41467-025-61321-z)
- Chen *et al.* (2024) *Nature Plants*. *Nicotiana benthamiana* genome assembly. [doi:10.1038/s41477-024-01849-y](https://doi.org/10.1038/s41477-024-01849-y)
