# TTISS On- and Off-Target Insertion Analysis

Pipeline for determining R2 retrotransposon insertion site specificity using tagmentation-based tag integration site sequencing (TTISS). Adapted from [Edmonds *et al.*, 2025](https://doi.org/10.1038/s41467-025-61321-z) for use in plant genomes. The pipeline exploits tagmentation-based enrichment for read pairs in which one end derives from the R2Tg transgene and the other maps to flanking genomic DNA, allowing insertion sites to be identified across the full *N. benthamiana* genome. Each pair is classified as on-target (25S rDNA), off-target, or plasmid-derived, and a payload-only negative control processed in parallel provides the false-positive rate for chimera correction.

## Set up

This workflow has been developed and tested on MacOS Apple Silicon. Of note is the use of `barrnap` to annotate 25S rDNA sites, which requires HMMER 3.1b2 — a tool that runs on x64 architecture only. On Apple Silicon, Rosetta 2 is required and barrnap must be installed in a separate x64 conda environment. Install Rosetta if not already present:

```bash
softwareupdate --install-rosetta
```

We provide two environment files, *ttiss_workflow.yaml* for all pipeline steps and *barrnap_workflow.yaml* for the barrnap annotation step only.

### ttiss-workflow Environment Install

```bash
conda env create -f ttiss_workflow.yaml
```

### barrnap Environment Install

```bash
CONDA_SUBDIR=osx-64 conda env create -f barrnap_workflow.yaml
```

Please use the `ttiss-workflow` environment for all steps except `Step 2: barrnap 25S annotation`.

## Step 1: File organization

The workflow requires the *N. benthamiana* genome assembly concatenated with the payload cassette sequence as a single FASTA file, and paired-end FASTQ files from TTISS sequencing. The genome assembly is available from [CNCB](https://ngdc.cncb.ac.cn/) under project number PRJCA022857, our FASTQ files are deposited at [SRA](https://www.ncbi.nlm.nih.gov/sra) under accession **PRJNA1463879**, and for reference we provide the payload cassette sequence here as plasmid.fasta. The payload contig must be named `plasmid` in the FASTA file. To be compatible with the workflow your directory should look like this:

```
- TTISS-analysis-main/ (name of this folder can be changed)
│   ├── ttiss_workflow.yaml
│   ├── barrnap_workflow.yaml
│   ├── data/
│   │   └── <sample>_R1.fastq.gz
│   │   └── <sample>_R2.fastq.gz
│   ├── genome/
│   │   └── genome_plasmid.fasta
│   ├── index/
│   ├── final/
```

Build the bowtie2 index once before running the pipeline:

```bash
bowtie2-build genome/genome_plasmid.fasta index/benth_plasmid
```

Set `SAMPLE` to the sample prefix (e.g. `rep1`, `neg_control`) at the start of each run. All commands below use this variable.

```bash
SAMPLE=rep1
```

## Step 2: barrnap 25S annotation

Activate the barrnap environment and run barrnap in eukaryotic mode to generate a GFF file of predicted rDNA positions. Run this step once per genome before processing any samples. For the 2.88 Gb *N. benthamiana* genome this takes approximately 30 minutes with 10 threads on M4 Pro. Adjust thread count as needed.

```bash
conda activate barrnap-workflow
conda config --env --set subdir osx-64
barrnap --kingdom euk --threads 10 genome/genome_plasmid.fasta > genome/genome.gff
```

Extract 25S rDNA intervals into a BED file for use in subsequent steps:

```bash
awk 'BEGIN{FS=OFS="\t"}
  $0 !~ /^#/ && $3 ~ /rRNA/ && tolower($9) ~ /28s/ {
    print $1, $4-1, $5, "25S_rDNA", ".", $7
  }' genome/genome.gff > genome/25S_features.bed
```

Switch back to the main environment before continuing:

```bash
conda activate ttiss-workflow
```

## Step 3: Run the analysis pipeline

Switch to the main conda environment for all remaining steps:

```bash
conda activate ttiss-workflow
```

Open `workflow.ipynb` and follow the in-built instructions. Steps 3–5 (read filtering, alignment, and read pair classification) and negative control correction are all completed within the notebook. All output files matching those described in our publication will be generated and saved to the `final/` directory.

## References

- Edmonds *et al.* (2025) *Nature Communications*. TTISS methodology. [doi:10.1038/s41467-025-61321-z](https://doi.org/10.1038/s41467-025-61321-z)
- Chen *et al.* (2024) *Nature Plants*. *Nicotiana benthamiana* genome assembly. [doi:10.1038/s41477-024-01849-y](https://doi.org/10.1038/s41477-024-01849-y)
