# TOMM40_WGS: Genotyping TOMM40'523 Poly-T Polymorphisms from Whole-Genome Sequencing

`TOMM40_WGS` is a Snakemake pipeline that enables automated genotyping of the TOMM40'523 poly-T repeat from WGS data. It integrates STR genotyping tools, k-mer analysis, and machine learning to provide accurate repeat length and genotype predictions.

---

## Background

The TOMM40'523 poly-T repeat polymorphism (`rs10524523`) is a variable-length poly-T repeat located in intron 6 of the *TOMM40* gene on chromosome 19. This polymorphism has been associated with age of onset of Alzheimer's disease (AD) and cognitive decline.

The repeat length is typically categorized into three classes:

* **Short (S)**: ≤19 T residues
* **Long (L)**: 20-29 T residues
* **Very Long (VL)**: ≥30 T residues

---

## 🚀 Quick Start

### Software Dependencies

- **Conda**: install via [Miniforge](https://conda-forge.org/download/) (recommended) 

### Create a conda environment
```bash
conda config --set channel_priority strict
conda create -y -n TOMM40_WGS_env -c conda-forge -c bioconda \
snakemake==9.3.0 samtools==1.21 pandas==2.2.3
conda activate TOMM40_WGS_env
```

### Clone the repository and install dependencies
```bash
git clone https://github.com/RushAlz/TOMM40_WGS.git
cd TOMM40_WGS
chmod u+x TOMM40_WGS

# Get help
./TOMM40_WGS --help
```

### (soon) Alternatively: Install with bioconda 
```bash
conda install -c conda-forge -c bioconda tomm40_wgs
```

Use the `TOMM40_WGS` launcher script to run the pipeline with minimal setup.

---

### Basic usage (recommended)

```bash
# Basic usage with a single BAM file (aligned to GRCh38 by default)
TOMM40_WGS --input_wgs sample.bam --ref_fasta GRCh38.fa --cores 2

# Basic usage with multiple BAM files using glob pattern (quotes are required)
TOMM40_WGS --input_wgs "*.bam" --ref_fasta GRCh38.fa --cores 2

# Basic usage with CRAM files 
TOMM40_WGS --input_wgs "*.cram" --ref_fasta GRCh38.fa --genome_build GRCh38 --cores 2

# Using a sample table (TSV file with sample_id and path columns)
TOMM40_WGS --input_table samples.tsv --ref_fasta GRCh38.fa --cores 2

# With a custom configuration file
TOMM40_WGS --configfile custom_config.yaml --cores 2
```

This script wraps Snakemake, applies default values when needed, and accepts any extra Snakemake flags (e.g., `--dryrun`).

---

### End-to-end testing with real data (data from 1000 Genomes Project)

[Reference for the data](https://doi.org/10.1016/j.cell.2022.08.004)

> Install Miniforge3

```bash
# Skip if conda is already installed
wget -O Miniforge3.sh "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3.sh -b -p "${HOME}/conda"
source "${HOME}/conda/etc/profile.d/conda.sh"
```

> Create TOMM40_WGS_env

```bash
# Skip if already created
conda config --set channel_priority strict
conda create -y -n TOMM40_WGS_env -c conda-forge -c bioconda \
snakemake==9.3.0 samtools==1.21 pandas==2.2.3
conda activate TOMM40_WGS_env
```

> Clone the repository and create a testing directory

```bash
git clone https://github.com/RushAlz/TOMM40_WGS.git
chmod u+x TOMM40_WGS/TOMM40_WGS

# Make a temporary working directory
mkdir -p tomm40_wgs_test
cd tomm40_wgs_test
```

> Test TOMM40_WGS with GRCh38 CRAMs

```bash
# Download the reference genome
REF_FASTA="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"
REF_FASTA_LOCAL=$PWD/"GRCh38_full_analysis_set_plus_decoy_hla.fa"
wget -O ${REF_FASTA_LOCAL} ${REF_FASTA}
wget -O ${REF_FASTA_LOCAL}.fai ${REF_FASTA}.fai

# Subset the CRAM file to the TOMM40 region
CRAMFILE_1="https://ftp.sra.ebi.ac.uk/vol1/run/ERR324/ERR3240114/HG00096.final.cram"
CRAM_1_LOCAL=$PWD/"HG00096.GRCh38.cram"
samtools view -h -T GRCh38_full_analysis_set_plus_decoy_hla.fa -C ${CRAMFILE_1} "chr19:44399792-45399826" > ${CRAM_1_LOCAL} 
samtools index ${CRAM_1_LOCAL}

# Download a second CRAM file for testing
CRAMFILE_2="https://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239334/NA12878.final.cram"
CRAM_2_LOCAL=$PWD/"NA12878.GRCh38.cram"
samtools view -h -T GRCh38_full_analysis_set_plus_decoy_hla.fa -C ${CRAMFILE_2} "chr19:44399792-45399826" > ${CRAM_2_LOCAL}
samtools index ${CRAM_2_LOCAL}

# Run the pipeline specifying the input files as a glob (quotes are required) 
../TOMM40_WGS/TOMM40_WGS --input_wgs "*.cram" \
--ref_fasta GRCh38_full_analysis_set_plus_decoy_hla.fa \
--output_dir results_GRCh38 \
--cores 2

# Alternatively, run the pipeline specifying the input files as a TSV file
echo -e "sample_id\tpath\nHG00096\t${CRAM_1_LOCAL}\nNA12878\t${CRAM_2_LOCAL}" > samples.tsv
../TOMM40_WGS/TOMM40_WGS --input_table samples.tsv \
--ref_fasta GRCh38_full_analysis_set_plus_decoy_hla.fa \
--output_dir results_GRCh38 \
--cores 2

# Check results
cat results_GRCh38/predictions/tomm40_predictions_summary.tsv | column -t
```

> Test TOMM40_WGS with b37 BAMs

```bash
# Download the reference genome
REF_FASTA="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz"
REF_FAI="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz.fai"
wget ${REF_FASTA}
wget ${REF_FAI}
REF_FASTA_LOCAL=$PWD/"hs37d5.fa.gz"

# Subset the BAM file to the TOMM40 region
BAMFILE_1="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00096/high_coverage_alignment/HG00096.wgs.ILLUMINA.bwa.GBR.high_cov_pcr_free.20140203.bam"
BAM_1_LOCAL=$PWD/"HG00096.hg37.bam"
samtools view -h -b ${BAMFILE_1} "19:44903049-45903083" > ${BAM_1_LOCAL}.tmp
samtools addreplacerg -r '@RG\tID:HG00096\tSM:HG00096\tLB:lib1' -o ${BAM_1_LOCAL} ${BAM_1_LOCAL}.tmp
samtools index ${BAM_1_LOCAL}

# Download a second CRAM file for testing
BAMFILE_2="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA12878/high_coverage_alignment/NA12878.mapped.ILLUMINA.bwa.CEU.high_coverage_pcr_free.20130906.bam"
BAM_2_LOCAL=$PWD/"NA12878.hg37.bam"
samtools view -h -b ${BAMFILE_2} "19:44903049-45903083" > ${BAM_2_LOCAL}.tmp
samtools addreplacerg -r '@RG\tID:NA12878\tSM:NA12878\tLB:lib1' -o ${BAM_2_LOCAL} ${BAM_2_LOCAL}.tmp

samtools index ${BAM_2_LOCAL}

# Run the pipeline specifying the input files as a glob (quotes are required) 
../TOMM40_WGS/TOMM40_WGS --input_wgs "*.bam" \
--ref_fasta hs37d5.fa.gz \
--genome_build b37 \
--output_dir results_b37 \
--cores 2 


# Check results
cat results_b37/predictions/tomm40_predictions_summary.tsv | column -t
```

---

## Pipeline Overview

The pipeline performs the following steps for each input sample:

1. **Subset CRAM/BAM** to the TOMM40 region
2. **STR genotyping** using:
   - HipSTR
   - GangSTR
   - ExpansionHunter
3. **K-mer extraction** from poly-T region using Jellyfish
4. **Feature processing** using R scripts
5. **ML-based prediction** of repeat lengths and TOMM40 genotype class
6. **Aggregate predictions** into a summary table

---

## Input Configuration

Input data can be specified in two ways:

1. **Using command line parameters**:
   - `--input_wgs`: A single file, a glob pattern (must be quoted, e.g., `"*.cram"`)
   - `--input_table`: Path to a TSV file with `sample_id` and `path` columns

2. **Using configuration file**:
   - `input_glob`: A glob pattern to CRAM/BAM files
   - `input_table`: Path to a TSV file with `sample_id` and `path` columns

When both `--input_wgs` and `--input_table` are specified, `--input_table` takes precedence.

`input_table` format:

```tsv
sample_id    path
sample1      /path/to/sample1.bam
sample2      /path/to/sample2.bam
```

- `genome_build`: Specify the genome build according to the reference genome used for alignment. The default is `GRCh38`. 

Supported builds include:

  - GRCh38
  - GRCh37
  - hg19
  - b37

```yaml
# TOMM40`523 (rs10524523) coordinates mapping to genome builds
"GRCh38": "chr19:44899792-44899826"
"GRCh37": "chr19:45403049-45403083"
"hg19": "chr19:45403049-45403083"
"b37": "19:45403049-45403083"
```

**Important** currently only GRCh38 has been extensively tested. 

---

##  Outputs

Each sample generates:

- `results/GangSTR/`, `HipSTR/`, `ExpansionHunter/`: Genotyping VCFs for each STR tool
- `results/jellyfish/{sample}.polyT_kmer.txt`: k-mer counts 
- `results/features/{sample}_features.txt`: Combined features for ML
- `results/predictions/{sample}.tomm40_prediction.txt`: Length/genotype predictions for {sample}
- `results/predictions/tomm40_predictions_summary.tsv`: Combined summary

---

## License

TOMM40_WGS is licensed under the GPL-2.0 License - see the [LICENSE](LICENSE) file for details. 
Important licensing information is available [here](docs/license_instructions.md).

---

## Citation
If you use TOMM40_WGS, please cite:

> Vialle RA et al. (2025). *Genotyping TOMM40'523 Poly-T Polymorphisms Using Whole-Genome Sequencing*. [medRxiv 2025.04.23.25326276; doi: https://doi.org/10.1101/2025.04.23.25326276](https://doi.org/10.1101/2025.04.23.25326276)
