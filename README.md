# TOMM40_WGS: Genotyping TOMM40 '523 Poly-T Polymorphisms from Whole-Genome Sequencing

TOMM40_WGS is a Snakemake pipeline that enables automated genotyping of the TOMM40 '523 poly-T repeat from WGS data. It integrates STR genotyping tools, k-mer analysis, and machine learning to provide accurate repeat length and genotype predictions.

---

## Background

The TOMM40 '523 poly-T repeat polymorphism (rs10524523) is a variable-length poly-T repeat located in intron 6 of the TOMM40 gene on chromosome 19. This polymorphism has been associated with age of onset of Alzheimer's disease (AD) and cognitive decline.

The repeat length is typically categorized into three classes:
- **Short (S)**: ≤19 T residues
- **Long (L)**: 20-29 T residues
- **Very Long (VL)**: ≥30 T residues

## 🚀 Quick Start

### Clone the repository and install dependencies
```bash
git clone https://github.com/RushAlz/TOMM40_WGS.git
cd TOMM40_WGS
chmod u+x TOMM40_WGS

# Pre-install dependencies (make sure to have snakemake and conda installed first)
./TOMM40_WGS --conda-create-envs-only
```

### (soon) Alternatively: Install with bioconda 
```bash
conda install -c bioconda tomm40_wgs
```

Use the `TOMM40_WGS` launcher script to run the pipeline with minimal setup.

### Software Dependencies

- **R**: 4.0.0 or higher
- **Conda** 

### Basic usage (recommended)

```bash
# Basic usage with a single BAM file (aligned to GRCh38 by default)
TOMM40_WGS --input_wgs sample.bam --ref_fasta GRCh38.fa --cores 8

# Basic usage with multiple BAM files using glob pattern (quotes are required)
TOMM40_WGS --input_wgs "*.bam" --ref_fasta GRCh38.fa --cores 8

# Basic usage with CRAM files 
TOMM40_WGS --input_wgs "*.cram" --ref_fasta GRCh38.fa --genome_build GRCh38 --cores 8

# Using a sample table (TSV file with sample_id and path columns)
TOMM40_WGS --input_table samples.tsv --ref_fasta GRCh38.fa --cores 8

# With a custom configuration file
TOMM40_WGS --configfile custom_config.yaml --cores 8

# To pre-create environments
TOMM40_WGS --conda-create-envs-only

# Get help
TOMM40_WGS --help
```

This script wraps Snakemake, applies default values when needed, and accepts any extra Snakemake flags (e.g., `--dryrun`, `--printshellcmds`, `--timestamp`).

### Testing with real data (data from 1000 Genomes Project)
```bash
# Make a temporary working directory
mkdir -p tomm40_wgs_test
cd tomm40_wgs_test

# Download the reference genome
REF_FASTA="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"
REF_FASTA_LOCAL=$PWD/"GRCh38_full_analysis_set_plus_decoy_hla.fa"
wget -O ${REF_FASTA_LOCAL} ${REF_FASTA}
wget -O ${REF_FASTA_LOCAL}.fai ${REF_FASTA}.fai

# Subset the CRAM file to the TOMM40 region
CRAMFILE_1="https://ftp.sra.ebi.ac.uk/vol1/run/ERR324/ERR3240157/HG00157.final.cram"
CRAM_1_LOCAL=$PWD/"HG00157.final.cram"
samtools view -h -T GRCh38_full_analysis_set_plus_decoy_hla.fa -C ${CRAMFILE_1} "chr19:44399792-45399826" > ${CRAM_1_LOCAL} 
samtools index ${CRAM_1_LOCAL}

# Download a second CRAM file for testing
CRAMFILE_2="https://ftp.sra.ebi.ac.uk/vol1/run/ERR324/ERR3240160/HG00171.final.cram"
CRAM_2_LOCAL=$PWD/"HG00171.final.cram"
samtools view -h -T GRCh38_full_analysis_set_plus_decoy_hla.fa -C ${CRAMFILE_2} "chr19:44399792-45399826" > ${CRAM_2_LOCAL}
samtools index ${CRAM_2_LOCAL}

# Run the pipeline specifying the input files as a glob (quotes are required) 
TOMM40_WGS --input_wgs "*.cram" \
--ref_fasta GRCh38_full_analysis_set_plus_decoy_hla.fa \
--cores 8

# Check results
cat results/predictions/tomm40_predictions_summary.tsv

# Test inputs as TSV file
echo -e "sample_id\tpath\nHG00157\t${CRAM_1_LOCAL}\nHG00171\t${CRAM_2_LOCAL}" > samples.tsv
# Run the pipeline specifying the input files as a TSV file
TOMM40_WGS --input_table samples.tsv \
--ref_fasta GRCh38_full_analysis_set_plus_decoy_hla.fa \
--cores 8

# Check results
cat results/predictions/tomm40_predictions_summary.tsv
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
   - `--input_wgs`: A single file, a glob pattern (must be quoted, e.g., `"*.cram"`), or a space-separated list of files
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

**Important** currently only GRCh38 is supported. Future configuration options:

- `genome_build`: Specify the genome build according to the reference genome used for alignment. The default is `GRCh38`. 
Supported builds include:

  - GRCh38
  - GRCh37
  - hg19
  - b37

```yaml
# TOMM40`523 region mapping to genome builds
"GRCh38": "chr19:44399792-45399826"
"GRCh37": "chr19:44399780-45399837"
"hg19": "chr19:44399780-45399837"
"b37": "19:44399780-45399837"
```

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
If you use this pipeline, please cite:

> Vialle RA et al. (2025). *Genotyping TOMM40'523 Poly-T Polymorphisms Using Whole-Genome Sequencing*. [medRxiv 2025.04.23.25326276; doi: https://doi.org/10.1101/2025.04.23.25326276](https://doi.org/10.1101/2025.04.23.25326276)
