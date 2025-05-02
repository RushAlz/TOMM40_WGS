import os
import glob
import pandas as pd

# Globals
SNAKEDIR = os.path.dirname(workflow.snakefile) + "/"
# if "output_dir" is not set =  $PWD/results
if not config.get("output_dir"):
    config["output_dir"] = os.getcwd() + "/results"


LOGDIR = config["output_dir"] + "/logs"
os.makedirs(LOGDIR, exist_ok=True)

# Load sample information from either input_glob, input_table
def load_samples():
    if config.get("input_table"):
        df = pd.read_csv(config["input_table"], sep="\t")
        return dict(zip(df["sample_id"], df["path"]))
    elif config.get("input_glob"):
        files = sorted(glob.glob(config["input_glob"]))
        return {os.path.basename(f).split(".")[0]: f for f in files}
    else: # Error message
        # Instead of returning a set with a print function call:
        print("No input files found. Please provide either input_glob or input_table in the config file.")
        # Return an empty dictionary instead:
        return {}

sample_dict = load_samples()
sample_ids = list(sample_dict.keys())
genome_build = config["genome_build"]

region_mapping = {
    "GRCh38": "chr19:44399792-45399826",
    "GRCh37": "chr19:44399780-45399837",
    "hg19": "chr19:44399780-45399837",
    "b37": "19:44399780-45399837"
}

def cram_lookup(wildcards):
    return sample_dict[wildcards.sample]

rule all:
    input:
        expand("results/{sample}/{sample}.{genome_build}.TOMM40.cram.crai", sample=sample_ids, genome_build=[genome_build]),
        expand("results/HipSTR/{sample}.TOMM40polyT.maxflankindel09.vcf.gz", sample=sample_ids),
        expand("results/jellyfish/{sample}/{sample}.polyT_kmer.txt", sample=sample_ids),
        expand("results/ExpansionHunter/{sample}/{sample}.expansionHunter.vcf", sample=sample_ids),
        expand("results/GangSTR/{sample}/{sample}.GangSTR.vcf", sample=sample_ids),
        expand("results/features/{sample}_features.txt", sample=sample_ids),
        expand("results/predictions/{sample}.tomm40_prediction.txt", sample=sample_ids),
        "results/predictions/tomm40_predictions_summary.tsv"

wildcard_constraints:
    genome_build=".+"
        
rule subset_cram:
    output:
        cram=temp("results/{sample}/{sample}." + genome_build + ".TOMM40.cram"),
        index=temp("results/{sample}/{sample}." + genome_build + ".TOMM40.cram.crai")
    input:
        cram=cram_lookup,
        ref_fasta=config["ref_fasta"]
    log:
        "results/logs/subset_cram_{sample}.log"
    params:
        region=region_mapping[genome_build]
    conda:
        "envs/samtools.yaml"
    shell:
        """
        mkdir -p results/{wildcards.sample}
        samtools view -h -T {input.ref_fasta} -C {input.cram} {params.region} > {output.cram} 2> {log}
        samtools index {output.cram} 2>> {log}
        """

rule hipstr:
    output:
        vcf="results/HipSTR/{sample}.TOMM40polyT.maxflankindel09.vcf.gz"
    input:
        cram="results/{sample}/{sample}." + genome_build + ".TOMM40.cram",
        fasta=config["ref_fasta"],
        bed=SNAKEDIR + "resources/HipSTR.tomm40.GRCh38.bed",
        hipstr_exec=SNAKEDIR + "resources/tools/HipSTR"
    log:
        "results/logs/hipstr_{sample}.log"
    params:
        min_reads=15
    conda:
        "envs/samtools.yaml"
    shell:
        """
        mkdir -p results/HipSTR
        {input.hipstr_exec} \
            --bams {input.cram} \
            --fasta {input.fasta} \
            --regions {input.bed} \
            --str-vcf {output.vcf} \
            --min-reads {params.min_reads} \
            --output-filters \
            --max-flank-indel 0.9  &> {log}
        """
        
rule strling:
    output:
        strling_cram="results/jellyfish/{sample}/jf/{sample}.strling.cram",
        string_fa="results/jellyfish/{sample}/jf/{sample}.strling.fa"
    input:
        cram="results/{sample}/{sample}." + genome_build + ".TOMM40.cram",
        ref_fasta=config["ref_fasta"]
    log:
        "results/logs/strling_{sample}.log"
    params:
        region_mapping=region_mapping[genome_build]
    conda:
        "envs/strling.yaml"
    shell:
        """
        strling pull_region -f {input.ref_fasta} -o {output.strling_cram} {input.cram} {params.region_mapping} &> {log}
        samtools fasta -f 4 {output.strling_cram} > {output.string_fa} 2>> {log}
        samtools fasta -F 4 {output.strling_cram} >> {output.string_fa} 2>> {log}
        """

rule jellyfish_kmer:
    output:
        kmer_count="results/jellyfish/{sample}/{sample}.polyT_kmer.txt"
    input:
        string_fa="results/jellyfish/{sample}/jf/{sample}.strling.fa"
    log:
        "results/logs/jellyfish_kmer_{sample}.log"
    params:
        region_mapping=region_mapping[genome_build]
    threads: 8
    conda:
        "envs/jellyfish.yaml"
    shell:
        """
        cat /dev/null > {output.kmer_count}
        for i in $(seq 3 50); do
            jellyfish count -m $i -s 10M -t {threads} -C -o results/jellyfish/{wildcards.sample}/jf/${{i}}mer.jf {input.string_fa} &>> {log}
            repeat_T=$(printf '%*s' "$i" | tr ' ' 'T')
            jellyfish query results/jellyfish/{wildcards.sample}/jf/${{i}}mer.jf "$repeat_T" | awk -v i=$i '{{print i"\t"$2}}' >> {output.kmer_count} 2>> {log}
        done
        """

rule expansion_hunter:
    output:
        "results/ExpansionHunter/{sample}/{sample}.expansionHunter.vcf"
    input:
        cram="results/{sample}/{sample}." + genome_build + ".TOMM40.cram",
        ref_fasta=config["ref_fasta"],
        expansion_hunter_locus=SNAKEDIR + "resources/expansionHunter.tomm40.GRCh38.json"
    log:
        "results/logs/expansion_hunter_{sample}.log"
    conda:
        "envs/expansionhunter.yaml"
    shell:
        """
        mkdir -p results/ExpansionHunter/{wildcards.sample}
        ExpansionHunter --reads {input.cram} \
            --reference {input.ref_fasta} \
            --variant-catalog {input.expansion_hunter_locus} \
            --output-prefix results/ExpansionHunter/{wildcards.sample}/{wildcards.sample}.expansionHunter \
            --analysis-mode streaming &> {log}
        """

rule gangstr:
    output:
        "results/GangSTR/{sample}/{sample}.GangSTR.vcf"
    input:
        cram="results/{sample}/{sample}." + genome_build + ".TOMM40.cram",
        ref_fasta=config["ref_fasta"],
        regions=SNAKEDIR + "resources/GangSTR.tomm40.GRCh38.bed"
    log:
        "results/logs/gangstr_{sample}.log"
    conda:
        "envs/gangstr.yaml"
    shell:
        """
        mkdir -p results/GangSTR/{wildcards.sample}
        GangSTR --bam {input.cram} \
            --ref {input.ref_fasta} \
            --targeted \
            --regions {input.regions} \
            --output-readinfo \
            --include-ggl \
            --out results/GangSTR/{wildcards.sample}/{wildcards.sample}.GangSTR &> {log}
        """

rule generate_features:
    input:
        expansionhunter="results/ExpansionHunter/{sample}/{sample}.expansionHunter.vcf",
        gangstr="results/GangSTR/{sample}/{sample}.GangSTR.vcf",
        jellyfish="results/jellyfish/{sample}/{sample}.polyT_kmer.txt"
    output:
        "results/features/{sample}_features.txt"
    log:
        "results/logs/generate_features_{sample}.log"
    params:
        script=SNAKEDIR + "scripts/process_features.R",
        script_lib=SNAKEDIR + "scripts/feature_parsing_functions.R",
        MLP_model1=SNAKEDIR + "resources/hap_classifier.RData",
        MLP_model2=SNAKEDIR + "resources/hap_classifier_3.RData"
    conda:
        "envs/R.yaml"
    shell:
        """
        mkdir -p results/features
        Rscript {params.script} \
            --script_lib {params.script_lib} \
            --run_id {wildcards.sample} \
            --expansionhunter {input.expansionhunter} \
            --gangstr {input.gangstr} \
            --jellyfish {input.jellyfish} \
            --MLP_model1 {params.MLP_model1} \
            --MLP_model2 {params.MLP_model2} \
            --output {output} &> {log}
        """

rule predict_tomm40:
    input:
        features="results/features/{sample}_features.txt"
    output:
        prediction="results/predictions/{sample}.tomm40_prediction.txt"
    params:
        script=SNAKEDIR + "scripts/predict_tomm40.R",
        model=SNAKEDIR + "resources/model_fit.RData"
    log:
        "results/logs/predict_tomm40_{sample}.log"
    conda:
        "envs/R.yaml"
    shell:
        """
        mkdir -p results/predictions
        Rscript {params.script} --features {input.features} \
            --model {params.model} \
            --sample_id {wildcards.sample} \
            --output {output.prediction} &> {log}
        """

rule summarize_predictions:
    input:
        predictions=expand("results/predictions/{sample}.tomm40_prediction.txt", sample=sample_ids)
    output:
        summary="results/predictions/tomm40_predictions_summary.tsv"
    log:
        f"{LOGDIR}/summarize_predictions.log"
    run:
        import pandas as pd
        import sys

        sys.stdout = open(log[0], "w")
        sys.stderr = sys.stdout

        print("[INFO] Generating formatted TOMM40 prediction summary...")

        combined_df = pd.DataFrame()
        for pred_file in input.predictions:
            df = pd.read_csv(pred_file, sep=None, engine='python')
            combined_df = pd.concat([combined_df, df], ignore_index=True)

        columns_to_keep = [
            "sample_id", "A1", "A2", "hap_T", "tomm40_cat_A1", "tomm40_cat_A2", "tomm40_genotype"
        ]
        combined_df = combined_df[columns_to_keep]
        combined_df.to_csv(output.summary, sep="\t", index=False)

        print("[INFO] Summary written to", output.summary)
        sys.stdout.close()
