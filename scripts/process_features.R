#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(vcfR)
  library(dplyr)
  library(purrr)
  library(data.table)
  library(tidyr)
  library(reshape2)
  library(readr)
  library(stringr)
  library(tibble)
  library(tidymodels)
})

option_list = list(
  make_option(c("--script_lib"), type="character", help="Path to the script library"),
  make_option(c("--run_id"), type="character", help="Run ID"),
  make_option(c("--expansionhunter"), type="character", help="ExpansionHunter VCF file"),
  make_option(c("--gangstr"), type="character", help="GangSTR VCF file"),
  make_option(c("--jellyfish"), type="character", help="Jellyfish k-mer file"),
  make_option(c("--MLP_model1"), type="character", help="MLP model Rfile (3 classes)"),
  make_option(c("--MLP_model2"), type="character", help="MLP model Rfile (6 classes)"),
  make_option(c("--output"), type="character", help="Output file path")
)

opt = if(!interactive()){
  parse_args(OptionParser(option_list=option_list))
}else{
  parse_args(OptionParser(option_list=option_list), args = c(
    "--script_lib", "/home/ricardo_a_vialle/TOMM40_WGS/scripts/feature_parsing_functions.R",
    "--run_id", "HG00171",
    "--expansionhunter", "results/ExpansionHunter/HG00171/HG00171.expansionHunter.vcf",
    "--gangstr", "results/GangSTR/HG00171/HG00171.GangSTR.vcf",
    "--jellyfish", "results/jellyfish/HG00171/HG00171.polyT_kmer.txt",
    "--MLP_model1", "/home/ricardo_a_vialle/TOMM40_WGS/resources/hap_classifier.RData",
    "--MLP_model2", "/home/ricardo_a_vialle/TOMM40_WGS/resources/hap_classifier_3.RData",
    "--output", "results/features/HG00171_features.txt"
  ))
}
source(opt$script_lib)

# Process ExpansionHunter results
expansion_res = parse_expansionHunter_wrap(vcfs = opt$expansionhunter, run_id = opt$run_id)

# Process GangSTR results
gangstr_res = parse_GangSTR_wrap(vcfs = opt$gangstr, run_id = opt$run_id)

# Collect STR tool features
str_features = collect_tool_feats(
  dat_tool_list = list(expansionHunter=expansion_res, GangSTR=gangstr_res),
  tools = c("expansionHunter", "GangSTR")
)

# Split alleles
str_features_A1 = str_features %>% select(contains("A1")) 
str_features_A2 = str_features %>% select(contains("A2")) 
colnames(str_features_A1) = gsub("A1_","",colnames(str_features_A1))
colnames(str_features_A2) = gsub("A2_","",colnames(str_features_A2))
str_features_A1 = bind_cols(str_features_A1, str_features %>% select(!contains(c("A1","A2"))))
str_features_A2 = bind_cols(str_features_A2, str_features %>% select(!contains(c("A1","A2"))))

str_features_A1 = str_features_A1 %>% 
  na.omit() %>%
  rownames_to_column("run_id") %>% mutate(allele = "A1") %>% as.data.frame()
str_features_A2 = str_features_A2 %>% 
  na.omit() %>%
  rownames_to_column("run_id") %>% mutate(allele = "A2") %>% as.data.frame()

str_features_combined = bind_rows(
  cbind(str_features_A1,
        otherAexpansionHunter = str_features_A2$expansionHunter,
        otherA_GangSTR = str_features_A2$GangSTR),
  cbind(str_features_A2,
        otherAexpansionHunter = str_features_A1$expansionHunter,
        otherA_GangSTR = str_features_A1$GangSTR)
)

# Collect Jellyfish k-mer features
kmer_features = fread(opt$jellyfish, col.names = c("kmer", "kmer_count")) %>%
  pivot_wider(names_from = kmer, values_from = kmer_count) %>%
  mutate(sample_id = opt$run_id) %>%
  column_to_rownames("sample_id") %>%
  as.data.frame()

colnames(kmer_features) = paste0("kmer_", colnames(kmer_features))

# Load the pre-trained MLP model (predicts 6 class haplotypes using k-mer counts)
if (!file.exists(opt$MLP_model1)) {
  stop("Model file not found.")
}
load(opt$MLP_model1)
dat_hap_predicted = bind_cols(
  .pred_class = predict(reg_fit$fit, newdata = kmer_features %>% select(contains("kmer")), type = "class"),
  # predict(reg_fit, kmer_features %>% select(contains("kmer"))),
  predict(reg_fit, kmer_features %>% select(contains("kmer")), type = "prob")) 

# Load the three pre-trained MLP models (predicts haplotypes S,L,VL using k-mer counts)
if (!file.exists(opt$MLP_model2)) {
  stop("Model file not found.")
}
load(opt$MLP_model2)
dat_hap_predicted2 = bind_cols(
  pred_S = predict(S_reg_fit, kmer_features %>% select(contains("kmer")), type = "prob")$.pred_1,
  pred_L = predict(L_reg_fit, kmer_features %>% select(contains("kmer")), type = "prob")$.pred_1,
  pred_VL = predict(VL_get_fit, kmer_features %>% select(contains("kmer")), type = "prob")$.pred_1) 

# Combine all features into a single dataframe
combined_features = str_features_combined %>% 
  left_join(bind_cols(kmer_features, dat_hap_predicted, dat_hap_predicted2) %>% rownames_to_column("run_id"), by="run_id")
  
# Write combined features to output file
write.table(combined_features, file=opt$output, sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)
