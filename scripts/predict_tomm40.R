#!/usr/bin/env Rscript
# predict_tomm40.R
# This script loads a feature file and uses a pre‐trained model to predict TOMM40 haplotypes.
# It writes the prediction output as a text file.

suppressPackageStartupMessages({
  library(optparse)
  library(mlr)
  library(dplyr)
  library(tidyr)
})

# Define command line options
option_list <- list(
  make_option(c("--features"), type = "character", help = "Path to input features file (tab-delimited)"),
  make_option(c("--model"), type = "character", default = "hap_classifier.RData",
              help = "Path to pre-trained model file [default %default]"),
  make_option(c("--sample_id"), type = "character", help = "Sample ID"),
  make_option(c("--output"), type = "character", help = "Path to output prediction text file")
)

opt <- if(!interactive()){
    parse_args(OptionParser(option_list = option_list))
  }else{
    setwd("/pastel/Github_scripts/TOMM40/TOMM40_WGS_testing")
    parse_args(OptionParser(option_list = option_list), args = c(
    "--features", "results/features/R1708627_00246264_features.txt",
    "--output", "results/predictions/R1708627_00246264.tomm40_prediction.txt",
    "--sample_id", "R1708627_00246264",
    "--model", "resources/model_fit.RData"
  ))
}

# Check for required arguments
if (is.null(opt$features) || is.null(opt$output)) {
  stop("Both --features and --output arguments must be provided.")
}

# Load the pre-trained model.
# The model file is expected to load an object 'model_fit' (and optionally 'xgboost_feat_names')
if (!file.exists(opt$model)) {
  stop("Model file not found.")
}
load(opt$model)
if (!exists("model_fit")) {
  stop("Pre-trained model object 'model_fit' not found in the model file.")
}

# Read the features file.
# Expecting a tab-delimited file with header (e.g. as produced by process_features.R)
features_df <- read.table(opt$features, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Ensure column names are valid R variable names.
colnames(features_df) <- make.names(colnames(features_df))

# If the training process saved expected predictor names, subset accordingly.
if (exists("xgboost_feat_names")) {
  predictors <- base::intersect(xgboost_feat_names, colnames(features_df))
  features_df <- features_df[, predictors, drop = FALSE]
  features_df$z = 0
  features_df = features_df %>% mutate_if(is.factor, as.numeric) 
}

predict_outcome <- function(model_fit, ds, param){
  colnames(ds)=make.names(colnames(ds))
  if (is.factor(ds[[param$outcome]])){
    model_fit$learner$predict.type="prob"
    traintask=makeClassifTask(data=ds, target=param$outcome)
  }else{
    traintask=makeRegrTask(data=ds, target=param$outcome)
  }
  model_prediction = stats::predict(model_fit, traintask)
  return(model_prediction)
}

# Perform prediction using the pre-trained mlr model.
pred <- predict_outcome(model_fit = model_fit$model_fit, ds = features_df, param = param)
predict_outcome(model_fit$model_fit,ds %>% mutate(z=0),param)

# Convert predictions to a data frame.
pred_df <- as.data.frame(pred$data)

pred_df_w = pred_df %>% 
  select(-truth) %>% mutate(sample_id = opt$sample_id) %>%
  group_by(sample_id) %>%
  reframe(A1 = round(min(response)),
          A2 = round(max(response))) %>%
  mutate(hap_T = paste0(A1,"/",A2)) %>%
  mutate(
    tomm40_cat_A1 = case_when(
      A1 <= 19 ~ "S",
      A1 >= 20 & A1 <= 29 ~ "L",
      A1 >= 30 ~ "VL"
    ),
    tomm40_cat_A2 = case_when(
      A2 <= 19 ~ "S",
      A2 >= 20 & A2 <= 29 ~ "L",
      A2 >= 30 ~ "VL"
    )
  ) %>%
  mutate(tomm40_genotype = paste0(tomm40_cat_A1,"/",tomm40_cat_A2)) 

# Write predictions to the output file.
write.table(pred_df_w, file = opt$output, sep = "\t", row.names = FALSE, quote = FALSE)

message("Prediction complete. Results written to: ", opt$output)
