parse_expansionHunter_tomm40 <- function(vcf_df){
  vcf_df = vcf_df %>% 
    separate(Alleles, into = c("REF","ALT"), sep = "\\/") %>% 
    mutate(REF_T = as.numeric(gsub("<STR(.*)>","\\1",REF)), ALT_T = as.numeric(gsub("<STR(.*)>","\\1",ALT))) %>%
    mutate(nT_Alleles = paste0(REF_T,"|",ALT_T))
  
  vcf_df$REF_polyT_cat = case_when(
    vcf_df$REF_T <= 19 ~ "S",
    vcf_df$REF_T >= 20 & vcf_df$REF_T <= 29 ~ "L",
    vcf_df$REF_T >= 30 ~ "VL"
  )
  vcf_df$ALT_polyT_cat = case_when(
    vcf_df$ALT_T <= 19 ~ "S",
    vcf_df$ALT_T >= 20 & vcf_df$ALT_T <= 29 ~ "L",
    vcf_df$ALT_T >= 30 ~ "VL"
  )
  
  vcf_df$TOMM40_Phased = paste0(vcf_df$REF_polyT_cat,"|",vcf_df$ALT_polyT_cat)
  vcf_df$tomm40_hap = NA
  vcf_df$tomm40_hap[vcf_df$TOMM40_Phased == "S|S"] = 1
  vcf_df$tomm40_hap[vcf_df$TOMM40_Phased %in% c("S|L","L|S")] = 2
  vcf_df$tomm40_hap[vcf_df$TOMM40_Phased %in% c("S|VL","VL|S")] = 3
  vcf_df$tomm40_hap[vcf_df$TOMM40_Phased == "L|L"] = 4
  vcf_df$tomm40_hap[vcf_df$TOMM40_Phased %in% c("L|VL","VL|L")] = 5
  vcf_df$tomm40_hap[vcf_df$TOMM40_Phased == "VL|VL" ] = 6
  return(vcf_df)
}
parse_expansionHunter_tomm40_features <- function(vcf){
  vcf.info = extract_info_tidy(vcf)
  ##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">
  ##INFO=<ID=REF,Number=1,Type=Integer,Description="Reference copy number">
  ##INFO=<ID=REPID,Number=1,Type=String,Description="Repeat identifier as specified in the variant catalog">
  ##INFO=<ID=RL,Number=1,Type=Integer,Description="Reference length in bp">
  ##INFO=<ID=RU,Number=1,Type=String,Description="Repeat unit in the reference orientation">
  ##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
  ##INFO=<ID=VARID,Number=1,Type=String,Description="Variant identifier as specified in the variant catalog">
  
  ##FORMAT=<ID=ADFL,Number=1,Type=String,Description="Number of flanking reads consistent with the allele">
  ##FORMAT=<ID=ADIR,Number=1,Type=String,Description="Number of in-repeat reads consistent with the allele">
  ##FORMAT=<ID=ADSP,Number=1,Type=String,Description="Number of spanning reads consistent with the allele">
  ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
  ##FORMAT=<ID=LC,Number=1,Type=Float,Description="Locus coverage">
  ##FORMAT=<ID=REPCI,Number=1,Type=String,Description="Confidence interval for REPCN">
  ##FORMAT=<ID=REPCN,Number=1,Type=String,Description="Number of repeat units spanned by the allele">
  ##FORMAT=<ID=SO,Number=1,Type=String,Description="Type of reads that support the allele; can be SPANNING, FLANKING, or INREPEAT meaning that the reads span, flank, or are fully contained in the repeat">
  vcf.format = extract_gt_tidy(vcf, verbose = F)
  
  gt_ADFL = str_split(vcf.format$gt_ADFL,"/") %>% unlist()
  gt_ADIR = str_split(vcf.format$gt_ADIR,"/") %>% unlist()
  gt_ADSP = str_split(vcf.format$gt_ADSP,"/") %>% unlist()
  gt_REPCI = str_split(vcf.format$gt_REPCI,"/") %>% unlist()
  gt_REPCN = str_split(vcf.format$gt_REPCN,"/") %>% unlist()
  gt_SO = str_split(vcf.format$gt_SO,"/") %>% unlist()
  
  whichA1 = which.min(str_split(vcf.format$gt_REPCN,"/") %>% unlist())
  A1_gt_ADFL = as.numeric(gt_ADFL[whichA1])
  A2_gt_ADFL = as.numeric(gt_ADFL[-whichA1])
  A1_gt_ADIR = as.numeric(gt_ADIR[whichA1])
  A2_gt_ADIR = as.numeric(gt_ADIR[-whichA1])
  A1_gt_ADSP = as.numeric(gt_ADSP[whichA1])
  A2_gt_ADSP = as.numeric(gt_ADSP[-whichA1])
  A1_gt_REPCI = as.numeric(str_split(gt_REPCI[whichA1],"-") %>% unlist())
  A1_gt_REPCI_1 = A1_gt_REPCI[1]
  A1_gt_REPCI_2 = A1_gt_REPCI[2]
  A2_gt_REPCI = as.numeric(str_split(gt_REPCI[-whichA1],"-") %>% unlist())
  A2_gt_REPCI_1 = A2_gt_REPCI[1]
  A2_gt_REPCI_2 = A2_gt_REPCI[2]
  A1_gt_REPCN = as.numeric(gt_REPCN[whichA1])
  A2_gt_REPCN = as.numeric(gt_REPCN[-whichA1])
  A1_gt_SO = as.character(gt_SO[whichA1])
  A2_gt_SO = as.character(gt_SO[-whichA1])
  
  return(data.frame(A1_gt_ADFL,A2_gt_ADFL,A1_gt_ADIR,A2_gt_ADIR,A1_gt_ADSP,A2_gt_ADSP,
                    A1_gt_REPCI_1,A1_gt_REPCI_2,A2_gt_REPCI_1,A2_gt_REPCI_2,
                    A1_gt_REPCN,A2_gt_REPCN,A1_gt_SO,A2_gt_SO))
}
parse_expansionHunter_wrap <- function(vcfs, run_id){
  res_list = list()
  i = 1
  vcf_i = vcfs[i]
  vcf = read.vcfR(vcf_i, verbose = F)
  
  vcf_df = na.omit(reshape2::melt(vcfR::extract.gt(vcf, return.alleles = T)))
  vcf_format = extract_gt_tidy(vcf, verbose = F)
  
  colnames(vcf_df) = c("ID", "Sample", "Alleles")
  vcf_df$run_id = run_id
  vcf_df$Alleles = vcf_format$gt_REPCN
  
  vcf_df = suppressWarnings(parse_expansionHunter_tomm40(vcf_df)) %>%
    mutate(run_id = run_id)
  
  vcf_feats = bind_cols(vcf_df[, c("run_id"),drop=F],
                        parse_expansionHunter_tomm40_features(vcf))

  res = list()
  res$vcf = vcf_df
  res$vcf_feats = vcf_feats
  res_list[[run_id]] = res
  
  res_features = purrr::map_df(res_list, ~bind_rows(.x$vcf_feats))
  res_expansionHunter = purrr::map_df(res_list, ~bind_rows(.x$vcf))
  
  return(list(res_list_expansionHunter = res_list,
              res_features_expansionHunter = res_features,
              res_expansionHunter = res_expansionHunter))
}

parse_GangSTR_tomm40 <- function(vcf_df){
  vcf_df = vcf_df %>% 
    separate(Alleles, into = c("REF","ALT"), sep = "/") %>% 
    mutate(REF_T = str_count(toupper(REF), pattern = "T"), ALT_T = str_count(toupper(ALT), pattern = "T")) %>%
    mutate(nT_Alleles = paste0(REF_T,"|",ALT_T))
  
  vcf_df$REF_polyT_cat = case_when(
    vcf_df$REF_T <= 19 ~ "S",
    vcf_df$REF_T >= 20 & vcf_df$REF_T <= 29 ~ "L",
    vcf_df$REF_T >= 30 ~ "VL"
  )
  vcf_df$ALT_polyT_cat = case_when(
    vcf_df$ALT_T <= 19 ~ "S",
    vcf_df$ALT_T >= 20 & vcf_df$ALT_T <= 29 ~ "L",
    vcf_df$ALT_T >= 30 ~ "VL"
  )
  
  vcf_df$TOMM40_Phased = paste0(vcf_df$REF_polyT_cat,"|",vcf_df$ALT_polyT_cat)
  vcf_df$tomm40_hap = NA
  vcf_df$tomm40_hap[vcf_df$TOMM40_Phased == "S|S"] = 1
  vcf_df$tomm40_hap[vcf_df$TOMM40_Phased %in% c("S|L","L|S")] = 2
  vcf_df$tomm40_hap[vcf_df$TOMM40_Phased %in% c("S|VL","VL|S")] = 3
  vcf_df$tomm40_hap[vcf_df$TOMM40_Phased == "L|L"] = 4
  vcf_df$tomm40_hap[vcf_df$TOMM40_Phased %in% c("L|VL","VL|L")] = 5
  vcf_df$tomm40_hap[vcf_df$TOMM40_Phased == "VL|VL" ] = 6
  return(vcf_df)
}
parse_GangSTR_tomm40_features <- function(vcf){
  vcf_df = na.omit(reshape2::melt(vcfR::extract.gt(vcf,return.alleles = T)))
  colnames(vcf_df) = c("ID","Sample","Alleles")
  
  vcf_df = vcf_df %>% 
    separate(Alleles, into = c("REF","ALT"), sep = "/") %>% 
    mutate(REF_T = str_count(toupper(REF), pattern = "T"), ALT_T = str_count(toupper(ALT), pattern = "T")) %>%
    mutate(nT_Alleles = paste0(REF_T,"|",ALT_T))
  
  ##INFO=<ID=END,Number=1,Type=Integer,Description="End position of variant">
  ##INFO=<ID=RU,Number=1,Type=String,Description="Repeat motif">
  ##INFO=<ID=PERIOD,Number=1,Type=Integer,Description="Repeat period (length of motif)">
  ##INFO=<ID=REF,Number=1,Type=Float,Description="Reference copy number">
  ##INFO=<ID=GRID,Number=2,Type=Integer,Description="Range of optimization grid">
  ##INFO=<ID=EXPTHRESH,Number=1,Type=Integer,Description="Threshold for calling expansions">
  ##INFO=<ID=STUTTERUP,Number=1,Type=Float,Description="Stutter model - up prob">
  ##INFO=<ID=STUTTERDOWN,Number=1,Type=Float,Description="Stutter model - down prob">
  ##INFO=<ID=STUTTERP,Number=1,Type=Float,Description="Stutter model - p">
  vcf.info = extract_info_tidy(vcf)
  info_GRID = str_split(vcf.info$GRID,",") %>% unlist() %>% as.numeric()
  info_STUTTERUP = vcf.info$STUTTERUP
  info_STUTTERDOWN = vcf.info$STUTTERDOWN
  info_STUTTERP = vcf.info$STUTTERP
  
  ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
  ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
  ##FORMAT=<ID=Q,Number=1,Type=Float,Description="Quality Score (posterior probability)">
  ##FORMAT=<ID=REPCN,Number=2,Type=Integer,Description="Genotype given in number of copies of the repeat motif">
  ##FORMAT=<ID=REPCI,Number=1,Type=String,Description="Confidence intervals">
  ##FORMAT=<ID=RC,Number=1,Type=String,Description="Number of reads in each class (enclosing, spanning, FRR, bounding)">
  ##FORMAT=<ID=ENCLREADS,Number=1,Type=String,Description="Summary of reads in enclosing class. Keys are number of copies and values show number of reads with that many copies.">
  ##FORMAT=<ID=FLNKREADS,Number=1,Type=String,Description="Summary of reads in flanking class. Keys are number of copies and values show number of reads with that many copies.">
  ##FORMAT=<ID=ML,Number=1,Type=Float,Description="Maximum likelihood">
  ##FORMAT=<ID=INS,Number=2,Type=Float,Description="Insert size mean and stddev">
  ##FORMAT=<ID=STDERR,Number=2,Type=Float,Description="Bootstrap standard error of each allele">
  ##FORMAT=<ID=QEXP,Number=3,Type=Float,Description="Prob. of no expansion, 1 expanded allele, both expanded alleles">
  vcf.format = extract_gt_tidy(vcf, verbose = F)
  gt_DP = vcf.format$gt_DP
  gt_Q = vcf.format$gt_Q
  gt_REPCN = str_split(vcf.format$gt_REPCN,",") %>% unlist()
  gt_REPCI = str_split(vcf.format$gt_REPCI,",") %>% unlist()
  gt_RC = str_split(vcf.format$gt_RC,",") %>% unlist()
  gt_ENCLREADS = str_split(vcf.format$gt_ENCLREADS,",") %>% unlist()
  gt_FLNKREADS = str_split(vcf.format$gt_FLNKREADS,",") %>% unlist()
  gt_ML = vcf.format$gt_ML
  gt_INS = str_split(vcf.format$gt_INS,",") %>% unlist()
  gt_STDERR = str_split(vcf.format$gt_STDERR,",") %>% unlist()
  gt_QEXP = str_split(vcf.format$gt_QEXP,",") %>% unlist()
  
  whichA1 = which.min(str_split(vcf_df$nT_Alleles,"\\|") %>% unlist())
  
  A1_info_GRID = as.numeric(info_GRID[whichA1])
  A2_info_GRID = as.numeric(info_GRID[-whichA1])
  gt_QEXP_1 = as.numeric(gt_QEXP[1])
  gt_QEXP_2 = as.numeric(gt_QEXP[2])
  gt_QEXP_3 = as.numeric(gt_QEXP[3])
  gt_RC_1 = as.numeric(gt_RC[1])
  gt_RC_2 = as.numeric(gt_RC[2])
  gt_RC_3 = as.numeric(gt_RC[3])
  gt_RC_4 = as.numeric(gt_RC[4])
  gt_QEXP_1 = as.numeric(gt_QEXP[1])
  gt_QEXP_2 = as.numeric(gt_QEXP[2])
  gt_QEXP_3 = as.numeric(gt_QEXP[3])
  A1_gt_REPCN = as.numeric(gt_REPCN[whichA1])
  A2_gt_REPCN = as.numeric(gt_REPCN[-whichA1])
  # if(gt_ENCLREADS[1] == "NULL"){
  #   A1_gt_ENCLREADS = NA
  #   A2_gt_ENCLREADS = NA
  # }else{
  #   A1_gt_ENCLREADS = as.numeric(gt_ENCLREADS[whichA1])
  #   A2_gt_ENCLREADS = as.numeric(gt_ENCLREADS[-whichA1])
  # }
  # if(gt_FLNKREADS[1] == "NULL"){
  #   A1_gt_FLNKREADS = NA
  #   A2_gt_FLNKREADS = NA
  # }else{
  #   A1_gt_FLNKREADS = as.numeric(gt_FLNKREADS[whichA1])
  #   A2_gt_FLNKREADS = as.numeric(gt_FLNKREADS[-whichA1])
  # }
  A1_gt_INS = as.numeric(gt_INS[whichA1])
  A2_gt_INS = as.numeric(gt_INS[-whichA1])
  A1_gt_STDERR = as.numeric(gt_STDERR[whichA1])
  A2_gt_STDERR = as.numeric(gt_STDERR[-whichA1])
  
  vcf.feats = data.frame(
    info_STUTTERUP,
    info_STUTTERDOWN,
    info_STUTTERP,
    gt_DP,
    gt_Q,
    gt_ML,
    gt_QEXP_1,
    gt_QEXP_2,
    gt_QEXP_3,
    gt_RC_1,
    gt_RC_2,
    gt_RC_3,
    gt_RC_4,
    A1_gt_REPCN,
    A2_gt_REPCN,
    # A1_gt_ENCLREADS,
    # A2_gt_ENCLREADS,
    # A1_gt_FLNKREADS,
    # A2_gt_FLNKREADS,
    A1_gt_INS,
    A2_gt_INS,
    A1_gt_STDERR,
    A2_gt_STDERR)
  
  return(vcf.feats)
}
parse_GangSTR_wrap <- function(vcfs, run_id){
  res_list = list()
  i = 1
  vcf_i = vcfs[i]
  vcf = read.vcfR(vcf_i, verbose = F)
    
  vcf_df = na.omit(reshape2::melt(vcfR::extract.gt(vcf,return.alleles = T)))
  colnames(vcf_df) = c("ID","Sample","Alleles")
  vcf_df$Sample = run_id

  vcf_df = suppressWarnings(parse_GangSTR_tomm40(vcf_df)) %>%
    mutate(run_id = run_id)
    
  vcf.feats = bind_cols(vcf_df[,c("run_id"),drop=F],
                        parse_GangSTR_tomm40_features(vcf))
  
  res = list()
  res$vcf = vcf_df
  res$vcf_feats = vcf.feats
  res_list[[run_id]] = res
  
  res_features = purrr::map_df(res_list,~{
    bind_rows(.x$vcf_feats)
  })

  res_list_GangSTR = res_list
  res_features_GangSTR = res_features
  
  res_GangSTR = purrr::map_df(res_list_GangSTR,~{
    bind_rows(.x$vcf)
  })
  
  return(list(res_list_GangSTR = res_list_GangSTR,
              res_features_GangSTR = res_features_GangSTR,
              res_GangSTR = res_GangSTR))
}

parse_hipstr_tomm40 <- function(vcf_df){
  vcf_df = vcf_df %>% 
    separate(Alleles, into = c("REF","ALT"), sep = "\\|") %>% 
    mutate(REF_T = str_count(REF, pattern = "T"), ALT_T = str_count(ALT, pattern = "T")) %>%
    mutate(nT_Alleles = paste0(REF_T,"|",ALT_T))
  
  vcf_df$REF_polyT_cat = case_when(
    vcf_df$REF_T <= 19 ~ "S",
    vcf_df$REF_T >= 20 & vcf_df$REF_T <= 29 ~ "L",
    vcf_df$REF_T >= 30 ~ "VL"
  )
  vcf_df$ALT_polyT_cat = case_when(
    vcf_df$ALT_T <= 19 ~ "S",
    vcf_df$ALT_T >= 20 & vcf_df$ALT_T <= 29 ~ "L",
    vcf_df$ALT_T >= 30 ~ "VL"
  )
  
  vcf_df$TOMM40_Phased = paste0(vcf_df$REF_polyT_cat,"|",vcf_df$ALT_polyT_cat)
  vcf_df$tomm40_hap = NA
  vcf_df$tomm40_hap[vcf_df$TOMM40_Phased == "S|S"] = 1
  vcf_df$tomm40_hap[vcf_df$TOMM40_Phased %in% c("S|L","L|S")] = 2
  vcf_df$tomm40_hap[vcf_df$TOMM40_Phased %in% c("S|VL","VL|S")] = 3
  vcf_df$tomm40_hap[vcf_df$TOMM40_Phased == "L|L"] = 4
  vcf_df$tomm40_hap[vcf_df$TOMM40_Phased %in% c("L|VL","VL|L")] = 5
  vcf_df$tomm40_hap[vcf_df$TOMM40_Phased == "VL|VL" ] = 6
  #table(vcf_df[,c("tomm40_hap","TOMM40_Phased")])
  
  return(vcf_df)
}
parse_hipstr_tomm40_features <- function(vcf){
  ##INFO=<ID=INFRAME_PGEOM,Number=1,Type=Float,Description="Parameter for in-frame geometric step size distribution">
  ##INFO=<ID=INFRAME_UP,Number=1,Type=Float,Description="Probability that stutter causes an in-frame increase in obs. STR size">
  ##INFO=<ID=INFRAME_DOWN,Number=1,Type=Float,Description="Probability that stutter causes an in-frame decrease in obs. STR size">
  ##INFO=<ID=OUTFRAME_PGEOM,Number=1,Type=Float,Description="Parameter for out-of-frame geometric step size distribution">
  ##INFO=<ID=OUTFRAME_UP,Number=1,Type=Float,Description="Probability that stutter causes an out-of-frame increase in read's STR size">
  ##INFO=<ID=OUTFRAME_DOWN,Number=1,Type=Float,Description="Probability that stutter causes an out-of-frame decrease in read's STR size">
  ##INFO=<ID=BPDIFFS,Number=A,Type=Integer,Description="Base pair difference of each alternate allele from the reference allele">
  ##INFO=<ID=START,Number=1,Type=Integer,Description="Inclusive start coodinate for the repetitive portion of the reference allele">
  ##INFO=<ID=END,Number=1,Type=Integer,Description="Inclusive end coordinate for the repetitive portion of the reference allele">
  ##INFO=<ID=PERIOD,Number=1,Type=Integer,Description="Length of STR motif">
  ##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
  ##INFO=<ID=REFAC,Number=1,Type=Integer,Description="Reference allele count">
  ##INFO=<ID=AC,Number=A,Type=Integer,Description="Alternate allele counts">
  ##INFO=<ID=NSKIP,Number=1,Type=Integer,Description="Number of samples not genotyped due to various issues">
  ##INFO=<ID=NFILT,Number=1,Type=Integer,Description="Number of samples whose genotypes were filtered due to various issues">
  ##INFO=<ID=DP,Number=1,Type=Integer,Description="Total number of valid reads used to genotype all samples">
  ##INFO=<ID=DSNP,Number=1,Type=Integer,Description="Total number of reads with SNP phasing information">
  ##INFO=<ID=DSTUTTER,Number=1,Type=Integer,Description="Total number of reads with a stutter indel in the STR region">
  ##INFO=<ID=DFLANKINDEL,Number=1,Type=Integer,Description="Total number of reads with an indel in the regions flanking the STR">
  vcf.info = extract_info_tidy(vcf)
  vcf.info_selected = vcf.info %>% select(-c(Key,START,END,PERIOD,AC,AN,NSKIP,NFILT,DSNP))
  # vcf.info_selected$A1_BPDIFFS = as.numeric(str_split(vcf.info$BPDIFFS, pattern = ",")[[1]][1])
  # vcf.info_selected$A2_BPDIFFS = as.numeric(str_split(vcf.info$BPDIFFS, pattern = ",")[[1]][2])
  
  ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
  ##FORMAT=<ID=GB,Number=1,Type=String,Description="Base pair differences of genotype from reference">
  ##FORMAT=<ID=Q,Number=1,Type=Float,Description="Posterior probability of unphased genotype">
  ##FORMAT=<ID=PQ,Number=1,Type=Float,Description="Posterior probability of phased genotype">
  ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Number of valid reads used for sample's genotype">
  ##FORMAT=<ID=DSNP,Number=1,Type=Integer,Description="Number of reads with SNP phasing information">
  ##FORMAT=<ID=PSNP,Number=1,Type=String,Description="Number of reads with SNPs supporting each haploid genotype">
  ##FORMAT=<ID=PDP,Number=1,Type=String,Description="Fractional reads supporting each haploid genotype">
  ##FORMAT=<ID=GLDIFF,Number=1,Type=Float,Description="Difference in likelihood between the reported and next best genotypes">
  ##FORMAT=<ID=DSTUTTER,Number=1,Type=Integer,Description="Number of reads with a stutter indel in the STR region">
  ##FORMAT=<ID=DFLANKINDEL,Number=1,Type=Integer,Description="Number of reads with an indel in the regions flanking the STR">
  ##FORMAT=<ID=AB,Number=1,Type=Float,Description="log10 of the allele bias pvalue, where 0 is no bias and more negative values are increasingly biased. 0 for all homozygous genotypes">
  ##FORMAT=<ID=FS,Number=1,Type=Float,Description="log10 of the strand bias pvalue from Fisher's exact test, where 0 is no bias and more negative values are increasingly biased. 0 for all homozygous genotypes">
  ##FORMAT=<ID=DAB,Number=1,Type=Integer,Description="Number of reads used in the AB and FS calculations">
  ##FORMAT=<ID=ALLREADS,Number=1,Type=String,Description="Base pair difference observed in each read's Needleman-Wunsch alignment">
  ##FORMAT=<ID=MALLREADS,Number=1,Type=String,Description="Maximum likelihood bp diff in each read based on haplotype alignments for reads that span the repeat region by at least 5 base pairs">
  ##FORMAT=<ID=FILTER,Number=1,Type=String,Description="Reason for filtering the current call, or PASS if the call was not filtered">
  vcf.format = extract_gt_tidy(vcf, verbose = F)
  
  vcf_df = na.omit(reshape2::melt(vcfR::extract.gt(vcf,return.alleles = T)))
  colnames(vcf_df) = c("ID","Sample","Alleles")
  
  if(nrow(vcf_df)>0){
    if(vcf_df$Alleles == "."){
      vcf_fix = getFIX(vcf)
      alt_i = str_split(vcf_fix[["ALT"]],",")[[1]]
      if(length(alt_i)>1){
        vcf_fix[["ALT"]] = str_flatten(rep("T",round(str_count(vcf_fix[["ALT"]],pattern = "T")/length(alt_i))))
      }
      vcf_df$Alleles = paste0(vcf_fix[["REF"]],"|",vcf_fix[["ALT"]])
    }
    vcf_df = suppressWarnings(parse_hipstr_tomm40(vcf_df)) 
    whichA1 = which.min(str_split(vcf_df$nT_Alleles,"\\|") %>% unlist())
    
    gt_GB = str_split(vcf.format$gt_GB,"\\|") %>% unlist() %>% as.numeric()
    A1_gt_GB = gt_GB[whichA1]
    A2_gt_GB = gt_GB[-whichA1]
    gt_PSNP = str_split(vcf.format$gt_PSNP,"\\|") %>% unlist() %>% as.numeric()
    A1_gt_PSNP = gt_PSNP[whichA1]
    A2_gt_PSNP = gt_PSNP[-whichA1]
    gt_PDP = str_split(vcf.format$gt_PDP,"\\|") %>% unlist() %>% as.numeric()
    A1_gt_PDP = gt_PDP[whichA1]
    A2_gt_PDP = gt_PDP[-whichA1]
    
    gt_ALLREADS = matrix(str_split(vcf.format$gt_ALLREADS,";") %>% unlist() %>% str_split("\\|") %>% unlist() %>% as.numeric(), ncol = 2)
    A1_gt_ALLREADS_sum = sum(gt_ALLREADS[,whichA1], na.rm = T)
    A1_gt_ALLREADS_mean = mean(gt_ALLREADS[,whichA1], na.rm = T)
    A2_gt_ALLREADS_sum = sum(gt_ALLREADS[,-whichA1], na.rm = T)
    A2_gt_ALLREADS_mean = mean(gt_ALLREADS[,-whichA1], na.rm = T)
    gt_MALLREADS = matrix(str_split(vcf.format$gt_MALLREADS,";") %>% unlist() %>% str_split("\\|") %>% unlist() %>% as.numeric(), ncol = 2)
    A1_gt_MALLREADS_sum = sum(gt_MALLREADS[,whichA1], na.rm = T)
    A1_gt_MALLREADS_mean = mean(gt_MALLREADS[,whichA1], na.rm = T)
    A2_gt_MALLREADS_sum = sum(gt_MALLREADS[,-whichA1], na.rm = T)
    A2_gt_MALLREADS_mean = mean(gt_MALLREADS[,-whichA1], na.rm = T)
    gt_GL_mean = mean(matrix(str_split(vcf.format$gt_GL,",") %>% unlist() %>% as.numeric()))
    gt_PL_mean = mean(matrix(str_split(vcf.format$gt_PL,",") %>% unlist() %>% as.numeric()))
    
    vcf.feats = data.frame(
      vcf.info_selected,
      vcf.format$gt_Q,
      vcf.format$gt_PQ,
      vcf.format$gt_DP,
      vcf.format$gt_DSNP,
      vcf.format$gt_GLDIFF,
      vcf.format$gt_DSTUTTER,
      vcf.format$gt_DFLANKINDEL,
      vcf.format$gt_AB,
      vcf.format$gt_FS,
      vcf.format$gt_DAB,
      A1_gt_ALLREADS_sum,
      A1_gt_ALLREADS_mean,
      A2_gt_ALLREADS_sum,
      A2_gt_ALLREADS_mean,
      A1_gt_MALLREADS_sum,
      A1_gt_MALLREADS_mean,
      A2_gt_MALLREADS_sum,
      A2_gt_MALLREADS_mean,
      gt_GL_mean,
      gt_PL_mean)
  }else{
    vcf.feats = unique(data.frame(
      t(setNames(rep(NA, length(colnames(vcf.info_selected))),colnames(vcf.info_selected))),
      gt_Q = NA,
      gt_PQ = NA,
      gt_DP = NA,
      gt_DSNP = NA,
      gt_GLDIFF = NA,
      gt_DSTUTTER = NA,
      gt_DFLANKINDEL = NA,
      gt_AB = NA,
      gt_FS = NA,
      gt_DAB = NA,
      A1_gt_ALLREADS_sum = NA,
      A1_gt_ALLREADS_mean = NA,
      A2_gt_ALLREADS_sum = NA,
      A2_gt_ALLREADS_mean = NA,
      A1_gt_MALLREADS_sum = NA,
      A1_gt_MALLREADS_mean = NA,
      A2_gt_MALLREADS_sum = NA,
      A2_gt_MALLREADS_mean = NA,
      gt_GL_mean = NA,
      gt_PL_mean = NA))
  }
  
  return(vcf.feats)
}
parse_hipstr_wrap <- function(vcfs, run_id){
  res_list = list()
  i = 1
  vcf_i = vcfs[i]
  vcf = read.vcfR(vcf_i, verbose = F)
  vcf_df = na.omit(reshape2::melt(vcfR::extract.gt(vcf, return.alleles = T)))
  if(nrow(vcf_df) > 0){
    colnames(vcf_df) = c("ID", "Sample", "Alleles")
    vcf_df$Sample = gsub("-DLPFC|\\.hg38\\.TOMM40", "", vcf_df$Sample)
    vcf_df$projid = gsub("(.*?)_(.*?)\\.(.*)", "\\2", run_id)
    vcf_df$projid = gsub("(.*?)_(.*?)", "\\2", run_id)
    
    vcf_df = suppressWarnings(parse_hipstr_tomm40(vcf_df)) %>%
      mutate(run_id = run_id)
    
    vcf_feats = bind_cols(vcf_df[, c("run_id"),drop=F],
                          parse_hipstr_tomm40_features(vcf))
    
    res = list()
    res$vcf = vcf_df
    res$vcf_feats = vcf_feats
    res_list[[run_id]] = res
  }

  res_features = purrr::map_df(res_list, ~bind_rows(.x$vcf_feats))
  res_HipSTR = purrr::map_df(res_list, ~bind_rows(.x$vcf))
  
  return(list(res_list_HipSTR = res_list,
              res_features_HipSTR = res_features,
              res_HipSTR = res_HipSTR))
}

collect_tool_feats <- function(dat_tool_list, tools){
  select = dplyr::select
  for(tool in tools){
    # tool = tools[1]
    tool_obj = dat_tool_list[[tool]]
    tool_res_features = tool_obj[[grep("res_features",names(tool_obj))]]
    tool_res = tool_obj[[paste0("res_",tool)]]
    
    tool_res %>% 
      group_by(run_id) %>%
      mutate(method = tool,
             A1 = min(REF_T,ALT_T),
             A2 = max(REF_T,ALT_T)) %>%
      select(run_id,method,A1,A2) %>% 
      left_join(tool_res_features, by = c("run_id")) %>% ungroup() %>%
      select(-c(method)) %>% as.data.frame() %>%
      column_to_rownames("run_id") %>% mutate_if(is.character, as.factor) %>% mutate_if(is.factor, as.numeric) %>%
      as.data.frame() -> dat_tool
    colnames(dat_tool) <- paste0(colnames(dat_tool),"_",tool)
    dat_tool = dat_tool %>% rownames_to_column("run_id")
    
    if(tool == tools[1]){
      tool_feats_df = dat_tool
    }else{
      tool_feats_df = full_join(tool_feats_df,dat_tool, by = "run_id")
    }
  }
  tool_feats_df %>% column_to_rownames("run_id") -> tool_feats_df
  return(tool_feats_df)
}

collect_kmer_feats <- function(datadir){
  kmer_files = list.files(path = datadir, pattern = ".polyT_kmer.txt", full.names = T, recursive = T)
  kmer_df = map_df(kmer_files, ~{
    datk = fread(.x, col.names = c("kmer","kmer_count")) %>% 
      mutate(projid = basename(.x) %>% gsub(".polyT_kmer.txt","",.))
    datk$projid = gsub("^R(.*)_(.*)","\\2",datk$projid)
    datk %>% pivot_wider(names_from = kmer, values_from = kmer_count) %>% as.data.frame()
  }) %>% column_to_rownames("projid") %>% as.data.frame()
  colnames(kmer_df) = paste0("kmer_",colnames(kmer_df))
  return(kmer_df)
}