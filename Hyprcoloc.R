#!/usr/bin/env Rscript
library(data.table)
library(hyprcoloc)
library(dplyr)  
gwas1_path <- "/scratch/users/s/h/shifang/data/GCST90728570_format.tsv"
gwas2_path <- "/scratch/users/s/h/shifang/data/GCST90627749_format.tsv"
gwas1 <- fread(gwas1_path, header = TRUE)
gwas2 <- fread(gwas2_path, header = TRUE)

gwas1 <- gwas1[, .(MarkerID, beta, standard_error)]
setnames(gwas1, c("beta", "standard_error"), c("beta_gwas1", "se_gwas1"))

gwas2 <- gwas2[, .(MarkerID, beta, sebeta)]
setnames(gwas2, c("beta", "sebeta"), c("beta_gwas2", "se_gwas2"))

gwas1 <- gwas1[se_gwas1 > 0 & se_gwas1 < 1]
gwas2 <- gwas2[se_gwas2 > 0 & se_gwas2 < 1]

eqtl_dir <- "/scratch/users/s/h/shifang/coloc/tenk10k/tenk10K_eQTL_sel"
eqtl_files <- list.files(path = eqtl_dir, pattern = "\\.txt$", full.names = TRUE)

if (length(eqtl_files) == 0) {
  stop("No find .txt：", eqtl_dir)
}

coloc_results <- data.frame()

# -----------------------------------------------------------------------------
# 3. In circle
# -----------------------------------------------------------------------------
for (eqtl_file in eqtl_files) {
  eqtl <- fread(eqtl_file, header = TRUE)
  eqtl <- eqtl[, .(MarkerID, BETA, SE)]
  setnames(eqtl, c("BETA", "SE"), c("beta_eqtl", "se_eqtl"))
  eqtl <- eqtl[se_eqtl > 0 & se_eqtl < 1]

  common_snps <- Reduce(intersect, list(gwas1$MarkerID, gwas2$MarkerID, eqtl$MarkerID))
  if (length(common_snps) == 0) {
    warning("file ", basename(eqtl_file), " No common SNP，skip")
    next
  }
  
  # 
  gwas1_sub <- gwas1[MarkerID %in% common_snps][order(MarkerID)]
  gwas2_sub <- gwas2[MarkerID %in% common_snps][order(MarkerID)]
  eqtl_sub  <- eqtl[MarkerID %in% common_snps][order(MarkerID)]
  
  # 
  betas <- cbind(
    GCST90728570 = gwas1_sub$beta_gwas1,
    GCST90627749 = gwas2_sub$beta_gwas2,
    !!basename(eqtl_file) := eqtl_sub$beta_eqtl  
  )
  
  ses <- cbind(
    GCST90728570 = gwas1_sub$se_gwas1,
    GCST90627749 = gwas2_sub$se_gwas2,
    !!basename(eqtl_file) := eqtl_sub$se_eqtl
  )
  
  betas <- as.matrix(betas)
  ses  <- as.matrix(ses)
  
  # run hyprcoloc
  result <- hyprcoloc(
    betas,
    ses,
    trait.names = colnames(betas),
    snp.id      = common_snps
  )
  
  # 
  if (!is.null(result$results) && nrow(result$results) > 0) {
    result_df <- as.data.frame(result$results)
    result_df$ID <- basename(eqtl_file)   
    coloc_results <- bind_rows(coloc_results, result_df)
  } else {
    warning("file ", basename(eqtl_file), "No coloc results")
  }
}

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
output_file <- "hycoloc_all.csv"
if (nrow(coloc_results) > 0) {
  fwrite(coloc_results, output_file)
} else {
  message("No results for colocalisation")
}
