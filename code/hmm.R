if (!require("futile.logger")) {
  install.packages("futile.logger")
  library(futile.logger)
}

if (!require("depmixS4")) {
  install.packages("depmixS4")
  library("depmixS4")
}

if (!require("PLIS")) {
  install.packages("PLIS")
  library(PLIS)
}

if(!exists("code/filtering.R", mode="function")) {
  source("code/filtering.R")
}


flog.threshold(INFO)
log_file_name <- sprintf("output/logs/%s/hmm.log", format(Sys.time(), "%Y-%m-%d"))
flog.appender(appender.file(), name="hmm")
phenotype_file <- "resources/adopted_TB_Phenotypes.csv"
phenotypes <- read.csv(phenotype_file, head=TRUE, row.names=1, sep="\t")

snps_csv <- "resources/finally-proper-snps-matrix.csv"
snps_data <- read.csv(snps_csv, head=TRUE, sep="\t", row.names=1)
filtered_snps_data <- minor.allele.frequency.filter(data = snps_data, minor_allele = 1)

get_snp_case_control_matrix <- function(snp_index, medicine) {  
  #flog.info("Calculate case-control matrix for SNP index %d", snp_index)
  result_matrix <- matrix(c(0, 0, 0, 0), nrow = 2, ncol = 2)
  strains <- rownames(phenotypes)
  i <- 1
  for (strain in strains) {                  
    if (strain %in% names(filtered_snps_data)) {            
      if (phenotypes[i, medicine] == 1) {
        if (filtered_snps_data[snp_index, strain] == 0) {
          result_matrix[1,2] <- result_matrix[1,2] + 1
        } else {
          result_matrix[2,2] <- result_matrix[2,2] + 1
        }
      } else {
        if (filtered_snps_data[snp_index, strain] == 0) {
          result_matrix[1,1] <- result_matrix[1,1] + 1
        } else {
          result_matrix[2,1] <- result_matrix[2,1] + 1
        }
      }     
    }
    i <- i + 1
  }
  return(result_matrix)
}

get_snp_stat_info <- function(medicine) {
  flog.info("Compute z-value for %s", medicine)
  snp_info.df <- data.frame()  
  for (snp_index in 1:nrow(filtered_snps_data)) {
    snp_matrix <- get_snp_case_control_matrix(snp_index, medicine)
    test_result <- fisher.test(snp_matrix)
    p_value <- test_result$p.value
    odds_ratio <- test_result$estimate  
    if (odds_ratio > 1) {
      z_value <- 1.0 / pnorm(1 - p_value * 0.5, mean = 0, sd = 1)
    } else {
      z_value <- 1.0 / pnorm(p_value * 0.5, mean = 0, sd = 1)
    }
    snp_info.df <- rbind(snp_info.df, c(filtered_snps_data[snp_index,1], p_value, odds_ratio, z_value))
  }
  names(snp_info.df) <- c("snp", "p_value", "odds_ratio", "z_value")
  return(snp_info.df)
}

get_hypothesis_states <- function(snp_info, fdr = 0.001) {
  flog.info("Analysing")
  BIC <- .Machine$double.xmax
  dimension <- 0
  final_result <- NULL
  for (i in 2:6) {
    results <- em.hmm(snp_info$z_value, L = i)
    flog.info("Converged - %s, BIC = %f",results$converged, results$BIC)
    flog.info("Converged - %s, BIC = %f",results$converged, results$BIC, name='hmm')
    if (results$converged && results$BIC < BIC) {
      BIC <- results$BIC
      dimension <- i
      final_result <- results
    }
  }  
  if (final_result$converged) {
    plis_result <- plis(final_result$LIS, fdr = fdr)
    return(plis_result$States)
  }  
}

out_dir <- sprintf("output/hmm_results/%s/", format(Sys.time(), "%Y-%m-%d"))

if(!file.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

flog.info("ETHA")
flog.info("ETHA", name='hmm')
etha_snp_info <- get_snp_stat_info("ETHA")
etha_states <- get_hypothesis_states(etha_snp_info, fdr = 0.0001)
if (!is.null(etha_states)) {
  write.csv(etha_states, file = paste(out_dir, "etha_hmm_plis.csv"),  row.names = etha_snp_info$snp)
}

flog.info("ISON")
flog.info("ISON", name='hmm')
ison_snp_info <- get_snp_stat_info("ISON")
ison_states <- get_hypothesis_states(ison_snp_info)
if (!is.null(ison_states)) {
  write.csv(ison_states, file = paste(out_dir, "ison_hmm_plis.csv"),  row.names = ison_snp_info$snp)
}

flog.info("STRE")
flog.info("STRE", name='hmm')
stre_snp_info <- get_snp_stat_info("STRE")
stre_states <- get_hypothesis_states(stre_snp_info)
if (!is.null(stre_states)) {
  write.csv(stre_states, file = paste(out_dir, "stre_hmm_plis.csv"),  row.names = stre_snp_info$snp)
}

flog.info("ETHI")
flog.info("ETHI", name='hmm')
ethi_snp_info <- get_snp_stat_info("ETHI")
ethi_states <- get_hypothesis_states(ethi_snp_info)
if (!is.null(ethi_states)) {
  write.csv(ethi_states, file = paste(out_dir, "ethi_hmm_plis.csv"),  row.names = ethi_snp_info$snp)
}

flog.info("PARA")
flog.info("PARA", name='hmm')
para_snp_info <- get_snp_stat_info("PARA")
para_states <- get_hypothesis_states(para_snp_info)
if (!is.null(para_states)) {
  write.csv(para_states, file = paste(out_dir, "para_hmm_plis.csv"),  row.names = para_snp_info$snp)
}

flog.info("AMIK")
flog.info("AMIK", name='hmm')
amik_snp_info <- get_snp_stat_info("AMIK")
amik_states <- get_hypothesis_states(amik_snp_info)
if (!is.null(amik_states)) {
  write.csv(amik_states, file = paste(out_dir, "amik_hmm_plis.csv"),  row.names = amik_snp_info$snp)
}

flog.info("OFLO")
flog.info("OFLO", name='hmm')
oflo_snp_info <- get_snp_stat_info("OFLO")
oflo_states <- get_hypothesis_states(oflo_snp_info)
if (!is.null(oflo_states)) {
  write.csv(oflo_states, file = paste(out_dir, "oflo_hmm_plis.csv"),  row.names = oflo_snp_info$snp)
}

