# Event SNPi happens if the snps at position SNPi is a mutation in the genome

#TODO; Optimise performance of cluster intersection cases calculation

setwd("~/Yandex.Disk/Курсовая 4 курс/")

if (!require("futile.logger")) {
  install.packages("futile.logger")
  library(futile.logger)
}

if(!exists("clustering.R", mode="function")) {
  source("clustering.R")
}

flog.threshold(INFO)
flog.appender(appender.file("distribution.log"), name="distribution")

snps_csv <- "finally-proper-snps-matrix.csv"
snps_data <- read.csv(snps_csv, head=TRUE, sep="\t", row.names=1)
snps_sums <- colSums(snps_data)
size <- length(snps_sums)


phenotype_file <- "adopted_TB_Phenotypes.csv"
phenotypes <- read.csv(phenotype_file, head=TRUE, row.names=1, sep="\t")

all_clusters <- get_all_clusters()
core_cluster <- get_core_cluster();

num_of_clusters <- length(all_clusters)
non_intersecting_clusters <- sapply(all_clusters, function(x){setdiff(x, core_cluster)})
non_intersecting_clusters[num_of_clusters+1] <- list(core_cluster)

cluster_statistics <- data.frame()

#update number of clusters
num_of_clusters <- length(non_intersecting_clusters)
cluster_probability <- rep(0, num_of_clusters)
#compute cluster probability
i <- 0
for (cluster in non_intersecting_clusters) {
  cluster_probability[i] <- length(cluster) / nrow(snps_data)
}

position_to_index <- function(snps_positions) {
  indexes <- sapply(snps_positions, function(x){match(x, snps_data$snp_pos)})
  return(indexes)
}

snps_set_probability <- function(snps) {
  genomes <- find_genomes_by_snp(snps[1]) # take any snp
  indexes <- position_to_index(snps)
  num_of_intersections <- 0
  for (genome in genomes) {
    intersect <- TRUE
    for (index in indexes) {
      if (genome[index] == 0) {
        intersect <- FALSE
      }
    }
    if (intersect) {
      num_of_intersections <- num_of_intersections + 1
    }
  }
  set_probability <- num_of_intersections / (ncol(snps_data) - 1)
  return(set_probability)
}

snp_probability <- function(snp) {
  genomes <- find_genomes_by_snp(snp)
  probability <- length(genomes) / (ncol(snps_data) - 1)
  return(probability)
}

cluster_occurences <- function(cluster) {
  genomes <- find_genomes_by_snp(cluster[1])
  indexes <- position_to_index(cluster)
  cluster_cases_number <- 0
  cluster_cases <- c()
  genome_index <- 0
  for (genome in genomes) {
    genome_index <- genome_index + 1
    cluster_presents <- TRUE
    for (index in indexes) {
      if (genome[index] == 0) {
        cluster_presents <- FALSE      
      }
    }
    if (cluster_presents) {
      cluster_cases_number <- cluster_cases_number + 1 
      cluster_cases[cluster_cases_number] <- list(genome)
      names(cluster_cases)[cluster_cases_number] <- names(genomes)[genome_index]
    }
  }
  return(cluster_cases)
}

find_mutation_cases <- function(resistance, genomes, medicine) {
  cases_number <- 0
  for (genome in genomes) {
    status <- phenotypes[genome, medicine]
    if (!is.na(status) & status == resistance) {
      cases_number <- cases_number + 1
    }
  }
  return(cases_number)
}

resistance_conditional_probability <- function(resistance_status, cluster, medicine) {
  cluster_cases <- cluster_occurences(cluster)
  mutation_cases_number <- find_mutation_cases(resistance_status, names(cluster_cases), medicine)
  if (is.null(cluster_cases) | length(cluster_cases) == 0) {
    probability <- 0
  } else {
    probability <- mutation_cases_number / length(cluster_cases)
  }
  return(probability)
}


medicine <- "ETHA"
resistance_probabilities <- rep(0, num_of_clusters)
cluster_index <- 0
for (cluster in non_intersecting_clusters) {
  cluster_index <- cluster_index + 1
  resistance_probabilities[cluster_index] <- resistance_conditional_probability(2, cluster, medicine)
}

plot(resistance_probabilities)
plot(resistance_probabilities, type = "b", xlab = "Clusters", ylab = "Resistance Probability")
