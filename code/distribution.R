if(!exists("code/clustering.R", mode="function")) {
  source("code/clustering.R")
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

resistance_probability <- function(medicine, resistance_status) {
  cases_number <- 0
  total_cases <- 0
  for (i in 1:nrow(phenotypes)) {
    status <- phenotypes[i, medicine]
    if (!is.na(status) & status == resistance_status) {
      cases_number <- cases_number + 1
    }
    if (!is.na(status)) {
      total_cases <- total_cases + 1
    }
  }
  probability <- cases_number / total_cases
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

snp_when_resistance_probability <- function(resistance_status, cluster, medicine) {
  cluster_cases <- cluster_occurences(cluster)
  if (is.null(cluster_cases) | length(cluster_cases) == 0) {
    probability <- 0
    
  } else {    
    genomes <- names(cluster_cases)
    total_cases <- 0  
    resistance_cases_number <- 0
    
    for (genome in genomes) {
      status <- phenotypes[genome, medicine]
      if (!is.na(status) & status == resistance_status) {
        resistance_cases_number <- resistance_cases_number + 1
      } 
      if (!is.na(status)) {
        total_cases <- total_cases + 1
      } 
    }  
    
    probability <- resistance_cases_number / total_cases
  }
  return(probability)
}

resistance_when_snp_probability <- function(resistance_status, cluster, medicine) {
  snp_when_resistance <- snp_when_resistance_probability(resistance_status, cluster, medicine)  
  
  cluster_probability <- snp_probability(cluster[1])  
  resistance_probability <- resistance_probability(medicine, resistance_status)
  
  probability <- 1.0 * snp_when_resistance * cluster_probability / resistance_probability

  return(probability)  
}

