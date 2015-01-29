if (!require("hash")) {
  install.packages("hash")
  library(hash)
}

if (!require("sets")) {
  install.packages("sets")
  library(sets)
}

or <- function(x,y) {return(x|y)}

and <- function(x, y) {return(x&y)}

get_all_snp_pos <- function() {
  snp_pos_vec <- rep(0,nrow(snps_data)) 
  for (genome in snps_data[,2:size]) {    
    snp_pos_vec <- mapply(or, snp_pos_vec, genome)
  }
  
  return(snp_pos_vec)
  #print(snp_pos_vec)
}

get_snp_position_labels <- function(snp_indicators) {
  num_of_pos <- 0
  for (indicator in snp_indicators) {
    if (indicator) {
      num_of_pos <- num_of_pos + 1
    }
  }
  positions <- numeric(num_of_pos)
  index <- 1
  label_index <- 1
  for (indicator in snp_indicators) {    
    if (indicator) {
      positions[index] <- snps_data$snp_pos[label_index]
      index <- index + 1
    }
    label_index <- label_index + 1
  }
  return(positions)
}

generate_next_base_vector <- function(index) {  
  for (column in 2:ncol(snps_data)) {
    if (snps_data[index, column] == 1) {      
      return(snps_data[,column])
    }
  }
  return(rep(0, nrow(snps_data)))
}

significant_snps <- as.vector(unlist(
  read.csv(file="resources/signif_snps_p_filter.csv",
           head=FALSE, sep=",")))

filter_significant_snps <- function(snps) {  
  filtered_snps <- snps
  snp_positions <- snps_data$snp_pos
  for (index in 1:length(snp_positions)) {
    if (!(snp_positions[index] %in% significant_snps)) {
      filtered_snps[index] = FALSE
    }
  }  
  return(filtered_snps)
}

get_core_cluster <- function() {
  cluster <- Reduce(function(x, y){mapply(and, x, y)}, snps_data[,2:size])                      
  return(get_snp_position_labels(cluster))      
}

find_mutated_genomes_indeces <- function(snp_index) {
  indeces <- c()
  for (column in 2:size) {
    if (snps_data[snp_index, column] == 1) {
      indeces <- c(indeces, column)
    } 
  }
  return(indeces)
}

find_genomes_by_snp <- function(snp) {
  snp_index <- match(snp, snps_data$snp_pos)
  return(find_genomes_by_snp_index(snp_index))
}

find_genomes_by_snp_index <- function(snp_index) {
  return(snps_data[find_mutated_genomes_indeces(snp_index)])
}

get_significant_core_snps <- function(core_cluster) {
  significant_core <- intersect(significant_snps, core_cluster)
  return(significant_core)
}

remove_core_snps <- function(cluster, core_cluster) {
  return(setdiff(cluster, core_cluster))
}

retrieve_cluster <- function(genomes) {
  cluster <- rep(1, nrow(genomes))
  for (i in 1:ncol(genomes)) {
    cluster <- mapply(and, genomes[,i], cluster)
  }
  return(cluster)
}

get_all_clusters <- function() {
  snp_positions <- get_all_snp_pos()
  snp_positions <- filter_significant_snps(snp_positions)
  core_cluster <- get_core_cluster()
  cluster_storage <- list()
  storage_size <- 0  
  len <- length(snp_positions)
  for (index in 1:len) {
    if (snp_positions[index]) {
      mutated_genomes <- find_genomes_by_snp_index(index)
      cluster <- retrieve_cluster(mutated_genomes)
      #      snp_positions <- mapply(xor, snp_positions, filter_significant_snps(cluster))
      if (length(remove_core_snps(get_snp_position_labels(cluster), core_cluster)) > 0) {
        storage_size <- storage_size + 1
        cluster_storage[[storage_size]] <- remove_core_snps(get_snp_position_labels(cluster),
                                                            core_cluster)
        names(cluster_storage)[[storage_size]] <- paste("Cluster of SNP", snps_data$snp_pos[index])
      }
    }
  }
  return(cluster_storage)
}

find.cluster <- function(target_snp_index, all_snps, processed_snps) {
  print(paste("Building cluster for ", target_snp_index))
  columns <- ncol(snps_data)
  cluster <- list(snps_data$snp_pos[target_snp_index])  
  for (index in 1:length(all_snps)) {
    snp <- all_snps[index]
    if (!has.key(snp, processed_snps)) {        
      equal = sum(mapply(xor, snps_data[index, 2:columns], snps_data[target_snp_index, 2:columns])) == 0
      if (equal) {
        cluster <- c(cluster, snps_data$snp_pos[index])
        processed_snps[snp] <- TRUE
      }
    }
  }
  return(cluster)
}

extract.clusters <- function() {
  # regard this hash as hash set
  processed_snps <- hash()
  clusters <- set()
  all_snps <- as.character(snps_data$snp_pos)  
  subsetting_mask <- rep(TRUE, length(all_snps))
  for (index in 1:length(all_snps)) {  
    snp <- all_snps[index]
    if (!has.key(snp, processed_snps)) {
      processed_snps[snp] <- TRUE      
      new_cluster <- find.cluster(index, all_snps, processed_snps)
      clusters <- set_union(clusters, new_cluster)
    } else {
      subsetting_mask[index] <- FALSE
    }
  }
  return(list(clusters = clusters, mask = subsetting_mask))
}

result <- extract.clusters()
clusters <- result$clusters
subsetting_mask <- result$mask
clustered_snps_data <- snps_data[subsetting_mask, ]
write.table(clustered_snps_data, file="resources/clustered_snps_data.csv", sep=",", 
            col.names=TRUE, row.names=FALSE)
