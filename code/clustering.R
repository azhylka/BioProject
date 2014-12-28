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
