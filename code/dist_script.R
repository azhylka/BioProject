snps_csv <- "/home/anjenson/Yandex.Disk/Курсовая 4 курс/finally-proper-snps-matrix.csv"
snps_data <- read.csv(snps_csv, head=TRUE, sep="\t", row.names=1)
snps_sums <- colSums(snps_data)
size <- length(snps_sums)


logged_sums <- sapply(snps_sums[2:size], log2)

plot(logged_sums)
title("Number of mutations per strain")

nearest_snp_distance <- function(snps_vector) {
  result <- rep(.Machine$integer.max, length(snps_vector))  
  prev_snp <- match(1, snps_vector)
  result[1:(prev_snp-1)] <- 0
  for (i in (prev_snp+1):length(snps_vector)) {
    if (snps_vector[i] == 1) {    
      result[i] <- snps_data$snp_pos[i] - snps_data$snp_pos[prev_snp]
      if (i != prev_snp) {
        result[prev_snp] <- pmin.int(result[prev_snp], result[i])
      } else {
        result[prev_snp] <- result[i]
      }
      prev_snp = i
    } else {
      result[i] <- 0
    }
  }
  return(result)
}

left_mutations_number <- function(snps_vector) {
  vec <- rep(0, length(snps_vector))
  sum <- 0
  for (i in 1:length(snps_vector)) {
    if (snps_vector[i] == 1) {
      sum <- sum + 1
    }
    vec[i] = sum
  }
  return(vec)
}

plot_mutation_distances_lines <- function() {
  for (vector in snps_data[,2:size]) {    
    dist_vec <- nearest_snp_distance(vector)
    plot(dist_vec, type="l")
  }
}

mutations_per_position <- function() {
  snps_per_position <- rowSums(snps_data[,2:size])
  # plot by clusters
  cluster_size <- 100
  right_bound <- cluster_size
  rows_num <- length(snps_per_position)
  cluster_number <- rows_num/cluster_size  
  for (i in 1:cluster_number) {
    left_bound <- right_bound - cluster_size + 1
    plot(snps_per_position[left_bound:right_bound],
         snps_data$snp_pos[left_bound:right_bound], type="l")
    
    right_bound <- right_bound + cluster_size    
    if (right_bound > rows_num) {
      right_bound <- rows_num
    }
  }
}

mutations_per_position()

#EMMA P-value filter
significant_snps_stat <- function() {
  significant_snps <- as.vector(unlist(
    read.csv(file="/home/anjenson/Yandex.Disk/Курсовая 4 курс/singif_snps_p_filter.csv",
                                         head=FALSE, sep=",")))
  snps_pos <- snps_data$snp_pos
  pos_index <- 1
  significant_snps_number <- c()
  for (pos in snps_pos) {
    if (pos %in% significant_snps) {    
      significant_snps_number <- c(significant_snps_number, sum(snps_data[pos_index, 2:size]))
    }
    pos_index <- pos_index + 1
  }
  hist(significant_snps_number, main = "Significant SNPs number histogram")
  print(significant_snps_number)
}

significant_snps_stat()

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

get_snp_frequency_in_cluster

significant_snps <- as.vector(unlist(
  read.csv(file="/home/anjenson/Yandex.Disk/Курсовая 4 курс/singif_snps_p_filter.csv",
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

find_genomes_by_snp <- function(snp_index) {
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
      mutated_genomes <- find_genomes_by_snp(index)
      cluster <- retrieve_cluster(mutated_genomes)
 #     snp_positions <- mapply(xor, snp_positions, filter_significant_snps(cluster))
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

all_clusters <- get_all_clusters()
length(all_clusters)
Reduce(intersect, all_clusters)

max_len <- 0
for (cluster in all_clusters) {
  cluster_length <- length(cluster)
  if (cluster_length > max_len) {    
    max_len <- cluster_length
  }
} 

all_cluster.df <- data.frame(row.names=1:max_len)
for (cluster in all_clusters) {  
  if (length(cluster) < max_len) {
    cluster[max_len] <- 0    
  }  
  all_cluster.df <- cbind(all_cluster.df, sort(cluster, na.last=TRUE))
}

names(all_cluster.df) <- paste("Cluster #", c(1:length(all_clusters)))

write.csv(all_cluster.df, file="clusters.csv")
