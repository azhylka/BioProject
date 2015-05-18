cluster_points <- rep(0, length(clusters))

for (i in 1:length(clusters)) {
 cluster_points[i] <- strtoi(clusters[[i]]@name)
}

snps_data <- original_snps_data[match(cluster_points, original_snps_data$snp_pos),]

phenotype_file <- "resources/adopted_TB_Phenotypes.csv"
phenotypes <- read.csv(phenotype_file, head=TRUE, row.names=1, sep="\t")
resist_status <- 2

drug <- "ETHA"

snps_when_resist_probabilities <- rep(0, nrow(snps_data))
names(snps_when_resist_probabilities) <- c(as.character(cluster_points))

for (snp in 1:nrow(snps_data)) {
  snp_when_resist_occur <- 0
  snp_occur <- 0
  for (genome in 2:ncol(snps_data)) {
    genome_name <- colnames(snps_data)[genome]
    if (snps_data[snp, genome] == 1 & !is.na(phenotypes[genome_name,drug])
        & phenotypes[genome_name,drug] == resist_status) {
      snp_when_resist_occur <- snp_when_resist_occur + 1      
    } 
    if (snps_data[snp, genome] == 1) {
      snp_occur <- snp_occur + 1
    }
  }
  snps_when_resist_probabilities[snp] <- snp_when_resist_occur / snp_occur  
}

plot(snps_when_resist_probabilities)
start_index <- 1
plot_step <- 40
clusters_number <- length(clusters)
while (start_index < length(clusters)) {
  next_index <- start_index + plot_step
  next_index <- ifelse(next_index > clusters_number, clusters_number, next_index)
  png(file = paste("output/graphics/18-05-2015/resistance_probability_", next_index, ".png", sep=""))
  plot(start_index:next_index, snps_when_resist_probabilities[start_index:next_index], type = "l", xlab = "Clusters", ylab = "Resistance Probability")
  dev.off()
  start_index <- next_index + 1 
}
