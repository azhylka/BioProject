
if(!exists("code/distribution.R", mode="function")) {
  source("code/distribution.R")
}

flog.threshold(INFO)
flog.appender(appender.file("output/logs/distribution.log"), name="distribution")

snps_csv <- "resources/full_proper_snps_matrix.csv"
original_snps_data <- read.csv(snps_csv, head=TRUE, sep=",")
snps_sums <- colSums(snps_data)
size <- length(snps_sums)


phenotype_file <- "resources/adopted_TB_Phenotypes.csv"
phenotypes <- read.csv(phenotype_file, head=TRUE, row.names=1, sep="\t")

clustering_result <- extract.clusters()
clusters <- clustering_result$clusters
subsetting_mask <- clustering_result$mask
snps_data <- original_snps_data[subsetting_mask, ]

medicine <- "ETHA"
clusters_number <- length(clusters)
resistance_probabilities <- rep(0, clusters_number)

cluster_index <- 0
for (cluster in clusters) {
  cluster_index <- cluster_index + 1
  # cluster name is always formed of the first element in cluster that is saved in filtered dataset  
  cluster_representative <- strtoi(cluster@name)
  # create list of representative in order to reuse code
  resistance_probabilities[cluster_index] <- 
    resistance_when_snp_probability(2, c(cluster_representative), medicine)  
}

plot(resistance_probabilities)
start_index <- 1
plot_step <- 40
while (start_index < clusters_number) {
  next_index <- start_index + plot_step
  next_index <- ifelse(next_index > clusters_number, clusters_number, next_index)
  png(file = paste("output/graphics/resistance_probability_", next_index, ".png", sep=""))
  plot(start_index:next_index, resistance_probabilities[start_index:next_index], type = "l", xlab = "Clusters", ylab = "Resistance Probability")
  dev.off()
  start_index <- next_index + 1 
}
