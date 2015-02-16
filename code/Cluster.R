if (!require("hash")) {
  install.packages("hash")
  library(hash)
}

if (!require("sets")) {
  install.packages("sets")
  library(sets)
}

cluster_storage <- hash()

Cluster <- setClass("Cluster", slots = c(name = "character", elements = "set"),
                    prototype = list(name="New cluster", elements = c()))

write <- function(Cluster) UseMethod("write")

write.Cluster <- function(cluster, file_name) {
  value <- c(cluster@name, as.list(cluster@elements))  
  write.table(value, file = file_name, sep = ",",
            append = TRUE, row.names = FALSE, col.names=FALSE)
}

# Cluster string format:
# cluster_name, cluster_el1[, cluster_el2, ...]
#
parse.cluster <- function(cluster_line) {
  cluster_info <- strsplit(cluster_line, "[,]")
  parsed_cluster <- Cluster()
  parsed_cluster@name <- cluster_info[[1]]
  elements <- as.set(cluster_info[2:length(cluster_info)])
  parsed_cluster@elements <- elements
  
  return(parsed_cluster)
}

read.clusters <- function(file_name) {  
  for (line in readLines(file_name)) {
    cluster <- parse.cluster(line)
    cluster_storage[cluster.name] <- cluster
  }
}

