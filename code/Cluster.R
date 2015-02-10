if (!require("hash")) {
  install.packages("hash")
  library(hash)
}

if (!require("sets")) {
  install.packages("sets")
  library(sets)
}

cluster_storage <- hash()

Cluster <- setClass("Cluster", slots = c(name = "character", elements = "list"),
                    prototype = list(name="New cluster", elements = set()))

write <- function(cluster) UseMethod("write")

write.Cluster <- function(cluster, file_name) {
  write.csv(c(cluster.name, cluster.elements), file = file_name,
            append = TRUE, row.names = FALSE)
}

parse.cluster <- function(cluster_line) {
  cluster_info <- strsplit(cluster_line, "[,]")
  parsed_cluster <- Cluster()
  parsed_cluster.name <- cluster_info[[1]]
  elements <- as.set(cluster_info[2:length(cluster_info)])
  parsed_cluster.elements <- elements
  
  return(parsed_cluster)
}

read.clusters <- function(file_name) {  
  for (line %in% readLines(file_name)) {
    cluster <- parse.cluster(line)
    cluster_storage[cluster.name] <- cluster
  }
}

