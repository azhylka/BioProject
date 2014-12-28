minor.allele.frequency.filter <- function(data, minor_allele) {
  filter_border <- 0.01
  to_remove <- c()
  genomes <-ncol(data)-1
  for (index in 1:nrow(data)) {
    allele_number <- Reduce(function(x, y){return(x + ifelse(y==minor_allele, 1, 0))}, data[index,], 0)
    frequency <- allele_number/genomes
    if (frequency < filter_border) {
      to_remove <- c(to_remove, index)
    }
  }
  return(data[-to_remove,])
}