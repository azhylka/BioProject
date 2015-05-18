start_index <- 0
finish_index <- 0

or <- function(x,y) {return(x|y)}

is.resistant <- function(status) {
  return (status == 2) # 2 indicates resistance
}


find.resistance.pairs <- function(medicine) {
  pairs <- hash() # key is a SNP with lower position
  pairs.df <- data.frame(c("key1", "key2"))
  
  snp_index <- 1
  while (snp_index < nrow(original_snps_data)-1) {
    vec1 <- original_snps_data[snp_index,]
    vec2 <- original_snps_data[snp_index+1,]
    
    result_resistant <- FALSE
    coincidence_border <- 10
    coincidences <- 0
    for (i in range(1, ncol(original_snps_data))) {
      pos1 <- vec1[i]
      pos2 <- vec2[i]
      if ((pos1 | pos2) & !xor(pos1, pos2)
          & is.resistant(phenotypes[i, medicine])) {
        result_resistance <- TRUE
      } else if (pos1 & pos2
                 & is.resistant(phenotypes[i, medicine])) {
        if (coincidences >= coincidence_border) {
          result_resistance <- FALSE
        }
        coincidences <- coincidences + 1
      }
    }
    if (result_resistant) {
      snp1 <- original_snps_data$snp_pos[snp_index]
      snp2 <- original_snps_data$snp_pos[snp_index+1]
      print("####")
      print(snp1)
      pairs[snp1] <- snp2
      pairs.df <- rbind(pairs.df, c(snp1, snp2))
      print(snp_index)
      snp_index <- snp_index + 2
    } else {
      snp_index <- snp_index + 1
    }
  }  
  return(pairs.df)
}

