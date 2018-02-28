## Functions defining the distribution of mu and 1/k of the Negative binomial distribution
MeanLoad <- function(L1, L2, alpha, hybridIndex){
  heterozygoty <- 2 * hybridIndex * (1 - hybridIndex)
  mean <- (L1 + (L2 - L1) * hybridIndex) * (1 - alpha * heterozygoty)
  mean <- sapply(mean, function(x) {
    return(max(x, 0.01))
  })
  return(mean)
}

Aggregation <- function(A1, A2, Z, hybridIndex){
  heterozygoty <- 2 * hybridIndex * (1 - hybridIndex)
  aggregation <- (A1 + (A2 - A1) * hybridIndex) + Z * heterozygoty 
  return(aggregation)
} 

SizeNegBin <- function(A1, A2, Z, hybridIndex){
  aggregation <- Aggregation(A1, A2, Z, hybridIndex)
  aggregation <- sapply(aggregation, function(x) {
    return(max(x, 0.01))
  })
  size <- 1/aggregation
  return(size)
} 