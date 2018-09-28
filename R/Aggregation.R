#' Aggregation parameter for negative binomial distribution
#'
#' @param A1 A number. Aggregation in one side of hybrid index
#' @param A2 A number. Aggregation in the other side of hybrid index
#' @param Z A number
#' @param hybridIndex A number
#' @return The aggregation parameter of negative binomial at one given point of \code{hybridIndex}
#' @export

Aggregation <- function(A1, A2, Z, hybridIndex){
  heterozygoty <- 2 * hybridIndex * (1 - hybridIndex)
  aggregation <- (A1 + (A2 - A1) * hybridIndex) + Z * heterozygoty
  return(aggregation)
}

