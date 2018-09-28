#' Size parameter for negative binomial distribution
#'
#' @param A1 A number. Aggregation in one side of hybrid index
#' @param A2 A number. Aggregation in the other side of hybrid index
#' @param Z A number
#' @param hybridIndex A number
#' @return The size parameter of negative binomial at one given point of \code{hybridIndex}
#' @export

SizeNegBin <- function(A1, A2, Z, hybridIndex){
  aggregation <- Aggregation(A1, A2, Z, hybridIndex)
  aggregation <- sapply(aggregation, function(x) {
    return(max(x, 0.01))
  })
  size <- 1/aggregation
  return(size)
}
