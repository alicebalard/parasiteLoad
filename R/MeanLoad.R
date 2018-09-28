#' Function defining the expected load at one given point of the index
#'
#' @param L1 A number. Aggregation in one side of hybrid index
#' @param L2 A number. Aggregation in the other side of hybrid index
#' @param alpha A number. Curvature of the expected load along hybrid index
#' @param hybridIndex A number
#' @return The expected load at one given point of \code{hybridIndex}
#' @export

##
MeanLoad <- function(L1, L2, alpha, hybridIndex){
  heterozygoty <- 2 * hybridIndex * (1 - hybridIndex)
  mean <- (L1 + (L2 - L1) * hybridIndex) * (1 - alpha * heterozygoty)
  mean <- sapply(mean, function(x) {
    return(max(x, 0.01))
  })
  return(mean)
}
