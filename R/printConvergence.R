#' Print if model converged or not
#'
#' @param fit A fit to be tested for convergence
#' @return A character string, "Did converge" or "Did not converge"
#' @export

printConvergence <- function(fit) {
  convergence <- fit@details[["convergence"]]
  print(ifelse(convergence == 0, "Did converge", "Did not converge"))
}
