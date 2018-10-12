#' Compare the hypotheses between each other : G-test. Test if the difference between 2 likelihood is significant
#'
#' @param model0 A model
#' @param model1 A model
#' @return A data.frame giving the difference of likelihood betweend the two models, the difference of degrees of freedom, and the pvalue of the G-test
#' @export

Gtest <- function(model0, model1){
  LL0 <- Reduce(x = c(model0), f = function(accum, model){
    accum + bbmle::logLik(model)}, init = 0)
  LL1 <- Reduce(x = c(model1), f = function(accum, model){
    accum + bbmle::logLik(model)}, init = 0)
  dLL <- LL1 - LL0
  N0 <- Reduce(x = c(model0), f = function(accum, model){
    accum + length(bbmle::coef(model))}, init = 0)
  N1 <- Reduce(x = c(model1), f = function(accum, model){
    accum + length(bbmle::coef(model))}, init = 0)
  dDF <- N1 - N0
  pvalue <- 1 - pchisq(2*dLL, df=dDF)
  out <- data.frame(dLL = round(dLL, 2),
                    dDF = dDF,
                    pvalue = round(pvalue, 6))
  print(out)
  return(out)
}
