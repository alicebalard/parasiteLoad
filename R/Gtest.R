## Compare the hypotheses between each other : G-test
## Test if the difference between 2 likelihood is significant
Gtest <- function(model0, model1){
  LL0 <- Reduce(x = c(model0), f = function(accum, model){
    accum + logLik(model)}, init = 0)
  LL1 <- Reduce(x = c(model1), f = function(accum, model){
    accum + logLik(model)}, init = 0)
  dLL <- LL1 - LL0
  N0 <- Reduce(x = c(model0), f = function(accum, model){
    accum + length(coef(model))}, init = 0)
  N1 <- Reduce(x = c(model1), f = function(accum, model){
    accum + length(coef(model))}, init = 0)
  dDF <- N1 - N0
  pvalue <- 1 - pchisq(2*dLL, df=dDF)
  out <- data.frame(dLL = round(dLL, 2),
                    dDF = dDF,
                    pvalue = 
                      if (pvalue < 0.01){
                      "< 0.01"
                    } else {
                      round(pvalue, 3)
                    }
  )
  print(out)
  return(out)
}
