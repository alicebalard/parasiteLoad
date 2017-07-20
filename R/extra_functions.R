

#name{anova.glm.hybrid}
#alias{anova.glm.hybrid}
#title
#  Likelihood Ratio Tests for Negative Binomial GLMs, in the specific case of glm.hybrid
#
#\description{
#  Method function to perform sequential likelihood ratio tests for Negative
#  Binomial generalized linear models, in the specific case of glm.hybrid

#Arguments

#object	
#Fitted model object of class "glm.hybrid", the output of glm.hybrid().

# formula
# give the formula used

anova.hybrid <- function(object, formula){
  ## Test if the difference between 2 likelihood is significant
  Gtest <- function(dLL, dDF){
    1 - pchisq(2*dLL, df = dDF) 
  }
  
  formula <- loads ~ group1 * HI * group2
  object <- glm.hybrid::glm.hybrid(formula = loads ~ group1 * HI * group2, data = simdata, "HI",
                                   alpha.start = 1)
  
  glm.hybrid::glm.hybrid(formula = loads ~ HI * group1, data = simdata, "HI",
                         alpha.start = 1)
  
  #  nullmodel <-
  
  
  2*dLL <- object$twologlik - nullmodel$twologlik
  dDF <- length(object$opt.param) - length(submodel$opt.param)
  p <- round(Gtest(dLL, dDF),4)
  ##
  out <- data.frame(Model = mds, theta = ths, Resid.df = dfs,
                    "2 x log-lik." = lls, Test = tss, df = df, LRtest = x2,
                    Prob = pr)
  names(out) <- c("Model", "theta", "Resid. df",
                  "   2 x log-lik.", "Test", "   df", "LR stat.", "Pr(Chi)")
  class(out) <- c("Anova", "data.frame")
  attr(out, "heading") <-
    c("Likelihood ratio tests of Negative Binomial Models\n",
      paste("Response:", rsp))
  out
}


print.anova.hybrid <- function(){}
