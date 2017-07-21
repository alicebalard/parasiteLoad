

#name{anova.glm.hybrid}
#alias{anova.glm.hybrid}
#title
#  Likelihood Ratio Tests for Negative Binomial GLMs, in the specific case of glm.hybrid
#
#\description{
#  Method function to perform sequential likelihood ratio tests for Negative
#  Binomial generalized linear models, in the specific case of glm.hybrid

#Arguments

#m1
#Fitted model object of class "glm.hybrid", the output of glm.hybrid().

#m2
#Fitted model object of class "glm.hybrid", the output of glm.hybrid().

anova.hybrid <- function(m1, m2){
  ## Test if the difference between 2 likelihood is significant
dLL = abs(m1$twologlik/2 - m2$twologlik/2)
dDF = length(m1$opt.param) - length(m2$opt.param)
p = 1 - pchisq(2*dLL, df = dDF) # G-test
}

## Get df from each model!


print.anova.hybrid <- function(){}
