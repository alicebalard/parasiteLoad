## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
collapse = TRUE,
comment = "#>"
)

## ----install-------------------------------------------------------------
library(parasiteLoad)
library(bbmle)
require(optimx) # for bbmle it needs to be required(?)
library(ggplot2)
library(MASS)
library(fitdistrplus) # evaluate distribution
library(epiR) # Sterne's exact method
library(simpleboot) # BS
library(boot) # BS

## ----prepData------------------------------------------------------------
# add 1 for worms/oocysts count!!
WATWMdata$`Aspiculuris.Syphacia+1` <- WATWMdata$Aspiculuris.Syphacia + 1
BALdata$`Aspiculuris.Syphacia+1` <- BALdata$Aspiculuris_Syphacia + 1
BALdata$`OPG+1` <- BALdata$OPG + 1

# and select variables for each analysis
pinworms_data <- BALdata[!is.na(BALdata$`Aspiculuris.Syphacia+1`) &
                           !is.na(BALdata$HI) &
                           !is.na(BALdata$Sex),]

WATWMdata <- WATWMdata[!is.na(WATWMdata$`Aspiculuris.Syphacia+1`) &
                         !is.na(WATWMdata$HI) &
                         !is.na(WATWMdata$Sex), ]

# all worms together for comparison
d1 <- WATWMdata[c("Sex", "Aspiculuris.Syphacia+1", "HI")]
d2 <- pinworms_data[c("Sex", "Aspiculuris.Syphacia+1", "HI")]
d1$batch <- "WATWM"
d2$batch <- "JENNY"
allWorms <- rbind(d1,d2)
allWorms$batch <- as.factor(allWorms$batch)

## ------------------------------------------------------------------------------------------------------
x <- BALdata$Aspiculuris_Syphacia[!is.na(BALdata$Aspiculuris_Syphacia)]

averageload <- round(mean(x, na.rm = T),1)
averageload 
Ni <- sum(x !=0, na.rm = T)
Ni 
# Confidence intervals calculated with Sterne's exact method
sternetest <- epiR::epi.prev(pos = length(x[x > 0]), tested = length(x),
                             se = 1, sp=1, conf.level = .95, method = "sterne")
cilow <- sternetest$ap["lower"]
cihigh <- sternetest$ap["upper"]
prevalence <- sternetest$ap["est"]

Ni
prevalence 
cilow
cihigh

Meanload <-  mean(x, na.rm = T)
max <- max(x, na.rm = T)
Meanload 
max

set.seed(8345) # to make it replicable
b1 <- one.boot(x, mean, R=10000)
b1$t0
boot.ci(b1, type=c("bca")) 

# Prevalence Ni (%)
# [CI 95%]
# 321 (52.9%)
# [48.85–56.84]

# Mean load (Max)
# [CI 95%]
# 19.67 (489)
# [16.22–24.56]

## MLE of parameter k negbin
hist(x, prob=T, breaks = max(x))
# fit the negative binomial distribution
fit <- fitdist(x, "nbinom")
# get the fitted densities. mu and size from fit.
fitD <- dnbinom(0:max(x), size=fit$estimate[1], mu=fit$estimate[2])
# add fitted line (blue) to histogram
lines(fitD, lwd="3", col="blue")
# Goodness of fit with the chi squared test  
# get the frequency table
t <- table(x)   
# convert to dataframe
df <- as.data.frame(t)
df$x <- as.numeric(as.character(df$x))
df$expp <- pnbinom(q = df$x, size=fit$estimate[1], mu = fit$estimate[2])
df$obsp <- cumsum(df$Freq / sum(df$Freq))
# perform the chi-squared test
chisq.test(x = df$expp, y = df$obsp)

# Pearson's Chi-squared test
# 
# data:  df$expp and df$obsp
# X-squared = 9120, df = 9025, p-value = 0.239

# Observed and expected frequencies do not differ significantly (P =
# 0.05) so there is no statistical evidence to reject H 0

size=fit$estimate[1]
size
# The negative binomial distribution seems to describe parasite load well for all parasites

## ----runfitJo, fig.width=7, fig.height=4---------------------------------------------------------------
fit_Joelle_negbin <- analyse(WATWMdata, "Aspiculuris.Syphacia+1",
                             model = "negbin", group = "Sex")
fit_Joelle_negbin

plot_Joelle_negbin <- bananaPlots(mod = fit_Joelle_negbin$H3,
                                  data = WATWMdata,
                                  response = "Aspiculuris.Syphacia+1",
                                  islog10 = TRUE, group = "Sex")
plot_Joelle_negbin

## ----runfitJe,  fig.width=7, fig.height=4--------------------------------------------------------------
fit_pinworms_negbin <- analyse(pinworms_data, response = "Aspiculuris.Syphacia+1",
                               model = "negbin", group = "Sex")
fit_pinworms_negbin

plot_pinworms_negbin <- bananaPlots(mod = fit_pinworms_negbin$H3,
                                    data = pinworms_data,
                                    response = "Aspiculuris.Syphacia+1",
                                    islog10 = TRUE, group = "Sex")
plot_pinworms_negbin

## ------------------------------------------------------------------------------------------------------
paramAlphaFixed <- getParamBounds(model = "negbin", data = pinworms_data, response= "Aspiculuris.Syphacia+1")

paramAlphaFixed[c("alphaStart", "alphaLB", "alphaUB")] <- c(1.4, 1.34, 1.44)

fit_pinworms_negbinB <- analyse(pinworms_data, response = "Aspiculuris.Syphacia+1",
                               model = "negbin", group = "Sex", 
                               myparamBounds = paramAlphaFixed)
fit_pinworms_negbinB$H3

Gtest(fit_pinworms_negbinB$H3, fit_pinworms_negbin$H3)

