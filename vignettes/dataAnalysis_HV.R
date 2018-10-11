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
library(fitdistrplus)

## ----prepData------------------------------------------------------------
# add 1 for worms/oocysts count!!
WATWMdata$`Aspiculuris.Syphacia+1` <- WATWMdata$Aspiculuris.Syphacia + 1
BALdata$`Aspiculuris.Syphacia+1` <- BALdata$Aspiculuris_Syphacia + 1
BALdata$`OPG+1` <- BALdata$OPG + 1

# and select variables for each analysis
qpcr_data <- BALdata[!is.na(BALdata$HI) &
                       !is.na(BALdata$Sex) &
                       !is.na(BALdata$delta_ct_MminusEAllPos), ]

flotation_data <- BALdata[!is.na(BALdata$`OPG+1`) &
                            !is.na(BALdata$HI) &
                            !is.na(BALdata$Sex), ]


pinworms_data <- BALdata[!is.na(BALdata$`Aspiculuris.Syphacia+1`) &
                           !is.na(BALdata$HI) &
                           !is.na(BALdata$Sex),]

body_data <- BALdata[!is.na(BALdata$Body_weight) &
                       !is.na(BALdata$Body_length) &
                       !is.na(BALdata$HI) &
                       !is.na(BALdata$Sex) &
                       !is.na(BALdata$qPCRstatus) ,]

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

## ----choosefitsqpcr, fig.width=6, fig.height=4---------------------------
qpcrPOS <- qpcr_data[qpcr_data$delta_ct_MminusEAllPos > 0, "delta_ct_MminusEAllPos"]

fitdistrplus::descdist(qpcrPOS, discrete = FALSE)

fits_test_qpcr <- list(
  normal = fitdistrplus::fitdist(qpcrPOS, "norm"),
  weibull = fitdistrplus::fitdist(qpcrPOS, "weibull"),
  student = fitdistrplus::fitdist(qpcrPOS, "t", start = list(df = 1))
)

plot(fits_test_qpcr$normal)
plot(fits_test_qpcr$weibull)
plot(fits_test_qpcr$student)

#Get the logliks for each model...
sapply(fits_test_qpcr, function(i) i$loglik)
sapply(fits_test_qpcr, function(i) i$aic)

## ----fitandplotweibull, fig.width=6, fig.height=4------------------------
fit_qpcr_weibull <- analyse(qpcr_data[qpcr_data$delta_ct_MminusEAllPos > 0, ],
        response = "delta_ct_MminusEAllPos", 
        model = "weibull", group = "Sex")

plot_qpcr_weibull <- bananaPlots(mod = fit_qpcr_weibull$H0, 
                                 data = qpcr_data[qpcr_data$delta_ct_MminusEAllPos > 0, ], 
                                 response = "delta_ct_MminusEAllPos", group = "Sex")
plot_qpcr_weibull

## ----runfitflot----------------------------------------------------------
fit_flotation_negbin <- analyse(flotation_data[flotation_data$OPG > 0, ], 
                                   response = "OPG+1", 
                                   model = "negbin", group = "Sex")
fit_flotation_negbin

## ----plotfitflot, fig.width=6, fig.height=4------------------------------
# bananaPlots(mod = fit_flotation_negbin$H0, data = flotation_data, 
#             response = "OPG+1")
# problems for the profiling...

## ----runfitJo------------------------------------------------------------
fit_Joelle_negbin <- analyse(WATWMdata, "Aspiculuris.Syphacia+1", 
                             model = "negbin", group = "Sex")
fit_Joelle_negbin

## ----plotfitJo, fig.width=7, fig.height=4--------------------------------
plot_Joelle_negbin <- bananaPlots(mod = fit_Joelle_negbin$H1, 
                                 data = WATWMdata, 
                                 response = "Aspiculuris.Syphacia+1", 
                                 islog10 = TRUE, group = "Sex") 
plot_Joelle_negbin

## ----runfitJe,  fig.width=7, fig.height=4--------------------------------
fit_pinworms_negbin <- analyse(pinworms_data, response = "Aspiculuris.Syphacia+1", 
                               model = "negbin", group = "Sex")
fit_pinworms_negbin

plot_pinworms_negbin <- bananaPlots(mod = fit_pinworms_negbin$H1, 
                                  data = pinworms_data, 
                                  response = "Aspiculuris.Syphacia+1", 
                                  islog10 = TRUE, group = "Sex") 
plot_pinworms_negbin

