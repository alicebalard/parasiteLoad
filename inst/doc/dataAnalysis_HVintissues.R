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

## ----dataprep------------------------------------------------------------
# Works if OPG are integers
BALdata$OPG <- round(BALdata$OPG)

# To pass positive I add 5 to all
BALdata$delta_ct_MminusE <- BALdata$delta_ct_MminusE + 5
# changes to do for avoiding plotting error
levels(BALdata$Sex) <- c(levels(BALdata$Sex), "female", "male")
BALdata$Sex[BALdata$Sex == "F"] <- "female"
BALdata$Sex[BALdata$Sex == "M"] <- "male"
BALdata$Sex <- droplevels(BALdata$Sex)

## First question: can we consider the zeros out? If our hypothesis is : all get infected the same way, 
# but hybrids got lower loads?
BALdata$qPCRstatus <- 0
BALdata$qPCRstatus[BALdata$delta_ct_MminusE > 0] <- 1

BALdata$flotStatus <- 0
BALdata$flotStatus[BALdata$OPG > 0] <- 1

BALdata$presenceAspiSypha[BALdata$Aspiculuris_Syphacia > 0] <- 1
BALdata$presenceAspiSypha[BALdata$Aspiculuris_Syphacia == 0] <- 0

## Prepare dataset for each analysis
qpcr_data <- BALdata[!is.na(BALdata$HI) &
                       !is.na(BALdata$Sex) &
                       !is.na(BALdata$delta_ct_MminusE), ]

flotation_data <- BALdata[!is.na(BALdata$OPG) &
                                        !is.na(BALdata$HI) &
                                        !is.na(BALdata$Sex), ]

pinworms_data <- BALdata[!is.na(BALdata$Aspiculuris_Syphacia) &
                                             !is.na(BALdata$HI) &
                                             !is.na(BALdata$Sex),]

body_data <- BALdata[!is.na(BALdata$Body_weight) &
                                   !is.na(BALdata$Body_length) &
                                   !is.na(BALdata$HI) &
                                   !is.na(BALdata$Sex) &
                                   !is.na(BALdata$delta_ct_MminusE) ,]
body_data$qPCRstatus <- as.factor(body_data$qPCRstatus)

## Prepare data WATWM
WATWMdata <- WATWMdata[complete.cases(WATWMdata$HI),]
levels(WATWMdata$Sex) <- c(levels(WATWMdata$Sex), "female", "male")
WATWMdata$Sex[WATWMdata$Sex == "F"] <- "female"
WATWMdata$Sex[WATWMdata$Sex == "M"] <- "male"
WATWMdata$Sex <- droplevels(WATWMdata$Sex)

WATWMdata <- WATWMdata[!is.na(WATWMdata$Aspiculuris.Syphacia),]
WATWMdata$presenceAspiSypha[WATWMdata$Aspiculuris.Syphacia > 0] <- 1
WATWMdata$presenceAspiSypha[WATWMdata$Aspiculuris.Syphacia == 0] <- 0

## ------------------------------------------------------------------------
library(MASS)
fits_test_qpcr <- list(
  normal = fitdistr(qpcr_data[qpcr_data$delta_ct_MminusE > 0, "delta_ct_MminusE"],"normal"),
  student = fitdistr(qpcr_data[qpcr_data$delta_ct_MminusE > 0, "delta_ct_MminusE"], "t", 
                     start = list(m = mean(qpcr_data[qpcr_data$delta_ct_MminusE > 0, "delta_ct_MminusE"]), 
                                  s = sd(qpcr_data[qpcr_data$delta_ct_MminusE > 0, "delta_ct_MminusE"]), 
                                  df=3), lower=c(-1, 0.001,1))
)

## ------------------------------------------------------------------------
sapply(fits_test_qpcr, function(i) i$loglik)

## ------------------------------------------------------------------------
fit_qpcr_student <- analyse(qpcr_data[qpcr_data$delta_ct_MminusE > 0, ],
                               response = "delta_ct_MminusE", 
                            model = "student", group = "Sex")

## ------------------------------------------------------------------------
plot_qpcr_student <- bananaPlots(mod = fit_qpcr_student$H0, 
                                 data = qpcr_data[qpcr_data$delta_ct_MminusE > 0, ], 
                                 response = "delta_ct_MminusE", group = "Sex")
plot_qpcr_student

## ------------------------------------------------------------------------
fit2016_qpcr_student <- analyse(qpcr_data[qpcr_data$delta_ct_MminusE > 0 & qpcr_data$Year %in% 2016, ],
                                   response = "delta_ct_MminusE", 
                                   model = "student", group = "Sex")
bananaPlots(mod = fit2016_qpcr_student$H1, 
            data = qpcr_data[qpcr_data$delta_ct_MminusE > 0 & qpcr_data$Year %in% 2016, ], 
            response = "delta_ct_MminusE", group = "Sex") # significant

fit2017_qpcr_student <- analyse(qpcr_data[qpcr_data$delta_ct_MminusE > 0 & qpcr_data$Year %in% 2017, ],
                                   response = "delta_ct_MminusE",
                                   model = "student", group = "Sex")
bananaPlots(mod = fit2017_qpcr_student$H1, 
            data = qpcr_data[qpcr_data$delta_ct_MminusE > 0 & qpcr_data$Year %in% 2017, ], 
            response = "delta_ct_MminusE", group = "Sex") # not significant

