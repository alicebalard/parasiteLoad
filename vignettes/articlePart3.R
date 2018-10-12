---
  title: "Data analysis : test of hybrid immune vigor in response to parasite infection"
subtitle: "Article part 3"
author: "Alice Balard"
date: "`r Sys.Date()`"
output:
  rmarkdown::html_vignette:
  vignette: >
  %\VignetteIndexEntry{Vignette Title}
%\VignetteEngine{knitr::rmarkdown}
%\VignetteEncoding{UTF-8}
---

  ```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

## Installation

```{r install}
library(parasiteLoad)
library(bbmle)
require(optimx) # for bbmle it needs to be required(?)
library(ggplot2)
library(MASS)
library(fitdistrplus) # evaluate distribution
library(epiR) # Sterne's exact method
library(simpleboot) # BS
library(boot) # BS
```

## Prepare dataset for each analysis

```{r prepData}
# add 1 for worms/oocysts count!!
WATWMdata$`Aspiculuris.Syphacia+1` <- WATWMdata$Aspiculuris.Syphacia + 1
BALdata$`Aspiculuris.Syphacia+1` <- BALdata$Aspiculuris_Syphacia + 1
BALdata$`OPG+1` <- BALdata$OPG + 1

# and select variables for each analysis
qpcr_data <- BALdata[!is.na(BALdata$HI) &
                       !is.na(BALdata$Sex) &
                       !is.na(BALdata$delta_ct_MminusEAllPos), ]

tocheck <- qpcr_data[qpcr_data$Mouse_ID %in%
                       qpcr_data[qpcr_data$delta_ct_ilwe_MminusE > -5 &
                                   qpcr_data$delta_ct_cewe_MminusE > -5,"Mouse_ID"],
                     c("Mouse_ID", "delta_ct_ilwe_MminusE", "delta_ct_cewe_MminusE", "qPCRsummary")]

speciesEimeria <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/raw_data/Eimeria_detection/Eimeria_species_assignment_14_17.csv")
speciesEimeria <- speciesEimeria[!names(speciesEimeria) %in% c("Code", "Transect")]

## NB. double infection for AA_0080 + cecum only for AA_0100
speciesEimeria$Mouse_ID[grep("i", speciesEimeria$Mouse_ID)]
speciesEimeria$Mouse_ID[grep("c", speciesEimeria$Mouse_ID)]

speciesEimeria[grep("AA_0080", speciesEimeria$Mouse_ID),]
  speciesEimeria[grep("Double", speciesEimeria$Species),]


speciesEimeria <- rbind(data.frame(Year = c(2016, 2016),
                                   Mouse_ID = c("AA_0080", "AA_0100"),
                                   n18S_Seq = "positive",
                                   COI_Seq = "positive",
                                   ORF470_Seq = "positive",
                                   Species = c("E_vermiformis_E_ferrisi", "Double")), speciesEimeria)

speciesEimeria$Mouse_ID <- gsub(" ", "", as.character(speciesEimeria$Mouse_ID))

qpcr_data_species <- merge(qpcr_data, speciesEimeria, by = "Mouse_ID")

qpcr_data_species$Mouse_ID[duplicated(qpcr_data_species$Mouse_ID)]

flotation_data <- BALdata[!is.na(BALdata$`OPG+1`) &
                            !is.na(BALdata$HI) &
                            !is.na(BALdata$Sex), ]

```

## Resistance in hybrids compared to pure strains mice : Eimeria

First, choose a correct distribution for our data

```{r}
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

averageload <- round(mean(qpcrPOS, na.rm = T),1)
averageload
Ni <- sum(qpcrPOS, na.rm = T)
Ni
# Confidence intervals calculated with Sterne's exact method
sternetest <- epiR::epi.prev(pos = length(qpcrPOS), tested = length(qpcrPOS),
                             se = 1, sp=1, conf.level = .95, method = "sterne")
cilow <- sternetest$ap["lower"]
cihigh <- sternetest$ap["upper"]
prevalence <- sternetest$ap["est"]

Ni
prevalence
cilow
cihigh

Meanload <-  mean(qpcrPOS, na.rm = T)
max <- max(qpcrPOS, na.rm = T)
Meanload
max

set.seed(8345) # to make it replicable
b1 <- one.boot(qpcrPOS, mean, R=10000)
b1$t0
boot.ci(b1, type=c("bca"))

## MLE of parameter k negbin
hist(qpcrPOS, prob=T, breaks = max(qpcrPOS))
# fit the negative binomial distribution
fit <- fitdist(qpcrPOS, "weibull")
# get the fitted densities. mu and size from fit.
fitD <- dnbinom(0:max(xPOS), size=fit$estimate[1], mu=fit$estimate[2])
# add fitted line (blue) to histogram
lines(fitD, lwd="3", col="blue")
# Goodness of fit with the chi squared test
# get the frequency table
t <- table(xPOS)
# convert to dataframe
df <- as.data.frame(t)
df$xPOS <- as.numeric(as.character(df$xPOS))
df$expp <- pnbinom(q = df$xPOS, size=fit$estimate[1], mu = fit$estimate[2])
df$obsp <- cumsum(df$Freq / sum(df$Freq))
# perform the chi-squared test
chisq.test(x = df$expp, y = df$obsp)

size=fit$estimate[1]
size
# The negative binomial distribution seems to describe parasite load well for all parasites
```


```{r,  fig.width=7, fig.height=4}
pinworms_dataPOS <- pinworms_data[pinworms_data$Aspiculuris_Syphacia >= 1, ]

fit_pinworms_negbinPOS <- analyse(pinworms_dataPOS, response = "Aspiculuris_Syphacia",
                                  model = "negbin", group = "Sex")
fit_pinworms_negbinPOS

plot_pinworms_negbinPOS <- bananaPlots(mod = fit_pinworms_negbinPOS$H3,
                                       data = pinworms_dataPOS,
                                       response = "Aspiculuris_Syphacia",
                                       islog10 = TRUE, group = "Sex")
plot_pinworms_negbinPOS
```
