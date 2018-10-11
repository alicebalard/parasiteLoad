---
  title: "Data analysis : test of hybrid immune vigor in response to parasite infection"
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
library(fitdistrplus)
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
```

## Worms count, 4 hypotheses, difference between sexes

### Choose a correct distribution for our data : negative binomial see WATWM

### pinworms (A. tetraptera and S. obvelata (WATWMdata$Aspiculuris.Syphacia))

```{r runfitJo, fig.width=7, fig.height=4}
fit_Joelle_negbin <- analyse(WATWMdata, "Aspiculuris.Syphacia+1",
                             model = "negbin", group = "Sex")
fit_Joelle_negbin

plot_Joelle_negbin <- bananaPlots(mod = fit_Joelle_negbin$H1,
                                  data = WATWMdata,
                                  response = "Aspiculuris.Syphacia+1",
                                  islog10 = TRUE, group = "Sex")
plot_Joelle_negbin
```

### And Jennys now:

```{r runfitJe,  fig.width=7, fig.height=4}
fit_pinworms_negbin <- analyse(pinworms_data, response = "Aspiculuris.Syphacia+1",
                               model = "negbin", group = "Sex")
fit_pinworms_negbin

plot_pinworms_negbin <- bananaPlots(mod = fit_pinworms_negbin$H3,
                                    data = pinworms_data,
                                    response = "Aspiculuris.Syphacia+1",
                                    islog10 = TRUE, group = "Sex")
plot_pinworms_negbin
```



## Part 1. qpcr results, 4 hypotheses, difference between sexes

### Choose a correct distribution for our data

Compute some fits...

```{r choosefitsqpcr, fig.width=6, fig.height=4}
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
```

WEIBULL IS FAR BETTER THAN STUDENT!!

### Run the fit with Weibull distribution

```{r fitandplotweibull, fig.width=6, fig.height=4}
fit_qpcr_weibull <- analyse(qpcr_data[qpcr_data$delta_ct_MminusEAllPos > 0, ],
response = "delta_ct_MminusEAllPos",
model = "weibull", group = "Sex")

plot_qpcr_weibull <- bananaPlots(mod = fit_qpcr_weibull$H0,
data = qpcr_data[qpcr_data$delta_ct_MminusEAllPos > 0, ],
response = "delta_ct_MminusEAllPos", group = "Sex")
plot_qpcr_weibull
```

## Part 2. OPG results, 4 hypotheses, difference between sexes

## Run the fit for load estimation along HI, negative binomial distribution

```{r runfitflot}
fit_flotation_negbin <- analyse(flotation_data[flotation_data$OPG > 0, ],
response = "OPG+1",
model = "negbin", group = "Sex")
fit_flotation_negbin
```

```{r plotfitflot, fig.width=6, fig.height=4}
# bananaPlots(mod = fit_flotation_negbin$H0, data = flotation_data,
#             response = "OPG+1")
# problems for the profiling...
```



### And Jennys now for only positive (RESISTANCE):

```{r runfitJeRES, fig.width=7, fig.height=4}
pinworms_dataPOS <- pinworms_data[pinworms_data$Aspiculuris.Syphacia >= 1,]

fit_pinworms_negbinRES <- analyse(pinworms_dataPOS, response = "Aspiculuris.Syphacia",
model = "negbin", group = "Sex")
fit_pinworms_negbinRES

plot_pinworms_negbinRES <- bananaPlots(mod = fit_pinworms_negbinRES$H1,
data = pinworms_dataPOS,
response = "Aspiculuris.Syphacia",
islog10 = TRUE, group = "Sex")
plot_pinworms_negbinRES
```

## qPCR vs worms (is there a correlation?)

```{r plotcorr, fig.width=7, fig.height=4}
ggplot2::ggplot(pinworms_data,
aes(delta_ct_MminusE, `Aspiculuris.Syphacia+1`)) +
geom_point(aes(fill = HI), pch = 21, size = 5) +
scale_fill_gradient(low = "blue", high = "red") +
scale_y_log10() +
geom_smooth(method = "lm") +
theme_classic()

summary(lm(`Aspiculuris.Syphacia+1` ~ delta_ct_MminusE, data = BALdata))
```

## Part 4. body data, 4 hypotheses, difference between infected/not infected
### TO DO : test coinfection data!! Are coinfected individuals more suffering?

Remove pregnant/post partum and juveniles

```{r removepregnant, fig.width=7, fig.height=4}
body_data$status[body_data$Sex %in% c("F")] <- "non pregnant/lactating female"
body_data$status[body_data$Status %in% c("post partum", "post partum (lactating)", "pregnant")] <- "pregnant/lactating female"
body_data$status[body_data$Sex %in% c("M")] <- "male"

# Test  our detection of pregnancy in females using BCI
ggplot2::ggplot(body_data,
ggplot2::aes(x = HI, y = BCI, fill = status, group = status)) +
geom_point(pch = 21, size = 3, alpha = .5)+
geom_smooth(aes(col = status)) +
theme_bw()

# Remove pregnant females
body_data <- body_data[!body_data$status %in% c("pregnant/lactating female"), ]
```

Regression of BM/BS for males and females (all together, then separate subsp.)
Advantage: independant of size!!

* Step 1: fit the model

```{r fitresmod}
fitRes <- lm(Body_weight ~ Body_length * Sex, data = body_data)
```

* Step 2: obtain predicted and residual values

```{r getpredval}
body_data$predicted <- predict(fitRes)   # Save the predicted values
body_data$residuals <- residuals(fitRes) # Save the residual values
```

* Step 3: plot the actual and predicted values

```{r plotpredval, fig.width=7, fig.height=4}
ggplot2::ggplot(body_data, ggplot2::aes(x = Body_length, y = Body_weight)) +
ggplot2::geom_smooth(method = "lm", se = FALSE, color = "lightgrey") +  # Plot regression slope
ggplot2::geom_segment(ggplot2::aes(xend = Body_length, yend = predicted), alpha = .2) +  # alpha to fade lines
ggplot2::geom_point(ggplot2::aes(col = delta_ct_MminusE), size = 3) +
ggplot2::scale_color_gradient(low = "lightgrey", high = "red") +
ggplot2::geom_point(ggplot2::aes(y = predicted), shape = 1) +
ggplot2::facet_grid(~ Sex, scales = "free_x") +  # Split panels here by `iv`
ggplot2::theme_bw()  # Add theme for cleaner look
```

* Step 4: use residuals as indice

```{r removeoutliers, fig.width=7, fig.height=4}
hist(body_data$residuals[body_data$Sex =="F"], breaks = 100) # remove outliers, keep [-5,5] interval
```

```{r, fig.width=7, fig.height=4}
body_data <- body_data[body_data$residuals <= 5,]
body_data$resBMBL <- body_data$residuals

# give positive values only
body_data$resBMBL <- body_data$resBMBL + 5

# Plot the actual and predicted values
ggplot2::ggplot(body_data, ggplot2::aes(x = Body_length, y = Body_weight)) +
ggplot2::geom_smooth(method = "lm", se = FALSE, color = "lightgrey") +  # Plot regression slope
ggplot2::geom_segment(ggplot2::aes(xend = Body_length, yend = predicted), alpha = .2) +  # alpha to fade lines
ggplot2::geom_point(ggplot2::aes(col = delta_ct_MminusE), size = 3) +
ggplot2::scale_color_gradient(low = "lightgrey", high = "red") +
ggplot2::geom_point(ggplot2::aes(y = predicted), shape = 1) +
ggplot2::facet_grid(~ Sex, scales = "free_x") +  # Split panels here by `iv`
ggplot2::theme_bw()  # Add theme for cleaner look
```

### Which distribution to choose?

Let's compute some fits...

```{r choosefitbody}
fits <- list(
  normal = MASS::fitdistr(body_data$resBMBL,"normal"),
  student = MASS::fitdistr(body_data$resBMBL, "t",
                           start = list(m = mean(body_data$resBMBL),
                                        s = sd(body_data$resBMBL), df = 3),
                           lower=c(-1, 0.001,1))
)
# get the logliks for each model...
sapply(fits, function(i) i$loglik)
```

STUDENT is the way to go!

  ```{r plotdistribbody, fig.width=7, fig.height=4}
ggplot2::ggplot(body_data, ggplot2::aes(resBMBL)) +
  ggplot2::geom_histogram(ggplot2::aes(y=..density..), bins = 100) +
  ggplot2::stat_function(fun = dnorm, n = 1e3, args = list(mean = fits$normal$estimate[1], sd = fits$normal$estimate[2]),
                         ggplot2::aes(color = "normal"), size = 2) +
  ggplot2::stat_function(fun = dt, n = 1e3, args = list(ncp = fits$student$estimate[1],
                                                        df = fits$student$estimate[3]),
                         ggplot2::aes(color = "student"), size = 2) +
  ggplot2::stat_function(fun = dt, n = 1e3, args = list(ncp = fits$student$estimate[1],
                                                        df = 10),
                         ggplot2::aes(color = "student10"), size = 2) +
  ggplot2::stat_function(fun = dt, n = 1e3, args = list(ncp = fits$student$estimate[1],
                                                        df = 400),
                         ggplot2::aes(color = "student500"), size = 2) +
  ggplot2::stat_function(fun = dt, n = 1e3, args = list(ncp = fits$student$estimate[1],
                                                        df = 1000),
                         ggplot2::aes(color = "student1000"), size = 2) +
  ggplot2::stat_function(fun = dt, n = 1e3, args = list(ncp = fits$student$estimate[1],
                                                        df = 5),
                         ggplot2::aes(color = "student5"), size = 2) +
  ggplot2::theme_bw(base_size = 24)
```

### Fit the model!

```{r runfitbody}
fit_body_student <- analyse(body_data, response = "resBMBL",
                            model = "student", group = "qPCRstatus")

fit_body_student
```

And plot

```{r plotfitbody, fig.width=6, fig.height=4}
# bananaPlots(mod = fit_body_student$H1, data = body_data, response = "resBMBL")
#
# bananaPlots(mod = fit_body_student$H3, data = body_data, response = "resBMBL", group = "qPCRstatus")
```


# Suggestion from Stuart

Sounds very cool. You may well have thought of this, but anyway (because likelihood is lovely)... a nice zonewise pinworm comparison is as follows:

  Do whole analysis on combined {WATWM, Jenny's} data -> Joint LLsurface
You already have the whole analyses on the zonewise-split of this joint data. -> {WATM LL surface,Jenny LL surface}
... and off we go: compare
Max(WATWM LL surface)+Max(Jenny LL surface) VS Max(Joint LLsurface)

On either side of the VS the same data is analysed - so LL comparison is valid. Any 'is there really a difference between the two zones re pinworms' style of question reduces to:
Force them to be the same (the joint analysis) VS allow them to be different (WATWM and Jenny's each have their own parameters).... once you have counted degrees of freedom you can look up the answer to any such  question in a ChiSquared table... its so sweeeet!
  (ok... there are lots of axes of heterogeneity... it not quite so simple as all that - but its basically what you do already re eg comparison of sexes)

```{r compare all worms, fig.width=7, fig.height=4}

fit_allWorms_negbin <- analyse(allWorms, "Aspiculuris.Syphacia+1",
                               model = "negbin", group = "batch")
fit_allWorms_negbin

bananaPlots(data = allWorms, response = "Aspiculuris.Syphacia+1",
            mod = fit_allWorms_negbin$H3,
            group = "batch", islog10 = T)


bananaPlots
```
