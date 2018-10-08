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

## ---- choosefitsqpcr-----------------------------------------------------

qpcr_data[qpcr_data$delta_ct_MminusEAllPos > 0,"delta_ct_MminusEAllPos"]

fits_test_qpcr <- list(
  normal = MASS::fitdistr(qpcr_data[qpcr_data$delta_ct_MminusEAllPos > 0,"delta_ct_MminusEAllPos"],
                          "normal"),
  student = MASS::fitdistr(qpcr_data[qpcr_data$delta_ct_MminusEAllPos > 0, "delta_ct_MminusEAllPos"], "t",
                     start = list(m = mean(qpcr_data[qpcr_data$delta_ct_MminusEAllPos > 0, "delta_ct_MminusEAllPos"]),
                                  s = sd(qpcr_data[qpcr_data$delta_ct_MminusEAllPos > 0, "delta_ct_MminusEAllPos"]),
                                  df=3), lower=c(-1, 0.001,1))
)

#Get the logliks for each model...

sapply(fits_test_qpcr, function(i) i$loglik)

## ----runfitqpcr----------------------------------------------------------
fit_qpcr_student <- analyse(qpcr_data[qpcr_data$delta_ct_MminusEAllPos > 0, ],
                               response = "delta_ct_MminusEAllPos", 
                            model = "student", group = "Sex")
fit_qpcr_student

## ----plotfitqpcr, fig.width=6, fig.height=4------------------------------
plot_qpcr_student <- bananaPlots(mod = fit_qpcr_student$H0, 
                                 data = qpcr_data[qpcr_data$delta_ct_MminusEAllPos > 0, ], 
                                 response = "delta_ct_MminusEAllPos", group = "Sex")
plot_qpcr_student

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

## ----runfitJe------------------------------------------------------------
fit_pinworms_negbin <- analyse(pinworms_data, "Aspiculuris.Syphacia+1", 
                               model = "negbin", group = "Sex")
fit_pinworms_negbin

## ----plotfitJe, fig.width=7, fig.height=4--------------------------------
plot_pinworms_negbin <- bananaPlots(mod = fit_pinworms_negbin$H1, 
                                  data = pinworms_data, 
                                  response = "Aspiculuris.Syphacia+1", 
                                  islog10 = TRUE, group = "Sex") 
plot_pinworms_negbin

## ----plotcorr, fig.width=7, fig.height=4---------------------------------
ggplot2::ggplot(pinworms_data, 
       aes(delta_ct_MminusE, `Aspiculuris.Syphacia+1`)) +
  geom_point(aes(fill = HI), pch = 21, size = 5) + 
  scale_fill_gradient(low = "blue", high = "red") +
  scale_y_log10() +
  geom_smooth(method = "lm") +
  theme_classic()

summary(lm(`Aspiculuris.Syphacia+1` ~ delta_ct_MminusE, data = BALdata))

## ----removepregnant, fig.width=7, fig.height=4---------------------------
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

## ----fitresmod-----------------------------------------------------------
fitRes <- lm(Body_weight ~ Body_length * Sex, data = body_data)

## ----getpredval----------------------------------------------------------
body_data$predicted <- predict(fitRes)   # Save the predicted values
body_data$residuals <- residuals(fitRes) # Save the residual values

## ----plotpredval, fig.width=7, fig.height=4------------------------------
ggplot2::ggplot(body_data, ggplot2::aes(x = Body_length, y = Body_weight)) +
  ggplot2::geom_smooth(method = "lm", se = FALSE, color = "lightgrey") +  # Plot regression slope
  ggplot2::geom_segment(ggplot2::aes(xend = Body_length, yend = predicted), alpha = .2) +  # alpha to fade lines
  ggplot2::geom_point(ggplot2::aes(col = delta_ct_MminusE), size = 3) +
  ggplot2::scale_color_gradient(low = "lightgrey", high = "red") +
  ggplot2::geom_point(ggplot2::aes(y = predicted), shape = 1) +
  ggplot2::facet_grid(~ Sex, scales = "free_x") +  # Split panels here by `iv`
  ggplot2::theme_bw()  # Add theme for cleaner look

## ----removeoutliers, fig.width=7, fig.height=4---------------------------
hist(body_data$residuals[body_data$Sex =="F"], breaks = 100) # remove outliers, keep [-5,5] interval

## ---- fig.width=7, fig.height=4------------------------------------------
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

## ----choosefitbody-------------------------------------------------------
fits <- list(
  normal = MASS::fitdistr(body_data$resBMBL,"normal"),
  student = MASS::fitdistr(body_data$resBMBL, "t", 
                     start = list(m = mean(body_data$resBMBL), 
                                  s = sd(body_data$resBMBL), df = 3), 
                     lower=c(-1, 0.001,1))
)
# get the logliks for each model...
sapply(fits, function(i) i$loglik)

## ----plotdistribbody, fig.width=7, fig.height=4--------------------------
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

## ----runfitbody----------------------------------------------------------
fit_body_student <- analyse(body_data, response = "resBMBL",
                            model = "student", group = "qPCRstatus")

fit_body_student

## ----plotfitbody, fig.width=6, fig.height=4------------------------------
# bananaPlots(mod = fit_body_student$H1, data = body_data, response = "resBMBL")
# 
# bananaPlots(mod = fit_body_student$H3, data = body_data, response = "resBMBL", group = "qPCRstatus")

## ----compare all worms, fig.width=7, fig.height=4------------------------

fit_allWorms_negbin <- analyse(allWorms, "Aspiculuris.Syphacia+1", 
                               model = "negbin", group = "batch")
fit_allWorms_negbin

bananaPlots(data = allWorms, response = "Aspiculuris.Syphacia+1",
            mod = fit_allWorms_negbin$H3,
            group = "batch", islog10 = T)


bananaPlots

