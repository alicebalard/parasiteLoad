## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----install-------------------------------------------------------------
library(parasiteLoad)
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

## ----runqpcr-------------------------------------------------------------
fit_qpcr_binom <- analyse(data = qpcr_data, response = "qPCRstatus", model = "binomial", group = "Sex")

## ----plot_qpcr_binom-----------------------------------------------------
plot_qpcr_binom <- bananaPlots(mod = fit_qpcr_binom$H0, data = qpcr_data, response = "qPCRstatus", group = "Sex") +
  ggplot2::coord_cartesian(ylim = c(0, 0.5)) # zoom in
plot_qpcr_binom

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

## ------------------------------------------------------------------------
fit_opg_binom <- analyse(flotation_data, "flotStatus", model = "binomial", group = "Sex")

## ------------------------------------------------------------------------
bananaPlots(mod = fit_opg_binom$H0, data = flotation_data, response = "flotStatus", group = "Sex") 
 
hist(flotation_data[flotation_data$OPG > 0, "OPG"], 100)

## ------------------------------------------------------------------------
fit_flotation_negbin <- analyse(flotation_data[flotation_data$OPG > 0, ], 
                                   response = "OPG", 
                                   model = "negbin", group = "Sex")

## ------------------------------------------------------------------------
# bananaPlots(mod = fit_flotation_negbin$H0, data = flotation_data, response = "OPG", group = "Sex")
# problems for the profiling...

## ------------------------------------------------------------------------
fit_Joelle_binom <- analyse(WATWMdata, "presenceAspiSypha", model = "binomial",
                            group = "Sex")
fit_Joelle_binom

## ------------------------------------------------------------------------
plot_Joelle_binom <- bananaPlots(mod = fit_Joelle_binom$H0, 
                                 data = WATWMdata, 
                                 response = "presenceAspiSypha", 
                                 group = "Sex")
plot_Joelle_binom 

## ------------------------------------------------------------------------
fit_pinworms_data <- analyse(pinworms_data, "presenceAspiSypha", model = "binomial",
               group = "Sex")

## ------------------------------------------------------------------------
plot_pinworms_data <- bananaPlots(mod = fit_pinworms_data$H0, 
                                  data = pinworms_data, 
                                  response = "presenceAspiSypha") 
plot_pinworms_data# + coord_cartesian(ylim = c(0.45, 0.9)) # zoom in

## ------------------------------------------------------------------------
fit_Joelle_negbin <- analyse(WATWMdata, "Aspiculuris.Syphacia", 
                             model = "negbin", group = "Sex")
fit_Joelle_negbin

## ------------------------------------------------------------------------
plot_Joelle_negbin <- bananaPlots(mod = fit_Joelle_negbin$H1, 
                                 data = WATWMdata, 
                                 response = "Aspiculuris.Syphacia", 
                                 islog10 = TRUE, group = "Sex") 
plot_Joelle_negbin

## ------------------------------------------------------------------------
fit_pinworms_negbin <- analyse(pinworms_data, "Aspiculuris_Syphacia", 
                               model = "negbin", group = "Sex")
fit_pinworms_negbin

## ------------------------------------------------------------------------
plot_pinworms_negbin <- bananaPlots(mod = fit_pinworms_negbin$H1, 
                                  data = pinworms_data, 
                                  response = "Aspiculuris_Syphacia", 
                                  islog10 = TRUE, group = "Sex") 
plot_pinworms_negbin

## ------------------------------------------------------------------------
ggplot2::ggplot(pinworms_data, 
       aes(delta_ct_MminusE, Aspiculuris_Syphacia + 1)) +
  geom_point(aes(fill = HI), pch = 21, size = 5) + 
  scale_fill_gradient(low = "blue", high = "red") +
  scale_y_log10() +
  geom_smooth(method = "lm") +
  theme_classic()

summary(lm(Aspiculuris_Syphacia ~ delta_ct_MminusE, data = BALdata))

## ------------------------------------------------------------------------
body_data$status[body_data$Sex %in% c("female")] <- "non pregnant/lactating female"
body_data$status[body_data$Status %in% c("post partum", "post partum (lactating)", "pregnant")] <- "pregnant/lactating female"
body_data$status[body_data$Sex %in% c("male")] <- "male"

# Test  our detection of pregnancy in females using BCI
ggplot2::ggplot(body_data, 
       ggplot2::aes(x = HI, y = BCI, fill = status, group = status)) +
  geom_point(pch = 21, size = 3, alpha = .5)+
  geom_smooth(aes(col = status)) +
  theme_bw()

# Remove pregnant females
body_data <- body_data[!body_data$status %in% c("pregnant/lactating female"), ]

## ------------------------------------------------------------------------
fitRes <- lm(Body_weight ~ Body_length * Sex, data = body_data)

## ------------------------------------------------------------------------
body_data$predicted <- predict(fitRes)   # Save the predicted values
body_data$residuals <- residuals(fitRes) # Save the residual values

## ------------------------------------------------------------------------
ggplot2::ggplot(body_data, ggplot2::aes(x = Body_length, y = Body_weight)) +
  ggplot2::geom_smooth(method = "lm", se = FALSE, color = "lightgrey") +  # Plot regression slope
  ggplot2::geom_segment(ggplot2::aes(xend = Body_length, yend = predicted), alpha = .2) +  # alpha to fade lines
  ggplot2::geom_point(ggplot2::aes(col = delta_ct_MminusE), size = 3) +
  ggplot2::scale_color_gradient(low = "lightgrey", high = "red") +
  ggplot2::geom_point(ggplot2::aes(y = predicted), shape = 1) +
  ggplot2::facet_grid(~ Sex, scales = "free_x") +  # Split panels here by `iv`
  ggplot2::theme_bw()  # Add theme for cleaner look

## ------------------------------------------------------------------------
hist(body_data$residuals[body_data$Sex =="female"], breaks = 100) # remove outliers, keep [-5,5] interval

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

## ------------------------------------------------------------------------
fits <- list(
  normal = MASS::fitdistr(body_data$resBMBL,"normal"),
  student = MASS::fitdistr(body_data$resBMBL, "t", 
                     start = list(m = mean(body_data$resBMBL), 
                                  s = sd(body_data$resBMBL), df = 3), 
                     lower=c(-1, 0.001,1))
)
# get the logliks for each model...
sapply(fits, function(i) i$loglik)

## ------------------------------------------------------------------------
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

## ------------------------------------------------------------------------
fit_body_student <- analyse(body_data, response = "resBMBL",
                            model = "student", group = "qPCRstatus")

## ------------------------------------------------------------------------
bananaPlots(mod = fit_body_student$H1, data = body_data, response = "resBMBL")

bananaPlots(mod = fit_body_student$H3, data = body_data, response = "resBMBL", group = "qPCRstatus")

