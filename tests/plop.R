# install.packages("bbmle")
library(bbmle)
library(ggplot2)

MeanLoad <- function(L1, L2, alpha, hybridIndex){
  heterozygoty <- 2 * hybridIndex * (1 - hybridIndex)
  mean <- (L1 + (L2 - L1) * hybridIndex) * (1 - alpha * heterozygoty)
  mean <- sapply(mean, function(x) {
    return(max(x, 0.01))
  })
  return(mean)
}

Aggregation <- function(A1, A2, Z, HI){
  He <- 2 * HI * (1 - HI)
  (A1 + (A2 - A1) * HI) + Z * He 
} 

Joelle_data <- read.csv("../examples/Reproduction_WATWM/EvolutionFinalData.csv")
Joelle_data[is.na(Joelle_data)] <- 0

simparams <- c(k = 2,
               "male:old.alpha" = 1.4,
               "male:young.alpha" = 1.2,
               "male:baby.alpha" = 1.0,
               "female:old.alpha" = 2.0,
               "female:young.alpha" = -1.8,
               "female:baby.alpha" = 1.1,
               "male:old.inter" = 14,
               "male:young.inter" = 12,
               "male:baby.inter" = 10,
               "female:old.inter" = 20,
               "female:young.inter" = 18,
               "female:baby.inter" = 11,
               "male:old.growth" = 2,
               "male:young.growth" = 1,
               "male:baby.growth" = -4,
               "female:old.growth" = 2,
               "female:young.growth" = 0,
               "female:baby.growth" = -1)

SimulatedData <- function(param, n){
  gdata <- data.frame(
    group1 = rep(c("male", "female"), each=n/2),
    group2 = sample(c("old", "young", "baby"), n, replace=TRUE)
  )
  gdata$HI<- round(runif(n), 2)
  xloads <- by(
    gdata, 
    gdata$group1:gdata$group2, 
    function (x) {
      pattern <- paste0("^", unique(x$group1), ":", unique(x$group2))
      this.param <- param[grepl(pattern, names(param))]
      loads <- rnbinom(
        n = nrow(x), 
        size = param["k"],
        mu = MeanLoad(
          L1 = this.param[grepl("\\.inter", names(this.param))],
          L2 = this.param[grepl("\\.growth", names(this.param))],
          alpha = this.param[grepl("\\.alpha", names(this.param))],
          hybridIndex = x$HI
        )
      )
      cbind(x, loads)
    }
  )
  as.data.frame(do.call("rbind", xloads))
}

set.seed(5)
simdata <- SimulatedData(simparams, 1000)

##### Input end #####

fit <- mle2(
  loads ~ dnbinom(mu = MeanLoad(L1, L2, alpha, hybridIndex = HI), size = abs(k)),
  parameters = list(L1~group1, L2~group1, alpha~group1, k~group1),
  data = simdata,
  start = list(L1 = 10, L2 = 10, alpha = 0, k = 1)
  # control = list(fnscale=-1)
  # method="L-BFGS-B",
  # lower = c(L1 = 0, L2 = 0, k = 0.01),
  # upper = c(L1 = 200, L2 = 200, alpha = 10, k = 10)
)
fit
mypred <- data.frame(prediction = predict(fit),
                     HI = simdata$HI,
                     group1 = simdata$group1)

ggplot() +
  geom_point(aes(x = HI, y = loads, color = group1), data = simdata) +
  geom_point(aes(x = HI, y = prediction, fill = group1), pch = 21, size = 3, data = mypred) +
  scale_fill_manual(values = c("red", "darkblue")) +
  theme_linedraw()

##
fitwormsA <- mle2(
  Aspiculuris.Syphacia ~ dnbinom(mu = MeanLoad(L1, L2, alpha, HI), size = 1/abs(Aggregation(A1, A2, Z, HI))),
  #  parameters = list(L1~Sex, L2~Sex, alpha~Sex, A1~Sex, A2~Sex, Z~Sex),
  data = Joelle_data,
  start = list(L1 = 10, L2 = 10, alpha = 0, A1 = 1, A2 = 1, Z = 0),
  method="L-BFGS-B",
  lower = c(L1 = 0, L2 = 0, A1 = 0, A2 = 0),
  upper = c(L1 = 200, L2 = 200, A1 = 1000, A2 = 1000)
)
summary(fitwormsA)

mypred <- data.frame(prediction = predict(fitwormsA),
                     HI = Joelle_data$HI,
                     Sex = Joelle_data$Sex)

ggplot() +
  scale_color_manual(values = c("red", "darkblue")) +
  geom_point(aes(x = HI, y = Aspiculuris.Syphacia, color = Sex), 
             data = Joelle_data) +
  geom_line(aes(x = HI, y = prediction),
            data = mypred) +
  scale_y_continuous(trans='log10')+
  theme_linedraw()

##
##
fitwormsT <- mle2(
  Taenia ~ dnbinom(mu = MeanLoad(L1, L2, alpha, HI), 
                   size = 1/abs(Aggregation(A1, A2, Z, HI))),
  #  parameters = list(L1~Sex, L2~Sex, alpha~Sex, A1~Sex, A2~Sex, Z~Sex),
  data = Joelle_data,
  start = list(L1 = 10, L2 = 10, alpha = 0, A1 = 1, A2 = 1, Z = 0),
  method="L-BFGS-B",
  lower = c(L1 = 0, L2 = 0, A1 = 0, A2 = 0),
  upper = c(L1 = 200, L2 = 200, A1 = 1000, A2 = 1000)
)
summary(fitwormsT)

mypred <- data.frame(prediction = predict(fitwormsT),
                     HI = Joelle_data$HI,
                     Sex = Joelle_data$Sex)

ggplot() +
  scale_color_manual(values = c("red", "darkblue")) +
  geom_point(aes(x = HI, y = Taenia, color = Sex), 
             data = Joelle_data) +
  geom_line(aes(x = HI, y = prediction),
            data = mypred) +
  scale_y_continuous(trans='log10')+
  theme_linedraw()



########### TO CLEAN AFTER



fitJ <- mle2(
  loads ~ dnbinom(mu = MeanLoad(L1, L2, alpha, hybridIndex = HI), size = abs(k)),
  parameters = list(L1~group1, L2~group1, alpha~group1, k~group1),
  data = simdata,
  start = list(L1 = 10, L2 = 10, alpha = 0, k = 1)
  # control = list(fnscale=-1)
  # method="L-BFGS-B",
  # lower = c(L1 = 0, L2 = 0, k = 0.01),
  # upper = c(L1 = 200, L2 = 200, alpha = 10, k = 10)
)
fit
mypred <- data.frame(prediction = predict(fit),
                     HI = simdata$HI,
                     group1 = simdata$group1)

ggplot() +
  geom_point(aes(x = HI, y = loads, color = group1), data = simdata) +
  geom_point(aes(x = HI, y = prediction, fill = group1), pch = 21, size = 3, data = mypred) +
  scale_fill_manual(values = c("red", "darkblue")) +
  theme_linedraw()



plot(predict(fit))
summary(fit)
# confint(fit, parm = c("alpha.(Intercept)"))
anova(fit)

y <- predict(fit, list(L1= seq(0,1,0.01)), type = "response")


fit


plot(simdata$HI, simdata$loads, pch = 16, xlab = "HI", ylab = "loads")
lines(y, seq(0,1,0.01))

y
##############

fit0 <- mle2(
  loads ~ dnbinom(mu = MeanLoad(L1, L2, alpha, HI), size = abs(k)),
  data = simdata,
  start = list(L1 = 10, L2 = 10, alpha = 0, k = 1)
)
summary(fit0)

fit1 <- mle2(
  loads ~ dnbinom(mu = MeanLoad(L1, L2, alpha, HI), size = abs(k)),
  parameters = list(L1~group1, L2~group1, alpha~group1, k~group1),
  data = simdata,
  start = list(L1 = 10, L2 = 10, alpha = 0, k = 1)
)
summary(fit1)

fit2 <- mle2(
  loads ~ dnbinom(mu = MeanLoad(L1, L2, alpha, HI), size = abs(k)),
  parameters = list(L1~group1*group2, L2~group1*group2, alpha~group1*group2, k~group1*group2),
  data = simdata,
  start = list(L1 = 10, L2 = 10, alpha = 0, k = 1)
)
summary(fit2)

anova(fit0, fit1, fit3)
anova(fit0, fit3)

range(simdata$loads)

summary(fit)
confint(fit)
anova(fit)

############# Real data


fitworms <- mle2(
  Trichuris ~ dnbinom(mu = MeanLoad(L1, L2, alpha, HI), size = abs(k)),
  parameters = list(L1~Sex, L2~Sex, alpha~Sex, k~Sex),
  data = Joelle_data,
  start = list(L1 = 10, L2 = 10, alpha = 0, k = 1)
)

plot(fitworms)
profileworms <- profile(fitworms)


class(fitworms)


fitwormsA <- mle2(
  Aspiculuris.Syphacia ~ dnbinom(mu = MeanLoad(L1, L2, alpha, HI), size = 1/abs(Aggregation(A1, A2, Z, HI))),
  parameters = list(L1~Sex, L2~Sex, alpha~Sex, A1~Sex, A2~Sex, Z~Sex),
  data = Joelle_data,
  start = list(L1 = 10, L2 = 10, alpha = 0, A1 = 1, A2 = 1, Z = 0),
  method="L-BFGS-B",
  lower = c(L1 = 0, L2 = 0, A1 = 0, A2 = 0),
  upper = c(L1 = 200, L2 = 200, A1 = 1000, A2 = 1000)
)
summary(fitwormsA)

fitwormsTaenia <- mle2(
  Taenia ~ dnbinom(mu = MeanLoad(L1, L2, alpha, HI), size = 1/abs(Aggregation(A1, A2, Z, HI))),
  data = Joelle_data,
  start = list(L1 = 10, L2 = 10, alpha = 0, A1 = 1, A2 = 1, Z = 0),
  method="L-BFGS-B",
  lower = c(L1 = 0, L2 = 0, A1 = 0, A2 = 0, Z = -5, alpha = -10),
  upper = c(L1 = 200, L2 = 200, A1 = 1000, A2 = 1000, Z = 5, alpha = 5)
)
summary(fitwormsTaenia)






fitwormsA <- mle2(
  Aspiculuris.Syphacia ~ dnbinom(mu = MeanLoad(L1, L2, alpha, HI), size = 1/abs(Aggregation(A1, A2, Z, HI))),
  parameters = list(L1~Sex, L2~Sex, alpha~Sex, A1~Sex, A2~Sex, Z~Sex),
  data = Joelle_data,
  start = list(L1 = 10, L2 = 10, alpha = 0, A1 = 1, A2 = 1, Z = 0),
  method="L-BFGS-B",
  lower = c(L1 = 0, L2 = 0, A1 = 0, A2 = 0),
  upper = c(L1 = 200, L2 = 200, A1 = 1000, A2 = 1000)
)
summary(fitwormsA)

fitwormsA <- mle2(
  Aspiculuris.Syphacia ~ dnbinom(mu = MeanLoad(L1, L2, alpha, HI), size = 1/abs(Aggregation(A1, A2, Z, HI))),
  parameters = list(L1~Sex, L2~Sex, alpha~Sex, A1~Sex, A2~Sex, Z~Sex),
  data = Joelle_data,
  start = list(L1 = 10, L2 = 10, alpha = 0, A1 = 1, A2 = 1, Z = 0),
  method="L-BFGS-B",
  lower = c(L1 = 0, L2 = 0, A1 = 0, A2 = 0),
  upper = c(L1 = 200, L2 = 200, A1 = 1000, A2 = 1000)
)
summary(fitwormsA)





