set.seed(1001)
lymax <- c(0,2)
lhalf <- 0

install.packages("bbmle")
library(bbmle)

x <- sort(runif(200))
g <- factor(sample(c("a","b"),200,replace=TRUE))
y <- rnbinom(200,mu=exp(lymax[g])/(1+x/exp(lhalf)),size=2)
d2 <- data.frame(x,g,y)

fit3 <- mle2(
  y~dnbinom(
    mu = exp(lymax) / (1+x/exp(lhalf)) ,
    size = exp(logk)
  ),
  parameters = list(lymax~g),
  data = d2,
  start = list(lymax=0,lhalf=0,logk=0)
)

plop <- function (lymax, x, lhalf){
  exp(lymax) / (1+x/exp(lhalf))
}

fit3 <- mle2(y~dnbinom(mu=plop(lymax, x, lhalf),size=exp(logk)),
             parameters=list(lymax~g),data=d2,
             start=list(lymax=0,lhalf=0,logk=0))

summary(fit3)


simpara <- c(k = 2,
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

################## input end ##################


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
simdata <- SimulatedData(simpara, 1000)

MeanLoad <- function(L1, L2, alpha, hybridIndex){
  heterozygoty <- 2 * hybridIndex * (1 - hybridIndex)
  mean <- (L1 + (L2 - L1) * hybridIndex) * (1 - alpha * heterozygoty)
  mean <- sapply(mean, function(x) {
    return(max(x, 0.01))
  })
  return(mean)
}

fit <- mle2(
  loads ~ dnbinom(mu = MeanLoad(L1, L2, alpha, HI), size = abs(k)),
  parameters = list(L1~group1*group2, L2~group1*group2, alpha~group1*group2, k~group1*group2),
  data = simdata,
  start = list(L1 = 10, L2 = 10, alpha = 0, k = 1)
  # control = list(fnscale=-1)
  # method="L-BFGS-B",
  # lower = c(L1 = 0, L2 = 0, k = 0.01),
  # upper = c(L1 = 200, L2 = 200, alpha = 10, k = 10)
)

range(simdata$loads)

summary(fit)
confint(fit)
anova(fit)


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






Aggregation <- function(A1, A2, Z, HI){
  He <- 2 * HI * (1 - HI)
  (A1 + (A2 - A1) * HI) + Z * He 
} 

fitworms <- mle2(
  Trichuris ~ dnbinom(mu = MeanLoad(L1, L2, alpha, HI), size = abs(k)),
  parameters = list(L1~Sex, L2~Sex, alpha~Sex, k~Sex),
  data = Joelle_data,
  start = list(L1 = 10, L2 = 10, alpha = 0, k = 1)
)


Joelle_data <- read.csv("../examples/Reproduction_WATWM/EvolutionFinalData.csv")
Joelle_data[is.na(Joelle_data)] <- 0
Joelle_data$Aspiculuris.Syphacia
#G3 <- glm.hybrid(Trichuris ~ HI * Sex, data = Joelle_data, alpha.along = "HI")

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





