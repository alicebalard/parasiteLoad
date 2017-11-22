# library(devtools)
# install_github("alicebalard/Parasite_Load")

source("../R/UserInput.R")
source("../R/ML_functions.R")

simpara <- c(k = 2, alpha = 1.92,
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

## Test on worms data (local alice) :
Joelle_data <- na.omit(read.csv("../../../EvolutionFinalData.csv"))

################## input end ################## 

SimulatedData <- function(param, n){
    gdata <- data.frame(group1 = rep(c("male", "female"), each=n/2),
                        group2 = sample(c("old", "young", "baby"),
                                        n, replace=TRUE))
    gdata$HI<- round(runif(n), 2)
    xloads <- by(gdata, gdata$group1:gdata$group2, function (x) {
      pattern <- paste0("^", unique(x$group1), ":", unique(x$group2))
      this.param <- param[grepl(pattern, names(param))]
      loads <- rnbinom(n = nrow(x), size = param["k"],
                       mu = glm.hybrid:::MeanLoad(intercept=this.param[grepl("\\.inter",
                                                                             names(this.param))],
                                                  growth=this.param[grepl("\\.growth",
                                                                          names(this.param))],
                                                  alpha=param["alpha"],
                                                  HI=x$HI))
      cbind(x, loads)
    })
    as.data.frame(do.call("rbind", xloads))
}

set.seed(5)
simdata <- SimulatedData(simpara, 1000)

################## simulation end ################## 

################## Test : one discrete group OK ##################

G1 <- glm.hybrid(loads ~ HI * group1, data = simdata, alpha.along = "HI", alpha.start = 1)

G2 <- glm.hybrid(loads ~ HI * group1, data = simdata, alpha.along = "HI")

G3 <- glm.hybrid(Trichuris ~ HI * Sex, data = Joelle_data, alpha.along = "HI")

################## Test : two discrete groups OK ##################

G4 <- glm.hybrid(loads ~ HI * group1 * group2, data = simdata, alpha.along = "HI")

################## Test : starting parameters OK ##################

# really bad when starting parameters just close to zero
opt.para <- glm.hybrid(formula=loads ~ HI*group1*group2, data=simdata, alpha.along = "HI",
                       alpha.start=1, start.values=simpara)

glm.h1 <- glm.hybrid(formula=loads~HI*group1*group2, data=simdata, alpha.along = "HI",
                     alpha.start=1)

glm.h1.5 <- glm.hybrid(formula=loads~ HI*group1*group2, data=simdata, alpha.along = "HI",
                       alpha.start=1.5)

glm.h1.9 <- glm.hybrid(formula=loads~ HI*group1*group2, data=simdata, alpha.along = "HI",
                       alpha.start=1.9)

glm.h2.5 <- glm.hybrid(formula=loads~HI*group1*group2, data=simdata, "HI",
                       alpha.start=2.5)

para.table <- cbind(simpara,
                    opt.sim = opt.para$opt.param[names(simpara)],
                    opt.nb1 = glm.h1$opt.param[names(simpara)],
                    opt.nb1.5 = glm.h1.5$opt.param[names(simpara)],
                    opt.nb1.9 = glm.h1.9$opt.param[names(simpara)],
                    opt.nb2.5 = glm.h2.5$opt.param[names(simpara)])

glm.h1$opt.param
opt.para
names(simpara)
para.table

opt.para$value
glm.h1$twologlik/2
glm.h1.5$twologlik/2
glm.h1.9$twologlik/2
glm.h2.5$twologlik/2

## replace some of the parameters to not come via glm.nb (start.mod)
## but being entered manually (in a named vector, names have to be
## same as the model would assign)
glm.try <- glm.hybrid(formula=loads~HI*group1, data=simdata, "HI",
                      alpha.start=2.5, start.values=simpara[5:8])

non.nb <- glm.hybrid(formula=loads~HI*group1, data=simdata, "HI",
                     alpha.start=1.5, start.mod = MASS::glm.nb)


