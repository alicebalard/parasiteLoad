### source or load package once done with packaging
source("ML_functions.R")

source("UserInput.R")

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

################### input end ###################################

#########################    
SimulatedData <- function(param, n){
    gdata <- data.frame(group1 = rep(c("male", "female"), each=n/2),
                        group2 = sample(c("old", "young", "baby"),
                                        n, replace=TRUE))
    gdata$HI<- round(runif(n), 2)
    xloads <- by(gdata, gdata$group1:gdata$group2, function (x) {
        pattern <- paste0("^", unique(x$group1), ":", unique(x$group2))
        this.param <- param[grepl(pattern, names(param))]
        loads <- rnbinom(n = nrow(x), size = param["k"],
                         mu = MeanLoad(intercept=this.param[grepl("\\.inter",
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

LogLik(simdata, simpara, c("group1", "group2"))

### really bad when starting paramters just close to zero
### parameters
opt.para <- glm.hybrid(formula=loads~group2*HI*group1, data=simdata, "HI",
                     alpha.start=1, start.values=simpara)

glm.h1 <- glm.hybrid(formula=loads~group2*HI*group1, data=simdata, "HI",
                     alpha.start=1)

glm.h1.5 <- glm.hybrid(formula=loads~group2*HI*group1, data=simdata, "HI",
                     alpha.start=1.5)

glm.h1.9 <- glm.hybrid(formula=loads~group2*HI*group1, data=simdata, "HI",
                     alpha.start=1.9)

glm.h2.5 <- glm.hybrid(formula=loads~group2*HI*group1, data=simdata, "HI",
                     alpha.start=2.5)

## replace some of the parameters to not come via glm.nb (start.mod)
## but being entered manually (in a named vector, names have to be
## same as the model would assign)
glm.try <- glm.hybrid(formula=loads~group2*HI*group1, data=simdata, "HI",
                      alpha.start=2.5, start.values=simpara[5:8])

non.nb <- glm.hybrid(formula=loads~group2*HI*group1, data=simdata, "HI",
                     alpha.start=1.5, start.mod=glm.nb)

para.table <- cbind(simpara,
                    opt.sim = opt.para$par[names(simpara)],
                    opt.nb1 = glm.h1$opt.param[names(simpara)],
                    opt.nb1.5 = glm.h1.5$opt.param[names(simpara)],
                    opt.nb1.9 = glm.h1.9$opt.param[names(simpara)],
                    opt.nb2.5 = glm.h2.5$opt.param[names(simpara)])

opt.para$value
glm.h1$twologlik/2
glm.h1.5$twologlik/2
glm.h1.9$twologlik/2
glm.h2.5$twologlik/2

## pairs(para.table)

glm.hybrid(formula=loads~group2*HI, data=simdata, "HI")$twologlik/2

glm.hybrid(formula=loads~group1*HI, data=simdata, "HI")$twologlik/2

NBglm <- glm.nb(formula=loads~group1*group2*HI, data=simdata)
NBglm$twologlik/2

## should give error because it is not implemented
## glm.hybrid(formula=loads~(group2+HI+group1)^2, data=simdata, "HI")


## Alice trying to fix MLbounds:

ML_bounds_2(data = simdata, maxL = glm.h1$twologlik / 2, 
            group.name = c("group1", "group2"),
            opt.param = glm.h1$opt.param, 
            threshold = qchisq(p = 0.95, df = 1) / 2)


