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

### looks like it find the starting paramters quite well when starting
### with good parameters
opt.para <- optim(par = simpara,
                  fn = LogLik, ## function to be maximized
                  control = list(fnscale=-1),
                  method = "L-BFGS-B",
                  data = simdata,
                  group.name=c("group1", "group2"))

### really bad when starting paramters just close to zero
### parameters
glm.h1 <- glm.hybrid(formula=loads~group2*HI*group1, data=simdata, "HI",
                     alpha.start=1)

glm.h1.5 <- glm.hybrid(formula=loads~group2*HI*group1, data=simdata, "HI",
                     alpha.start=1.5)

glm.h1.9 <- glm.hybrid(formula=loads~group2*HI*group1, data=simdata, "HI",
                     alpha.start=1.9)

glm.h2.5 <- glm.hybrid(formula=loads~group2*HI*group1, data=simdata, "HI",
                     alpha.start=2.5)

para.table <- cbind(simpara,
                    opt.sim = opt.para$par[names(simpara)],
                    opt.nb1 = glm.h1$opt.param[names(simpara)],
                    opt.nb1.5 = glm.h1.5$opt.param[names(simpara)],
                    opt.nb1.9 = glm.h1.9[[1]]$opt.param[names(simpara)],
                    opt.nb2.5 = glm.h2.5[[1]]$opt.param[names(simpara)])

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
glm.hybrid(formula=loads~(group2+HI+group1)^2, data=simdata, "HI")


