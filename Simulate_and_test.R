### source or load package once done with packaging
source("ML_functions.R")

### This should come from input script... ###########################
param <- vector(mode="list", length=3)
names(param) <- c("group1", "group2", "group1:group2")

## parameters as starting values
param[[1]] <- c(k=0.5, alpha=0.1,
                male.inter=0.1,
                female.inter=0.2,
                male.growth=0.1,
                female.growth=0.2)

## parameters as starting values
param[[2]] <- c(k=0.5, alpha=0.1,
                baby.inter=0.2,
                young.inter=0.02,
                old.inter=0.02,
                baby.growth=0.1,
                young.growth=0.1,
                old.growth=0.2)

## parameters as starting values
param[[3]] <- c(k = 0.01, alpha = 0.01,
                "male:old.inter" = 0.01,
                "male:young.inter" = 0.01,
                "male:baby.inter" = 0.01,
                "female:old.inter" = 0.01,
                "female:young.inter" = 0.01,
                "female:baby.inter" = 0.01,
                "male:old.growth" = 0.01,
                "male:young.growth" = 0.01,
                "male:baby.growth" = 0.01,
                "female:old.growth" = 0.01,
                "female:young.growth" = 0.01,
                "female:baby.growth" = 0.01)

## parameters for simulation
simpara <- c(k = 2, alpha = 0.64,
             "male:old.inter" = 14,
             "male:young.inter" = 12,
             "male:baby.inter" = 22,
             "female:old.inter" = 12,
             "female:young.inter" = 23,
             "female:baby.inter" = 32,
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

library(MASS)
## nb <- glm.nb(loads~group2*HI+group1*HI, simdata)
nb <- glm.nb(loads~group2*HI*group1, simdata)

library(effects)

nb.e  <- allEffects(nb, xlevels=2)

### plot the effects if needed ## for now steal code from it
## effects:::plot.efflist(nb.effects)

get.params.from.glm.nb <- function (x){
    var.names <- sort(names(x[["variables"]]))
    original.inverse <- x$transformation$inverse
    y <- as.data.frame(x, transform = original.inverse)
    nb.param <- y[, "fit"]
    names(nb.param) <- apply(y[,var.names], 1,
                             paste, collapse=":")
    nb.param[grepl("\\:1", names(nb.param))] <-
        nb.param[grepl("\\:1", names(nb.param))] - 
        nb.param[grepl("\\:0", names(nb.param))]
    names(nb.param) <- gsub("\\:1", "\\.growth", names(nb.param))
    names(nb.param) <- gsub("\\:0", "\\.inter", names(nb.param))
    nb.param
}

nb.effects <- get.params.from.glm.nb(nb.e[[1]])

nb.param <- c(k=nb$theta, alpha=0.1, nb.effects)

## the glm paramter estimates
merge(nb.param, simpara, by=0)

LogLik(simdata, param=nb.param,
       group.name="group1:group2")
                                  
LogLik(simdata, param=simpara,
       group.name="group2:group1")
                                  

LogLik(simdata, param=param[["group1:group2"]],
       group.name=names(param["group1:group2"]))

LogLik(simdata, param=param[["group1"]],
       group.name=names(param["group1"]))

LogLik(simdata, param=param[["group2"]],
       group.name=names(param["group2"]))

all.optim <- lapply(names(param), function (x){
    optim(par = c(param[[x]]),
          fn = LogLik, ## function to be maximized
          control = list(fnscale=-1),
          data = simdata,
          group.name=x)
})

names(all.optim) <- names(param)

all.optim

### looks like it find the starting paramters quite well when starting
### with good parameters
opt.para <- optim(par = c(simpara),
                  fn = LogLik, ## function to be maximized
                  control = list(fnscale=-1),
                  data = simdata,
                  group.name="group1:group2")$par

### really bad when starting paramters just close to zero
### parameters
opt.para.nb <- optim(par = nb.param,
                  fn = LogLik, ## function to be maximized
                  control = list(fnscale=-1),
                  data = simdata,
                  group.name="group2:group1")$par

cbind(simpara,
      opt.sim = opt.para[names(simpara)],
      opt.nb = opt.para.nb[names(simpara)])
