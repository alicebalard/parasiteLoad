### source or load package once done with packaging
source("ML_functions.R")


### This should come from input script... ###########################
param <- vector(mode="list", length=3)
names(param) <- c("group1", "group2", "group1:group2")

## parameters as starting values
param[[1]] <- c(k=0.5, alpha=0.1,
                int.male=0.1,
                int.female=0.2,
                growth.male=0.1,
                growth.female=0.2)

## parameters as starting values
param[[2]] <- c(k=0.5, alpha=0.1,
                int.baby=0.2,
                int.young=0.02,
                int.old=0.02,
                growth.baby=0.1,
                growth.young=0.1,
                growth.old=0.2)

## parameters as starting values
param[[3]] <- c(k = 0.01, alpha = 0.01,
                "int.male:old" = 0.01,
                "int.male:young" = 0.01,
                "int.male:baby" = 0.01,
                "int.female:old" = 0.01,
                "int.female:young" = 0.01,
                "int.female:baby" = 0.01,
                "growth.male:old" = 0.01,
                "growth.male:young" = 0.01,
                "growth.male:baby" = 0.01,
                "growth.female:old" = 0.01,
                "growth.female:young" = 0.01,
                "growth.female:baby" = 0.01)

## parameters for simulation
simpara <- c(k = 2, alpha = 1.64,
             "int.male:old" = 1,
             "int.male:young" = 12,
             "int.male:baby" = 22,
             "int.female:old" = 2,
             "int.female:young" = 23,
             "int.female:baby" = 32,
             "growth.male:old" = 2,
             "growth.male:young" = 1,
             "growth.male:baby" = -4,
             "growth.female:old" = 6,
             "growth.female:young" = 0,
             "growth.female:baby" = -8)

################### input end ###################################

#########################    
SimulatedData <- function(param, n){
    gdata <- data.frame(group1 = rep(c("male", "female"), each=n/2),
                        group2 = sample(c("old", "young", "baby"),
                                        n, replace=TRUE))
    gdata$HI<- round(runif(n), 2)
    xloads <- by(gdata, gdata$group1:gdata$group2, function (x) {
        pattern <- paste0("\\.", unique(x$group1), ":", unique(x$group2))
        this.param <- param[grepl(pattern, names(param))]
        loads <- rnbinom(n = nrow(x), size = param["k"],
                         mu = MeanLoad(intercept=this.param[grepl("int",
                                                                  names(this.param))],
                                       growth=this.param[grepl("growth",
                                                               names(this.param))],
                                       alpha=param["alpha"],
                                       HI=x$HI))
        cbind(x, loads)
    })
    as.data.frame(do.call("rbind", xloads))
}

set.seed(5)
simdata <- SimulatedData(simpara, 1000)

LogLik(simdata, param=param[["group1"]],
       group.name=names(param["group1"]))

LogLik(simdata, param=param[["group2"]],
       group.name=names(param["group2"]))

LogLik(simdata, param=param[["group1:group2"]],
       group.name=names(param["group1:group2"]))

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

cbind(simpara, opt.para)

### really bad when starting paramters just close to zero
### parameters
opt.para <- optim(par = c(param[[3]]),
                  fn = LogLik, ## function to be maximized
                  control = list(fnscale=-1),
                  data = simdata,
                  group.name="group1:group2")$par

cbind(simpara, opt.para)
