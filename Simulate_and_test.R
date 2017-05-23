

### source or load package once done with packaging
source("ML_functions.R")

#########################    
## Simulate data for 1 independant variable (male/female)
## NB: add later "age", and generalise...
SimulatedData <- function(NindsM, NindsF, interceptM, interceptF, growthM, growthF, alpha, k){
    ##    set.seed(5)
    IndHIsM<- round(runif(NindsM), 2)
    IndHIsF <- round(runif(NindsF), 2)
    MalLoads  <- data.frame(HI = IndHIsM,
                            loads= rnbinom(NindsM, size=k,
                                           mu= MeanLoad(interceptM, growthM, alpha, IndHIsM)),
                            group1 = "male")
    FemLoads  <- data.frame(HI = IndHIsF,
                            loads= rnbinom(NindsF, size=k,
                                           mu= MeanLoad(interceptF, growthF, alpha, IndHIsF)),
                            group1 = "female")
    ## Bind the dataframes:
    return(rbind(MalLoads, FemLoads))
}

## Example:
## 1.Choose our parameters for the simulation:
NindsF_exp <- 80; NindsM_exp <- 85; interceptF_exp <- 15; interceptM_exp <- 10
growthF_exp <- 2; growthM_exp <- 0; alpha_exp <- 2; k_exp <- 4
## And simulate data:
alicedata <- SimulatedData(NindsM_exp, NindsF_exp, interceptM_exp, interceptF_exp,
                              growthM_exp, growthF_exp, alpha_exp, k_exp)


## Trial
fakedata <- cbind(rbind(head(alicedata),tail(alicedata)),
                  group2 = c(rep("old", 3), rep("young", 9)))


param <- vector(mode="list", length=3)
names(param) <- c("group1", "group2", "group1:group2")


param[[1]] <- c(k=0.5, alpha=0.1, int.male=0.1,
                       int.female=0.2, growth.male=0.1, growth.female=0.2)

param[[2]] <- c(k=0.5, alpha=0.1, int.group2.young=0.02,
                       int.old=0.02, growth.young=0.1, growth.old=0.2)

param[[3]] <- c(k=0.5, alpha=0.1,
                rep(0.01, length(with(fakedata, levels(group1:group2)))*2))

names(param[[3]]) <- c("k", "alpha", 
                       paste0("int.", with(fakedata, levels(group1:group2))),
                       paste0("growth.", with(fakedata, levels(group1:group2))))

PrLoadbyRow(fakedata, sub.param=param[["group1"]])

PrLoadbyRow(fakedata, sub.param=param[["group1:group2"]])

## works for non-ineraction effects
LikelihoodFunction(fakedata,
                   ugly.param=param[[2]], name=names(param[2]))

## works for non-ineraction effects
LikelihoodFunction(fakedata,
                   ugly.param=param[[1]], name=names(param[1]))

## does not work for ineraction effects
LikelihoodFunction(fakedata,
                   ugly.param=param[[3]], name=names(param[3]))

optim(par = c(param[["group2"]]),
      fn = LikelihoodFunction, ## function to be maximized
##      upper=c(male=list(1, 1, 1, 1),
##              female=list(1, 1, 1, 1)),
##      lower=c(male=list(-0.1, -0.1, 0.001, 0.0001),
##              female=list(-0.1, -0.1, 0.001, 0.0001)),
      control = list(fnscale=-1), ##turn the default minimizer into
                                  ##maximizer
      ##      method = "L-BFGS-B",
      data = fakedata,
      name="group2")
