

### source or load package once done with packaging
source("ParasiteLoadCoverageTest.R")

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
