########################## On reprend... 8 Mai
rm(list=ls()) 
library(ggplot2)
library(reshape)
library(gridExtra)
library(tableHTML)

#########################    
## Generate an hypotheses matrix:
## 4 parameters (with various number for each):
## k, alpha, intercept, growth ## 
## All hypotheses (so far...):
## H0: k, alpha, intercept
## H1: k, alpha, intercept, growth
## H2: k, alpha, intercept1, intercept2
## H3: k, alpha, intercept1, intercept2, growth1, growth2
HypMat <- t(matrix(c(1,1,1,0,0,
                     1,1,1,1,0,
                     1,1,1,0,1,
                     1,1,1,1,1),
                   nrow=5))

start <- c(k=1, alpha=0, intercept=20, growth=2)
lower <- c(k=1, alpha=-3, intercept=0, growth=0)
upper <- c(k=8, alpha=3, intercept=30, growth=10)

control <- list(HypMat, start, lower, upper)

class(control) <- "MLcontrol"


#########################
## MeanLoad model: How does the mean load vary depending on the parameters?
MeanLoad <- function(intercept, growth, alpha, HI){
    (intercept + growth*HI)*(1 - alpha*2*HI*(1 - HI))
}

#########################    
## Simulate data for 1 independant variable (male/female)
## NB: add later "age", and generalise...
SimulatedData <- function(NindsM, NindsF, interceptM, interceptF, growthM, growthF, alpha, k){
    ##    set.seed(5)
    IndHIsM<- round(runif(NindsM), 2)
    IndHIsF <- round(runif(NindsF), 2)
    MalLoads  <- data.frame(HI = IndHIsM, loads= rnbinom(NindsM, size=k, mu= MeanLoad(interceptM, growthM, alpha, IndHIsM)), group1 = "male")
    FemLoads  <- data.frame(HI = IndHIsF, loads= rnbinom(NindsF, size=k, mu= MeanLoad(interceptF, growthF, alpha, IndHIsF)), group1 = "female")
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

#########################    
## PrLoad: Probability an Ind has some
## observed Load, given its HI (and the model)        
## Generalised, ie the MeanLoad model is passed as an argument
PrLoad <- function(data, k, intercept, growth, alpha, ind){
    HI <- data$HI[ind]
    Load <- data$loads[ind]
    ## Always a good idea to make contributions to the likelihood function Bombproof, and never return zero
    max(10^-20, dnbinom(Load, size=abs(k), mu=abs(MeanLoad(intercept, growth, alpha, HI))))
}

## Test:
PrLoad(alicedata, k=2, intercept=10, growth=2, alpha=1, ind=2)

### By row?
PrLoadRow <- function(datarow, k, intercept, growth, alpha){
    HI <- datarow[[1]]  ## e.g. alicedata[1,][[1]]
    Load <- datarow[[2]]
    ## Always a good idea to make contributions to the likelihood function Bombproof, and never return zero
    max(10^-20, dnbinom(Load, size=abs(k), mu=abs(MeanLoad(intercept, growth, alpha, HI))))
}

## Test:
PrLoadRow(alicedata[2,], k=2, intercept=10, growth=2, alpha=1)

#########################
## Likelihood of a given subgroup of the data (aka male)

## nLev : number of levels to take into consideration
nLev <- length(levels(alicedata[,3:ncol(alicedata)]))



## Trial
fakedata <- cbind(rbind(head(alicedata),tail(alicedata)), group2 = c(rep("old",3),rep("young",9)))




lk <- function(data, k, intercept, growth, alpha){
    L <- 0
    L <- L + log(PrLoad(datarow, k, intercept, growth, alpha))
by(fakedata[,1], fakedata[,3:], summary)







L <- L + log(PrLoad(data, k, intercept, growth, alpha, ind))








#########################
## The likelihood function over a set of inds
## param: k, alpha, mM, mF, slapeM, slapeF
LikelihoodFunction <- function(param, data, multi) { 
    ## Likelihood initialisation:
    S <- 0
    ## Parameters X HypMat line:
    mypar <- param*multi
    ## Assign each parameter to its value:
    for (i in 1:length(names(mypar))) assign(names(mypar[i]),as.numeric(mypar[i]))
    ## Group the m & slape of each level of Group1:
    m <- c(mM,mF)
    slape <- c(slapeM, slapeF)
    ## In the case of H0 and H1, there are NO independant variables considered:
    if (length(m[m!=0]) == 1){
        mB <- m[m!=0]
        ## slape is either 0 or a given value:
        if (length(slape[slape!=0]) == 1) {slapeB <- slape[slape!=0] ## H1 case
        } else { slapeB <- 0 } ## H0 case
        ## For all individuals in each level:
        for (ind in 1:nrow(data)){
            S <- S + log(PrLoad(data, k, mB, slapeB, alpha, ind))
        }
    } else {
        ## We consider so far 1 independant variable: sex, with 2 levels (TO IMPROVE)
        by(data, data$group1, function(x)  {

          
            S <- S + log(PrLoad(subdata, k, m[lev], slape[lev], alpha, indiv))


        ## Run over all its levels:
        for (lev in 1:2) {
            ## subset the dataframe:
            subdata <- data[levels(data[ ,3]) %in% levels(data[ ,3])[lev],]
            ## For all individuals of this level:
            for (indiv in 1:nrow(subdata)){
                S <- S + log(PrLoad(subdata, k, m[lev], slape[lev], alpha, indiv))
            }
        }
    }
    return(S)
}








##########################################









#########################
## The likelihood function over a set of inds
## param: k, alpha, mM, mF, DmM, DmF
LikelihoodFunction <- function(param, data, multi) { 
    ## Likelihood initialisation:
    S <- 0
    ## Parameters X HypMat line:
    mypar <- param*multi
    ## Assign each parameter to its value:
    for (i in 1:length(names(mypar))) assign(names(mypar[i]),as.numeric(mypar[i]))
    ## Group the m & Dm of each level of Group1:
    m <- c(mM,mF)
    Dm <- c(DmM, DmF)
    ## In the case of H0 and H1, there are NO independant variables considered:
    if (length(m[m!=0]) == 1){
        mB <- m[m!=0]
        ## Dm is either 0 or a given value:
        if (length(Dm[Dm!=0]) == 1) {DmB <- Dm[Dm!=0]
        } else { DmB <- 0 }
        ## For all individuals in each level:
        for (ind in 1:nrow(data)){
            S <- S + log(PrLoad(data, k, mB, DmB, alpha, ind))
        }
    } else {
        ## We consider so far 1 independant variable: sex, with 2 levels (TO IMPROVE)


        ## Run over all its levels:
        for (lev in 1:2) {
            ## subset the dataframe:
            subdata <- data[levels(data[ ,3]) %in% levels(data[ ,3])[lev],]
            ## For all individuals of this level:
            for (indiv in 1:nrow(subdata)){
                S <- S + log(PrLoad(subdata, k, m[lev], Dm[lev], alpha, indiv))
            }
        }
    }
    return(S)
}





LikelihoodFunctionOLD <- function(param, data, multi) { 
    ## Likelihood initialisation:
    S <- 0
    ## Parameters X HypMat line:
    mypar <- param*multi
    ## Assign each parameter to its value:
    for (i in 1:length(names(mypar))) assign(names(mypar[i]),as.numeric(mypar[i]))
    ## Group the m & Dm of each level of Group1:
    m <- c(mM,mF)
    Dm <- c(DmM, DmF)
    ## In the case of H0 and H1, there are NO independant variables considered:
    if (length(m[m!=0]) == 1){
        mB <- m[m!=0]
        ## Dm is either 0 or a given value:
        if (length(Dm[Dm!=0]) == 1) {DmB <- Dm[Dm!=0]
        } else { DmB <- 0 }
        ## For all individuals in each level:
        for (ind in 1:nrow(data)){
            S <- S + log(PrLoad(data, k, mB, DmB, alpha, ind))
        }
    } else {
        ## We consider so far 1 independant variable: sex, with 2 levels (TO IMPROVE)
        ## Run over all its levels:
        for (lev in 1:2) {
            ## subset the dataframe:
            subdata <- data[levels(data[ ,3]) %in% levels(data[ ,3])[lev],]
            ## For all individuals of this level:
            for (indiv in 1:nrow(subdata)){
                S <- S + log(PrLoad(subdata, k, m[lev], Dm[lev], alpha, indiv))
            }
        }
    }
    return(S)
}

## Test
k <- 5; kmin <- 0; kmax <- 8
alpha <- 0; alphamin <- -2; alphamax <- 2
m1 <- 0; m1min <- 0; m1max <- 30
Dm <- 0; Dmmin <- 0; Dmmax <- 1

LikelihoodFunctionOLD(c(k=k, alpha=alpha, mM=m1, mF=m2, DmM=Dm1, DmF=Dm2), alicedata, HypMat[1,])












#########################
## The Maximum likelihood analysis

## 1. Function to maximize the likelihood:
mymle <- function(data, lower, upper, start, multi){
    ## Set limits to the lower and upper parameters
    ## in the case where they must stay NULL (or close to):
    lower <- lower*multi - 0.1 
    upper <- upper*multi + 0.1
    start <- start*multi
    ## And run the optimisation:    
    optim(par = start, ## initial values for the parameters
          fn = LikelihoodFunction, ## function to be maximized
          lower = lower, ## lower bounds for the parameters
          upper = upper, ## upper bounds for the parameters
          control = list(fnscale=-1), ##turn the default minimizer into maximizer
          method = "L-BFGS-B",
          data = alicedata,
          multi = multi) ## extra param for LikelihoodFunction
}

## 2. Function to maximize the parameters in MLE+/- 2 :

## MLEbounds : optim function with constraint on the search
UpperBoundsFunction <- function(data, lower, upper, start,MyMLE, multi){
    ## Set limits to the lower and upper parameters
    ## in the case where they must stay NULL (or close to):
    lower <- lower*multi - 0.01 
    upper <- upper*multi + 0.01
    start <- start*multi
    ## empty dataframe to store the data
    summaryVec <- vector()
    max <- length(MyMLE$par)
    for (rankpar in 1:max){ ## run over the parameters
        ## Functional constraint on search (L > MLE-2) + max&min each param
        OptimBounds <- function(param){
            LK <- LikelihoodFunction(param, data, multi)
            if (LK <= MyMLE$value - 2){ 
                param[rankpar] <- -100
            }
            return(param[rankpar])
        }
        MyBoundsMax <- optim(par = MyMLE$par, ## initial values for the parameters
                             fn = OptimBounds, ## function to be maximized
                             lower = lower, ## lower bounds for the parameters
                             upper = upper, ## upper bounds for the parameters
                             method = "L-BFGS-B",
                             control = list(fnscale=-1)) ##turn the default minimizer into maximizer
        summaryVec <- c(summaryVec, MyBoundsMax$par[rankpar])
    }
    return(summaryVec)
}

## 3. Function to minimize the parameters in MLE+/- 2 :

## MLEbounds : optim function with constraint on the search
LowerBoundsFunction <- function(data, lower, upper, start,MyMLE, multi){
    ## Set limits to the lower and upper parameters
    ## in the case where they must stay NULL (or close to):
    lower <- lower*multi - 0.01 
    upper <- upper*multi + 0.01
    start <- start*multi
    ## empty dataframe to store the data
    summaryVec <- vector()
    max <- length(MyMLE$par)
    for (rankpar in 1:max){ ## run over the parameters
        ## Functional constraint on search (L > MLE-2) + max&min each param
        OptimBounds <- function(param){
            LK <- LikelihoodFunction(param, data, multi)
            if (LK <= MyMLE$value - 2){ 
                param[rankpar] <- 100
            }
            return(param[rankpar])
        }
        MyBoundsMax <- optim(par = MyMLE$par, ## initial values for the parameters
                             fn = OptimBounds, ## function to be maximized
                             lower = lower, ## lower bounds for the parameters
                             upper = upper, ## upper bounds for the parameters
                             method = "L-BFGS-B")
        summaryVec <- c(summaryVec, MyBoundsMax$par[rankpar])
    }
    return(summaryVec)
}

## 4. Function to run the optimisations over the 4 hypotheses :
RunOptim <- function(data, lower, upper, start, multi) {
    MyMLE <- mymle(data, lower, upper, start, multi)
    LB <- LowerBoundsFunction(data, lower, upper, start, MyMLE, multi)
    UB <- UpperBoundsFunction(data, lower, upper, start, MyMLE, multi)
    return(list(MyMLE, LB, UB))
}
    
## Test:
data <- alicedata
start <- c(k=1, alpha=0, mM=20, mF=20, DmM=2, DmF=2)
lower <- c(k=1, alpha=-3, mM=0, mF=0, DmM=0, DmF=0)
upper <- c(k=8, alpha=3, mM=30, mF=30, DmM=10, DmF=10)

before <- Sys.time()
MySupaResult <- lapply(1:4, function(x) RunOptim(data, lower, upper, start, HypMat[x,]))
after <-  Sys.time()
print("The optimisation took that long:")
after - before



## Better visualisation:
MyResults <- function(x){
    DF <- t(data.frame(MySupaResult[[x]][[2]],MySupaResult[[x]][[1]]$par,MySupaResult[[x]][[3]]))
    rownames(DF) <- c("low", "est", "high")
    DF
}

lapply(1:4, function(x) MyResults(x))

## To compare with the simulation:
simpar <- c(k_exp, alpha_exp, mM_exp, mF_exp, DmM_exp, DmF_exp)

list(lapply(1:4, function(x) MyResults(x)), simpar)









############################
## Test if the difference between 2 likelihood is significant
    Gtest <- function(dLL, dDF){
        1 - pchisq(2*dLL, df=dDF) 
    }
    GtestOnNestedMLEs <- function(a){
        MyMLE1 <- Summary$MLE[1]
        MyMLE2 <- Summary$MLE[a+1]
        dLL <- round(abs(MyMLE1 - MyMLE2),4)
        dDF <- Summary$'Ddf(H0)'[a+1]
        p <- round(Gtest(dLL, dDF),4)
        Gresults <- c(dLL, dDF, p)
        return(Gresults)
    }
    dLL <- c(NA, GtestOnNestedMLEs(1)[1],GtestOnNestedMLEs(2)[1],GtestOnNestedMLEs(3)[1])
    p_Gtest <- c(NA, GtestOnNestedMLEs(1)[3],GtestOnNestedMLEs(2)[3],GtestOnNestedMLEs(3)[3])
    Summary <- cbind(Summary[1], dLL, p_Gtest, Summary[-1])
    return(Summary) 
}

## example:
k <- 5; kmin <- 0; kmax <- 8
alpha <- 0; alphamin <- -2; alphamax <- 2
m1 <- 0; m1min <- 0; m1max <- 30
Dm <- 0; Dmmin <- 0; Dmmax <- 1
RunOptim(alicedata, "L-BFGS-B",
                        k , kmin, kmax,
                        alpha, alphamin, alphamax,
                        m1, m1min, m1max,
                        Dm, Dmmin, Dmmax)


######## Let's calculate the coverage of the parameters
## (which % of time are the parameters taken to simulate data
## found in the interval estimated)
MyBS <- function(){
    ## 1.Choose our parameters for the simulation:
    NindsF_exp <- 80; NindsM_exp <- 85; mF_exp <- 15; mM_exp <- 22
    DmF_exp <- 5; DmM_exp <- 0; alpha_exp <- 1.1; k_exp <- 4
    ## And simulate data:
    TestData <- SimulatedData(NindsF_exp, NindsM_exp, mF_exp, mM_exp,
                               DmF_exp, DmM_exp, alpha_exp, k_exp)
    ## 2.Choose our parameters for start of search of MaxLik (should be different):
    k <- 1; kmin <- 0; kmax <- 8;
    alpha <- 1; alphamin <- -2; alphamax <- 2;
    m1 <- 1; m1min <- 0; m1max <- 30;
    Dm <- 1; Dmmin <- 0; Dmmax <- 10
    ## 3.Run the full optimisation:
    return(RunOptim(TestData, "L-BFGS-B",
                        k , kmin, kmax,
                        alpha, alphamin, alphamax,
                        m1, m1min, m1max,
                        Dm, Dmmin, Dmmax))
}

before <- Sys.time()
before
result1 <- MyBS()
after <-  Sys.time()
after
after - before

## Check if param_exp belongs to [parammin-parammax] for H3: 
myCoverageFunction <- function(myResult){ 
    myMins <- c("k_min", "alpha_min", "m1_min", "Dm1_min", "m2_min", "Dm2_min")
    myMaxs <- c("k_max", "alpha_max", "m1_max", "Dm1_max", "m2_max", "Dm2_max")
    myExps <- c(k_exp, alpha_exp, mF_exp, DmF_exp, mM_exp, DmM_exp)
    ## Check if param_exp belongs to [parammin-parammax] for H3:
    cov <- 0
    for (j in 1:length(myExps)){
        A <- vector()
        for (i in 1:N_BS){
            A <- c(A,
                   data.table::between(myExps[j],
                                       as.numeric(as.character(myResult[which(rownames(myResult)==myMins[j]),][[i]][4])),
                                       as.numeric(as.character(myResult[which(rownames(myResult)==myMaxs[j]),][[i]][4]))))
            cov[j] <- prop.table(table(A))[2]*100
        }
    }
    names(cov) <- c("k", "alpha", "m1", "Dm1", "m2", "Dm2")
    return(cov)
}

## NB: find the parameters chosen for simulation on top of the code:
MyCoverageFunction(1, result1)

result1



result1[which(rownames(result1)==myMins[1]),]

,
                                       as.numeric(as.character(myResult[which(rownames(myResult)==myMaxs[j]),][[i]][4]))))
            cov[j] <- prop.table(table(A))[2]*100
  
################ NOT WORKING ???
MyParallel <- function(N_BS){ ## N_BS <- number of repeats
    ##create cluster
    library(parallel)
    cl <- makeCluster(detectCores()-2)
    ##put objects in place that might be needed for the code
    clusterExport(cl,c("MeanLoad", "PrLoad", "LikelihoodFunction", "mymle","UpperBoundsFunction", "LowerBoundsFunction", "RunOptim", "SimulatedData", "melt", "MyBS"))
                                        #... then parallel replicate...
    myResult <- parSapply(cl, 1:N_BS, function(i,...) { MyBS() } )
    ##stop the cluster
    stopCluster(cl)
    return(myResult)
}

before <- Sys.time()
before
TestAlice <- MyParallel(1000)
after <-  Sys.time()
after
after - before

## Check if param_exp belongs to [parammin-parammax] for H3: 
MyCoverageFunction <- function(N_BS, myResult){     ## N_BS <- number of repeats
    myMins <- c("k_min", "alpha_min", "m1_min", "Dm1_min", "m2_min", "Dm2_min")
    myMaxs <- c("k_max", "alpha_max", "m1_max", "Dm1_max", "m2_max", "Dm2_max")
    myExps <- c(k_exp, alpha_exp, mF_exp, DmF_exp, mM_exp, DmM_exp)
    ## Check if param_exp belongs to [parammin-parammax] for H3:
    cov <- 0
    for (j in 1:length(myExps)){
        A <- vector()
        for (i in 1:N_BS){
            A <- c(A,
                   data.table::between(myExps[j],
                                       as.numeric(as.character(myResult[which(rownames(myResult)==myMins[j]),][[i]][4])),
                                       as.numeric(as.character(myResult[which(rownames(myResult)==myMaxs[j]),][[i]][4]))))
            cov[j] <- prop.table(table(A))[2]*100
        }
    }
    names(cov) <- c("k", "alpha", "m1", "Dm1", "m2", "Dm2")
    return(cov)
}
