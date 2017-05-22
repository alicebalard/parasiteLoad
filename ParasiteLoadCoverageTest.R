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
## H2: k, alpha, intercept1*factor1
## H3: k, alpha, intercept1*factor1, 
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
## PrLoad: Probability an Ind has some
## observed Load, given its HI (and the model)        
## Generalised, ie the MeanLoad model is passed as an argument


## Aim for such an interface... 

## INPUT
## alice.ml.nb(bristles ~ latitude + sex + strain, data = drosophila)

## OUTPUT
## alpha                          -   estimate, lower, upper
## k                              -   estimate, lower, upper
## factor1..N * intercept         -   estimate, lower, upper
## factor1..N * growth            -   estimate, lower, upper
##                                    ML estimate

## in your input script use e.g.

contained <- all.vars(HI~ loads + group1*group2)

## same slope but different intercepts
try1 <- lm(loads ~ HI + group1 +  group2, data=fakedata)


## different slopes and different intercepts
try2 <- lm(loads ~ HI * group1 +  HI*group2, data=fakedata)


#########################
## MeanLoad model: How does the mean load vary depending on the parameters?
MeanLoad <- function(intercept, growth, alpha, HI){
    (intercept + growth*HI)*(1 - alpha*2*HI*(1 - HI))
}

## By row
PrLoadbyRow <- function (data, sub.param) { 
    PrLoad <- function(data, k=k, alpha=alpha, intercept=intercept, growth=growth){
        with(data, 
             ## Always a good idea to make contributions to the likelihood
             ## function Bombproof, and never return zero
             max(10^-20, dnbinom(loads, size=abs(k),
                                 mu=abs(MeanLoad(intercept=intercept,
                                                 growth=growth,
                                                 alpha=alpha,
                                                 HI=HI))))
             )
    }
    sapply(1:nrow(data), function (i){
        PrLoad(data=data[i, ],
               k = sub.param[[1]],
               alpha = sub.param[[2]],
               intercept = sub.param[[3]],
               growth = sub.param[[4]])})
}


PrLoadbyRow(fakedata, sub.param=param[[1]])
PrLoadbyRow(fakedata, sub.param=param[[2]])

#########################
## Likelihood of a given subgroup of the data (aka male)

## The likelihood function over a set of inds
## param: k, alpha, mM, mF, slapeM, slapeF
LikelihoodFunction <- function(data, ugly.param) { 
    param <- list()
    param[["male"]] <-  ugly.param[1:4]
    param[["female"]] <-  ugly.param[c(1,2,5,6)]
    ## We consider so far 1 independant variable: sex, with 2 levels (TO IMPROVE)
    ## levels(group1)
    split.L<- by(data, data$group1, function(x)  {
        PrLoadbyRow(x, param[[unique(x$group1)]])}) ### now hardcoded for group1
##    print(split.L)
    sum(log(unlist(split.L)))
}

LikelihoodFunction(fakedata,
                   ugly.param=unlist(param)[c(1:4, 7, 8)])


optim(par = unlist(param)[c(1:4, 7, 8)],
      fn = LikelihoodFunction, ## function to be maximized
##      upper=c(male=list(1, 1, 1, 1),
##              female=list(1, 1, 1, 1)),
##      lower=c(male=list(-0.1, -0.1, 0.001, 0.0001),
##              female=list(-0.1, -0.1, 0.001, 0.0001)),
      control = list(fnscale=-1), ##turn the default minimizer into
                                  ##maximizer
##      method = "L-BFGS-B",
      data = fakedata)


#########################
## The Maximum likelihood analysis
## 1. Function to maximize the likelihood:
mymle <- function(data, fn, lower, upper, start){
    ## Set limits to the lower and upper parameters
    ## in the case where they must stay NULL (or close to):
    ##    lower <- lower*multi - 0.1 
    ##    upper <- upper*multi + 0.1
    ##    start <- start*multi
    ## And run the optimisation:    
    optim(par = start, ## initial values for the parameters
          fn = fn, ## function to be maximized
          lower = lower, ## lower bounds for the parameters
          upper = upper, ## upper bounds for the parameters
          control = list(fnscale=-1), ##turn the default minimizer into maximizer
          method = "L-BFGS-B",
          data = data)
}

## 2. Function to maximize the parameters in MLE+/- 2 :

## MLEbounds : optim function with constraint on the search
UpperBoundsFunction <- function(data, lower, upper, start, MyMLE, multi){
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
