#########################################
### Translation of code from Stuart Baird
## http://onlinelibrary.wiley.com/doi/10.1111/j.1558-5646.2012.01633.x/pdf
## WHERE ARE THE WORMY MICE? A REEXAMINATION OF HYBRID PARASITISM IN THE
## EUROPEAN HOUSE MOUSE HYBRID ZONE. Stuart J. E. Baird, Alexis Ribas,
## Milos Macholan, Tomas Albrecht, Jaroslav Pialek, Joelle Gouy de Bellocq
rm(list=ls()) 
library(ggplot2)
library(reshape)

#####################
### MeanLoad model### How does the mean load vary depending on the parameters?
#####################

## We define a function where m1 is the load of Mmmusculus,
## Dm the additional load of Mmdomesticus
## alpha the hybridization effect (also called V in the article,
## We will investigate it later through MaxLikelihood)
## HI the hybrid index ( He = 2HI(1-HI) )
MeanLoad <- function(m1,Dm,alpha,HI){
    (m1 + Dm*HI)*(1 - alpha*2*HI*(1 - HI))
}

## First numerical example
HI <- seq(from = 0, to = 1, by = 0.025)
ALPHA <- c(-0.4, -0.2, 0, 0.2, 0.4)
m <- as.data.frame(matrix(0, ncol = 5, nrow = 41))
i <- 1
for (alpha in c(-0.4, -0.2, 0, 0.2, 0.4)){
    m[i] <- MeanLoad(40,20,alpha,HI)
    i <- i+1
}

names(m) <- as.character(ALPHA)
m <- melt(m)
names(m) <- c("alpha", "load")
m$HI <- rep(HI,5)

ggplot(m, aes(x=HI, y=load, group=alpha, color=alpha))+
    geom_line()+
    ggtitle("MeanLoad model (alpha coloured)")+
    annotate("text", x = 0.05, y = 40, label = "m1")+
    annotate("text", x = 0.95, y = 40+20, label = "m1+Dm")

######################################################
### NegativeBinomial parameterisation for MeanLoad ###
######################################################
## Probability of getting a mouse with a certain load, knowing:
## meanload : MeanLoad
## aggregation parameter : k   We will investigate it later through MaxLik

## For understanding: E(X)=k*(1-prob/prob) for negbino
## E(X)=MeanLoad <=>k*(1-prob/prob)=MeanLoad <=>prob=k/(MeanLoad+k)

## First, numerical application to check that the MeanLoad model
## has "the right" mean, being the input MeanLoad
myvec <- vector()
for (k in 1:3){
    for (meanload in seq(from=20, to=40, by=10)){
        myvec <- c(myvec,(mean(rnbinom(10000, size=k, mu= meanload))))
        ## (n=10000 number of observation)
    }
}
t(matrix(myvec,ncol=3))

## Secondly, let's have a look on the plot for different values:
## Notes : q is a quantile vector, meaning the values of probabilities
## for which we want to get an x. We set here 0 to 0.1 with a 0.001 step.
q <- seq(from=0,to=0.1, by=0.001)
myvec <- vector()
mynames <- vector()
for (x in 0:200){
    for (k in 1:3){
        for (meanload in seq(from=20, to=40, by=10)){
            myvec <- c(myvec,dnbinom(x, size=k, mu=meanload)) #alternative:prob=c/k+MeanLoad
            mynames <- c(mynames,paste(k, meanload, sep=", "))
        }
    }
}

x <- vector()
for (i in 0:200) {
    x <- c(x,rep(i,9))
}

mydf <- data.frame(myvec,mynames,x)

ggplot(mydf, aes(x=x, y=myvec, group=mynames, color=mynames))+
    geom_line()+
    scale_x_discrete(name="Number of worms found in 1 individual")+
    scale_y_continuous(name="Probability")
## When k increases, A=1/k decreases, less aggregation, tends towards Poisson

######################
### Simulated data ###
######################
SimulatedData <- function(Ninds, m1, Dm, alpha, k){
    set.seed(5)
    IndHIs <- runif(Ninds)
    TheMeanLoadModel <- function(HI){
        MeanLoad(m1, Dm, alpha, HI)
    }
    TheMeanLoads <- sapply(IndHIs, TheMeanLoadModel)
    IndLoads  <- rnbinom(Ninds, size=k, mu= TheMeanLoads)
    data.frame(IndHIs, IndLoads)}
TestData <- list(males = SimulatedData(50, 20, 10, 1, 2), 
                 females = SimulatedData(50, 10, 30, 1, 2))
TestData[[1]][1,] <- c(0.5, 31) ## Stick to Stuart example
testvec <- c(TestData$males$IndLoads, TestData$females$IndLoads)
summary(testvec)

###########################################################
### PrNegativeBinomialLoad: Probability an Ind has some ###
### observed Load, given its HI (and the model)         ###
###########################################################

## Function to calculate the density of probability for 1 ind,
## given his ID (and his HI and Load associated)
PrNegativeBinomialLoad <- function(data, sex, ind, k, m1, Dm, alpha){
    if (sex=="male"){
        HI <- data[[1]][ind,1]
        Load <- data[[1]][ind,2]
    } else {
        HI <- data[[2]][ind,1]
        Load <- data[[2]][ind,2]
    }
    dnbinom(Load, size=k, mu=MeanLoad(m1, Dm, alpha, HI))
}

## Let's calculate an example natural Log
log(PrNegativeBinomialLoad(TestData, "male", 5, 2, 0.4, 0.2, 0.0))
## Knowing that the mean load of musculus is 0.4, the mean load of domesticus
## is 0.6 (0.4+0.2), k is 2 (aggregation), alpha is 0 (hybrid effect), the load of our
## toy sample is 31 for a HI of 0.5, then the probability of getting this load is
## reaaally low.

## Generalised ... ie the MeanLoad model is passed as an argument (PDF4cMeanLoad)
PrLoad <- function(data, sex, ind, k, m1, Dm, alpha, PDF4cMeanLoad){
    if (sex=="male"){
        HI <- data[[1]][ind,1]
        Load <- data[[1]][ind,2]
    } else {
        HI <- data[[2]][ind,1]
        Load <- data[[2]][ind,2]
    }
    PDF4cMeanLoad(k, MeanLoad(m1, Dm, alpha, HI), Load)
}

SimplePDFNegBin <- function(k, meanLoad, Load){
    dnbinom(Load, size=k, mu=meanLoad)
}
## Always a good idea to make contributions to the likelihood function Bombproof, and never return zero.
PDFNegBin <- function(k, meanLoad, Load){
    max(10^-20, SimplePDFNegBin(abs(k), abs(meanLoad), Load))
}
## This previous function will be our MeanLoad model function, to be tested for arguments
## k and alpha through maximum likelihood. 

## Test: 
PrLoad(TestData, "male", 1, 2, 40, 20, 0.2, PDFNegBin)

##################################################
### The likelihood function over a set of inds ###
##################################################

## 'data' is two sets/lists of inds (e.g. male and female, hence m1M vs m1F)
## PDF4cMeanLoad will be a chosen MeanLoad model function (e.g. PDFNegBin)
## the first argument will be optimised with "optim" later (k, alpha)

## All hypotheses:
## H0: k, alpha, mB
## H1: k, alpha, mB, DmB
## H2: k, alpha, m1M, m1F
## H3: k, alpha, m1M, m1F, DmM, DmF

LikelihoodFunction <- function(param, PDF4cMeanLoad, data, hyp) {
    Sm <- 0 ; Sf <- 0 ## initialisation
    if (hyp== "H0"){ ## H0: k, alpha, mB
        k <- param[1]; alpha <- param[2]
        mB <- param[3]
        ## males
        for (ind in (1:nrow(data[[1]]))){
            Sm <- Sm + log(PrLoad(data, "male", ind, k, mB, 0, alpha, PDF4cMeanLoad))
        }
        ## females
        for (ind in (1:nrow(data[[2]]))){
            Sf <- Sf + log(PrLoad(data, "female", ind, k, mB, 0, alpha, PDF4cMeanLoad))
        }
    } else if (hyp== "H1"){ ## H1: k, alpha, mB, DmB
        k <- param[1]; alpha <- param[2]
        mB <- param[3]; DmB <- param[4]
        ## males
        for (ind in (1:nrow(data[[1]]))){
            Sm <- Sm + log(PrLoad(data, "male", ind, k, mB, DmB, alpha, PDF4cMeanLoad))
        }
        ## females
        for (ind in (1:nrow(data[[2]]))){
            Sf <- Sf + log(PrLoad(data, "female", ind, k, mB, DmB, alpha, PDF4cMeanLoad))
        }
    } else if (hyp== "H2"){ ## H2: k, alpha, m1M, m1F
        k <- param[1]; alpha <- param[2]
        m1M <- param[3]; m1F <- param[4]
        ## males
        for (ind in (1:nrow(data[[1]]))){
            Sm <- Sm + log(PrLoad(data, "male", ind, k, m1M, 0, alpha, PDF4cMeanLoad))
        }
        ## females
        for (ind in (1:nrow(data[[2]]))){
            Sf <- Sf + log(PrLoad(data, "female", ind, k, m1F, 0, alpha, PDF4cMeanLoad))
        }
    } else if (hyp== "H3"){ ## H3: k, alpha, m1M, m1F, DmM, DmF
        k <- param[1]; alpha <- param[2]
        m1M <- param[3]; m1F <- param[4]
        DmM <- param[5]; DmF <- param[6]
        ## males
        for (ind in (1:nrow(data[[1]]))){
            Sm <- Sm + log(PrLoad(data, "male", ind, k, m1M, DmM, alpha, PDF4cMeanLoad))
        }
        ## females
        for (ind in (1:nrow(data[[2]]))){
            Sf <- Sf + log(PrLoad(data, "female", ind, k, m1F, DmF, alpha, PDF4cMeanLoad))
        }
    }
    ## all
    S <- Sm + Sf #sum of the logs (same as product of the row values) of males and females
    return(S)
}

## Test
LikelihoodFunction(c(2, 0.3, 10, 12, 5, 10), PDFNegBin, TestData, "H3")
## -388.5384

#######################################
### The Maximum likelihood analysis ###
#######################################

## 1. Function to maximize the likelihood:
mymle <- function(data, PDF4cMeanLoad, hyp, lower, upper, start, method){
    mymle <- optim(par = start, ## initial values for the parameters
          fn = LikelihoodFunction, ## function to be maximized
          lower = lower, ## lower bounds for the parameters
          upper = upper, ## upper bounds for the parameters
          method = method, ## set the optimization method
          control = list(fnscale=-1), ##turn the default minimizer into maximizer
          PDF4cMeanLoad = PDF4cMeanLoad,
          data = data,
          hyp = hyp)
}

## 2. Function to maximize the parameters in MLE+/- 2 :

## MLEbounds : optim function with constraint on the search
UpperBoundsFunction <- function(data, PDF4cMeanLoad, hyp, lower, upper, start, method,
                                MyMLE){
    ## empty dataframe to store the data
    summaryVec <- vector()
    max <- length(MyMLE$par)
    for (rankpar in 1:max){ ## run over the parameters
        ## Functional constraint on search (L > MLE-2) + max&min each param
        OptimBounds <- function(param){
            LK <- LikelihoodFunction(param, PDF4cMeanLoad, data, hyp)
            if (LK <= (MyMLE$value - 2)){
                param[rankpar] <- -100
            }
            return(param[rankpar])
        }
        MyBoundsMax <- optim(par = MyMLE$par, ## initial values for the parameters
                             fn = OptimBounds, ## function to be maximized
                             lower = lower, ## lower bounds for the parameters
                             upper = upper, ## upper bounds for the parameters
                             method = method, ## set the method
                             control = list(fnscale=-1)) ##turn the default minimizer into maximizer
        summaryVec <- c(summaryVec, MyBoundsMax$par[rankpar])
    }
    return(summaryVec)
}

## 3. Function to minimize the parameters in MLE+/- 2 :

## MLEbounds : optim function with constraint on the search
LowerBoundsFunction <- function(data, PDF4cMeanLoad, hyp, lower, upper, start, method,
                                MyMLE){
    ## empty dataframe to store the data
    summaryVec <- vector()
    max <- length(MyMLE$par)
    for (rankpar in 1:max){ ## run over the parameters
        ## Functional constraint on search (L > MLE-2) + max&min each param
        OptimBounds <- function(param){
            LK <- LikelihoodFunction(param, PDF4cMeanLoad, data, hyp)
            if (LK <= (MyMLE$value - 2)){
                param[rankpar] <- 100
            }
            return(param[rankpar])
        }
        MyBoundsMax <- optim(par = MyMLE$par, ## initial values for the parameters
                             fn = OptimBounds, ## function to be maximized
                             lower = lower, ## lower bounds for the parameters
                             upper = upper, ## upper bounds for the parameters
                             method = method) ## set the method
        summaryVec <- c(summaryVec, MyBoundsMax$par[rankpar])
    }
    return(summaryVec)
}

## Test ## ******************

## ## Constraints:
## lower <- c(kmin = 0, Alphamin = -5, Meanloadmin1 = 0, Meanloadmin2 = 0, Dmmin1 = 0, Dmmin2 = 0)
## upper <- c(kmax = 8, Alphamax = 5, Meanloadmax1 = 50, Meanloadmax2 = 50, Dmmax1 = 10, Dmmax2 = 10)
## ## Starter (likely values) :
## start <- c(k = 2, Alpha = 0, Meanload1 = 10, Meanload2 = 10, Dm1 = 2, Dm2 = 2)

## ## H0:
## lower_H0 <- lower[c(1,2,3)]; upper_H0 <- upper[c(1,2,3)]; start_H0 <- start[c(1,2,3)]
## ## H1:
## lower_H1 <- lower[c(1,2,3,5)]; upper_H1 <- upper[c(1,2,3,5)]; start_H1 <- start[c(1,2,3,5)]
## ## H2:
## lower_H2 <- lower[c(1,2,3,4)]; upper_H2 <- upper[c(1,2,3,4)]; start_H2 <- start[c(1,2,3,4)]
## ## H3:
## lower_H3 <- lower; upper_H3 <- upper; start_H3 <- start

## testRun <- function(hyp){
##     lower <- get(paste0("lower_",hyp))
##     upper <- get(paste0("upper_",hyp))
##     start <- get(paste0("start_",hyp))
##     ## 1
##     MyMLE <- mymle(TestData, PDFNegBin, hyp, lower, upper, start, "L-BFGS-B")
##     MLEvalue <- MyMLE$value
##     MyparamMLE <- MyMLE$par
##     ## 2
##     MyU <- UpperBoundsFunction(TestData, PDFNegBin, hyp, lower, upper, start, method =  "L-BFGS-B", MyMLE)
##     ##3
##     MyL <- LowerBoundsFunction(TestData, PDFNegBin, hyp, lower, upper, start, method =  "L-BFGS-B", MyMLE)
##     return(list(MLEvalue, MyparamMLE, MyU, MyL))
## }

## testRun("H3") ## ok
## testRun("H2") ## ok
## testRun("H1") ## ok
testRun("H0") ## ok

## Test ## ******************

## 4. Function to run the optimisations over the 4 hypotheses :
RunOptim <- function(data, PDF4cMeanLoad, method,
                     k, kmin, kmax,
                     Alpha, Alphamin, Alphamax,
                     Meanload, Meanloadmin, Meanloadmax,
                     Dm, Dmmin, Dmmax) {
    ## Create a summary table for all hypotheses:
    SumTab <- data.frame(hyp=c("H0", "H1", "H2", "H3"),
                         LL = NA, deltaLL_H0 = NA, deltaDF_H0 = NA, p = NA,
                         alpha.min = NA, alpha.MLE = NA, alpha.max = NA,
                         M1m.min = NA, M1m.MLE = NA, M1m.max = NA,
                         F1m.min = NA, F1m.MLE = NA, F1m.max = NA,
                         DmM.min = NA, DmM.MLE = NA, DmM.max = NA,
                         DmF.min = NA, DmF.MLE = NA, DmF.max = NA,
                         k.min = NA, k.MLE = NA, k.max = NA)
    ## Constraints:
    lower <- c(kmin = kmin, Alphamin = Alphamin, MeanLoadmin1 = Meanloadmin,
               MeanLoadmin2 = Meanloadmin, Dmmin = Dmmin, Dmmin = Dmmin)
    upper <- c(kmax = kmax, Alphamax = Alphamax, MeanLoadmax1 = Meanloadmax,
               MeanLoadmax2 = Meanloadmax, Dmmax1 = Dmmax, Dmmax2 = Dmmax)
    ## Starter (likely values) :
    start <- c(k = k, Alpha = Alpha, MeanLoad1 = Meanload,
               MeanLoad2 = Meanload, Dm1 = Dm, Dm2 = Dm)
    ## For each hypotese:
    ## H0:
    lower_H0 <- lower[c(1,2,3)]; upper_H0 <- upper[c(1,2,3)]; start_H0 <- start[c(1,2,3)]
    ## H1:
    lower_H1 <- lower[c(1,2,3,5)]; upper_H1 <- upper[c(1,2,3,5)]; start_H1 <- start[c(1,2,3,5)]
    ## H2:
    lower_H2 <- lower[c(1,2,3,4)]; upper_H2 <- upper[c(1,2,3,4)]; start_H2 <- start[c(1,2,3,4)]
    ## H3:
    lower_H3 <- lower; upper_H3 <- upper; start_H3 <- start
    ## Run the function mymle over the 4 hypotheses:
    for (i in 0:3){
        assign(paste0("MyMLE_H",i), mymle(data, PDF4cMeanLoad, paste0("H",i),
                                          get(paste0("lower_H",i)),
                                          get(paste0("upper_H",i)),
                                          get(paste0("start_H",i)), method))
        assign(paste0("MLEvalue_H",i), get(paste0("MyMLE_H",i))$value)
        assign(paste0("MyparamMLE_H",i), get(paste0("MyMLE_H",i))$par)
    }
    ## Run the upperbounds over the 4 hypotheses:
    for (i in 0:3){
        assign(paste0("MyUp_H", i), UpperBoundsFunction(data, PDF4cMeanLoad,
                                                        paste0("H",i),
                                                        get(paste0("lower_H",i)),
                                                        get(paste0("upper_H",i)),
                                                        get(paste0("start_H",i)),
                                                        method,
                                                        get(paste0("MyMLE_H",i))))
    }
    ## Run the lowerbounds over the 4 hypotheses:
    for (i in 0:3){
        assign(paste0("MyLow_H", i), LowerBoundsFunction(data, PDF4cMeanLoad,
                                                        paste0("H",i),
                                                        get(paste0("lower_H",i)),
                                                        get(paste0("upper_H",i)),
                                                        get(paste0("start_H",i)),
                                                        method,
                                                        get(paste0("MyMLE_H",i))))
    }
    ## Store in a list over the 4 hypotheses :
    for (i in 0:3){
        assign(paste0("Result_H",i), list(get(paste0("MLEvalue_H",i)),
                                   get(paste0("MyparamMLE_H",i)),
                                   get(paste0("MyUp_H",i)),
                                   get(paste0("MyLow_H",i))))
    }
    return(list(H0 = Result_H0, H1 = Result_H1, H2 = Result_H2, H3 = Result_H3))
}

##All_MLE <- list(HO = MyMLE_H0, H1 = MyMLE_H1, H2 = MyMLE_H2, H3 = MyMLE_H3)
## Test ## ******************
TEST <- RunOptim(TestData, PDFNegBin, "L-BFGS-B",
                     k = 2, kmin = 0, kmax = 8,
                     Alpha = 0, Alphamin = -5, Alphamax = 5,
                     Meanload = 2, Meanloadmin = 0, Meanloadmax = 50,
                 Dm = 2, Dmmin = 0, Dmmax = 10)

TEST
## okaaay :D
## Test ## ******************

## Storage and G-tests

for (i in 0:3){
    HYP <- TEST[[paste0("H",i)]]
    df <- data.frame(HYP[-1])
    df <- t(df)
    row.names(df) <- c("min", "max", "est")
    assign(paste0("HYP",i), df)
}

as.data.frame(HYP0, nrows = 1)

reshape(HYP0, direction = "long", varying=list(names(HYP0)[1:2]))

##############################################################
    ## Test if the difference between 2 likelihood is significant
    Gtest <- function(dLL, dDF){
        1 - pchisq(2*dLL, df=dDF) 
    }
    GtestOnNestedMLEs <- function(hyp1, hyp2){
        MyMLE1 <- get(paste0("MyMLE_",hyp1))
        MyMLE2 <- get(paste0("MyMLE_",hyp2))
        dLL <- round(abs(MyMLE1$value - MyMLE2$value),4)
        dDF <- length(MyMLE1$par) - length(MyMLE2$par)
        p <- round(Gtest(dLL, dDF),4)
        Gresults <- c(dLL, dDF, p)
        return(Gresults)
    }




    ## H0:
    hyp <- "H0"    
    ## Optimisation of the parameters:
    MyMLE_H0 <- mymle(data, PDF4cMeanLoad, hyp, lower_H0, upper_H0, start_H0, method)
    SumTab$k.MLE[1] <- MyMLE_H0$par[1]
    SumTab$alpha.MLE[1] <- MyMLE_H0$par[2]
    SumTab$M1m.MLE[1] = SumTab$F1m.MLE[1] = SumTab$DmM.MLE[1] = SumTab$DmF.MLE[1] = MyMLE_H0$par[3]

##store upper parameters
    up <- UpperBoundsFunction(data, PDF4cMeanLoad, hyp, lower_H0, upper_H0, start_H0, method, MyMLE_H0)
    SumTab$k.max[1] <- up[1]
    SumTab$alpha.max[1] <- up[2]
    SumTab$M1m.max[1] = SumTab$F1m.max[1] = SumTab$DmM.max[1] = SumTab$DmF.max[1] = up[3]
    ##store lower parameters
    low <- LowerBoundsFunction(data, PDF4cMeanLoad, hyp, lower_H0, upper_H0, start_H0, method, MyMLE_H0)
    SumTab$k.min[1] <- low[1]
    SumTab$alpha.min[1] <- low[2]
    SumTab$M1m.min[1] = SumTab$F1m.min[1] = SumTab$DmM.min[1] = SumTab$DmF.min[1] = low[3]

    ## H1:
    hyp <- "H1"
    ## Optimisation of the parameters:
    MyMLE_H1 <- mymle(data, PDF4cMeanLoad, hyp, lower_H1, upper_H1, start_H1, method)
    SumTab$k.MLE[2] <- MyMLE_H1$par[1]
    SumTab$alpha.MLE[2] <- MyMLE_H1$par[2]
    SumTab$M1m.MLE[2] = SumTab$F1m.MLE[2] <- MyMLE_H1$par[3]
    SumTab$DmM.MLE[2] = SumTab$DmF.MLE[2] <- MyMLE_H1$par[4]
    ##store upper parameters
    up <- UpperBoundsFunction(data, PDF4cMeanLoad, hyp, lower_H1, upper_H1, start_H1, method, MyMLE_H1)
    SumTab$k.max[2] <- up[1]
    SumTab$alpha.max[2] <- up[2]
    SumTab$M1m.max[2] = SumTab$F1m.max[2] = up[3]
    SumTab$DmM.max[2] = SumTab$DmF.max[2] = up[4]
    ##store lower parameters
    low <- LowerBoundsFunction(data, PDF4cMeanLoad, hyp, lower_H1, upper_H1, start_H1, method, MyMLE_H1)
    SumTab$k.min[2] <- low[1]
    SumTab$alpha.min[2] <- low[2]
    SumTab$M1m.min[2] = SumTab$F1m.min[2] = low[3]
    SumTab$DmM.min[2] = SumTab$DmF.min[2] = low[4]

    ## H2:
    hyp <- "H2"
    ## Optimisation of the parameters:
    MyMLE_H2 <- mymle(data, PDF4cMeanLoad, hyp, lower_H2, upper_H2, start_H2, method)
    SumTab$k.MLE[3] <- MyMLE_H2$par[1]
    SumTab$alpha.MLE[3] <- MyMLE_H2$par[2]
    SumTab$M1m.MLE[3] = SumTab$DmM.MLE[3] = MyMLE_H2$par[3]
    SumTab$F1m.MLE[3] = SumTab$DmF.MLE[3] = MyMLE_H2$par[4]
    ##store upper parameters
    up <- UpperBoundsFunction(data, PDF4cMeanLoad, hyp, lower_H2, upper_H2, start_H2, method, MyMLE_H2)
    SumTab$k.max[3] <- up[1]
    SumTab$alpha.max[3] <- up[2]
    SumTab$M1m.max[3] = SumTab$DmM.max[3] = up[3]
    SumTab$F1m.max[3] = SumTab$DmF.max[3] = up[4]
    ##store lower parameters
    low <- LowerBoundsFunction(data, PDF4cMeanLoad, hyp, lower_H2, upper_H2, start_H2, method, MyMLE_H2)
    SumTab$k.min[3] <- low[1]
    SumTab$alpha.min[3] <- low[2]
    SumTab$M1m.min[3] = SumTab$DmM.min[3] = low[3]
    SumTab$F1m.min[3] = SumTab$DmF.min[3] = low[4]

    ## H3:
    hyp <- "H3"
    ## Optimisation of the parameters:
    MyMLE_H3 <- mymle(data, PDF4cMeanLoad, hyp, lower_H3, upper_H3, start_H3, method)
    SumTab$k.MLE[4] <- MyMLE_H3$par[1]
    SumTab$alpha.MLE[4] <- MyMLE_H3$par[2]
    SumTab$M1m.MLE[4] <- MyMLE_H3$par[3]
    SumTab$F1m.MLE[4] <- MyMLE_H3$par[4]
    SumTab$DmM.MLE[4] <- MyMLE_H3$par[5]
    SumTab$DmF.MLE[4] <- MyMLE_H3$par[6]
    ##store upper parameters
    up <- UpperBoundsFunction(data, PDF4cMeanLoad, hyp, lower_H3, upper_H3, start_H3, method, MyMLE_H3)
    SumTab$k.max[4] <- up[1]
    SumTab$alpha.max[4] <- up[2]
    SumTab$M1m.max[4] <- up[3]
    SumTab$F1m.max[4] <- up[4]
    SumTab$DmM.max[4] <- up[5]
    SumTab$DmF.max[4] <- up[6]
    ##store lower parameters
    low <- LowerBoundsFunction(data, PDF4cMeanLoad, hyp, lower_H3, upper_H3, start_H3, method, MyMLE_H3)
    SumTab$k.min[4] <- low[1]
    SumTab$alpha.min[4] <- low[2]
    SumTab$M1m.min[4] <- low[3]
    SumTab$F1m.min[4] <- low[4]
    SumTab$DmM.min[4] <- low[5]
    SumTab$DmF.min[4] <- low[6]
    ## Add the max likelihoods
    SumTab$LL[1] <- MyMLE_H0$value
    SumTab$LL[2] <- MyMLE_H1$value
    SumTab$LL[3] <- MyMLE_H2$value
    SumTab$LL[4] <- MyMLE_H3$value

    ##############################################################
    ## Test if the difference between 2 likelihood is significant
    Gtest <- function(dLL, dDF){
        1 - pchisq(2*dLL, df=dDF) 
    }
    GtestOnNestedMLEs <- function(hyp1, hyp2){
        MyMLE1 <- get(paste0("MyMLE_",hyp1))
        MyMLE2 <- get(paste0("MyMLE_",hyp2))
        dLL <- round(abs(MyMLE1$value - MyMLE2$value),4)
        dDF <- length(MyMLE1$par) - length(MyMLE2$par)
        p <- round(Gtest(dLL, dDF),4)
        Gresults <- c(dLL, dDF, p)
        return(Gresults)
    }
    ## Test the differences of likelihood between each hypotheses and the H0
    SumTab$deltaLL_H0[1] <- " "
    SumTab$deltaLL_H0[2] <- MyMLE_H1$value - MyMLE_H0$value
    SumTab$deltaLL_H0[3] <- MyMLE_H2$value - MyMLE_H0$value
    SumTab$deltaLL_H0[4] <- MyMLE_H3$value - MyMLE_H0$value 
    ##
    SumTab$deltaDF_H0[1] <- " "
    SumTab$deltaDF_H0[2] <- round(GtestOnNestedMLEs("H1", "H0")[2],4)
    SumTab$deltaDF_H0[3] <- round(GtestOnNestedMLEs("H2", "H0")[2],4)
    SumTab$deltaDF_H0[4] <- round(GtestOnNestedMLEs("H3", "H0")[2],4)
    ##
    SumTab$p[1] <- " "
    SumTab$p[2] <- round(GtestOnNestedMLEs("H1", "H0")[3],4)
    SumTab$p[3] <- round(GtestOnNestedMLEs("H2", "H0")[3],4)
    SumTab$p[4] <- round(GtestOnNestedMLEs("H3", "H0")[3],4)
    return(SumTab)
}

## Test :
RunOptim(data = TestData, PDF4cMeanLoad = PDFNegBin, method = "L-BFGS-B",
         kstart=2, kmin=0, kmax=8,
         Alphastart=0, Alphamin=-5, Alphamax=5,
         Meanloadstart=10, Meanloadmin=1, Meanloadmax=50,
         Dmstart=2, Dmmin=1, Dmmax=10) 

## Crash. Ideas : small grid with movable steps OR no constraints...
RunOptim(data = TestData, PDF4cMeanLoad = PDFNegBin, method = "BFGS",
         kstart=2, kmin=NA, kmax=NA,
         Alphastart=0, Alphamin=NA, Alphamax=NA,
         Meanloadstart=10, Meanloadmin=NA, Meanloadmax=NA,
         Dmstart=2, Dmmin=NA, Dmmax=NA) 

RunOptim(data = TestData, PDF4cMeanLoad = PDFNegBin, method = "L-BFGS-B",
         kstart=1, kmin=0, kmax=8,
         Alphastart=0, Alphamin=-5, Alphamax=5,
         Meanloadstart=15, Meanloadmin=0, Meanloadmax=85,
         Dmstart=2, Dmmin=0, Dmmax=50) 

## Grid for the start
k <- c(4,5,6)
alpha <- c(-3, 0, 3)
MeanLoad <- c(1, 5, 10)
Dm <- c(1, 5, 10)

gridA <- expand.grid(k, alpha, MeanLoad, Dm)

##Create list to store updated models
mod.list=list()

for (i in 1:nrow(gridA)){
    mod2 <- try(RunOptim(data = TestData, PDF4cMeanLoad = PDFNegBin, method = "L-BFGS-B",
                         gridA[i,1], kmin=0, kmax=8,
                         gridA[i,2], Alphamin=-5, Alphamax=5,
                         gridA[i,3], Meanloadmin=10, Meanloadmax=85,
                         gridA[i,4], Dmmin=10, Dmmax=50), TRUE)
    if(isTRUE(class(mod2)=="try-error")) {next} else {mod.list[i] <- mod2}
}


mod.list





################################################################
## It's plotting time!!
## Plot the mean load model along the hybrid zone, according to the hyp.

## H0:
HI <- seq(0,1,0.001)
df <- data.frame(HI,
                 LB=MeanLoad(m1=SumTab$M1m.min[1],
                          Dm=0,
                          alpha=SumTab$alpha.min[1],
                          HI),
                 est=MeanLoad(m1=SumTab$M1m.MLE[1],
                          Dm=0,
                          alpha=SumTab$alpha.MLE[1],
                          HI),
                 UB=MeanLoad(m1=SumTab$M1m.max[1],
                          Dm=0,
                          alpha=SumTab$alpha.max[1],
                          HI))
plot0 <- ggplot(df, aes(x=HI, y=est))+
    geom_line()+
    geom_ribbon(aes(x=HI, ymax=UB, ymin=LB), fill="grey", alpha=.5)+
    ggtitle("H0")+
    theme_bw()
plot0

## H1:
df <- data.frame(HI,
                 LB=MeanLoad(m1=summaryDF_H1[1,3],
                          Dm=summaryDF_H1[1,4],
                          alpha=summaryDF_H1[1,2],
                          HI),
                 est=MeanLoad(m1=summaryDF_H1[2,3],
                          Dm=summaryDF_H1[2,4],                          
                          alpha=summaryDF_H1[2,2],
                          HI),
                 UB=MeanLoad(m1=summaryDF_H1[3,3],
                          Dm=summaryDF_H1[3,4],                         
                          alpha=summaryDF_H1[3,2],
                          HI))
plot1 <- ggplot(df, aes(x=HI, y=est))+
    geom_line()+
    geom_ribbon(aes(x=HI, ymax=UB, ymin=LB), fill="grey", alpha=.5)+
    ggtitle("H1")+
    theme_bw()
plot1

## H2:
df <- data.frame(HI=HI,
                 LBF=MeanLoad(m1=summaryDF_H2[1,4],
                          Dm=0,
                          alpha=summaryDF_H2[1,2],
                          HI),
                 UBF=MeanLoad(m1=summaryDF_H2[3,4],
                          Dm=0,                          
                          alpha=summaryDF_H2[3,2],
                          HI),
                 estF=MeanLoad(m1=summaryDF_H2[2,4],
                          Dm=0,                         
                          alpha=summaryDF_H2[2,2],
                          HI),                 
                 LBM=MeanLoad(m1=summaryDF_H2[1,3],
                          Dm=0,
                          alpha=summaryDF_H2[1,2],
                          HI),
                 UBM=MeanLoad(m1=summaryDF_H2[3,3],
                          Dm=0,                          
                          alpha=summaryDF_H2[3,2],
                          HI),
                 estM=MeanLoad(m1=summaryDF_H2[2,3],
                          Dm=0,                         
                          alpha=summaryDF_H2[2,2],
                          HI)
                 )
plot2 <-ggplot(df, aes(x=HI, y=estM))+
    geom_line(color="darkblue")+
    geom_line(aes(x=HI, y=estF), color="red")+
    geom_ribbon(aes(x=HI, ymax=UBF, ymin=LBF), fill="pink", alpha=.5)+
    geom_ribbon(aes(x=HI, ymax=UBM, ymin=LBM), fill="blue", alpha=.5)+
    ggtitle("H2")+
    theme_bw()
plot2

## H3:
summaryDF_H3

df <- data.frame(HI=HI,
                 LBF=MeanLoad(m1=summaryDF_H3[1,4],
                          Dm=summaryDF_H3[1,6],
                          alpha=summaryDF_H3[1,2],
                          HI),
                 UBF=MeanLoad(m1=summaryDF_H3[3,4],
                          Dm=summaryDF_H3[3,6],           
                          alpha=summaryDF_H3[3,2],
                          HI),
                 estF=MeanLoad(m1=summaryDF_H3[2,4],
                          Dm=summaryDF_H3[2,6],                         
                          alpha=summaryDF_H3[2,2],
                          HI),                 
                 LBM=MeanLoad(m1=summaryDF_H3[1,3],
                          Dm=summaryDF_H3[1,5],
                          alpha=summaryDF_H3[1,2],
                          HI),
                 UBM=MeanLoad(m1=summaryDF_H3[3,3],
                          Dm=summaryDF_H3[3,5],            
                          alpha=summaryDF_H3[3,2],
                          HI),
                 estM=MeanLoad(m1=summaryDF_H3[2,3],
                          Dm=summaryDF_H3[2,5],           
                          alpha=summaryDF_H3[2,2],
                          HI)
                 )
plot3 <-ggplot(df, aes(x=HI, y=estM))+
    geom_line(color="darkblue")+
    geom_line(aes(x=HI, y=estF), color="red")+
    geom_ribbon(aes(x=HI, ymax=UBF, ymin=LBF), fill="pink", alpha=.5)+
    geom_ribbon(aes(x=HI, ymax=UBM, ymin=LBM), fill="blue", alpha=.5)+
    ggtitle("H3")+
    theme_bw()
plot3

###########################################
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }
 if (numPlots==1) {
    print(plots[[1]])
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }}


multiplot(plot0,plot1,plot2,plot3, cols=2)
