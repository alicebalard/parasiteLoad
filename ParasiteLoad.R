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
    IndHIs <- runif(Ninds)
    TheMeanLoadModel <- function(HI){
        MeanLoad(m1, Dm, alpha, HI)
    }
    TheMeanLoads <- sapply(IndHIs, TheMeanLoadModel)
    IndLoads  <- rnbinom(Ninds, size=k, mu= TheMeanLoads)
    data.frame(IndHIs, IndLoads)}

TestData <- list(males = SimulatedData(20, 20, 10, 1, 2), 
                 females = SimulatedData(20, 10, 30, 1, 2))
TestData[[1]][1,] <- c(0.5, 31) ## Stick to Stuart example

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
LikelihoodFunction <- function(param, PDF4cMeanLoad, data) {
    k <- param[1]; alpha <- param[2]
    m1M <- param[3]; m1F <- param[4]
    DmM <- param[5]; DmF <- param[6]
    Sm <- 0 ; Sf <- 0 ## initialisation
    ## males
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
LikelihoodFunction <- function(param, PDF4cMeanLoad, data) {
    k <- param[1]; alpha <- param[2]
    m1M <- param[3]; m1F <- param[4]
    DmM <- param[5]; DmF <- param[6]
    Sm <- 0 ; Sf <- 0 ## initialisation
        ## males                       for (ind in (1:nrow(data[[1]]))){
        Sm <- Sm + log(PrLoad(data, "male", ind, k, m1M, DmM, alpha, PDF4cMeanLoad))
    }
    ## females
    for (ind in (1:nrow(data[[2]]))){
        Sf <- Sf + log(PrLoad(data, "female", ind, k, m1F, DmF, alpha, PDF4cMeanLoad))
    }  
    ## all
    S <- Sm + Sf #sum of the logs (same as product of the row values) of males and females
    return(S)
}
   
## Test
LikelihoodFunction(c(2, 0.3, 10, 12, 5, 10), PDFNegBin, TestData)

############### --> turns until there perfectly

## Divide the 4 hypotheses :

## H3 Full model (reduce later on):
LikelihoodFunction <- function(param, PDF4cMeanLoad, data) {
    k <- param[1]; alpha <- param[2]
    m1M <- param[3]; m1F <- param[4]
    DmM <- param[5]; DmF <- param[6]
    Sm <- 0 ; Sf <- 0 ## initialisation
    ## males
    for (ind in (1:nrow(data[[1]]))){
        Sm <- Sm + log(PrLoad(data, "male", ind, k, m1M, DmM, alpha, PDF4cMeanLoad))
    }
    ## females
    for (ind in (1:nrow(data[[2]]))){
        Sf <- Sf + log(PrLoad(data, "female", ind, k, m1F, DmF, alpha, PDF4cMeanLoad))
    }  
    ## all
    S <- Sm + Sf #sum of the logs (same as product of the row values) of males and females
    return(S)
}

## Define the function for which we will optimize the parameters:
OptimLikelihood <- function(param){
    return(LikelihoodFunction(param, PDFNegBin, TestData))}

## Constraints:
lower <- c(kmin=0, AlphaLB=-5, min_m1M=0, min_m1F=0, min_DmM=0, min_DmF=0)
upper <- c(kmax=8, AlphaUB=5, max_m1M=50, max_m1F=50, max_DmM=50, max_DmF=50)

## Starter (likely values) :
start <- c(k=2, Alpha=0, m1M=10, m1F=12, DmM=2, DmF=3)
   
## Optimisation of the parameters:
O <- optim(par = start, ## initial values for the parameters
      fn = OptimLikelihood, ## function to be maximized
      lower = lower, ## lower bounds for the parameters
      upper = upper, ## upper bounds for the parameters
      method = "L-BFGS-B", ## set the method (Method "L-BFGS-B" from Byrd et. al. (1995))
      control = list(fnscale=-1)) ##turn the default minimizer into maximizer

O

## Optimisation of the errors :
## For each parameters (here N=5) find max and min as MLE > O$value -2

## Define the ranges of parameters
k_range <- seq(0, 8, 0.1)
alpha_range <- seq(-5, 5, 0.1)
m1M_range <- seq(0, 50, 0.1)
m1F_range <- seq(0, 50, 0.1)
m2M_range <- seq(0, 50, 0.1)
m2F_range <- seq(0, 50, 0.1)


## Run the likelihood function accross all these values,
## and find the minimum and maximum of k if MLE > MLE_O - 2
LikelihoodFunction(c(2, 0.3, 10, 12, 5, 10), PDFNegBin, TestData)


combn(letters[1:4], 2)

## expand.grid(k_range, alpha_range, m1M_range, m1F_range, m2M_range, m2F_range)

sapply(k_range, LikelihoodFunction, arg1=c(0.3, 10, 12, 5, 10, PDFNegBin, TestData))


mapply(LikelihoodFunction, k_range, alpha_range, m1M_range, m1F_range, m2M_range, m2F_range, PDFNegBin, TestData)
                                          

a <- 1:5
b <- 2:6
combn(a,b,1)

start[1]

Minimum <- function(MLE){
    myoptim <- optim(par = start, ## initial values for the parameters
    fn = OptimLikelihood, ## function to be maximized
    lower = lower, ## lower bounds for the parameters
#   upper = upper, ## upper bounds for the parameters
    method = "L-BFGS-B", ## set the method (Method "L-BFGS-B" from Byrd et. al. (1995))
    control = list(fnscale=-1)) ##turn the default minimizer into maximizer
    return(myoptim$par[1])

    

    return(param)
    

optim(par = O$value - 2, ## initial values for the parameters
      fn = , ## function to be maximized
      lower = lower[1], ## lower bounds for the parameters
#      upper = upper, ## upper bounds for the parameters
      method = "L-BFGS-B", ## set the method (Method "L-BFGS-B" from Byrd et. al. (1995))
      control = list(fnscale=-1)) ##turn the default minimizer into maximizer



## Define the function for which we will find the lower & upper bounds :
OptimBounds <- function(param){
    MLE <- LikelihoodFunction(param, PDFNegBin, TestData)
    return(O$value))
    }

## Lower bound:
optim(par = start, ## initial values for the parameters
      fn = OptimBounds, ## function to be maximized
#      lower = lower, ## lower bounds for the parameters
      upper = upper, ## upper bounds for the parameters
      method = "L-BFGS-B") ## set the method (Method "L-BFGS-B" from Byrd et. al. (1995))

## Upper bound:
optim(par = start, ## initial values for the parameters
      fn = OptimBounds, ## function to be maximized
      lower = lower, ## lower bounds for the parameters
#      upper = upper, ## upper bounds for the parameters
      method = "L-BFGS-B", ## set the method (Method "L-BFGS-B" from Byrd et. al. (1995))
      control = list(fnscale=-1)) ##turn the default minimizer into maximizer












######################################################################

## Test if the difference between 2 likelihood is significant
Gtest <- function(dLL, dDF){
    1 - CDF[ChiSquareDistribution[dDF], 2*dLL]
}

##
GtestOnNestedMLEs <- function(MLEformat1, MLEformat2){
    dLL <- abs(MLEformat1[[1]] - MLEformat2[[1]]) ## format in R
    dDF <- abs(length(MLEformat1[[2]]] - Length[MLEformat2[[2]]]) ## format in R
    p <- Gtest(dLL, dDF)
    print(paste("dLL = ", dLL, " dDF = ", dDF, " Gtest p = ", p))
}

        
### Equivalent of FindOptim fucntion in R is optim
## http://stat.ethz.ch/R-manual/R-devel/library/stats/html/optim.html
