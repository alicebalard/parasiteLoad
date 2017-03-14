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


#######################################
### The Maximum likelihood analysis ###
#######################################

#################################################################
## Define the function for which we will optimize the parameters:
## MLE : optim function with box constraints
OptimLikelihood <- function(param){
    return(LikelihoodFunction(param, PDFNegBin, TestData, hyp))}

#####################################################
## Define the functions to find the parameters bounds:
## MLEbounds : optim function with constraint on the search
UpperBoundsFunction <- function(hyp, MyMLE){
    if (hyp == "H0"){
        max <- 3
    } else if (hyp =="H1") {
        max <- 4
    } else if (hyp =="H2") {
        max <- 4
    } else if (hyp =="H3") {
        max <- 6
    }
    ## empty dataframe to store the data
    summaryVec <- vector()
    for (rankpar in 1:max){ ## run over the parameters
        ## Functional constraint on search (L > MLE-2) + max&min each param
        OptimBounds <- function(param){
            LK <- LikelihoodFunction(param, PDFNegBin, TestData, hyp)
            if (LK <= (MyMLE$value - 2)){
                if (which_side== "lower"){
                    param[rankpar] <-100
                } else {
                    param[rankpar] <- -100
                }
            }
            return(param[rankpar])
        }
        ##  Upper bound
        which_side <- "upper"
        MyBoundsMax <- optim(par = MyMLE$par, ## initial values for the parameters
                             fn = OptimBounds, ## function to be maximized
                             lower = lower, ## lower bounds for the parameters
                             upper = upper, ## upper bounds for the parameters
                             method = "L-BFGS-B", ## set the method
                             control = list(fnscale=-1)) ##turn the default minimizer into maximizer
        summaryVec[rankpar] <- MyBoundsMax$par[rankpar]
    }
    return(summaryVec)
}

LowerBoundsFunction <- function(hyp, MyMLE){
    if (hyp == "H0"){
        max <- 3
    } else if (hyp =="H1") {
        max <- 4
    } else if (hyp =="H2") {
        max <- 4
    } else if (hyp =="H3") {
        max <- 6
    }
    ## empty dataframe to store the data
    summaryVec <- vector()
    for (rankpar in 1:max){ ## run over the parameters
        ## Functional constraint on search (L > MLE-2) + max&min each param
        OptimBounds <- function(param){
            LK <- LikelihoodFunction(param, PDFNegBin, TestData, hyp)
            if (LK <= (MyMLE$value - 2)){
                if (which_side== "lower"){
                    param[rankpar] <-100
                } else {
                    param[rankpar] <- -100
                }
            }
            return(param[rankpar])
        }
        ##  Lower bound
        which_side <- "lower"
        MyBoundsMin <- optim(par = MyMLE$par, ## initial values for the parameters
                             fn = OptimBounds, ## function to be maximized
                             lower = lower, ## lower bounds for the parameters
                             upper = upper, ## upper bounds for the parameters
                             method = "L-BFGS-B") ## set the method (Method "L-BFGS-B" from Byrd et. al. (1995))
        summaryVec[rankpar] <- MyBoundsMin$par[rankpar]
    }
    return(summaryVec)
}

###################### H0 ###################### 
hyp <- "H0"

## Constraints:
lower <- c(kmin=0, AlphaLB=-5, min_mB=10)
upper <- c(kmax=8, AlphaUB=5, max_mB=50)
## Starter (likely values) :
start <- c(k=2, Alpha=0, mB=10)

## Optimisation of the parameters:
MyMLE_H0 <- optim(par = start, ## initial values for the parameters
                  fn = OptimLikelihood, ## function to be maximized
                  lower = lower, ## lower bounds for the parameters
                  upper = upper, ## upper bounds for the parameters
                  method = "L-BFGS-B", ## set the method (Method "L-BFGS-B" from Byrd et. al. (1995))
                  control = list(fnscale=-1)) ##turn the default minimizer into maximizer

## Create the dataframe where I will store the results
summaryDF_H0 <- data.frame(k = numeric(),
                           alpha = numeric(),
                           mB = numeric())
summaryDF_H0[2, ] <- MyMLE_H0$par

##store upper parameter
summaryDF_H0[3, ] <- UpperBoundsFunction("H0", MyMLE_H0)

##store lower parameter
summaryDF_H0[1, ] <- LowerBoundsFunction("H0", MyMLE_H0)

## Results:
print("My Maximum Likelihood is")
MyMLE_H0$value
summaryDF_H0

###################### H1 ###################### 
hyp <- "H1"

## Constraints:
lower <- c(kmin=0, AlphaLB=-5, min_mB=10, min_DmB=10)
upper <- c(kmax=8, AlphaUB=5, max_mB=50, max_DmB=50)
## Starter (likely values) :
start <- c(k=2, Alpha=0, mB=10, DmB=2)

## Optimisation of the parameters:
MyMLE_H1 <- optim(par = start, ## initial values for the parameters
                  fn = OptimLikelihood, ## function to be maximized
                  lower = lower, ## lower bounds for the parameters
                  upper = upper, ## upper bounds for the parameters
                  method = "L-BFGS-B", ## set the method (Method "L-BFGS-B" from Byrd et. al. (1995))
                  control = list(fnscale=-1)) ##turn the default minimizer into maximizer

## Create the dataframe where I will store the results
summaryDF_H1 <- data.frame(k = numeric(),
                           alpha = numeric(),
                           mB = numeric(),
                           DmB = numeric())
summaryDF_H1[2, ] <- MyMLE_H1$par

##store upper parameter
summaryDF_H1[3, ] <- UpperBoundsFunction("H1", MyMLE_H1)

##store lower parameter
summaryDF_H1[1, ] <- LowerBoundsFunction("H1", MyMLE_H1)

## Results:
print("My Maximum Likelihood is")
MyMLE_H1$value
summaryDF_H1

###################### H2 ###################### 
hyp <- "H2"

## Constraints:
lower <- c(kmin=0, AlphaLB=-5, min_m1M=10, min_m1F=10)
upper <- c(kmax=8, AlphaUB=5, max_m1M=50, max_m1F=50)
## Starter (likely values) :
start <- c(k=2, Alpha=0, m1M=10, m1F=12)

## Optimisation of the parameters:
MyMLE_H2 <- optim(par = start, ## initial values for the parameters
                  fn = OptimLikelihood, ## function to be maximized
                  lower = lower, ## lower bounds for the parameters
                  upper = upper, ## upper bounds for the parameters
                  method = "L-BFGS-B", ## set the method (Method "L-BFGS-B" from Byrd et. al. (1995))
                  control = list(fnscale=-1)) ##turn the default minimizer into maximizer

## Create the dataframe where I will store the results
summaryDF_H2 <- data.frame(k = numeric(),
                           alpha = numeric(),
                           m1M = numeric(),
                           m1F = numeric())
summaryDF_H2[2, ] <- MyMLE_H2$par

##store upper parameter
summaryDF_H2[3, ] <- UpperBoundsFunction("H2", MyMLE_H2)

##store lower parameter
summaryDF_H2[1, ] <- LowerBoundsFunction("H2", MyMLE_H2)

## Results:
print("My Maximum Likelihood is")
MyMLE_H2$value
summaryDF_H2

###################### H3 ###################### 
hyp <- "H3"

## Constraints:
lower <- c(kmin=0, AlphaLB=-5, min_m1M=10, min_m1F=10, min_DmM=0, min_DmF=0)
upper <- c(kmax=8, AlphaUB=5, max_m1M=50, max_m1F=50, max_DmM=50, max_DmF=50)
## Starter (likely values) :
start <- c(k=2, Alpha=0, m1M=10, m1F=12, DmM=2, DmF=3)

## Optimisation of the parameters:
MyMLE_H3 <- optim(par = start, ## initial values for the parameters
                  fn = OptimLikelihood, ## function to be maximized
                  lower = lower, ## lower bounds for the parameters
                  upper = upper, ## upper bounds for the parameters
                  method = "L-BFGS-B", ## set the method (Method "L-BFGS-B" from Byrd et. al. (1995))
                  control = list(fnscale=-1)) ##turn the default minimizer into maximizer

## Create the dataframe where I will store the results
summaryDF_H3 <- data.frame(k = numeric(),
                           alpha = numeric(),
                           m1M = numeric(),
                           m1F = numeric(),
                           DmM = numeric(),
                           DmF = numeric())
summaryDF_H3[2, ] <- MyMLE_H3$par

##store upper parameter
summaryDF_H3[3, ] <- UpperBoundsFunction("H3", MyMLE_H3)

##store lower parameter
summaryDF_H3[1, ] <- LowerBoundsFunction("H3", MyMLE_H3)

## Results:
print("My Maximum Likelihood is")
MyMLE_H3$value
summaryDF_H3

##############################################################
## Test if the difference between 2 likelihood is significant
Gtest <- function(dLL, dDF){
    1 - pchisq(2*dLL, df=dDF) 
}

GtestOnNestedMLEs <- function(hyp1, hyp2){
    MyMLE1 <- get(paste0("MyMLE_",hyp1))$value
    MyMLE2 <- get(paste0("MyMLE_",hyp2))$value
    summaryDF1 <- get(paste0("summaryDF_",hyp1))
    summaryDF2 <- get(paste0("summaryDF_",hyp2))
    dLL <- round(abs(MyMLE1 - MyMLE2),4)
    dDF <- abs(length(summaryDF1) - length(summaryDF2))
    p <- round(Gtest(dLL, dDF),4)
    print(paste("Likelihood test between", hyp1, "and", hyp2, ":"))
    print(paste("dLL = ", dLL, " dDF = ", dDF, " Gtest p = ", p))
}

## Test the differences of likelihood between each hypotheses and the H0
GtestOnNestedMLEs("H1", "H0")
GtestOnNestedMLEs("H2", "H0")
GtestOnNestedMLEs("H3", "H0")

################################################################
## It's plotting time!!
## Plot the mean load model along the hybrid zone, according to the hyp.

## H0:
HI <- seq(0,1,0.001)
df <- data.frame(HI,
                 LB=MeanLoad(m1=summaryDF_H0[1,3],
                          Dm=0,
                          alpha=summaryDF_H0[1,2],
                          HI),
                 est=MeanLoad(m1=summaryDF_H0[2,3],
                          Dm=0,
                          alpha=summaryDF_H0[2,2],
                          HI),
                 UB=MeanLoad(m1=summaryDF_H0[3,3],
                          Dm=0,
                          alpha=summaryDF_H0[3,2],
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


### Multiplot:
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
    }
}


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
