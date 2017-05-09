#########################################
### Translation of code from Stuart Baird
## http://onlinelibrary.wiley.com/doi/10.1111/j.1558-5646.2012.01633.x/pdf
## WHERE ARE THE WORMY MICE? A REEXAMINATION OF HYBRID PARASITISM IN THE
## EUROPEAN HOUSE MOUSE HYBRID ZONE. Stuart J. E. Baird, Alexis Ribas,
## Milos Macholan, Tomas Albrecht, Jaroslav Pialek, Joelle Gouy de Bellocq
rm(list=ls()) 
library(ggplot2)
library(reshape)
library(gridExtra)
library(tableHTML)
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
TestData <- list(males = SimulatedData(500, 20, 10, 1, 2), 
                 females = SimulatedData(500, 10, 30, 1, 2))
TestData[[1]][1,] <- c(0.5, 31) ## Stick to Stuart example
testvec <- c(TestData$males$IndLoads, TestData$females$IndLoads)

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
## testRun("H0") ## ok

## Test ## ******************
## 4 Test ## ******************
    ## data= TestData; PDF4cMeanLoad= PDFNegBin
    ## method= "L-BFGS-B";
    ## k = 2; kmin = 0; kmax = 8;
    ## Alpha = 0; Alphamin = -5; Alphamax = 5;
    ## Meanload = 2; Meanloadmin = 0; Meanloadmax = 50;
    ## Dm = 2; Dmmin = 0; Dmmax = 10

## 4. Function to run the optimisations over the 4 hypotheses :
RunOptim <- function(data, PDF4cMeanLoad, method,
                     k, kmin, kmax,
                     Alpha, Alphamin, Alphamax,
                     Meanload, Meanloadmin, Meanloadmax,
                     Dm, Dmmin, Dmmax) {
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
        assign(paste0("df_H",i), length(get(paste0("MyMLE_H",i))$par))
        assign(paste0("Ddf_H",i), get(paste0("df_H",i)) - df_H0)
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
                                          get(paste0("MyLow_H",i)),
                                          get(paste0("MyparamMLE_H",i)),
                                          get(paste0("MyUp_H",i)),
                                          get(paste0("Ddf_H",i))))
    }
    ## Storage:
    Tot <- list(H0 = Result_H0, H1 = Result_H1, H2 = Result_H2, H3 = Result_H3)
    ## Loop over hypotheses:
    for (i in 0:3){
        HYP <- Tot[[paste0("H",i)]]
        ## Store maximum likelihood value:
        df0 <- data.frame(unlist(Tot[[paste0("H",i)]][1]))
        names(df0) <- "MLE"
        ## Store delta degree of freedom:
        df1 <- data.frame(unlist(Tot[[paste0("H",i)]][5]))
        names(df1) <- "Ddf(H0)"
        ## Store parameters:
        df2 <- data.frame(HYP[-c(1,5)])
        df2 <- t(df2)
        row.names(df2) <- c("min","est", "max")
        df2 <- as.data.frame(df2)
        df2$names <- row.names(df2)
        df2 <- melt(df2, id.vars = "names")
        df2$param <- paste(df2$variable, df2$names, sep = "_")
        df2 <- df2[c(4,3)]
        df2 <- as.data.frame(t(df2))
        Names <- NULL
        for (j in 1:ncol(df2)){
            Names <- c(Names,(as.character(df2[1,j])))
        }
        colnames(df2) <- Names
        df2 <- df2[-1,]
        ## Collapse the dataframes:
        df <- merge(merge(df0,df1, all=TRUE), df2, all= TRUE)
        Hyp <- paste0("H", i)
        df <- merge(Hyp, df, all=TRUE)
        names(df)[1] <- "hyp"
        ## Give name of the hypothese:
        assign(paste0("HYP",i), df)
    }
    Summary <- merge(merge(HYP0,HYP1, all=TRUE), merge(HYP2,HYP3, all=TRUE), all= TRUE)
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

## Test ## ******************
Summary <- RunOptim(TestData, PDFNegBin, "L-BFGS-B",
                      k = 2, kmin = 0, kmax = 8,
                      Alpha = 0, Alphamin = -5, Alphamax = 5,
                      Meanload = 2, Meanloadmin = 0, Meanloadmax = 50,
                      Dm = 2, Dmmin = 0, Dmmax = 10)


### create an html table as string
tableHTML(Summary)
###and to export in a file
## write_tableHTML(tableHTML(TESTFinal), file = 'TestFinal.html')
## write.csv(TESTFinal, file = "TESTFinal.csv")
## Summary <- read.csv("TESTFinal.csv")
## okaaay :D
## Test ## ******************

################################################################
## It's plotting time!!
## Plot the mean load model along the hybrid zone, according to the hyp.

ParaPlot <- function(data, PDF4cMeanLoad, method,
                     k, kmin, kmax,
                     Alpha, Alphamin, Alphamax,
                     Meanload, Meanloadmin, Meanloadmax,
                     Dm, Dmmin, Dmmax, Summary){
    ## Define HI
    HI <- seq(0,1,0.01)
    anar <- function(x){
        as.numeric(as.character(x))
    }
    ## Case 1 : group 1/no groups, no Dm (H0, H2)
    case1 <- function(i){
        data.frame(LB=MeanLoad(m1 = anar(Summary$MeanLoad1_min[i]),
                               Dm = 0,
                               alpha = anar(Summary$Alpha_min[i]),
                               HI),
                   UB=MeanLoad(m1 = anar(Summary$MeanLoad1_max[i]),
                               Dm = 0,
                               alpha = anar(Summary$Alpha_max[i]),
                               HI),
                   est=MeanLoad(m1 = anar(Summary$MeanLoad1_est[i]),
                                Dm = 0,
                                alpha = anar(Summary$Alpha_est[i]),
                                HI),
                   HI)
    }
    ## Case 2 : group 1/no group, Dm (H1, H3)
    case2 <- function(i){
        data.frame(LB=MeanLoad(m1 = anar(Summary$MeanLoad1_min[i]),
                               Dm = anar(Summary$Dm1_min[i]),
                               alpha = anar(Summary$Alpha_min[i]),
                               HI),
                   UB=MeanLoad(m1 = anar(Summary$MeanLoad1_max[i]),
                               Dm = anar(Summary$Dm1_max[i]),
                               alpha = anar(Summary$Alpha_max[i]),
                               HI),
                   est=MeanLoad(m1 = anar(Summary$MeanLoad1_est[i]),
                                Dm = anar(Summary$Dm1_max[i]),
                                alpha = anar(Summary$Alpha_est[i]),
                                HI),
                   HI)
    }
    ## Case 3 : group 2, no Dm (H2)
    case3 <- function(i){
        data.frame(LB=MeanLoad(m1 = anar(Summary$MeanLoad2_min[i]),
                               Dm = 0,
                               alpha = anar(Summary$Alpha_min[i]),
                               HI),
                   UB=MeanLoad(m1 = anar(Summary$MeanLoad2_max[i]),
                               Dm = 0,
                               alpha = anar(Summary$Alpha_max[i]),
                               HI),
                   est=MeanLoad(m1 = anar(Summary$MeanLoad2_est[i]),
                                Dm = 0,
                                alpha = anar(Summary$Alpha_est[i]),
                                HI),
                   HI)
    }
    ## Case 4 : group 2, Dm (H3)
    case4 <- function(i){
        data.frame(LB=MeanLoad(m1 = anar(Summary$MeanLoad2_min[i]),
                               Dm = anar(Summary$Dm2_min[i]),
                               alpha = anar(Summary$Alpha_min[i]),
                               HI),
                   UB=MeanLoad(m1 = anar(Summary$MeanLoad2_max[i]),
                               Dm = anar(Summary$Dm2_max[i]),
                               alpha = anar(Summary$Alpha_max[i]),
                               HI),
                   est=MeanLoad(m1 = anar(Summary$MeanLoad2_est[i]),
                                Dm = anar(Summary$Dm2_max[i]),
                                alpha = anar(Summary$Alpha_est[i]),
                                HI),
                   HI)
    }
    ##Formatting data for gg_point
    dataDF_1 <- as.data.frame(data[1])
    dataDF_1$group <- names(data[1])
    names(dataDF_1) <- c("IndHIs", "IndLoads", "group")
    dataDF_2 <- as.data.frame(data[2])
    dataDF_2$group <- names(data[2])
    names(dataDF_2) <- c("IndHIs", "IndLoads", "group")
    dataDF <- rbind(dataDF_1, dataDF_2)
    ## H0:
    df <- case1(1)
    plot0 <- ggplot(df, aes(x=HI, y=est))+
        ggtitle("H0")+
        geom_point(data = dataDF, aes(x=IndHIs, y=IndLoads), color="darkgreen", size=0.5)+
        geom_line()+
        geom_ribbon(aes(x=HI, ymax=UB, ymin=LB), fill="grey", alpha=.3)+
        theme_bw()+
        theme(legend.position="none")
    ## H1:
    df <- case2(2)
    plot1 <- ggplot(df, aes(x=HI, y=est))+
        ggtitle("H1")+
        geom_point(data = dataDF, aes(x=IndHIs, y=IndLoads), color="darkgreen", size=0.5)+
        geom_line()+
        geom_ribbon(aes(x=HI, ymax=UB, ymin=LB), fill="grey", alpha=.3)+
        theme_bw()+
        theme(legend.position="none")
    ## H2:
    df1 <- case1(3)
    df2 <- case3(3)
    df1$group <- "group1"
    df2$group <- "group2"
    df <- rbind(df1, df2)
    plot2 <-ggplot(df, aes(x=HI, y=est, color=group))+
        geom_line()+
        geom_point(data = dataDF, aes(x=IndHIs, y=IndLoads, color=group), size=0.5)+
        geom_ribbon(aes(x=HI, ymax=UB, ymin=LB), alpha=.3)+
        scale_fill_manual(values=c("red", "red", "blue", "blue"))+
        scale_color_manual(values=c("red", "red", "blue", "blue"))+
        ggtitle("H2")+
        theme_bw()+
        theme(legend.position="none")+
        annotate("text", x = .1, y = 100, label = "Group 1", color = "red", size=5)+
        annotate("text", x = .1, y = 75, label = "Group 2", color = "blue", size=5)
    ## H3:
    df1 <- case2(4)
    df2 <- case4(4)
    df1$group <- "group1"
    df2$group <- "group2"
    df <- rbind(df1, df2)
    plot3 <-ggplot(df, aes(x=HI, y=est, color=group))+
        geom_line()+
        geom_point(data = dataDF, aes(x=IndHIs, y=IndLoads), size=0.5)+
        geom_ribbon(aes(x=HI, ymax=UB, ymin=LB), alpha=.3)+
        ggtitle("H3")+
        theme_bw()+
        theme(legend.position="none")+
        scale_fill_manual(values=c("red", "red", "blue", "blue"))+
        scale_color_manual(values=c("red", "red", "blue", "blue"))+
        annotate("text", x = .1, y = 100, label = "Group 1", color = "red", size=5)+
        annotate("text", x = .1, y = 75, label = "Group 2", color = "blue", size=5)
    ## Final plot:
    plotfinal <- arrangeGrob(plot0, plot1, plot2, plot3)
    ## To save the plot:
    ## ggsave(file="FinalPlot.pdf", plotfinal)
    ## To read the plot:
    ## plot(plotfinal)
    return(plotfinal)
}
## Test: 
PLOTFinal <- ParaPlot(TestData, PDFNegBin, "L-BFGS-B",
                      k = 2, kmin = 0, kmax = 8,
                      Alpha = 0, Alphamin = -5, Alphamax = 5,
                      Meanload = 2, Meanloadmin = 0, Meanloadmax = 50,
                      Dm = 2, Dmmin = 0, Dmmax = 10,
                      Summary)

plot(PLOTFinal)




