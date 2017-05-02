########################## On reprend... 29 april
rm(list=ls()) 

library(ggplot2)
library(reshape)
library(gridExtra)
library(tableHTML)

#########################    
## Simulate data for 1 independant variable (male/female)
## NB: add later "age", and generalise...

SimulatedData <- function(NindsF, NindsM, mF, mM, DmF, DmM, alpha, k){
    ##    set.seed(5)
    IndHIsF <- round(runif(NindsF), 2)
    IndHIsM<- round(runif(NindsM), 2)
    FemLoads  <- data.frame(HI = IndHIsF, loads= rnbinom(NindsF, size=k, mu= MeanLoad(mF, DmF, alpha, IndHIsF)), group1 = "female")
    MalLoads  <- data.frame(HI = IndHIsM, loads= rnbinom(NindsM, size=k, mu= MeanLoad(mM, DmM, alpha, IndHIsM)), group1 = "male")
    ## Bind the dataframes:
    return(rbind(FemLoads, MalLoads))
}

## example:
## 1.Choose our parameters for the simulation:
     NindsF_exp <- 80; NindsM_exp <- 85; mF_exp <- 15; mM_exp <- 22
     DmF_exp <- 5; DmM_exp <- 0; alpha_exp <- 1.1; k_exp <- 4
## And simulate data:
alicedata <- SimulatedData(NindsF_exp, NindsM_exp, mF_exp, mM_exp,
                              DmF_exp, DmM_exp, alpha_exp, k_exp)

#####################
### MeanLoad model### How does the mean load vary depending on the parameters?
#####################
MeanLoad <- function(m,Dm,alpha,HI){
    (m + Dm*HI)*(1 - alpha*2*HI*(1 - HI))
}

###########################################################
### PrNegativeBinomialLoad: Probability an Ind has some ###
### observed Load, given its HI (and the model)         ###
###########################################################

## Generalised, ie the MeanLoad model is passed as an argument
PrLoad <- function(data, k, m, Dm, alpha, ind){
    HI <- data$HI[ind]
    Load <- data$loads[ind]
    ## Always a good idea to make contributions to the likelihood function Bombproof, and never return zero (???)
    max(10^-20, dnbinom(Load, size=abs(k), mu=abs(MeanLoad(m, Dm, alpha, HI))))
}

## example:
## PrLoad(alicedata, k=2, m=10, Dm=2, alpha=1, ind=2)

##################################################
### The likelihood function over a set of inds ###
##################################################

## 'data' is a dataframe with HI, loads, group1 (so far...)
## PDF4cMeanLoad will be a chosen MeanLoad model function (e.g. PDFNegBin)
## the first argument will be optimised with "optim" later (k, alpha)

## All hypotheses (so far...):
## H0: k, alpha, mB
## H1: k, alpha, mB, DmB
## H2: k, alpha, m1M, m1F
## H3: k, alpha, m1M, m1F, DmM, DmF

## parametre given for optimisation: k, alpha, m, Dm
LikelihoodFunction <- function(param, data, hyp) {
    S <- 0 ## initialisation
    if (hyp== "H0"){ ## H0: k, alpha, m
        k <- param[1]; alpha <- param[2];m <- param[3]; Dm <- 0
        ## For all levels of group1:
        G1_levels <- levels(data$group1)
        for (lev in 1:length(G1_levels)){
            ## Subset the dataframe for each level:
            subdata <- data[which(data$group1 %in% G1_levels[lev]),]
            ## For all individuals in each level:
            for (ind in 1:nrow(subdata)){
                S <- S + log(PrLoad(subdata, k, m, Dm, alpha, ind))
            }
        }
    } else if (hyp== "H1"){ ## H1: k, alpha, m, Dm
        k <- param[1]; alpha <- param[2];m <- param[3]; Dm <- param[4]
        ## For all levels of group1:
        G1_levels <- levels(data$group1)
        for (lev in 1:length(G1_levels)){
            ## Subset the dataframe for each level:
            subdata <- data[which(data$group1 %in% G1_levels[lev]),]
            ## For all individuals in each level:
            for (ind in 1:nrow(subdata)){
                S <- S + log(PrLoad(subdata, k, m, Dm, alpha, ind))
            }
        }
    } else if (hyp== "H2"){ ## H2: k, alpha, m1, m2
        k <- param[1]; alpha <- param[2]; m1 <- param[3]; m2 <- param[4]
        m <- c(m1,m2) ## group the m of each level of Group1
        ## For all levels of group1:
        G1_levels <- levels(data$group1)
        for (lev in 1:length(G1_levels)){
            ## Subset the dataframe for each level:
            subdata <- data[which(data$group1 %in% G1_levels[lev]),]
            ## For all individuals in each level:
            for (ind in 1:nrow(subdata)){
                S <- S + log(PrLoad(subdata, k, m[lev], 0, alpha, ind))
            }
        }
    } else if (hyp== "H3"){ ## H3: k, alpha, m1, Dm1, m2, Dm2
        k <- param[1]; alpha <- param[2]
        m1 <- param[3]; m2 <- param[4]; Dm1 <- param[5]; Dm2 <- param[6]
        m <- c(m1,m2) ## group the m of each level of Group1
        Dm <- c(Dm1, Dm2) ## same for Dm
        ## For all levels of group1:
        G1_levels <- levels(data$group1)
        for (lev in 1:length(G1_levels)){
            ## Subset the dataframe for each level:
            subdata <- data[which(data$group1 %in% G1_levels[lev]),]
            ## For all individuals in each level:
            for (ind in 1:nrow(subdata)){
                S <- S + log(PrLoad(subdata, k, m[lev], Dm[lev], alpha, ind))
            }
        }
    }
    return(S)
}

## example:
LikelihoodFunction(c(2,1,10), alicedata, "H0")
LikelihoodFunction(c(2,1,10,2), alicedata, "H1")
LikelihoodFunction(c(2,1,10,15), alicedata, "H2")
LikelihoodFunction(c(k=2, alpha=1, m1=10, Dm1=2,
                      m2=15, Dm2=3), alicedata, "H3")
        
#######################################
### The Maximum likelihood analysis ###
#######################################

## 1. Function to maximize the likelihood:
mymle <- function(data, hyp, lower, upper, start, method){
    mymle <- optim(par = start, ## initial values for the parameters
                   fn = LikelihoodFunction, ## function to be maximized
                   lower = lower, ## lower bounds for the parameters
                   upper = upper, ## upper bounds for the parameters
                   method = method, ## set the optimization method
                   control = list(fnscale=-1), ##turn the default minimizer into maximizer
                   data = data, ## extra param for LikelihoodFunction
                   hyp = hyp)## extra param for LikelihoodFunction
}

## 2. Function to maximize the parameters in MLE+/- 2 :

## MLEbounds : optim function with constraint on the search
UpperBoundsFunction <- function(data, hyp, lower, upper, start, method,
                                MyMLE){
    ## empty dataframe to store the data
    summaryVec <- vector()
    max <- length(MyMLE$par)
    for (rankpar in 1:max){ ## run over the parameters
        ## eliminate errors:
        if (!is.na(MyMLE$value)){
            ## Functional constraint on search (L > MLE-2) + max&min each param
            OptimBounds <- function(param){
                LK <- LikelihoodFunction(param, data, hyp)
                if (LK <= !is.na(MyMLE$value - 2)){ ## why !is.na works better? unclear
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
        } else {
            summaryVec <- c(summaryVec, "MLE optimisation error")
        }
    }
    return(summaryVec)
}

## 3. Function to minimize the parameters in MLE+/- 2 :

## MLEbounds : optim function with constraint on the search
LowerBoundsFunction <- function(data, hyp, lower, upper, start, method,
                                MyMLE){
    ## empty dataframe to store the data
    summaryVec <- vector()
    max <- length(MyMLE$par)
    for (rankpar in 1:max){ ## run over the parameters
        ## eliminate errors:
        if (!is.na(MyMLE$value)){
            ## Functional constraint on search (L > MLE-2) + max&min each param
            OptimBounds <- function(param){
                LK <- LikelihoodFunction(param, data, hyp)
                if (LK <= MyMLE$value - 2){
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
        } else {
            summaryVec <- c(summaryVec, "MLE optimisation error")
        }
    }
    return(summaryVec)
}

## 4. Function to run the optimisations over the 4 hypotheses :
RunOptim <- function(data, method,
                     kstart, kmin, kmax,
                     alphastart, alphamin, alphamax,
                     mstart, mmin, mmax,
                     Dmstart, Dmmin, Dmmax) {
    ## Constraints:
    lower <- c(kmin = kmin, alphamin = alphamin, mmin1 = mmin,
               mmin2 = mmin, Dmmin = Dmmin, Dmmin = Dmmin)
    upper <- c(kmax = kmax, alphamax = alphamax, mmax1 = mmax,
               mmax2 = mmax, Dmmax1 = Dmmax, Dmmax2 = Dmmax)
    ## Starter (likely values) :
    start <- c(k = kstart, alpha = alphastart, m1 = mstart,
               m2 = mstart, Dm1 = Dmstart, Dm2 = Dmstart)
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
        assign(paste0("MyMLE_H",i), mymle(data, paste0("H",i),
                                          get(paste0("lower_H",i)),
                                          get(paste0("upper_H",i)),
                                          get(paste0("start_H",i)), method))
        assign(paste0("MLEvalue_H",i), get(paste0("MyMLE_H",i))$value)
        assign(paste0("MyparamMLE_H",i), get(paste0("MyMLE_H",i))$par)
        assign(paste0("df_H",i), length(get(paste0("MyMLE_H",i))$par))
        assign(paste0("Ddf_H",i), get(paste0("df_H",i)) - df_H0)
    }
    ## Run the upperbounds over the 4 hypotheses:
    for (j in 0:3){
        assign(paste0("MyUp_H", j), UpperBoundsFunction(data, 
                                                        paste0("H",j),
                                                        get(paste0("lower_H",j)),
                                                        get(paste0("upper_H",j)),
                                                        get(paste0("start_H",j)),
                                                        method,
                                                        get(paste0("MyMLE_H",j))))
    }
    ## Run the lowerbounds over the 4 hypotheses:
    for (h in 0:3){
        assign(paste0("MyLow_H", h), LowerBoundsFunction(data,
                                                         paste0("H",h),
                                                         get(paste0("lower_H",h)),
                                                         get(paste0("upper_H",h)),
                                                         get(paste0("start_H",h)),
                                                         method,
                                                         get(paste0("MyMLE_H",h))))
    }
    ## Store in a list over the 4 hypotheses :
    for (l in 0:3){
        assign(paste0("Result_H",l), list(get(paste0("MLEvalue_H",l)),
                                          get(paste0("MyLow_H",l)),
                                          get(paste0("MyparamMLE_H",l)),
                                          get(paste0("MyUp_H",l)),
                                          get(paste0("Ddf_H",l))))
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

## example:
## k <- 5; kmin <- 0; kmax <- 8
## alpha <- 0; alphamin <- -2; alphamax <- 2
## m1 <- 0; m1min <- 0; m1max <- 30
## Dm <- 0; Dmmin <- 0; Dmmax <- 1
## RunOptim(alicedata, "L-BFGS-B",
##                         k , kmin, kmax,
##                         alpha, alphamin, alphamax,
##                         m1, m1min, m1max,
##                         Dm, Dmmin, Dmmax)


######## Let's calculate the coverage of the parameters
## (which % of time are the parameters taken to simulate data
## found in the interval estimated)
MyBS <- function(){
    ## 1.Choose our parameters for the simulation:
    NindsF_exp <- 80; NindsM_exp <- 85; mF_exp <- 20; mM_exp <- 22
    DmF_exp <- 5; DmM_exp <- 1; alpha_exp <- 1; k_exp <- 3
    ## And simulate data:
    alicedata <- SimulatedData(NindsF_exp, NindsM_exp, mF_exp, mM_exp,
                               DmF_exp, DmM_exp, alpha_exp, k_exp)
    ## 2.Choose our parameters for start of search of MaxLik (should be different):
    k <- 1; kmin <- 0; kmax <- 8;
    alpha <- 1; alphamin <- -2; alphamax <- 2;
    m1 <- 1; m1min <- 0; m1max <- 30;
    Dm <- 1; Dmmin <- 0; Dmmax <- 10
    ## 3.Run the full optimisation:
    return(RunOptim(alicedata, "L-BFGS-B",
                        k , kmin, kmax,
                        alpha, alphamin, alphamax,
                        m1, m1min, m1max,
                        Dm, Dmmin, Dmmax))
}


################
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
TestAlice <- MyParallel(1000)
after <-  Sys.time()
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

## NB: find the parameters chosen for simulation on top of the code:
MyCoverageFunction(100, myResult)
