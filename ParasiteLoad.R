##########################################
### NOT gitted yet (private or public?)###
##########################################
### Translation of code from Stuart Baird
## http://onlinelibrary.wiley.com/doi/10.1111/j.1558-5646.2012.01633.x/pdf
## WHERE ARE THE WORMY MICE? A REEXAMINATION OF HYBRID PARASITISM IN THE
## EUROPEAN HOUSE MOUSE HYBRID ZONE. Stuart J. E. Baird, Alexis Ribas,
## Milos Macholan, Tomas Albrecht, Jaroslav Pialek, Joelle Gouy de Bellocq
rm(list=ls()) 
library(ggplot2)
library(reshape)

#############################
## How to format your data?##
#############################

############# Toy data
males <- list(HIS=c(0.5, 0.2),LOADS=c(31,35),INDS =c("A","B"), SEX=rep("male",2))
females <- list(HIS=c(0.3,0.5,0.4),LOADS=c(20,10,80),INDS =c("C","D","E"), SEX=rep("female",3)) 
keys <- unique(c(names(males), names(females)))
data <- data.frame(setNames(mapply(c, males[keys], females[keys]), keys))
data$HIS <- as.numeric(as.character(data$HIS))
data$LOADS <- as.numeric(as.character(data$LOADS))
#############

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
    for (MeanLoad in seq(from=20, to=40, by=10)){
        myvec <- c(myvec,(mean(rnbinom(10000, size=k, mu= MeanLoad))))
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
        for (MeanLoad in seq(from=20, to=40, by=10)){
            myvec <- c(myvec,dnbinom(x, size=k, mu=MeanLoad)) #alternative:prob=c/k+MeanLoad
            mynames <- c(mynames,paste(k, MeanLoad, sep=", "))
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

###########################################################
### PrNegativeBinomialLoad: Probability an Ind has some ###
### observed Load, given its HI (and the model)         ###
###########################################################

## As we stated before (need to reput it there) :
MeanLoad <- function(m1,Dm,alpha,HI){
    (m1 + Dm*HI)*(1 - alpha*2*HI*(1 - HI))
}
## Function to calculate the density of probability for 1 ind,
## given his ID (and his HI and Load associated)
PrNegativeBinomialLoad <- function(ind, k, m1, Dm, alpha){
    HI <- data$HIS[ind]
    Load <- data$LOADS[ind]
    dnbinom(Load, size=k, mu=MeanLoad(m1, Dm, alpha, HI))
}

## Let's calculate an example natural Log
log(PrNegativeBinomialLoad(1, 2, 0.4, 0.2, 0.0))
## Knowing that the mean load of musculus is 0.4, the mean load of domesticus
## is 0.6 (0.4+0.2), k is 2 (aggregation), alpha is 0 (hybrid effect), the load of our
## toy sample is 31 for a HI of 0.5, then the probability of getting this load is
## reaaally low.

## Generalised ... ie the MeanLoad model is passed as an argument (PDF4cMeanLoad)
PrLoad <- function(ind, k, m1, Dm, alpha, PDF4cMeanLoad){
    HI <- data$HIS[ind]
    Load <- data$LOADS[ind]
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

##################################################
### The likelihood function over a set of inds ###
##################################################
MeanLoad <- function(m1,Dm,alpha,HI){
    (m1 + Dm*HI)*(1 - alpha*2*HI*(1 - HI))
}

## 'data' is two sets/lists of inds (eg male and female, hence m1M vs m1F)
## PDF4cMeanLoad will be a chosen MeanLoad model function (e.g PDFNegBin)
LikelihoodFunction <- function(PDF4cMeanLoad, k, m1M, m1F, DmM, DmF, alpha, data) {
    Sm <- 0 ; Sf <- 0 ## initialisation
    for (ind in 1:nrow(data)) {
        if (data$SEX[ind] == "male") {
            Sm <- Sm + log(PrLoad(ind,k, m1M,DmM, alpha,PDF4cMeanLoad))
        } else if ( data$SEX[ind] == "female") {
            Sf <- Sf + log(PrLoad(ind,k, m1F,DmF, alpha,PDF4cMeanLoad))        
            print(data$INDS[ind])
        } else
            print("error")
    }
    S <- Sm + Sf #sum of the logs (same as product of the row values) of males and females
    return(S)
}

## Test
LikelihoodFunction(PDFNegBin, 2, 10, 12, 5, 10, 3, data)
data
############## Turns so far :)

#############################
###The likelihood analysis###
#############################
## CDF = cumulative distribution function

## Test if the difference between 2 likelihood is significant
Gtest <- function(dLL, dDF){
    1 - CDF[ChiSquareDistribution[dDF], 2*dLL]
}

## A wrapper to clean optima search format:
MLEformat <- function(l){
                                        # what the hell is "l"
    l[1]

"MLEformat[l_] := {l[[1]], Table[l[[2, i, 2]], {i, 1, Length[l[[2]]]}]}"}


##
GtestOnNestedMLEs <- function(MLEformat1, MLEformat2){
    dLL <- abs(MLEformat1[[1]] - MLEformat2[[1]]) ## format in R
    dDF <- abs(length(MLEformat1[[2]]] - Length[MLEformat2[[2]]]) ## format in R
    p <- Gtest(dLL, dDF)
    print(paste("dLL = ", dLL, " dDF = ", dDF, " Gtest p = ", p))
}


## Let's go go go!
LikelihoodAnalysis <- function(PDF4cMeanLoad, Data, HIpos, LOADpos, start, ## HIpos LOADpos not used
                               AlphaLB , AlphaUB , ALconstrainedToStart , ALstarts ,
                               verboseLikelihoodAnalysis){
    kMax <- 10 ## init maximum aggregation parameter (higher would be a Poisson)
    mMax <- 200 ## init maximum parasite load m (to check??or increase of load??) 

    ## vPrint will print stuff if we decided to be verbose
    vPrint <- function(s){
        if verboseLikelihoodAnalysis==TRUE
        print(s)
    }

    ## Core function :
    MLEest = MLEformat[
        FindMaximum[{Apply[LikelihoodFunction,
                           Join[{PDF4cMeanLoad}, psFl[[i]], {Data, HIpos, LOADpos}]],
                                cons[[i]]}, Transpose[{pars[[i]], starter[MLEest, i]}]]]

    LikelihoodFunction(PDF4cMeanLoad, k, m1M, m1F, DmM, DmF, alpha, data) {
}

######################
### Simulated data ###
######################
    SimulatedData <- function(Ninds, m1, Dm, alpha, k){
        IndHIs <- runif(Ninds)
        TheMeanLoadModel <- function(HI){
            MeanLoad(m1, Dm, alpha, HI)
        }
        TheMeanLoads <- sapply(IndHIs, TheMeanLoadModel)
        IndLoads  <- rnbinom(Ninds, size=k, mu= TheMeanLoads[i])
        data.frame(IndHIs, IndLoads)}

    TestData <- rbind(SimulatedData(500, 20, 10, 1, 2), ## males
        SimulatedData(500, 10, 30, 1, 2)) ## females

    ## NB, in Mathematica TestData[[1, 1]] is in R :
    TestData[1,]

    
    

        
    ## Creation of pars and psFL

## 
    
    ## pars is a list of  --------- 
    pars <- list(c(k, mB, alpha),
                 c(k, mB, DmB, alpha),
                 c(k, m1M, m1F, alpha),
                 c(k, m1M, m1F, DmM, DmF, alpha))
    
k <- numeric()

    ## psFL is a dataframe of -------
    psFl <- data.frame((k, mB, mB,0,0,alpha),(k, mB, mB,DmB, DmB, alpha),
    (k, m1M, m1F, 0, 0, alpha),(k, m1M, m1F, DmM, DmF, alpha))
    

    if(ALconstrainedToStart==TRUE){## to cheeeck what is alconstrainedToStart
        for (i in 1:4) ## for the 4 values in the list pars
            ## in original code, delete the last element of the list each step (???)
#            "TRpsFl = Transpose[psFl];
#       TRpsFl[[6]] = ALstarts;
#      psFl = Transpose[TRpsFl];" ## Ca n'a aucun sens
        ## Set up the list of contraints
        cons <- list(c(
           (0 < k < kMax, 0 < mB < mMax),
           (0 < k < kMax, 0 < mB < mMax, 0 < mB + DmB < mMax),
           (0 < k < kMax, 0 < m1M < mMax, 0 < m1F < mMax),
           (0 < k < kMax, 0 < m1M < mMax, 0 < m1M + DmM < mMax,
               0 < m1F < mMax, 0 < m1F + DmF < mMax)
        ))
        } else{  ## If we constrain the hybridization effect alpha          
       cons <- list(c(
           (0 < k < kMax, 0 < mB < mMax, AlphaLB < alpha < AlphaUB),
           (0 < k < kMax, 0 < mB < mMax, 0 < mB + DmB < mMax,
               AlphaLB < alpha < AlphaUB),
           (0 < k < kMax, 0 < m1M < mMax, 0 < m1F < mMax,
               AlphaLB < alpha < AlphaUB),
           (0 < k < kMax, 0 < m1M < mMax, 0 < m1M + DmM < mMax,
               0 < m1F < mMax, 0 < m1F + DmF < mMax, AlphaLB < alpha < AlphaUB)
       ))
        }

"
pars = {
   {c, mB, alpha},
   {c, mB, DmB, alpha},
   {c, m1M, m1F, alpha},
   {c, m1M, m1F, DmM, DmF, alpha}
   };
psFl = {
   {c, mB, mB,    0,     0,    alpha},
   {c, mB, mB,    DmB, DmB, alpha},
   {c, m1M, m1F, 0,     0,    alpha},
   {c, m1M, m1F, DmM, DmF, alpha}
   };

If[ALconstrainedToStart,
  For[k = 1, k <= 4, pars[[k]] = Delete[pars[[k]], Length[pars[[k]]]];
    k++];
  TRpsFl = Transpose[psFl];
  TRpsFl[[6]] = ALstarts;
  psFl = Transpose[TRpsFl];
  cons = {
    {0 < c < cMax, 0 < mB < mMax},
    {0 < c < cMax, 0 < mB < mMax, 0 < mB + DmB < mMax},
    {0 < c < cMax, 0 < m1M < mMax, 0 < m1F < mMax},
    {0 < c < cMax, 0 < m1M < mMax, 0 < m1M + DmM < mMax,
     0 < m1F < mMax, 0 < m1F + DmF < mMax}
    },(* otherwise *)
  cons = {
    {0 < c < cMax, 0 < mB < mMax, AlphaLB < alpha < AlphaUB},
    {0 < c < cMax, 0 < mB < mMax, 0 < mB + DmB < mMax,
     AlphaLB < alpha < AlphaUB},
    {0 < c < cMax, 0 < m1M < mMax, 0 < m1F < mMax,
     AlphaLB < alpha < AlphaUB},
    {0 < c < cMax, 0 < m1M < mMax, 0 < m1M + DmM < mMax,
     0 < m1F < mMax, 0 < m1F + DmF < mMax, AlphaLB < alpha < AlphaUB}
    }
  ];"
    
testList <- split(1:10,1:10)
append(testList, list(x=42),3)    

    ## Compare j to 1,2,3,4 and return the matching outputs 
    starter <- function(prevMLEest, j){
        if (j==1) {
            output <- start
        } else if (j==2) {
            ## Add 0 after 2nd element of prevMLEest[[2]]
            output <- append(prevMLEest[[2]], list(0), 2] 
        } else if (j==3) {
             ## Add ?? after 2nd element of prevMLEest[[2]]
            output <- append(prevMLEest[[2]], prevMLEest[[2,2]], 2)
        } else
            append(prevMLEest[[2]], prevMLEest[[2,2]], 2] ## Add ?? after 2nd element of prevMLEest[[2]]
            
"   1, start,
   2, Insert[prevMLEest[[2]], 0, 3],   
   3, Insert[Delete[prevMLEest[[2]], 3], prevMLEest[[2, 2]], 3],
   4, Insert[Insert[prevMLEest[[2]], prevMLEest[[2, 2]], 3],
    prevMLEest[[2, 3]], 5]
   ];
"
        ## What is MLEest??? preMLEest???

#### (*Print[Apply[LikelihoodFunction,Join[{PDF4cMeanLoad},psFl[[1]],{Data, HIpos,LOADpos}]]];*)
      
    MLEest, MLEbounds, FindExtremum, Extremum,
            TRpsFl,
            starter, starterCons, Temp, Hierarchy,
    vPrint, Mu, StdDev, k, Ans

## If verboselikeanalysis, Hierarchy takes H0, H1, H2 or H3
        Hierarchy <- vPrint(paste0("H", (i - 1))) 

        ## Define MLEest (enfin!)
     MLEest <- MLEformat[
      FindMaximum[{Apply[LikelihoodFunction,
         Join[{PDF4cMeanLoad}, psFl[[i]], {Data, HIpos, LOADpos}]],
         cons[[i]]}, Transpose[{pars[[i]], starter[MLEest, i]}]]];






        
    vPrint[MLEest[[1]]];
    vPrint[MatrixForm[{pars[[i]], MLEest[[2]]}]];
    MLEbounds = Table[
      FindExtremum = If[Extremum == -1, FindMinimum, FindMaximum];
      vPrint[
       StringJoin["...Finding ", If[Extremum > 0, "+", "-"],
        ToString[pars[[i, parID]]], " bound..."]];
      Ans =
       MLEformat[
        FindExtremum[{pars[[i, parID]],
          Apply[LikelihoodFunction,
            Join[{PDF4cMeanLoad},
             psFl[[i]], {Data, HIpos, LOADpos}]] >= MLEest[[1]] - 2,
          cons[[i]]}, Transpose[{pars[[i]], MLEest[[2]]}]]];
      vPrint[MatrixForm[{pars[[i]], Ans[[2]]}]];
      Ans,
      {parID, 1, Length[pars[[i]]]}, {Extremum, -1, 1, 2}];
    {MLEest, MLEbounds},
    {i, 1, Length[pars]}];
 Flatten[Join[
   Hierarchy, {Table[
     GtestOnNestedMLEs[Hierarchy[[i, 1]], Hierarchy[[1, 1]]], {i, 2,
      Length[Hierarchy]}]}], 1]
 ]

        "
LikelihoodAnalysis[PDF4cMeanLoad_, Data_, HIpos_, LOADpos_, start_,
  AlphaLB_, AlphaUB_, verboseLikelihoodAnalysis_] :=
 LikelihoodAnalysis[PDF4cMeanLoad, Data, HIpos, LOADpos, start,
  AlphaLB, AlphaUB, False, {}, verboseLikelihoodAnalysis]
LikelihoodAnalysis[PDF4cMeanLoad_, Data_, HIpos_, LOADpos_, start_,
  verboseLikelihoodAnalysis_] :=
 LikelihoodAnalysis[PDF4cMeanLoad, Data, HIpos, LOADpos, start, 2, 5,
  verboseLikelihoodAnalysis]
"

 ##????
LikelihoodAnalysis <- function(PDF4cMeanLoad, Data, HIpos, LOADpos, start,
                               AlphaLB, AlphaUB, verboseLikelihoodAnalysis){
    LikelihoodAnalysis(PDF4cMeanLoad, Data, HIpos, LOADpos, start,
  AlphaLB, AlphaUB, False, {}, verboseLikelihoodAnalysis]
LikelihoodAnalysis[PDF4cMeanLoad_, Data_, HIpos_, LOADpos_, start_,
  verboseLikelihoodAnalysis_] :=
 LikelihoodAnalysis[PDF4cMeanLoad, Data, HIpos, LOADpos, start, 2, 5,
  verboseLikelihoodAnalysis]


### Equivalent of FindOptim fucntion in R is optim
## http://stat.ethz.ch/R-manual/R-devel/library/stats/html/optim.html
