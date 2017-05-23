## Aim for such an interface... in input functions in a seperately
## tested script before this!

## INPUT
## alice.ml.nb(bristles ~ latitude + sex + strain, data = drosophila)

## OUTPUT
## alpha                          -   estimate, lower, upper
## k                              -   estimate, lower, upper
## factor1..N * intercept         -   estimate, lower, upper
## factor1..N * growth            -   estimate, lower, upper
##                                    ML estimate


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

#########################
## Likelihood of a given subgroup of the data (aka male)

## The likelihood function over a set of inds
## param: k, alpha, mM, mF, slapeM, slapeF
LikelihoodFunction <- function(data, ugly.param, name="group2") { 
    param <- list()
    param[["old"]] <-  ugly.param[grepl("^k$|alpha|\\.old", names(ugly.param))]
    param[["young"]] <-  ugly.param[grepl("^k$|alpha|\\.young", names(ugly.param))]
    ## We consider so far 1 independant variable: sex, with 2 levels (TO IMPROVE)
    ## levels(group1)
    split.L<- by(data, data[, name], function(x)  {
        ## selecting the right parameters for this row
        PrLoadbyRow(x, param[[unique(x[, name])]])}) ### now hardcoded for group1
##    print(split.L)
    sum(log(unlist(split.L)))}
