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
LikelihoodFunction <- function(data, ugly.param, name) { 
    ## ugly: get the name back as two strings, but at least it works
    ## for one factor or many in interactions
    gname <- strsplit(name, ":")[[1]] 
    split.L<- by(data, data[, gname], function(x)  {
        ## get the name of the paramter from the values within the by
        ## "loop"
        param.pattern <- unique(interaction(x[, gname], sep=":"))
        ## construct a regex with it
        par.string <- paste0("^k$|alpha|\\.", param.pattern)
        ## select from our ugly paramter collection
        this.param <- ugly.param[grepl(par.string, names(ugly.param))]
        PrLoadbyRow(x, this.param) # hit and run
    })
    sum(log(unlist(split.L)))
}
