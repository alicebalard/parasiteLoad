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
LogLik <- function(data, param, group.name){
    ## split the name into two
    gname <- sort(group.name)
    split.L<- by(data, data[, gname], function(x)  {
        ## by makes sure we get all levels: get the name of the paramter
        ## from the values within the by "loop"
        param.pattern <- unique(interaction(x[, gname], sep=":"))
        ## construct a regex with it
        par.regex <- paste0("^k$|alpha|^", param.pattern)
        ## select from our ugly paramter collection
        sub.param <- param[grepl(par.regex, names(param))]
        l.lik <- dnbinom(x$loads,
                         size=abs(sub.param[names(sub.param) %in% "k"]),
                         mu=abs(MeanLoad(alpha=sub.param[names(sub.param) %in% "alpha"],
                                         intercept=sub.param[grepl("inter",
                                                                   names(sub.param))],
                                         growth=sub.param[grepl("growth", names(sub.param))],
                                         HI=x$HI)),
                         log = TRUE)
       l.lik
    })
    all.l.lik <- unlist(split.L)
    if(length(all.l.lik)!=nrow(data)){
        stop("Not all likelihoods considered, group/paramter matching problem")
    } else{
        sum(all.l.lik)
    }
}

hybrid.maxim <- function (param, data, group.name){
    optim(par = param, 
          fn = LogLik, ## function to be maximized
          control = list(fnscale=-1),
          method = "L-BFGS-B",
          data = data,
          group.name = group.name)
}


library(parallel)

ML_bounds <- function(data, maxL,group.name, opt.param,
                      threshold = qchisq(p = 0.95, df = 1) / 2){
    ## maxL given by twologlik/2 from the model
    myfun <- function(par.name){
        sinP = opt.param[par.name]
        allbut1P = opt.param[!opt.param%in%par.name]
        ## Distance between likelihood and max likelihood:
        distance.L.ML.oneparfixed <- function(sinP, allbut1P, data, group.name, maxL){
            L <- LogLik(data = data, param = c(sinP, allbut1P), group.name = group.name)
            abs(L - (maxL - threshold))
        }
        ## Find the parameters (but one) that minimize the distance L - maxL:
        optim_at_threshold <- function(sinP, allbut1P, data, group.name, maxL){
            ## consider changing an algorithm that actually works on our problem
            myOpt <- optim(par = allbut1P,
                           fn = distance.L.ML.oneparfixed,
                           sinP = sinP,
                           data = data,
                           group.name = group.name,
                           maxL = maxL)
            myOpt$par[sinP]
        }
    ## Upper bound
    ## optim instead of optimize to avoid censoring??
        UB <- optim(f = optim_at_threshold,
                    interval = c(opt.param[[par.name]],
                                 opt.param[[par.name]] + 
                                 abs(opt.param[[par.name]]*100)),
                    maximum=TRUE)
        ## Lower bounds
        LB <- optimise(f = optim_at_threshold,
                       interval = c(opt.param[[par.name]]-
                                    abs(opt.param[[par.name]]*100),
                                    opt.param[[par.name]]),
                       maximum=FALSE)
        result <- list(c(LB = LB[1], UB = UB[1]))
        return(result)
    }
    ## Run over parameters:
    mclapply(names(opt.param), myfun, mc.cores=20)
}

