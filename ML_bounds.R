##Calculate upper and lower bounds of each parameters
source("ML_functions.R")

ML_bounds <- function(data, group.name, startpara, threshold){
    maxL <- optim(par = startpara,
                  fn = LogLik, ## function to be maximized
                  control = list(fnscale=-1),
                  data = data,
                  group.name = group.name)
    ## Define the function for ONE parameter:
    myfun <- function(i){
        sinP = startpara[i]
        allbut1P = startpara[-i]
        ## Distance between likelihood and max likelihood:
        distance.L.ML.oneparfixed <- function(sinP, allbut1P, i, data, group.name, maxL){
            L <- LogLik(data = data, param = c(sinP, allbut1P), group.name = group.name)
            abs(L - (maxL$value - threshold))
        }
        ## Find the parameters (but one) that minimize the distance L - maxL:
        Singlepar_ML_bounds <- function(sinP, allbut1P, i, data, group.name, maxL){
            myOpt <- optim(par = allbut1P,
                           fn = distance.L.ML.oneparfixed,
                           sinP = sinP,
                           i= i,
                           data = data,
                           group.name = group.name,
                           maxL = maxL)
            myOpt$value
        }
        ## For one chosen parameter (rank i)
        chosenP <- function(sinP, i){
            names(sinP) <- names(startpara)[i]
            Singlepar_ML_bounds(sinP = sinP,
                                allbut1P = allbut1P,
                                i = i,
                                data = data,
                                group.name = group.name,
                                maxL = maxL)
        }
        ## Upper bound
        UB <- optimise(f = chosenP, i = i, interval = c(maxL$par[[i]], maxL$par[[i]] + 10))
        ## Lower bounds
        LB <- optimise(f = chosenP, i = i, interval = c(maxL$par[[i]] -10, maxL$par[[i]]))
        result <- list(c(LB = LB[1], est = maxL$par[[i]], UB = UB[1]))
        names(result) <- names(startpara[i])
        return(result)
    }
    ## Run over parameters:
    sapply(1:length(startpara), myfun)
}


## Example (NB takes +/- 15min to run):
source("Simulate_and_test.R")

AliceTest <- ML_bounds(data = simdata,
                       group.name = "group1:group2",
                       startpara = simpara,
                       threshold = qchisq(p = 0.95, df = 1) / 2)

## Compare with the actual parameters used for simulation:
M <- lapply(1:length(simpara), function(i) {
    M <- data.frame(c(simpara[i], unlist(AliceTest[i])))
    names(M) <- names(simpara[i])
    M
})

M <- as.data.frame(M)
rownames(M) <- c("simu", "LB", "est", "UB")
M





ML_bounds <- function(data, maxL, group.name, opt.param,
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
                   maximum=F83;40500;0cALSE)
    result <- list(c(LB = LB[1], UB = UB[1]))
    return(result)
  }
  ## Run over parameters:
  mclapply(names(opt.param), myfun, mc.cores=20)
}

####### TESTS ALICE:

ML_bounds_2 <- function(data, maxL, group.name, opt.param,
                        threshold = qchisq(p = 0.95, df = 1) / 2){
  ## maxL given by twologlik/2 from the model

  # set initial parameters:
  param <- opt.param  
  
  # for one given parameter:
  f1 <- function(param){
    ## Distance between likelihood and max likelihood:
      L <- LogLik(data, param, group.name)
      distance <- abs(L - (maxL - threshold))
      return(distance)
  }
  # optimise the distance, output ONE parameter
  f2 <- function(param, f1, par.name){
    opt_f1 <- optim(fn = f1,
                    par = param,
                    method = "L-BFGS-B")
   onepar <- opt_f1$par[par.name]
   return(onepar)
  }
  # optimisation of the optimisation:
  opt_f2 <- function(f1, f2, par.name, param){
    optim(fn = f2,
          par = param,
          f1 = f1,
          par.name = par.name)
  }
  ## Run over parameters:
  mclapply("k", opt_f2(f1, f2, "k", param), mc.cores=3)
}

## FAIL
#system.time(result <- ML_bounds_2(data = simdata, maxL = glm.h1$twologlik / 2, 
 #                                 group.name = c("group1", "group2"),
  #                                opt.param = glm.h1$opt.param, 
   #                               threshold = qchisq(p = 0.95, df = 1) / 2))


### Genius or lame?
# Wald test (cf "Max Lik estimation and Inference book) p46:


### PLEASE DON'T DUPLICATE CODE! TRY it, REALLY!

# use start values inferred from glm.nb:
parasite.fit <- optim(par = simpara, 
                      fn = nLogLik, ## function to be maximized
                      method = "L-BFGS-B",
                      data = simdata,
                      group.name =  c("group1", "group2"),
                      hessian = TRUE)

MLE <- parasite.fit$par
ObsInfo <- parasite.fit$hessian # observed Fisher information matrix
Vhat <- solve(ObsInfo) # inverse of observed Fisher information matrix
Std.errors <- sqrt(diag(Vhat))
# obtain the MLEs, estimated std errors, and approx Wald 95% CIs
Wald.table <- cbind(MLE,
                    Std.errors,
                    LowerBounds = MLE - qnorm(0.975)*Std.errors,
                    UpperBounds = MLE + qnorm(0.975)*Std.errors)
parnames <- names(opt.para$opt.param)
rownames(Wald.table) <- parnames
round(Wald.table, 4)
