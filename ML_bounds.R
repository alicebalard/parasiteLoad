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




