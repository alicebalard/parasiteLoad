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


# We calculate the expected load as follow (see WATWM paper)
# He(HI) = 2 * HI * (1 - HI)
# L(HI) = L1 + (L2 - L1) * HI
# L(HI, alpha) = max(L(HI) * (1 - alpha * He), delta) // delta fixed at 0.01
MeanLoad <- function(intercept, growth, alpha, hybridIndex){
  # print(intercept)
  # print(growth)
  print(alpha)
  # print(hybridIndex)
  heterozygoty <- 2 * hybridIndex * (1 - hybridIndex)
  mean <- (intercept + growth * hybridIndex) * (1 - alpha * heterozygoty)
  sapply(mean, function(x) {
    return(max(x, 0.01))
  })
  # print(mean)
}

## The likelihood function over a set of inds
LogLikelihood <- function(dataForGroup, paramsForGroup, response, gradient){
  # TODO : we could return -Infinity if the paramseters are outside boundaries
  # Here data is already splitted by groups
  hybridIndexesVector <- dataForGroup[, gradient]
  mu <- MeanLoad(
    alpha = paramsForGroup[names(paramsForGroup) %in% "alpha"],
    intercept = paramsForGroup[grepl("inter", names(paramsForGroup))],
    growth = paramsForGroup[grepl("growth", names(paramsForGroup))],
    hybridIndex = hybridIndexesVector
  )
  logLik <- stats::dnbinom(
    dataForGroup[, response],
    size = abs(paramsForGroup[names(paramsForGroup) %in% "k"]),
    mu = mu,
    log = TRUE
  )
  return(sum(logLik))
}

GetParamsForGroup <- function(params, dataForGroup, sortedGroupNames) {
  ## by makes sure we get all levels: get the name of the paramseter
  ## from the values within the by "loop"
  paramsSelectionPattern <- unique(
    interaction(
      dataForGroup[, sortedGroupNames], 
      sep=":"
      )
    )
  ## construct a regex with it
  paramsSelectionRegex <- paste0("^k$|alpha|^", paramsSelectionPattern)
  ## select from our ugly paramster collection
  paramsForGroup <- params[grepl(paramsSelectionRegex, names(params))]
  return(paramsForGroup)
}

# The likelihood analysis
MaximumLikelihood <- function (params, data, groupNames, response,
                          gradient, hessian = FALSE, control = list(fnscale=-1)){
  ## split the name into two
  sortedGroupNames <- sort(groupNames)
  OptimizationResultsByGroups <- by(
    data,
    data[, sortedGroupNames],
    function(dataForGroup)  {
      ## we now have the data splitted by groups
      ## we are still missing the paramss splitted by groups
      paramsForGroup <- GetParamsForGroup(params, dataForGroup, sortedGroupNames)
      # TODO test nloptr optimzer (bobyK ?)
      stats::optim(par = paramsForGroup,
                   fn = LogLikelihood, ## function to be maximized
                   control = control, ## maximise by default
                   hessian = hessian,
                  # method = "L-BFGS-B",
                   dataForGroup = dataForGroup,
                   response = response,
                   gradient = gradient)
    }
  )
  print(OptimizationResultsByGroups)
  return(OptimizationResultsByGroups)
}
