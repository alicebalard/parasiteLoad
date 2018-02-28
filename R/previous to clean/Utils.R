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

## The likelihood function over a set of inds
LogLik <- function(data, param, group.name, response, alpha.along, individualLogLik) {
  ## split the name into two _ TODO if needed only!
  gname <- sort(group.name)
  split.L<- by(data, data[, gname], individualLogLik(???))
  all.l.lik <- unlist(split.L)
  if(length(all.l.lik)!=nrow(data)){
    stop("Not all likelihoods considered, group/parameter matching problem")
  } else{
    sum(all.l.lik)
  }
}

# The likelihood analysis
hybrid.maxim <- function (param, data, group.name, response = response,
                          alpha.along, hessian=FALSE, control = list(fnscale=-1),
                          individualLogLik){
  stats::optim(par = param,
               fn = LogLik, ## function to be maximized
               control = control, ## maximise by default
               method = "L-BFGS-B",
               data = data,
               group.name = group.name,
               response = response,
               alpha.along = alpha.along,
               individualLogLik = individualLogLik, 
               hessian = hessian)
}
