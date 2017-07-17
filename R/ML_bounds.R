##Approximation of the CI by hessian matrix
# Wald test (cf "Max Lik estimation and Inference book) p46:
nloglik <- function(data, param, group.name){
  -LogLik(data, param, group.name)
}

hybrid.minim <- function (param, data, group.name, hessian=FALSE){
  optim(par = param, 
        fn = nloglik, ## function to be maximized
        method = "L-BFGS-B",
        data = data,
        group.name = group.name,
        hessian = hessian)
}

ML_bounds_Wald <- function(param, data, group.name){
  # use start values inferred from glm.nb:
  fit.include.hessian <- hybrid.minim(param, data, group.name, hessian=TRUE)
  MLE <- fit.include.hessian$par
  ObsInfo <- fit.include.hessian$hessian # observed Fisher information matrix
  Vhat <- solve(ObsInfo) # inverse of observed Fisher information matrix
  Std.errors <- sqrt(diag(Vhat))
  # obtain the MLEs, estimated std errors, and approx Wald 95% CIs
  Wald.table <- cbind(MLE,
                      Std.errors,
                      LowerBounds = MLE - qnorm(0.975)*Std.errors,
                      UpperBounds = MLE + qnorm(0.975)*Std.errors)
  round(Wald.table, 4)
}