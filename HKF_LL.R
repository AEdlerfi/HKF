#-----------------------------------------------------------------------------------------
# Function name: HKF.LL (Hamilton Kalman Filter Log liklihood)
# Description: This function takes as inputs the coefficient matrices of the given ss model and starting values for initial state vector and CV and returns a vector with the log likelihood as well as the cumulative sum
#-----------------------------------------------------------------------------------------


HKF.LL <- function(params,build) {
  
  mod <- build(params)
  
  tt <- dim(mod$y)[1]
  n <- dim(mod$y)[2]
  
  ll.vec <- matrix(NA,tt,1)
  ll.cum <- 0
  xi.tt <- mod$X0
  P.tt  <- mod$P0
  
  for (t in 1:tt){
    
    xi.ttm1 <- mod$FF %*% xi.tt + mod$cons
    P.ttm1 <- mod$FF %*% P.tt %*% t(mod$FF) + mod$Q
    
    prediction.error <- (as.vector(mod$y.data[t,]) - as.vector(t(mod$A) %*% as.vector(mod$x.data[t,])) - as.vector(t(mod$H) %*% xi.ttm1))
    
    HPHR <- t(mod$H) %*% P.ttm1 %*% mod$H + mod$R
    
    ll.vec[t] <- drop(-(n / 2) * log(2 * atan(1) * 4) - 0.5 * log(det(HPHR))
                      -0.5 * prediction.error %*% solve(HPHR, prediction.error))
    
    ll.cum <- ll.cum + ll.vec[t]
    
    xi.tt <- xi.ttm1 + P.ttm1 %*% mod$H %*% solve(HPHR, prediction.error)
    P.tt  <- P.ttm1 - P.ttm1 %*% mod$H %*% solve(HPHR, t(mod$H) %*% P.ttm1)
  }
  
  return(list("ll.vec"=ll.vec,"ll.cum"=ll.cum))

  }

#-----------------------------------------------------------------------------------------
# Function name: HKF.MLE (Hamilton Kalman Filter Maximum Liklihood Estimation)
# Description: This function takes as inputs the "build" model which is the SS model with parameters to be estimated in it. It computes the MLE using the function above
#-----------------------------------------------------------------------------------------


HKF.MLE <- function(build, params, Param_LB = NA, Param_UB = NA){
  
  
  
    f <- function(params) {
    
    return(-HKF.LL(params, build)$ll.cum)
    
  }
  
  
  gradient <- function(f, x, delta = x * 0 + 1.0e-5) {
     
    g <- x * 0
    for (i in 1:length(x)) {
      x1 <- x
      x1[i] <- x1[i] + delta[i]
      f1 <- f(x1)
      x2 <- x
      x2[i] <- x2[i] - delta[i]
      f2 <- f(x2)
      g[i] <- (f1 - f2) / delta[i] / 2
    }
    return(g)
  }
  
  if(is.na(Param_LB) && is.na(Param_UB)) # no constraints
  nloptr.out <- nloptr(params, f, eval_grad_f=function(x) {gradient(f, x)}, opts=list("algorithm"="NLOPT_LD_LBFGS","xtol_rel"=1.0e-8,"maxeval"=5000))
  if(is.na(Param_LB) && !is.na(Param_UB)) # Upper constraint
    nloptr.out <- nloptr(params, f, eval_grad_f=function(x) {gradient(f, x)}, ub=theta.ub, opts=list("algorithm"="NLOPT_LD_LBFGS","xtol_rel"=1.0e-8,"maxeval"=5000))
  if(!is.na(Param_LB) && is.na(Param_UB)) # lower constraint
    nloptr.out <- nloptr(params, f, eval_grad_f=function(x) {gradient(f, x)},lb=theta.lb, opts=list("algorithm"="NLOPT_LD_LBFGS","xtol_rel"=1.0e-8,"maxeval"=5000))
  if(!is.na(Param_LB) && !is.na(Param_UB)) # both upper and lower constraint
    nloptr.out <- nloptr(params, f, eval_grad_f=function(x) {gradient(f, x)},lb=theta.lb, ub=theta.ub, opts=list("algorithm"="NLOPT_LD_LBFGS","xtol_rel"=1.0e-8,"maxeval"=5000))
  
  
  
  
  return(nloptr.out)
  
}

#-----------------------------------------------------------------------------------------
# TESTING
#----------------------------------------------------------------------------------------
# Test this returns a object of class HKF and does all the checks

