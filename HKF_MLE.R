#-----------------------------------------------------------------------------------------
# Function name: HKF.MLE (Hamilton Kalman Filter Maximum Liklihood Estimation)
# Description: This function takes as inputs the "build" model which is the SS model with parameters to be estimated in it. It computes the MLE using the function above
#-----------------------------------------------------------------------------------------
# @build a list containing the state space model to be estimated.
# @params a vector of initial parameters
# @Param_LB lower bound constraints on the parameter vector
# @Param_UB upper bound constraint on the parameter vector 

HKF.MLE <- function(build, params, Param_LB = NA, Param_UB = NA, useoptim = FALSE, ...){
  
  
  
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
  
  
  f <- function(params) {
    
    mod <- build(params)
    
    return(-HKF.LL(params, mod)$ll.cum)
    
  }
  
  
  if(useoptim ==TRUE){
    
    nloptr.out <- optim(params, f, ...)
    
  } else{
  
    if(is.na(Param_LB) && is.na(Param_UB)) # no mod$constraints
      
      nloptr.out <- nloptr(params, f, eval_grad_f= function(x) {gradient(f, x)},
                           opts=list("algorithm"="NLOPT_LD_LBFGS","xtol_rel"=1.0e-8,"maxeval"=5000))
    
    if(is.na(Param_LB) && !is.na(Param_UB)) # Upper mod$constraint
      nloptr.out <- nloptr(params, f, eval_grad_f=function(x) {gradient(f, x)}, ub=Param_UB, opts=list("algorithm"="NLOPT_LD_LBFGS","xtol_rel"=1.0e-8,"maxeval"=5000))
    
    if(!is.na(Param_LB) && is.na(Param_UB)) # lower mod$constraint
      nloptr.out <- nloptr(params, f, eval_grad_f=function(x) {gradient(f, x)},lb=Param_LB, opts=list("algorithm"="NLOPT_LD_LBFGS","xtol_rel"=1.0e-8,"maxeval"=5000))
    
    if(!is.na(Param_LB) && !is.na(Param_UB)) # both upper and lower mod$constraint
      nloptr.out <- nloptr(params, f, eval_grad_f=function(x) {gradient(f, x)},lb=Param_UB, ub=Param_LB, opts=list("algorithm"="NLOPT_LD_LBFGS","xtol_rel"=1.0e-8,"maxeval"=5000))
    
    
      
  }
    
    
  
  return(nloptr.out)
  
}

