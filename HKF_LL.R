#-----------------------------------------------------------------------------------------
# Function name: HKF.LL (Hamilton Kalman Filter Log liklihood)
# Description: This function takes as inputs the coefficient matrices of the given ss model and starting values for initial state vector and CV and returns a vector with the log likelihood as well as the cumulative sum
#-----------------------------------------------------------------------------------------


HKF.LL <- function(params,build) {
  
  mod <- build(params)
  
  tt <- dim(mod$y)[1]
  n <- dim(mod$y)[2]
  
  ll.vec <- matrix(mod$A,tt,1)
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
  nloptr.out <- nloptr(params, f, eval_grad_f=function(x) {gradient(f, x)}, opts=list("algorithm"="NLOPT_LD_LBFGS","xtol_rel"=1.0e-8,"maxeval"=200))
  if(is.na(Param_LB) && !is.na(Param_UB)) # Upper constraint
    nloptr.out <- nloptr(params, f, eval_grad_f=function(x) {gradient(f, x)}, ub=theta.ub, opts=list("algorithm"="NLOPT_LD_LBFGS","xtol_rel"=1.0e-8,"maxeval"=200))
  if(!is.na(Param_LB) && is.na(Param_UB)) # lower constraint
    nloptr.out <- nloptr(params, f, eval_grad_f=function(x) {gradient(f, x)},lb=theta.lb, opts=list("algorithm"="NLOPT_LD_LBFGS","xtol_rel"=1.0e-8,"maxeval"=200))
  if(!is.na(Param_LB) && !is.na(Param_UB)) # both upper and lower constraint
    nloptr.out <- nloptr(params, f, eval_grad_f=function(x) {gradient(f, x)},lb=theta.lb, ub=theta.ub, opts=list("algorithm"="NLOPT_LD_LBFGS","xtol_rel"=1.0e-8,"maxeval"=200))
  
  
  params <- nloptr.out$solution
  
  return(params)
  
}

#-----------------------------------------------------------------------------------------
# TESTING
#----------------------------------------------------------------------------------------
# Test this returns a object of class HKF and does all the checks


stage1 <-  HKF(
  
  A = matrix(0, 7, 2),
  H = matrix(0, 3, 2),
  R = diag(2),
  Q = matrix(0, 3, 3),
  FF = matrix(0,3,3),
  cons = matrix(0,3,1),
  
  y.data = matrix(rnorm(25),
                  rnorm(25), nrow = 25, ncol = 2),
  x.data = matrix(1,25,7),
  
  X0 = matrix(1,3,1),
  P0 = diag(0,3)
  
)  



# Test this returns parameters of MLE 

#Create "build" which is passed MLE

stage1Est <- function(params){
  
  stage1$A[1:2, 1] <- params[1:2]
  stage1$A[1, 2]   <- params[5]
  stage1$A[3:4, 2] <- params[3:4]
  stage1$A[5, 2]   <- 1-sum(stage1$A[3:4, 2])
  stage1$A[6:7, 2] <- params[6:7]
  
  stage1$H[1, 1]   <- 1
  stage1$H[2:3, 1] <- -params[1:2]
  stage1$H[2, 2]   <- -params[5]
  
  stage1$R         <- diag(c(params[9]^2, params[10]^2))
  stage1$Q         <- matrix(0, 3, 3)
  stage1$Q[1, 1]   <- params[11]^2
  
  stage1$FF[1, 1] <- stage1$FF[2, 1] <- stage1$FF[3, 2] <- 1
  
  stage1$cons[1, 1] <- params[8]
  
  return(stage1)
  
}

#MLE

HKF.MLE(build = stage1Est, params = c(rep(3,11)))

