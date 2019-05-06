#-----------------------------------------------------------------------------------------
# Function name: HKF.LL (Hamilton Kalman Filter Log liklihood)
# Description: This function takes as inputs the coefficient matrices of the given ss model and starting values for initial state vector and CV and returns a vector with the log likelihood as well as the cumulative sum
#-----------------------------------------------------------------------------------------
# This function is fine - it is MLE which needs looking at

HKF.LL <- function(params, mod) {
  
  #mod <- build(params)
  
#  mod <- stage1Est(initial.parameters)
  
  
  tt <- dim(mod$y)[1]
  n <- dim(mod$y)[2]
  
  ll.vec <- matrix(NA,tt,1)
  ll.cum <- 0
  
  xi.tt <- mod$X0
  P.tt  <- mod$P0
  
  for (t in 1:tt){

    xi.ttm1 <- mod$FF %*% xi.tt + mod$cons # State t+1
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
  
  
  if(is.na(Param_LB) && is.na(Param_UB)) # no mod$constraints
  
    nloptr.out <- nloptr(params, f, eval_grad_f= function(x) {gradient(f, x)},
                         opts=list("algorithm"="NLOPT_LD_LBFGS","xtol_rel"=1.0e-8,"maxeval"=5000))
  
  if(is.na(Param_LB) && !is.na(Param_UB)) # Upper mod$constraint
    nloptr.out <- nloptr(params, f, eval_grad_f=function(x) {gradient(f, x)}, ub=theta.ub, opts=list("algorithm"="NLOPT_LD_LBFGS","xtol_rel"=1.0e-8,"maxeval"=5000))
  
  if(!is.na(Param_LB) && is.na(Param_UB)) # lower mod$constraint
    nloptr.out <- nloptr(params, f, eval_grad_f=function(x) {gradient(f, x)},lb=theta.lb, opts=list("algorithm"="NLOPT_LD_LBFGS","xtol_rel"=1.0e-8,"maxeval"=5000))
  
  if(!is.na(Param_LB) && !is.na(Param_UB)) # both upper and lower mod$constraint
    nloptr.out <- nloptr(params, f, eval_grad_f=function(x) {gradient(f, x)},lb=theta.lb, ub=theta.ub, opts=list("algorithm"="NLOPT_LD_LBFGS","xtol_rel"=1.0e-8,"maxeval"=5000))
  
  
  
  
  return(nloptr.out)
  
}


#-----------------------------------------------------------------------------------------
# Function name: HKF.filter
#----------------------------------------------------------------------------------------

HKF.filter <-   function(mod, t=1, xi.tt1 = NA, P.tt1 = NA) {
  
  
   if(!is.na(xi.tt1) && !is.na(P.tt1)){
  #   
     xi.ttm1 <- as.vector(mod$FF %*% xi.tt1 + mod$cons)
     P.ttm1 <- mod$FF %*% P.tt1 %*% t(mod$FF) + mod$Q
     
   }else{
    
    # State updating eq
    xi.ttm1 <- as.vector(mod$FF %*% mod$X0 + mod$cons)
    
    # MSE updating eq
    P.ttm1 <- mod$FF %*% mod$P0 %*% t(mod$FF) + mod$Q
    
  }
  
  
  # Equations (13.2.9 - 13.2.10) from Hamilton Chapter 13, page 379
  prediction.error <- (as.vector(mod$y.data[t,]) - as.vector(t(mod$A) %*% as.vector(mod$x.data[t,])) - as.vector(t(mod$H) %*% xi.ttm1))
  
  # Mean Squared Error
  HPHR <- t(mod$H) %*% P.ttm1 %*% mod$H + mod$R
  
  # State updating equation
  xi.tt <- xi.ttm1 + as.vector(P.ttm1 %*% mod$H %*% solve(HPHR, prediction.error))
  
  # Updated MSE 
  P.tt <- P.ttm1 - P.ttm1 %*% mod$H %*% solve(HPHR, t(mod$H) %*% P.ttm1)
  
  if (t == dim(mod$y.data)[1]) {
    
    print("Kalman filter has run")
    return(list("xi.ttm1"=xi.ttm1, "P.ttm1"=P.ttm1, "xi.tt"=xi.tt, "P.tt"=P.tt))
    
    
  } else {
    
    tmp <- HKF.filter(mod, t+1, xi.tt1 = xi.tt , P.tt1 = P.tt)
    
    return(list("xi.ttm1"=rbind(xi.ttm1, tmp$xi.ttm1),
                "P.ttm1"=rbind(P.ttm1, tmp$P.ttm1),
                "xi.tt"=rbind(xi.tt, tmp$xi.tt),
                "P.tt"=rbind(P.tt, tmp$P.tt),
                "SS model" = mod
                )
           
           )
    
  }
}

#-----------------------------------------------------------------------------------------
# Function name: HKF.smoother
#----------------------------------------------------------------------------------------

HKF.smooth <- function(mod, t = 0,  xi.tp1T=NA, P.tp1T=NA){
  
  #if(t = 0){
   # stop("You must enter the number of time steps in Y")
  #}
  
  
  n <- dim(mod$xi.ttm1)[2]
  
  if (t == dim(mod$`SS model`$y.data)[1]) {
    
    xi.tT <- mod$xi.tt[t,]
    
    P.tT <- mod$P.tt[((t-1)*n+1):(t*n),]
    
    tmp <- HKF.smooth(mod, t = t-1, xi.tT, P.tT)
    
    print("Kalman smoother has run")
    
    return(list("xi.tT"=rbind(tmp$xi.tT, xi.tT),
                "P.tT" =rbind(tmp$P.tT, P.tT)))
    
    
                                                } else {
    P.tt <- mod$P.tt[((t-1)*n+1):(t*n),]
    P.tp1t <- mod$P.ttm1[(t*n+1):((t+1)*n),]
    
    J.t <- P.tt %*% t(mod$`SS model`$FF) %*% solve(P.tp1t)
    
    xi.tt <- mod$xi.tt[t,]
    xi.tp1t <- mod$xi.ttm1[t+1,]
    
    
    xi.tT <- xi.tt + as.vector(J.t %*% (xi.tp1T - xi.tp1t))
    P.tT <- P.tt + J.t %*% (P.tp1T - P.tp1t) %*% t(J.t)
    
    if (t > 1) {
      
      tmp <- HKF.smooth(mod, t = t-1 ,xi.tp1T = xi.tT, P.tp1T =  P.tT)
      
      return(list("xi.tT"=rbind(tmp$xi.tT, xi.tT),
                  "P.tT" =rbind(tmp$P.tT, P.tT)))
      
    } else {
      
      return(list("xi.tT"=xi.tT, "P.tT"=P.tT))
    }
  }
  
}



#-----------------------------------------------------------------------------------------
# Function name: HKF.initalstates
#----------------------------------------------------------------------------------------
# Test this returns a object of class HKF and does all the checks

