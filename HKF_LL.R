#-----------------------------------------------------------------------------------------
# Function name: HKF.LL (Hamilton Kalman Filter Log liklihood)
# Description: This function takes as inputs the coefficient matrices of the given ss model and starting values for initial state vector and CV and returns a vector with the log likelihood as well as the cumulative sum
#-----------------------------------------------------------------------------------------
# Change all xi.ttm1 to xi.p (prediction)
# Change all xi.tt to xi.f (filtered)
# Change all xi.tT to xi.s (smoothed)


HKF.LL <- function(params, mod, type) {
  
    
  
  tt <- dim(mod$y)[1]
  n <- dim(mod$y)[2]
  
  ll.vec <- matrix(NA,tt,1)
  ll.cum <- 0
  
  xi.f <- mod$X0
  P.f  <- mod$P0

  if(type == 1){
    
  for (t in 1:tt){

    # Forcast state
    xi.p <- mod$FF %*% xi.f + mod$cons# State t+1
    # Forecast error covariance
    P.p <- mod$FF %*% P.f %*% t(mod$FF) + mod$Q 
    
    # Prediction error (exogenous variables in measurement eq)
    prediction.error <- (as.vector(mod$y.data[t,]) - as.vector(t(mod$A) %*% as.vector(mod$x.data[t,])) - as.vector(t(mod$H) %*% xi.p))
    
    # Prediction error decompostion
    HPHR <- t(mod$H) %*% P.p %*% mod$H + mod$R
    
    # log likelihood
    ll.vec[t] <- drop(-(n / 2) * log(2 * atan(1) * 4) - 0.5 * log(det(HPHR))
                      -0.5 * prediction.error %*% solve(HPHR, prediction.error))
    
    ll.cum <- ll.cum + ll.vec[t]
    
    # Filtered state (previous state + Kalman gain)
    xi.f <- xi.p + P.p %*% mod$H %*% solve(HPHR, prediction.error)
    P.f  <- P.p - P.p %*% mod$H %*% solve(HPHR, t(mod$H) %*% P.p)
  }
  
  }
  
  if(type == 2){
    
    for (t in 1:tt){
      
      # Forecast state 
      xi.p <- mod$FF %*% xi.f + mod$cons %*% as.vector(mod$x.data[t,]) # State t+1
      # Forecast error covariance 
      P.p <- mod$FF %*% P.f %*% t(mod$FF) + mod$Q 
      
      # Forecast error
      prediction.error <- (as.vector(mod$y.data[t,])-as.vector(t(mod$H) %*% xi.p))
      
      # Forecast error decomposition
      HPHR <- t(mod$H) %*% P.p %*% mod$H + mod$R
      
      # log likelihood
      ll.vec[t] <- drop(-(n / 2) * log(2 * atan(1) * 4) - 0.5 * log(det(HPHR))
                        -0.5 * prediction.error %*% solve(HPHR, prediction.error))
      
      ll.cum <- ll.cum + ll.vec[t]
      
      # Filter (state in time t + Kalman gain)
      xi.f <- xi.p + P.p %*% mod$H %*% solve(HPHR, prediction.error)
      P.f  <- P.p - P.p %*% mod$H %*% solve(HPHR, t(mod$H) %*% P.p)
    }
    
    
    
  }
    
  return(list("ll.vec"=ll.vec,"ll.cum"=ll.cum))

  }


