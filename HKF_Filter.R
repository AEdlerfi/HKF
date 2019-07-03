#-----------------------------------------------------------------------------------------
# Function name: HKF.filter
#----------------------------------------------------------------------------------------

# @mod a list containing a State Space model
# @t the number of time steps
# Change all xi.p to xi.p (prediction)
# Change all xi.f to xi.f (filtered)
# Change all xi.f to xi.s (smoothed)


HKF.filter <-   function(mod, t=1, xi.f1 = NA, P.f1 = NA, type = 1) {
  
  
  if(!is.na(xi.f1) && !is.na(P.f1)){
    
    # State prediction    
    xi.p <- as.vector(mod$FF %*% xi.f1 + mod$cons)
    # State prediction error variance - covariance matrix
    P.p <- mod$FF %*% P.f1 %*% t(mod$FF) + mod$Q
    
  }else{ # Starts recursions 
    
    # State prediction
    xi.p <- as.vector(mod$FF %*% mod$X0 + mod$cons)
    
    # State prediction error variance - covariance matrix
    P.p <- mod$FF %*% mod$P0 %*% t(mod$FF) + mod$Q
    
  }
  
  
  # Prediction error
  prediction.error <- (as.vector(mod$y.data[t,]) - as.vector(t(mod$A) %*% as.vector(mod$x.data[t,])) - as.vector(t(mod$H) %*% xi.p))
  
  # Conditional variance
  HPHR <- t(mod$H) %*% P.p %*% mod$H + mod$R
  
  # Filtered state (Predicted state + Kalman gain)
  xi.f <- xi.p + as.vector(P.p %*% mod$H %*% solve(HPHR, prediction.error))
  
  # Filtered error variance covariance 
  P.f <- P.p - P.p %*% mod$H %*% solve(HPHR, t(mod$H) %*% P.p)
  
  
  # If we have reached the end of the data stop!
  if (t == dim(mod$y.data)[1]) {
    
    print("Kalman filter has run")
    return(list("xi.p"=xi.p, "P.p"=P.p, "xi.f"=xi.f, "P.f"=P.f))
    
    
  } else { # create next prediction
    
    tmp <- HKF.filter(mod, t+1, xi.f1 = xi.f , P.f1 = P.f)
    
    return(list("xi.p"=rbind(xi.p, tmp$xi.p),
                "P.p"=rbind(P.p, tmp$P.p),
                "xi.f"=rbind(xi.f, tmp$xi.f),
                "P.f"=rbind(P.f, tmp$P.f),
                "SS model" = mod
    )
    
    )
    
  }
}
