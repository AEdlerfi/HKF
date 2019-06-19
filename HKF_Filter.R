#-----------------------------------------------------------------------------------------
# Function name: HKF.filter
#----------------------------------------------------------------------------------------

# @mod a list containing a State Space model
# @t the number of time steps
# @xi.tt1 the state vector at time =1
# @xi.tt1 the state variance-covariance matrix at time =1


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
