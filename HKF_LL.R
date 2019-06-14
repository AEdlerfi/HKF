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

    if(!is.na(mod$x.data.s)){
      
      mod$x.data.s[dim(mod$H)[1],1] <- 1 
      
        }
    
    xi.ttm1 <- mod$FF %*% xi.tt + mod$cons%*%as.vector(mod$x.data.s[t,])# State t+1
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


