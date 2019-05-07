#-----------------------------------------------------------------------------------------
# Function name: HKF.smoother
#----------------------------------------------------------------------------------------

# @mod a model of the class HKF.filtered 
# @t the number of time steps
# @xi.tpT the final smoothed state
# @P.tpT the final smoothed variance - covariance


HKF.smooth <- function(mod, t = 0,  xi.tp1T=NA, P.tp1T=NA){
  
  if(t == 0){
    stop("t = 0, You must enter the number of time steps in Y")
  }
  
  
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

