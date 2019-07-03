#-----------------------------------------------------------------------------------------
# Function name: HKF.smoother
#----------------------------------------------------------------------------------------

# @mod a model of the class HKF.filtered 
# @t the number of time steps
# Change all xi.p to xi.p (prediction)
# Change all xi.s to xi.f (filtered)
# Change all xi.s to xi.s (smoothed)



HKF.smooth <- function(mod, t = 0,  xi.sT=NA, P.sT=NA){
  
  if(t == 0){
    stop("t = 0, You must enter the number of time steps in Y")
  }
  
  
  n <- dim(mod$xi.p)[2]
  
  if (t == dim(mod$`SS model`$y.data)[1]) {
    
    xi.s <- mod$xi.f[t,]
    
    P.s <- mod$P.f[((t-1)*n+1):(t*n),]
    
    tmp <- HKF.smooth(mod, t = t-1, xi.s, P.s)
    
    print("Kalman smoother has run")
    
    return(list("xi.s"=rbind(tmp$xi.s, xi.s),
                "P.s" =rbind(tmp$P.s, P.s)))
    
    
  } else {
    P.f <- mod$P.f[((t-1)*n+1):(t*n),]
    P.tp1t <- mod$P.p[(t*n+1):((t+1)*n),]
    
    J.t <- P.f %*% t(mod$`SS model`$FF) %*% solve(P.tp1t)
    
    xi.f <- mod$xi.f[t,]
    xi.tp1t <- mod$xi.p[t+1,]
    
    
    xi.s <- xi.f + as.vector(J.t %*% (xi.tp1T - xi.tp1t))
    P.s <- P.f + J.t %*% (P.tp1T - P.tp1t) %*% t(J.t)
    
    if (t > 1) {
      
      tmp <- HKF.smooth(mod, t = t-1 ,xi.sT = xi.s, P.sT =  P.s)
      
      return(list("xi.s"=rbind(tmp$xi.s, xi.s),
                  "P.s" =rbind(tmp$P.s, P.s)))
      
    } else {
      
      return(list("xi.s"=xi.s, "P.s"=P.s))
    }
  }
  
}