#-----------------------------------------------------------------------------------------
# Function name: HKF.startingVals 
# Description: This function takes a given HKF model and creates starting values for the state mean and variance
#-----------------------------------------------------------------------------------------


HKF.startingVals <- function(y, mod)
 
n.state.vars <- ncol(mod$FF)

if (any(is.na(xi.00))) {
  xi.00 <- rep(0, n.state.vars)
  
  ##  Starting values for xi.00 and P.00
  x  <- rbind(t(mod$H), t(mod$H) %*% mod$FF, t(mod$H) %*% mod$FF %*% mod$FF, t(mod$H) %*% mod$FF %*% mod$FF %*% mod$FF)
  
  # Check if exogenous variable design matrix is compatible
  nmXdata <- c("cons")
  nmXdataind <- match(nmXdata, names(x))
  
  if(!is.na(nmXdataind)){
  
    m <- nrow(mod$A)+nrow(mod$cons)  
    
  }else{
    
    m <- nrow(mod$A) 
    
  }
  
  om <-  matrix(0, m, m)
  om[5:6, 3:4] <- t(mod$H) %*% mod$FF %*% mod$mod$Q %*% mod$H
  om[7:8, 3:4] <- t(mod$H) %*% mod$FF %*% mod$FF %*% mod$Q %*% mod$H
  om[7:8, 5:6] <- t(mod$H) %*% (mod$FF %*% mod$FF %*% mod$Q %*% t(mod$FF) + mod$FF %*% mod$Q) %*% mod$H
  om           <- om + t(om)
  om[1:2, 1:2] <- mod$R
  om[3:4, 3:4] <- t(mod$H) %*% mod$Q %*% mod$H + mod$R
  om[5:6, 5:6] <- t(mod$H) %*% (mod$FF %*% mod$Q %*% t(mod$FF) + mod$Q) %*% mod$H + mod$R
  om[7:8, 7:8] <- t(mod$H) %*% (mod$FF %*% mod$FF %*% mod$Q %*% t(mod$FF) %*% t(mod$FF) + mod$FF %*% mod$Q %*% t(mod$FF) + mod$Q) %*% mod$H + mod$R
  
  p1 <- t(x) %*% solve(om, x) ## x' * inv(om) * x
  yy <- c(y.data[1,], y.data[2,], y.data[3,], y.data[4,])
  tmp <- c(t(mod$FF) %*% x.data[1,],
           (t(mod$A) %*% x.data[2,] + t(mod$H) %*% mod$cons),
           (t(mod$A) %*% x.data[3,] + t(mod$H) %*% mod$cons + t(mod$H) %*% mod$FF %*% cons),
           (t(mod$A) %*% x.data[4,] + t(mod$H) %*% (diag(n.state.vars)+mod$FF+mod$FF%*%mod$FF) %*% mod$cons))
  
  xi.00 <- solve(p1, t(x)) %*% solve(om, yy - tmp)
  tmp <- yy - tmp - x %*% xi.00
  P.00 <- solve(p1, (diag(nrow=3) * sum(tmp^2) / 3))
}

return(list("xi.00"=xi.00, "P.00"=P.00, "F"=F, "Q"=Q, "mod$A"=A, "H"=H, "R"=R, "cons"=cons, "x.data"=x.data, "y.data"=y.data))



