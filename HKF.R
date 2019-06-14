#-----------------------------------------------------------------------------------------
# Function name: HKF (Hamilton Kalman Filter)
# Description: This function builds the required matricies that can then be passed to the KF algorithms. Also does some sense checks on the size and dimensions of the matricies
#-----------------------------------------------------------------------------------------

HKF <- function(...)
  {
  if (nargs()== 1){
    
    x <- as.list(...)
      
  }else{
    
    x <- list(...)  
    
    }
    
  # Create new model
  
  nm <- c("A","H","R","Q","FF","cons","X0","P0","y.data","x.data", "x.data.s")
  nmInd <- match(nm,names(x))
  # Check for na's
  if(any(is.na(nmInd)))
    stop(paste("Components(s)", paste(nm[is.na(nmInd)], collapse = ", "),"is (are) missing" ))
  
  x[nmInd[-1]] <- lapply(x[nmInd[-1]], as.matrix)
  
  if(!is.numeric(x$H)){
    stop("Component H must be numeric")
  }
  # Get dimensions of FF matrix
  m <- nrow(x$H)
  p <- ncol(x$H)
  if(!is.numeric(x$R)){
    stop("Component R must be numeric")
  }
  # Check if measurement variance matrix is diagonal
  if(!(nrow(x$R) == p && ncol(x$R) == p)){
    stop("Incompatiable dimensions of matricies - maybe R is not diagonal?")
  }
  
  if(!is.numeric(x$FF)){
    stop("Component FF must be numeric")
  }
  # Check if state transition matrix is compatiable (remember H is transposed)
  if(!(nrow(x$FF) == m && ncol(x$FF) == m)){
    stop("Incompatiable dimensions of matricies - maybe FF is not right?")
  }
  if (!is.numeric(x$Q)){
    stop("Component Q must be numeric")
  } 
  # Check if state variance is compatiable
  if (!(nrow(x$Q) == m && ncol(x$Q) == m)){
    stop("Incompatible dimensions of matrices")
  } 
  
  if (!is.numeric(x$A)){
    stop("Component A must be numeric")
  } 
  # Check if exogenous variable design matrix is compatible
  nmXdata <- c("x.data")
  nmXdataind <- match(nmXdata, names(x))
  if(!is.na(nmXdataind)){
    
    if (!(nrow(x$A) == ncol(x$x.data) && ncol(x$A) == p)){
      
      stop("Incompatible dimensions of matrices check A and x.data")
    } 
    
  } else{ # Not doing the correct thing here
    
    x.data <- matrix(0,dim(y.data[1]),dim(y.data)[2])
    
    
  }
  
  if (!is.numeric(x$cons)){
    stop("Component cons must be numeric")
  } 
    
  if (!(nrow(x$cons) == m)){
    stop("Incompatible dimensions of matrices - cons needs to be same rows as H")
  } 
    
  
  # Check if starting values for state variables are compatible
  if (!(nrow(x$P0) == m && ncol(x$P0) == m)) 
    stop("Incompatible dimensions of matrices")
  if (!(is.numeric(x$X0) && NCOL(x$X0) == 1 && NROW(x$X0) == 
        m)) 
    stop(paste("Component X0 must be a numeric vector of length", 
               "\n equal to ncol of component H', or a matrix with one column and", 
               "\n number of rows equal to ncol of component H'"))
  
  if (!(all.equal(x$P0, t(x$P0)) && all(eigen(x$P0)$values >= 
                                        0))) 
    stop("C0 is not a valid variance matrix")
  
  if (any(c(is.na(x$P0), is.na(x$X0)))) 
    stop("Missing values are not allowed in components x0 and P0")
  
  if (!(all.equal(x$R, t(x$R)) && all(eigen(x$R)$values >= 
                                        0))) 
      stop("V is not a valid variance matrix")
    if (!(all.equal(x$Q, t(x$Q)) && all(eigen(x$Q)$values >= 
                                        0))) 
      stop("W is not a valid variance matrix")
    
  
  mod <- x[nmInd]
  class(mod) <- "HKF"
    return(mod)
  }  
  

