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
  
  nm <- c("A","H","R","Q","FF","cons")
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
  nmXdata <- c("xdata")
  nmXdataind <- match(nmXdata, names(x))
  if(!is.na(nmXdataind)){
    
    if (!(nrow(x$A) == ncol(x$xdata) && ncol(x$A) == p)){
      
      stop("Incompatible dimensions of matrices")
    } 
    
  }
  
  if (!is.numeric(x$cons)){
    stop("Component C0 must be numeric")
  } 
    
  if (!(nrow(x$cons) == m && NCOL(x$cons) == 1)){
    stop("Incompatible dimensions of matrices")
  } 
    
  mod <- x[nmInd]
  class(mod) <- "HKF"
    return(mod)
  }  
  

# Test this returns a object of class HKF and does all the checks

  
HKF(
  
  A = matrix(0, 7, 2),
  H = matrix(0, 3, 2),
  R = diag(2),
  Q = matrix(0, 3, 3),
  FF = matrix(0,3,3),
  cons = matrix(0,3,1),
  xdata = matrix(cbind(1,1,1,1,1,1,1),
                 cbind(1,1,1,1,1,1,1), nrow = 2, ncol = 7)
  
  
)  

