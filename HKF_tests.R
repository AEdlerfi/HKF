#-----------------------------------------------------------------------------------------
# Testing HKF package - using the NILE data set
#----------------------------------------------------------------------------------------

# Set up model

HKFnile <- HKF(
  

    A = matrix(0, 1, 1),
    H = matrix(1, 1, 1),
    R = diag(1),
    Q = matrix(0,1,1),
    FF = matrix(1,1,1),
    cons = matrix(0,1,1),
    
    y.data =  NA, 
    x.data =  NA ,
    
    X0 = 0,
    P0 = 1000000
    
   
   
)

# Add data to model - need to create the functionality to create x.data if and hide it if univariate model is specified 

HKFnile$y.data <- as.matrix(Nile)
HKFnile$x.data <- matrix(0,dim(as.matrix(Nile))[1],dim(as.matrix(Nile))[2])

# Create build function to be passsed to MLE

HKFnileEst <- function(params){
  
  HKFnile$R <-  params[1] 
  HKFnile$Q <-  params[2] 
  
  return(HKFnile)
  
}

# Estimate parameters


NileParms <- HKF.MLE(HKFnileEst, c(100,100), useoptim = TRUE, hessian = TRUE)

# Estimate parameters

HKFnile <- HKFnileEst(NileParms$par)  # note the change if using optim

# Filter

HKFnilefilt <- HKF.filter(HKFnile)

# Smooth

timestep <- dim(HKFnilefilt$`SS model`$y.data)[1]

HKFnilesmooth <- HKF.smooth(HKFnilefilt, t = timestep)

#-----------------------------------------------------------------------------------------
# Compare with example from DLM vingette
#----------------------------------------------------------------------------------------

buildNile <- function(x) dlmModPoly(1, dV = x[1], dW = x[2])

fitNile <- dlmMLE(Nile, parm = rep(100, 2), build = buildNile, lower = rep(1e-8, 2), hessian = TRUE)

fitNile$par

NileParms$par

DLMNilefilt <- dlmFilter(y = Nile, mod= buildNile( fitNile$par ) ) 

tibble(`Filtered level HKF` = HKFnilefilt$xi.tt,
       `Filtered level DLM` = DLMNilefilt$m[-1],
       `Actual data` = Nile,
       Date = 1871:1970) %>% 
  gather(Model, Value, -Date) %>% 
  ggplot()+
  geom_line(aes(Date, Value, colour = Model))+
  tst_theme()+
  scale_colour_tst()
                         