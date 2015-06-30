ur.ls <- function(y, model = c("crash", "break"), breaks = 1, lags = NULL, method = c("GTOS","Fixed"), pi = 0.1){
  

  #Check sanity of the function call
  if (any(is.na(y))) 
    stop("\nNAs in y.\n")
  
  if(pi >= 1 || pi <= 0){
    stop("\n pi has to be between 0 and 1.")
  }
  if(method == "Fixed" && is.null(lags) == TRUE){
    stop("\n If fixed lag length should be estimated, the number 
          \n of lags to be included should be defined explicitely.")
  }
  
  #Add lagmatrix function
  lagmatrix <- function(x,max.lag){
    embed(c(rep(NA,max.lag),x),max.lag+1)
                                  }
  #Add diffmatrix function
  diffmatrix <- function(x,max.diff,max.lag){
    #Add if condition to make it possible to differentiate between matrix and vector                  
    if(is.vector(x) == TRUE ){
      embed(c(rep(NA,max.lag),diff(x,max.lag,max.diff)),max.diff)
    }
    
    else if(is.matrix(x) == TRUE){
      rbind(matrix(rep(NA,max.lag), ncol = ncol(x)), matrix(diff(x,max.lag,max.diff), ncol = ncol(x)))
    }
    #if matrix the result should be 0, if only a vector it should be 1
    else if(as.integer(is.null(ncol(x))) == 0 ){
      rbind(matrix(rep(NA,max.lag), ncol = ncol(x)), matrix(diff(x,max.lag,max.diff), ncol = ncol(x)))
      
    }
  }
  
  #Number of observations
  n <- length(y)
  model <- match.arg(model)
  lags <- as.integer(lags)
  method <- match.arg(method)
  breaks <- as.integer(breaks)
  #Percentage to eliminate endpoints in the lag calculation
  pi <- 0.1
  #Critical Values for the one break test
  model.one.crash.cval <- c(-4,239, -3.566, -3.211)
  
  model.one.break.cval <- matrix(c(.1 , -5.11, -4.5, -4.21,
                               .2 , -5.07, -4.47, -4.20,
                               .3 , -5.15, -4.45, -4.18,
                               .4 , -5.05, -4.50, -4.18,
                               .5 , -5.11, -4.51, -4.17), nrow = 5, ncol = 4, byrow = TRUE)
  colnames(model.one.break.cval) <- c("lambda","1%","5%","10%")
#   All critical values were derived in samples of size T = 100. Critical values
#   in Model C (intercept and trend break) depend (somewhat) on the location of the
#   break (λ = T_B /T) and are symmetric around λ and (1-λ). Model C critical values
#   at additional break locations can be interpolated.
#
#   Critical values for the endogenous two break test  
#   Model A - "crash" model  
#   invariant to the location of the crash
 model.two.crash.cval <- matrix(c("LM_thau", -4.545, -3.842, -3.504,
                                  "LM_rho", -35.726, -26.894, -22.892), nrow = 2, ncol = 4, byrow = TRUE ) 

 colnames(model.two.crash.cval) <-  c("Statistic","1%","5%","10%")
  

# Model C (i) - "break" model, breaks in the data generating process
# Model C (i) - "break" model invariant to the location of the crash
 
 model.two.break.dgp.cval <- matrix(c("LM_thau", -5.823, -5.286, -4.989,
                                      "LM_rho", -52.550, -45.531, -41.663), nrow = 2, ncol = 4, byrow = TRUE ) 
 
 
 # Model C (ii) - "break" model, breaks are not considered in the data generating process
 # Model C (ii) - "break" model depends on the location of the crash
 ## highest level of list is the location of the second breakpoint - so the share inside 
 ## the matrix refers to the first breakpoint
 model.two.break.thau.cval <- list("0.4" = list( matrix(c("LM_thau - 0.2", -6.16, -5.59, -5.27)
                                                                  , nrow = 1, ncol = 4, byrow = TRUE )),
                                  "0.6" =   list(matrix(c("LM_thau - 0.2", -6.41, -5.74, -5.32,
                                                                    "LM_thau - 0.4", -6.45, -5.67, -5.31)
                                                                  , nrow = 2, ncol = 4, byrow = TRUE )),
                                  "0.8" = list( matrix(c( "LM_thau - 0.2", -6.33, -5.71, -5.33,
                                                          "LM_thau - 0.4", -6.42, -5.65, -5.32,
                                                          "LM_thau - 0.6", -6.32, -5.73, -5.32)
                                                                  , nrow = 3, ncol = 4, byrow = TRUE )
 ))
 model.two.break.rho.cval <- list("0.4" = list( matrix(c("LM_rho - 0.2", -55.4, -47.9, -44.0)
                                                       , nrow = 1, ncol = 4, byrow = TRUE )),
                                  "0.6" =   list(matrix(c("LM_rho - 0.2", -58.6, -49.9, -44.4,
                                                          "LM_rho - 0.4", -59.3, -49.0, -44.3)
                                                        , nrow = 2, ncol = 4, byrow = TRUE )),
                                  "0.8" = list( matrix(c( "LM_rho - 0.2", -57.6, -49.6, -44.6,
                                                          "LM_rho - 0.4", -58.8, -48.7, -44.5,
                                                          "LM_rho - 0.6", -57.4, -49.8, -44.4)
                                                       , nrow = 3, ncol = 4, byrow = TRUE )
                                  ))
 
  #Number of observations to eliminate in relation to the sample length
  pinobs <- round(pi*n)
  
  
  #Define the start values
  startl <- 0
  myBreakStart <- startl + 1 + pinobs
  myBreakEnd <- n - pinobs
  
  #Calculate Dy
  y.diff <- diffmatrix(y, max.diff = 1, max.lag = 1)
  
  #Calculation
  #trend for 1:n like in ur.sp 
  trend <- 1:n
  
  #Define minimum gap between the two possible break dates.
  #the gap is 2 in the crash case and 3 periods in the case of a break model
  gap <- 2 + as.numeric(model == "break")
  
  
  myBreaks <- matrix(NA, nrow = n - 2 * pinobs, ncol =  breaks)
  if(breaks == 1){
    
    myBreaks [,1] <- myBreakStart:myBreakEnd
  } else if (breaks == 2){
          myBreaks[, 1:breaks] <- cbind(myBreakStart:myBreakEnd,(myBreakStart:myBreakEnd)+gap)
    }
  #Define the variables to hold the minimum t-stat
  tstat <- NA
  mint <- 1000
  tstat.matrix <- matrix(NA, nrow = n, ncol = n )
  tstat.result <- matrix()
  #Create lists to store the results
  #Lists for the one break case
  result.reg.coef <- list()
  result.reg.resid <- list()
  result.matrix <- list()
  
  #Function to analyze the optimal lags to remove autocorrelation from the residuals
  #Lag selection with general to specific procedure based on Ng,Perron (1995)
  #ToDo make it possible, that zero lags are included
  myLagSelection <- function(y.diff, S.tilde, datmat, pmax, Dummy.diff){
    n <- length(y.diff)
    
    #               General to specific approach to determine the lag which removes autocorrelation from the residuals
    #               Ng, Perron 1995
    qlag <- pmax
    while (qlag >= 0) {
      
      #               Define optimal lags to include to remove autocorrelation from the residuals
      #               select p.value of the highest lag order and check if significant
      test.coef <- coef(summary(lm(y.diff ~ 0 + lagmatrix(S.tilde,1)[,-1] + datmat[,-1][, 1:(qlag + 1)]  + Dummy.diff)))
      
      #                print(test.coef)
      #                print(test.coef[qlag + 1 , 4])
      #               print(c("Number of qlag",qlag)) 
      if(test.coef[qlag +1 , 4] <= 0.1){
        slag <- qlag
        #                  print("break")
        break
      }
      qlag <- qlag - 1
      slag <- qlag
      
    }
    #            print(slag)
    return(slag)
  }
  
  # Function to calculate the test statistic and make the code shorter, because the function can be used in both cases for 
  # the one break as well as the two break case
  myLSfunc <- function(Dt, DTt, y.diff){
    
    Dt.diff <- diffmatrix(Dt, max.diff = 1, max.lag = 1)
    
    DTt.diff <- diffmatrix(DTt, max.diff = 1, max.lag = 1)
    
    S.tilde <- 0
    S.tilde <- c(0, cumsum(residuals((summary(lm(y.diff ~ 0 + Dt.diff[,]))))))
    S.tilde.diff <-  diffmatrix(S.tilde,max.diff = 1, max.lag = 1)
    #       Define optimal lags to include to remove autocorrelation
    #       max lag length pmax to include is based on Schwert (1989) 
    pmax <- min(round((12*(n/100)^0.25)),lags)
    
    #      Create matrix of lagged values of S.tilde.diff from 0 to pmax
    #      and check if there is autocorrelation if all these lags are included and iterate down to 1          
    datmat <- matrix(NA,n, pmax + 2)
    datmat[ , 1] <- S.tilde.diff
    #      Add column of 0 to allow the easy inclusion of no lags into the test estimation
    datmat[ , 2] <- rep(0, n)
    
    if(pmax > 0){
      datmat[, 3:(pmax + 2) ] <- lagmatrix(S.tilde.diff, pmax)[,-1]
      colnames(datmat) <- c("S.tilde.diff", "NoLags",  paste("S.tilde.diff.l",1:pmax, sep = ""))
    } else if(lags == 0){
      colnames(datmat) <- c("S.tilde.diff", "NoLags")
    }
    
    
    if(method == "Fixed"){
      slag <- lags
    } else if(method == "GTOS"){
      
      slag <- NA
      
      if(model == "crash"){
        slag <- myLagSelection(y.diff, S.tilde, datmat, pmax, Dt.diff)
        
      } else if(model == "break"){
        slag <- myLagSelection(y.diff, S.tilde, datmat, pmax, DTt.diff)
      }
    }
    
    S.tilde <- NA
    

    
    if(model == "crash"){
      S.tilde <- c(0, cumsum(residuals((summary(lm(y.diff ~ 0 + Dt.diff))))))
      S.tilde.diff <-  diff(S.tilde,1)
      
      #Add lag of S.tilde.diff to control for autocorrelation in the residuals
      roll.reg <- summary(lm(y.diff ~ 0 + lagmatrix(S.tilde, 1)[,-1] + datmat[,2:(slag+2)]  + Dt.diff))
      
      
     #print(roll.reg)
     
      return(roll.reg)
      
    } else if(model =="break"){
      S.tilde <- c( 0,cumsum(residuals((summary(lm(y.diff ~ 0 +  DTt.diff))))))
      S.tilde.diff <-  diff(S.tilde)
      
      
      #Add lag of S.tilde.diff to control for autocorrelation in the residuals
      roll.reg <- summary(lm(y.diff ~ 0 + lagmatrix(S.tilde,1)[,-1] + datmat[,2:(slag+2)] + DTt.diff))
      
      #Return Value
    # print(roll.reg)
      #print("break")
      return(roll.reg)
      
    }
    
    print(roll.reg)
    return(roll.reg)
    
  }

# Start of the actual function call
  
if(breaks == 1)
{
# Function to calculate the rolling t-stat
# One Break Case
  for(myBreak1 in myBreaks[,1]){
    #Dummies for one break case
    Dt1 <-  as.matrix(cbind(trend, trend >= (myBreak1 + 1)))
    
    #       Dummy with break in intercept and in trend
    DTt1 <- as.matrix(cbind(Dt1, c(rep(0, myBreak1), 1:(n - myBreak1))))
    colnames(Dt1) <- c("Trend","D")
    colnames(DTt1) <- c("Trend","D","DTt")
    #print(paste("Break1: ",myBreak1, sep = ""))

      #Combine all Dummies into one big matrix to make it easier to include in the regressions
      
      Dt <- cbind(Dt1)
      DTt <- cbind(DTt1)
      
      result.reg <- myLSfunc(Dt,DTt,y.diff)
      
      
      
      #Extract the t-statistic and if it is smaller than all previous 
      #t-statistics replace it and store the values of all break variables
      #Extract residuals and coefficients and store them in a list
      
      result.matrix[[myBreak1]] <- result.reg
      result.reg.resid[[myBreak1]] <- residuals(result.reg)
      result.reg.coef[[myBreak1]] <- coefficients(result.reg)
      
      
      tstat <- coef(result.reg)[1,3]
      tstat.result[myBreak1] <- coef(result.reg)[1,3]
      #print(tstat)
      if(tstat < mint){
        mint <- tstat
        mybestbreak1 <- myBreak1
      }
      
   
  }#End of first for loop
  
} else if(breaks == 2) {
  
  
## Two Break Case
  #First for loop for the two break case
  for(myBreak1 in myBreaks[,1]){
    #Dummies for one break case
    Dt1 <-  as.matrix(cbind(trend, trend >= (myBreak1 + 1)))
    
    #       Dummy with break in intercept and in trend
    DTt1 <- as.matrix(cbind(Dt1, c(rep(0, myBreak1), 1:(n - myBreak1))))
    colnames(Dt1) <- c("Trend","D")
    colnames(DTt1) <- c("Trend","D","DTt")
    #print(paste("Break1: ",myBreak1, sep = ""))
    
    #Second for loop for the two break case
    for(myBreak2 in (myBreak1+gap):myBreakEnd){
      
      #Dummies for two break case
      Dt2 <-  as.matrix(trend >= (myBreak2 + 1))
      DTt2 <- as.matrix(cbind(Dt2, c(rep(0, myBreak2), 1:(n - myBreak2))))
      colnames(Dt2) <- c("D2")
      colnames(DTt2) <- c("D2","DTt2")
      #print(paste("Break2: ",myBreak2, sep = ""))
      
      #Combine all Dummies into one big matrix to make it easier to include in the regressions
      
      Dt <- cbind(Dt1, Dt2)
      DTt <- cbind(DTt1, DTt2)
      
      result.reg <- myLSfunc(Dt,DTt,y.diff)
      
      #Extract the t-statistic and if it is smaller than all previous 
      #t-statistics replace it and store the values of all break variables
      #Extract residuals and coefficients and store them in a list
      
      
      #matrix to hold all the tstats
      tstat.matrix[myBreak1, myBreak2] <- coef(result.reg)[1,3]
      tstat <- coef(result.reg)[1,3]
      
      #print(tstat)
      if(tstat < mint){
              mint <- tstat
              mybestbreak1 <- myBreak1
              mybestbreak2 <- myBreak2
                  }
     
      #Break the second for loop,
      if(myBreak2 >= myBreakEnd){
        break
      }
      }#End of second for loop
    }#End of first for loop
} else if(breaks > 2){
  
  print("Currently more than two possible structural breaks are not implemented.")
}
  
  #Estimate regression results, based on the determined breaks and the selected lag 
  # to obtain all necessary information
  Dt1 <-  as.matrix(cbind(trend, trend >= (mybestbreak1 + 1)))
  
  #       Dummy with break in intercept and in trend
  DTt1 <- as.matrix(cbind(Dt1, c(rep(0, mybestbreak1), 1:(n - mybestbreak1))))
  colnames(Dt1) <- c("Trend","D")
  colnames(DTt1) <- c("Trend","D","DTt")
  #print(paste("Break1: ",myBreak1, sep = ""))
  
    if(breaks == 2){
        #Dummies for two break case
        Dt2 <-  as.matrix(trend >= (mybestbreak2 + 1))
        DTt2 <- as.matrix(cbind(Dt2, c(rep(0, mybestbreak2), 1:(n - mybestbreak2))))
        colnames(Dt2) <- c("D2")
        colnames(DTt2) <- c("D2","DTt2")
        #print(paste("Break2: ",myBreak2, sep = ""))
        
        #Combine all Dummies into one big matrix to make it easier to include in the regressions
        
        Dt <- cbind(Dt1, Dt2)
        DTt <- cbind(DTt1, DTt2)
    } else if (breaks == 1){
      Dt <- Dt1
      DTt <- DTt1
      
    }
  
    
    break.reg <- myLSfunc(Dt,DTt,y.diff) 

  print(mint)
  print(paste("First possible structural break at position:", mybestbreak1))
  print(paste("The location of the first break - lambda_1:", round(mybestbreak1/n, digits = 1),", with the number of total observations:", n))  
  
  if(breaks == 2){
    print(paste("Second possible structural brgiteak at position:", mybestbreak2))
    print(paste("The location of the second break - lambda_2:", round(mybestbreak2/n, digits = 1),", with the number of total observations:", n))  
  }
  if(method == "Fixed"){
    print(paste("Number of lags used:",lags))
  } else if(method == "GTOS"){
    print(paste("Number of lags determined by general-to-specific lag selection:" 
                ,as.integer(substring(unlist(attr(delete.response(terms(break.reg)), "dataClasses")[3]),9))-1))
  } 
  #print(break.reg)
  return(break.reg)
  
}#End of ur.ls function
myLS_test <- ur.ls(y=y , model = "crash", breaks = 1, lags = 1, method = "Fixed",pi = 0.1 )