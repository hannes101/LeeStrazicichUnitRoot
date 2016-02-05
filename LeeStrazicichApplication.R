# Application of the basic Lee-Strazizich test without any parallelization
# Definition of values to be used in the application 
# Number of possible structural breaks
myBreaks <- 2
# Assumed break in the series, "crash" - break in intercept; "break" - break in intercept and trend
myModel <- "break"
# Number of lags to be used in fixed specification or maximum number of lags, when using the GTOS method
myLags <- 5

myLS_test <- ur.ls(y=y , model = myModel, breaks = myBreaks, lags = myLags, method = "GTOS",pn = 0.1, print.results = "print" )

#Application of the parallelized Lee-Strazizich test
library(foreach)
library(doSNOW)
library(parallel)
#Define number of cores to use. By default the maximum available number minus one core is used
cl <- makeCluster(max(1, detectCores() - 1))
registerDoSNOW(cl)

myParallel_LS <- ur.ls.bootstrap(y=y , model = myModel, breaks = myBreaks, lags = myLags, method = "Fixed",pn = 0.1, critval = "bootstrap", print.results = "print")
