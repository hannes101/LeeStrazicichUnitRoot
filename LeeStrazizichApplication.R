# Application of the basic Lee-Strazizich test without any parallelization
# Definition of values to be used in the application 
# Data series, which should be examined
y <- myData
# Number of breaks assumed in the series
myBreaks = 2
# Number of lags, if method "Fixed" and maximum number of lags from which the GTOS is used.
myLags <- 5
# Nature of the break, which is assumed in the series
# "crash" assumes breaks in intercept
# "break" assumes breaks in intercept and trend
myBreak <- "crash"

myLS_test <- ur.ls(y=y , model = myBreak, breaks = myBreaks, lags = myLags, method = "GTOS",pn = 0.1, print.results = "print" )

#Application of the parallelized Lee-Strazizich test
library(foreach)
library(doSNOW)
library(parallel)
#Define number of cores to use. By default the maximum available number minus one core is used
cl <- makeCluster(max(1, detectCores() - 1))
registerDoSNOW(cl)

myParallel_LS <- ur.ls.bootstrap(y=y , model = myBreak, breaks = myBreaks, lags = myLags, method = "Fixed",pn = 0.1, critval = "bootstrap", print.results = "print")
