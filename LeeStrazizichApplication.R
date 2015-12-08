# Application of the basic Lee-Strazizich test without any parallelization
# Definition of values to be used in the application 

myLS_test <- ur.ls(y=y , model = "crash", breaks = 1, lags = NULL, method = "GTOS",pn = 0.1, print.results = "print" )

#Application of the parallelized Lee-Strazizich test
library(foreach)
library(doSNOW)
library(parallel)
#Define number of cores to use. By default the maximum available number minus one core is used
cl <- makeCluster(max(1, detectCores() - 1))
registerDoSNOW(cl)

myParallel_LS <- ur.ls.bootstrap(y=y , model = "crash", breaks = 1, lags = 5, method = "Fixed",pn = 0.1, critval = "bootstrap", print.results = "print")
