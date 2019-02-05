##############################  preamble  #############################
# simulations for grVBEM                                              #
# version: 01                                                         #
# author: Magnus M?nch                                                #
# created: 15-03-2017                                                 #
# last edited: 15-03-2017                                             #
#######################################################################

###############################  notes  ###############################
#                                                                     #
#######################################################################

### paths
path.code <- as.character(ifelse(Sys.info()[1]=="Darwin","/Users/magnusmunch/Documents/PhD/EBEN/code/" ,"~/EBEN/code/"))
path.res <- as.character(ifelse(Sys.info()[1]=="Darwin","/Users/magnusmunch/Documents/PhD/EBEN/results/" ,"~/EBEN/results/"))
path.data <- as.character(ifelse(Sys.info()[1]=="Darwin","/Users/magnusmunch/Documents/PhD/EBEN/data/hepB/" ,"~/EBEN/data/hepB/"))
path.graph <- "/Users/magnusmunch/Documents/PhD/EBEN/graphs/"

### libraries
library(Rcpp)
library(glmnet)
library(penalized)
library(GRridge)

### source grENVB functions
source(paste(path.code, "grVBEM.R", sep=""))

### a robust grridge function, that replaces lambda=Inf with lambda=10e5
rob.grridge <- function(highdimdata, response, partitions, unpenal, trace=FALSE) {
  tryCatch(grridge(highdimdata=highdimdata, response=response, partitions=partitions, unpenal=unpenal, trace=trace),
           error=function(e) {print("Infinite lambda");
             grridge(highdimdata=highdimdata, response=response, partitions=partitions, unpenal=unpenal, trace=trace,
                     optl=10e5)})
}

### the simulation
set.seed(1003)
p <- 400
n <- 100
G <- 5
groups <- rep(1:G, each=p/G)
beta <- rep(0.1, p)
nreps <- 100
m <- rep(1, n)

# mat.mult1 <- matrix(NA, nrow=nreps, ncol=2*G)
# for(r in 1:nreps) {
#   
#   set.seed(r + 100)
#   print(paste("Iteration: ", r, sep=""))
#   x <- matrix(rnorm(n*p), ncol=p, nrow=n)
#   y <- rbinom(n, 1, as.numeric(exp(x %*% as.matrix(beta))/(1 + exp(x %*% as.matrix(beta)))))
#   
#   fit.grVBEM <- grVBEM(x, y, m, groups, lambda1=NULL, lambda2=NULL, intercept=TRUE, 
#                        eps=0.001, maxiter=500, trace=TRUE)
#   fit.grridge <- rob.grridge(t(x), y, list(CreatePartition(as.factor(groups))), unpenal=~1, trace=TRUE)
#   
#   mat.mult1[r, ] <- c(fit.grVBEM$lambdag[[1]][, fit.grVBEM$nouteriter + 1], 
#                      fit.grridge$lambdamults[[1]])
#   
# }
# 
# colnames(mat.mult1) <- c(paste(rep("grEBEN", 5), 1:5, sep=""),
#                         paste(rep("GRridge", 5), 1:5, sep=""))
# save(mat.mult1, file=paste(path.res, "grVBEM_sim3_res1.Rdata", sep=""))

load(paste(path.res, "grVBEM_sim3_res1.Rdata", sep=""))
boxplot(mat.mult1)

opar <- par()
par(mfrow=c(2, 1), mar=c(2.1,4.1,2.1,2.1))
boxplot(mat.mult1[, 1:5], ylim=c(0, 3))
boxplot(mat.mult1[, 6:10], ylim=c(0, 3))
apply(mat.mult1, 2, median)
par(opar)


