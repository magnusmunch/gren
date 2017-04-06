##############################  preamble  #############################
# simulations for grVBEM                                              #
# version: 01                                                         #
# author: Magnus Münch                                                #
# created: 15-03-2017                                                 #
# last edited: 15-03-2017                                             #
#######################################################################

###############################  notes  ###############################
#                                                                     #
#######################################################################

### paths
path.rcode <- "~/EBEN/code/"
path.res <- "~/EBEN/results/"
path.tab <- "C:/Users/Magnus/Documents/phd/ENVB/tables/"

### libraries
library(Rcpp)
library(glmnet)
library(penalized)
library(GRridge)

### source grENVB functions
source(paste(path.rcode, "grVBEM.R", sep=""))

### a robust grridge function, that replaces lambda=Inf with lambda=10e5
rob.grridge <- function(highdimdata, response, partitions, unpenal, trace=FALSE) {
  tryCatch(grridge(highdimdata=highdimdata, response=response, partitions=partitions, unpenal=unpenal, trace=trace),
           error=function(e) {print("Infinite lambda");
             grridge(highdimdata=highdimdata, response=response, partitions=partitions, unpenal=unpenal, trace=trace,
                     optl=10e5)})
}

### the simulation
p <- 300
n <- 100
G <- 5
groups <- rep(1:G, each=p/G)
beta <- rep(0.1, p)
nreps <- 100
m <- rep(1, n)

mat.mult <- matrix(NA, nrow=nreps, ncol=2*G)
for(r in 1:nreps) {
  
  set.seed(r + 100)
  print(paste("Iteration: ", r, sep=""))
  x <- matrix(rnorm(n*p), ncol=p, nrow=n)
  y <- rbinom(n, 1, as.numeric(exp(x %*% as.matrix(beta))/(1 + exp(x %*% as.matrix(beta)))))
  
  fit.grVBEM <- grVBEM(x, y, m, groups, lambda1=NULL, lambda2=NULL, intercept=TRUE, 
                       eps=0.001, maxiter=500, trace=TRUE)
  fit.grridge <- rob.grridge(t(x), y, list(CreatePartition(as.factor(groups))), unpenal=~1, trace=TRUE)
  
  mat.mult[r, ] <- c(fit.grVBEM$lambdag[[1]][, fit.grVBEM$nouteriter + 1], 
                     fit.grridge$lambdamults[[1]])
  
}

save(mat.mult, file=paste(path.res, "grVBEM_res3.Rdata", sep=""))

load("C:/Users/Magnus/Documents/phd/ENVB/results/grVBEM_res3.Rdata")
boxplot(mat.mult)

par(mfrow=c(2, 1), mar=c(2.1,4.1,2.1,2.1))
boxplot(mat.mult[, 1:5], ylim=c(0, 3))
boxplot(mat.mult[, 6:10], ylim=c(0, 3))
apply(mat.mult, 2, median)
par(opar)



