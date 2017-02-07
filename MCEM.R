##############################  preamble  #############################
# MCEM implemented in R                                               #
# version: 01                                                         #
# author: Magnus Münch                                                #
# created: 03-02-2016                                                 #
# last edited: 07-02-2016                                             #
#######################################################################

###############################  notes  ###############################
#                                                                     #
#######################################################################

### paths
path.res <- "C:/Users/Magnus/Documents/phd/review/results/"
path.code <- "C:/Users/Magnus/Documents/phd/ENVB/code/"
path.graph <- "C:/Users/Magnus/Documents/phd/ENVB/graphs/"

### libraries
library(Rcpp)
library(penalized)
library(glmnet)

### functions
# source function to sample from posterior
sourceCpp(paste(path.code, "ENGibbs.cpp", sep=""))

# functions to simulate elastic net prior model parameters
qtrunc <- function(p, spec, a = -Inf, b = Inf, ...) {
  
  tt <- p
  G <- get(paste("p", spec, sep = ""), mode = "function")
  Gin <- get(paste("q", spec, sep = ""), mode = "function")
  tt <- Gin(G(a, ...) + p*(G(b, ...) - G(a, ...)), ...)
  return(tt)
  
}

rtrunc <- function(n, spec, a = -Inf, b = Inf, ...){
  
  x <- u <- runif(n, min = 0, max = 1)
  x <- qtrunc(u, spec, a = a, b = b,...)
  return(x)
  
}

renbeta <- function(p, lambda1, lambda2) {
  
  if(lambda1==0) {
    beta <- rnorm(p, 0, sqrt(1/lambda2))
  } else if(lambda2==0) {
    u <- runif(p) - 0.5
    beta <- -2*sign(u)*log(1 - 2*abs(u))/lambda1
  } else {
    shape <- 0.5
    scale <- 8*lambda2/lambda1^2
    tau <- rtrunc(n=p, spec="gamma", a=1, b=Inf, shape=shape, scale=scale)
    tau <- ifelse(tau < 1, 1, tau)
    sigma <- ifelse(is.infinite(tau), 1/lambda2, 
                    (tau - 1)/(lambda2*tau))
    beta <- rnorm(p, 0, sqrt(sigma))
  }
  return(beta)
  
}

mll <- function(lambda2, p, K, lambda1, sum1, sum2) {
  return(-p*pnorm(-lambda1/sqrt(4*lambda2), log.p=TRUE) - 0.5*lambda2*sum1/K - lambda1^2*sum2/(8*lambda2*K))
}

MCEM <- function(x, y, m, n, p, lambda1, lambda2, b0, intercept, K, eps, maxiter) {
  
  sam <- gibbsC(x, y, m, n, p, lambda1, rep(lambda2, p), b0, TRUE, 5000)
  beta <- apply(sam$beta, 1, function(b) {d <- density(b); return(d$x[which.max(d$y)])})
  if(intercept) {
    sum1 <- sum(sam$beta[-1, ]^2*sam$tau/(sam$tau - 1))
  } else {
    sum1 <- sum(sam$beta^2*sam$tau/(sam$tau - 1))
  }
  sum2 <- sum(sam$tau)
  
  mllseq <- mll(lambda2, p, K, lambda1, sum1, sum2)
  lambda2old <- lambda2seq <- lambda2
  
  conv <- FALSE
  iter <- 0
  while(!conv & (iter < maxiter)) {
    iter <- iter + 1
    print(paste("iteration ", iter, ", current lambda2 is ", round(lambda2, 2), sep=""))
    
    if(intercept) {
      sum1 <- sum(sam$beta[-1, ]^2*sam$tau/(sam$tau - 1))
    } else {
      sum1 <- sum(sam$beta^2*sam$tau/(sam$tau - 1))
    }
    sum2 <- sum(sam$tau)
    
    opt <- optim(par=lambda2old, fn=mll, p=p, K=K, lambda1=lambda1, sum1=sum1, sum2=sum2, control=list(fnscale=-1), 
                 method="Brent", lower=0.00001, upper=10000)
    lambda2 <- opt$par
    mllseq <- c(mllseq, opt$value)
    lambda2seq <- c(lambda2seq, lambda2)
    
    conv <- abs(lambda2 - lambda2old) < eps
    
    sam <- gibbsC(x, y, m, n, p, lambda1, rep(lambda2, p), beta, TRUE, 5000)
    beta <- apply(sam$beta, 1, function(b) {d <- density(b); return(d$x[which.max(d$y)])})
    
  }
  
  out <- list(beta=sam$beta, tau=sam$tau, omega=sam$omega, lambda2=lambda2seq, mll=mllseq)
  return(out)
  
}

### simulation
set.seed(123)
n <- 100
p <- 50
x <- matrix(rnorm(n*p), ncol=p, nrow=n)
lambda1 <- 1
lambda2 <- 1
#beta <- renbeta(p, lambda1=lambda1, lambda2=lambda2)
beta <- seq(0.1, 0.5, length.out=p)
m <- rep(1, n)
y <- rbinom(n, m, exp(x %*% beta)/(1 + exp(x %*% beta)))
b0 <- rnorm(p + 1)

test1 <- MCEM(x, y, m, n, p, lambda1, 10, b0, intercept=TRUE, K=5000, eps=0.00001, maxiter=100)
test2 <- optL2(y, x, unpenalized=~0, lambda1=1, model="logistic", fold=n)

plot(test1$lambda2, type="l")
  
test <- gibbsC(x, y, m, n, p, lambda1, rep(lambda2, p), b0, TRUE, 5000)
hist(test$beta[1, ], breaks=80)
abline(v=b0[1])
hist(test$beta[2, ], breaks=80)




