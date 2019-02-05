##############################  preamble  #############################
# test cross-validation of lambda1 and lambda2                        #
# version: 01                                                         #
# author: Magnus M?nch                                                #
# created: 14-07-2017                                                 #
# last edited: 14-07-2017                                             #
#######################################################################

###############################  notes  ###############################
#                                                                     #
#######################################################################

### paths
path.code <- as.character(ifelse(Sys.info()[1]=="Darwin","/Users/magnusmunch/Documents/PhD/EBEN/code/" ,"~/EBEN/code/"))
path.res <- as.character(ifelse(Sys.info()[1]=="Darwin","/Users/magnusmunch/Documents/PhD/EBEN/results/" ,"~/EBEN/results/"))
path.data <- as.character(ifelse(Sys.info()[1]=="Darwin","/Users/magnusmunch/Documents/PhD/EBEN/data/" ,"~/EBEN/data/"))
path.graph <- "/Users/magnusmunch/Documents/PhD/EBEN/graphs/"

### libraries
library(glmnet)
library(mvtnorm)

### functions
# cross validate penalty parameters
cv.pen <- function(x, y, intercept) {
  n <- nrow(x)
  seq.alpha <- seq(0.01, 0.99, length.out=50)
  seq.lam <- seq.df <- seq.cvll <- numeric(length(seq.alpha))
  for(a in 1:length(seq.alpha)) {
    cv.fit <- cv.glmnet(x, y, family="binomial", alpha=seq.alpha[a], standardize=FALSE,
                        intercept=intercept)
    ind <- which(cv.fit$lambda==cv.fit$lambda.min)
    seq.lam[a] <- cv.fit$lambda.min
    seq.df[a] <- cv.fit$nzero[ind]
    seq.cvll[a] <- cv.fit$cvm[ind]
  }
  lambda1 <- seq.lam*seq.alpha*n*2
  #lambda1 <- n*seq.alpha*seq.lam
  lambda2 <- seq.lam*(1 - seq.alpha)*n
  
  out <- list(lambda1=lambda1, lambda2=lambda2, alpha=seq.alpha, lambda=seq.lam, cvll=seq.cvll)
  return(out)
}

# functions to simulate elastic net prior model parameters
qtrunc <- function(p, spec, a = -Inf, b = Inf, ...) {
  tt <- p
  G <- get(paste("p", spec, sep = ""), mode = "function")
  Gin <- get(paste("q", spec, sep = ""), mode = "function")
  tt <- Gin(G(a, ...) + p*(G(b, ...) - G(a, ...)), ...)
  return(tt)
}

rtrunc <- function(n, spec, a = -Inf, b = Inf, ...) {
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

denbeta <- function(x, lambda1, lambda2) {
  d <- sapply(x, function(x) {0.5*sqrt(lambda2)*dnorm(0.5*lambda1/sqrt(lambda2))*exp(-0.5*(lambda1*abs(x) + lambda2*x^2))/
                pnorm(-0.5*lambda1/sqrt(lambda2))})
  return(d)
}

varbeta <- function(lambda1, lambda2) {
  varb <- 0.25*lambda1^2/lambda2^2 + 1/lambda2 - 0.5*lambda1*dnorm(0.5*lambda1/sqrt(lambda2))/
    (lambda2^(3/2)*pnorm(-0.5*lambda1/sqrt(lambda2)))
  return(varb)
}

### testing parameterizations of prior
p <- 10000
lambda1 <- 2
lambda2 <- 6
test1 <- renbeta(p, lambda1, lambda2)
test2 <- sapply(seq(min(test1), max(test1), 0.01), function(x) {
  0.5*sqrt(lambda2)*dnorm(0.5*lambda1/sqrt(lambda2))*exp(-0.5*(lambda1*abs(x) + lambda2*x^2))/
    pnorm(-0.5*lambda1/sqrt(lambda2))})
hist(test1, freq=FALSE, breaks=80)
lines(seq(min(test1), max(test1), 0.01), test2)
varbeta(lambda1, lambda2)
var(test1)

### simulation
set.seed(2017)

# model parameter settings
lambda1 <- 1
lambda2 <- 1
p <- 1000

# predictor data settings
n <- 100
sigma2 <- 1
pblock <- 10
sigma <- matrix(rho, ncol=pblock, nrow=pblock)
diag(sigma) <- sigma2
sigma <- kronecker(diag(p/pblock), sigma)

# dependent data settings
m <- rep(1, n)

nreps <- 100
vec.lambda1 <- vec.lambda2 <- vec.varbeta <- numeric(nreps)
for(r in 1:nreps) {

  # create model parameters
  beta <- renbeta(p, lambda1, lambda2)

  # create predictor data
  x <- rmvnorm(n, mean=rep(0, p), sigma=sigma)

  # create dependent data
  prob <- as.numeric(exp(x %*% beta)/(1 + exp(x %*% beta)))
  y <- rbinom(n, m, prob)

  # cross-validate the lambdas and calculate variance of betas
  est.pen <- cv.pen(x, y, intercept=TRUE)
  vec.lambda1[r] <- est.pen$lambda1[which.min(est.pen$cvll)]
  vec.lambda2[r] <- est.pen$lambda2[which.min(est.pen$cvll)]
  vec.varbeta[r] <- var(beta)
}
true.varbeta <- varbeta(lambda1, lambda2)
emp.varbeta <- median(vec.varbeta)
cv.varbeta <- median(varbeta(vec.lambda1, vec.lambda2))
true.varbeta
emp.varbeta
cv.varbeta
hist(vec.cvvarbeta)

ndens <- 1000
seq.beta <- seq(min(beta), max(beta), length.out=ndens)
est.den <- denbeta(seq.beta, est.lambda1, est.lambda2)
h <- hist(beta, breaks=50, plot=FALSE)
plot(h, freq=FALSE, ylim=range(c(h$density, est.den)), xlim=range(c(seq.beta, h$breaks)))
lines(seq.beta, est.den, col=2)
text(x=c(1.65, 0.85), y=c(0.3, 1.0), labels=c(as.expression(bquote(sigma^2==.(round(var(beta), 2)))), bquote(sigma^2==.(round(varbeta(est.lambda1, est.lambda2), 2)))), col=c(1:2))

varbeta(est.lambda1, est.lambda2)
var(beta)
log(exp(1))

seq.varbeta <- varbeta(seq.lambda1, seq.lambda2)
plot(seq.lambda1, seq.varbeta, type="l")
plot(seq.lambda2, seq.varbeta, type="l")

  