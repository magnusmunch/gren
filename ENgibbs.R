library(SuppDists)
library(HyperbolicDist)
lambda1 <- 5
lambda2 <- 5
beta <- 3
p <- 0.5
a <- lambda1^2/(4*lambda2)
b <- lambda2*beta^2

test1 <- rgig(100000, c(p, a, b))
test2 <- rgig(100000, c(-p, a, b))
hist(test1, breaks=60, freq=FALSE)
hist(1/((b/a)*test2), breaks=60, freq=FALSE, add=TRUE)

test3 <- rinvGauss(100000, sqrt(a/b), a)
hist(test1, breaks=60, freq=FALSE)
hist(1/((b/a)*test3), breaks=60, freq=FALSE, add=TRUE)

test4 <- rinvGauss(100000, sqrt(b/a), b)
hist(test1, breaks=60, freq=FALSE)
hist(1/test4, breaks=60, freq=FALSE, add=TRUE)

test5 <- rgig(100000, c(0.5, lambda1^2/(4*lambda2), lambda2*beta^2))
test6 <- rinvGauss(100000, 2*lambda2*abs(beta)/lambda1, lambda2*beta^2)
hist(test5, breaks=60, freq=FALSE)
hist(1/test6, breaks=60, freq=FALSE, add=TRUE)


rtau <- function(n, beta, lambda1, lambda2) {
  
  psi = lambda1^2/(4.0*lambda2)
  chi = lambda2 * beta
  
  Y <- rnorm(n)^2
  U <- runif(n)
  
  z <- sqrt(chi/psi) + Y/(2.0*psi) - sqrt(sqrt(chi)*Y/psi^1.5 + Y^2/(4.0*psi^2))
  prob <- 1.0/(1.0 + z*sqrt(psi/chi))
  tau <- ifelse(U <= prob, 1.0/z + 1.0, psi*z/chi + 1.0)
  
  return(tau)
  
}

test7 <- rgig(100000, c(p, a, b)) + 1
test8 <- rtau(100000, 1, 1, 1)
hist(test7, breaks=60, freq=FALSE)
hist(test8, breaks=60, freq=FALSE, add=TRUE)



library(Rcpp)
library(penalized)
sourceCpp("C:/Users/Magnus/Documents/phd/ENVB/code/ENgibbs_V01.cpp")
ENgibbs <- function(x, y, lambda1, lambda2, b0, K) {
  
  n <- nrow(x)
  p <- ncol(x)
  if(length(lambda2)==1) {
    lambda2 <- rep(lambda2, p)
  }
  
  sam <- gibbsC(x, y, n, p, lambda1, lambda2, b0, K)
  
  out <- list(beta=sam[1:p, (K + 1):(2*K)], tau=sam[1:p, (2*K + 1):(3*K)], omega=sam[1:n, 1:K])
  return(out)
  
}
set.seed(123)
n <- 50
p <- 1000
x <- matrix(rnorm(n*p), ncol=p, nrow=n)
beta <- matrix(rnorm(p))
y <- rbinom(n, 1, exp(x %*% beta)/(1 + exp(x %*% beta)))

test1 <- ENgibbs(x, y, 1, 1, rnorm(p), 10000)
test2 <- penalized(y, x, unpenalized=~0, lambda1=1, lambda2=1, model="logistic")
estgibbs <- apply(test1$beta, 1, function(b) {d <- density(b); return(d$x[which.max(d$y)])})
estpen <- test2@penalized
plot(estgibbs, estpen)


