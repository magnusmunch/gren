##############################  preamble  #############################
# comparison of marginal likelihood VB and Gibbs                      #
# version: 01                                                         #
# author: Magnus Münch                                                #
# created: 01-11-2016                                                 #
# last edited: 01-12-2016                                             #
#######################################################################

###############################  notes  ###############################
# 01-01-2016: Uses data from review graph                             #
#######################################################################

### paths
path.res <- "C:/Users/Magnus/Documents/phd/review/results/"
path.code <- "C:/Users/Magnus/Documents/phd/ENVB/code/"

### libraries
library(Rcpp)
library(penalized)
library(mvtnorm)
library(GeneralizedHyperbolic)

### functions
# import the Cpp Gibbs sampler
sourceCpp(paste(path.code, "ENgibbs.cpp", sep=""))

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

# simple ENVB function
ENVB <- function(y, x, lambda1, lambda2, b0, S0, eps, maxiter) {
  
  p <- ncol(x)
  kappa <- y - 0.5
  phi <- lambda1^2/(4*lambda2)
  
  mu <- muold <- b0
  sigma <- sigmaold <- S0
  
  ci <- ciold <- as.numeric(sqrt(colSums(t(x) * (sigma %*% t(x))) + (colSums(t(x)*mu))^2))
  chi <- chiold <- as.numeric(lambda2*(diag(sigma) + mu^2))
  
  niter <- 0
  conv <- FALSE
  while(!(conv | (niter==maxiter))) {
    
    niter <- niter + 1
    
    Om <- diag((0.5/ciold)*tanh(ciold/2))
    Z <- diag(sqrt(phi/chiold))
    
    sigma <- solve(t(x) %*% Om %*% x + lambda2*diag(p) + lambda2*Z)
    mu <- as.numeric(sigma %*% (t(x) %*% kappa))
    ci <- as.numeric(sqrt(colSums(t(x) * (sigma %*% t(x))) + (colSums(t(x)*mu))^2))
    chi <- as.numeric(lambda2*(diag(sigma) + mu^2))
    
    conv <- max(abs(c(sigma - sigmaold, mu - muold, ci - ciold, chi - chiold))) < eps
    
    sigmaold <- sigma
    muold <- mu
    ciold <- ci
    chiold <- chi
    
  }
  
  mll <- mllapprox(lambda2, lambda1, mu, sigma, phi, chi)
  
  out <- list(mu=mu, sigma=sigma, ci=ci, chi=chi, mll=mll, conv=conv, niter=niter)
  return(out)

}

mllapprox <- function(lambda2, lambda1, mu, sigma, phi, chi) {
  
  mll <- p*log(lambda1) - 0.5*lambda2*sum((diag(sigma) + mu^2)*(1 + sqrt(phi/chi))) - 
    p*pnorm(-lambda1/sqrt(4*lambda2), log.p=TRUE) - 
    lambda1^2*sum(1/phi + sqrt(chi/phi) + 1)/(8*lambda2)
  return(mll)
  
}

# Gibbs sampler
ENgibbs <- function(x, y, m, lambda1, lambda2, b0, intercept, K, mll=FALSE) {
  
  n <- nrow(x)
  p <- ncol(x)
  if(length(lambda2)==1) {
    lambda2 <- rep(lambda2, p)
  }

  sam <- gibbsC(x, y, m, n, p, lambda1, lambda2, b0, intercept, K)
  beta <- sam$beta
  tau <- sam$tau
  omega <- sam$omega
  
  if(mll) {
    mll <- mllGibbs(y, x, beta, omega, tau)
    out <- list(beta=beta, tau=tau, omega=omega, mll=mll)
  } else {
    out <- list(beta=beta, tau=tau, omega=omega)
  }
  
  return(out)
  
}

# marginal likelihood calculation on Gibbs samples
mllGibbs <- function(y, x, beta, omega, tau) {
  n <- nrow(x)
  K <- nrow(beta)
  bstar <- apply(beta, 2, function(b) {d <- density(b); return(d$x[which.max(d$y)])})
  ostar <- apply(omega, 2, function(o) {d <- density(o); return(d$x[which.max(d$y)])})
  
  ll <- -n*log(2) - 0.5*t(bstar) %*% t(x) %*% diag(ostar) %*% x %*% bstar + 
    t(y) %*% x %*% bstar - 0.5*sum(x %*% bstar)
  pom <- sum(dpg(ostar, 1, 0, log.d=TRUE)) - 
    sum(dpg(ostar, 1, as.numeric(x %*% bstar), log.d=TRUE))
  pbe <- sum(denbeta(bstar, lambda1, lambda2, log.d=TRUE))
  pgi <- log(mean(sapply(1:K, function(k) {
    S <- solve(t(x) %*% diag(omega[k, ]) %*% x + lambda2*diag(tau[k, ]/(tau[k, ] - 1))); 
    m <- S %*% t(x) %*% (y - 0.5); return(dmvnorm(bstar, mean=S %*% m, sigma=S))})))
  mll <- ll + pom + pbe - pgi
  return(mll)
  
}

# Polya-gamma density
dpg <- function(x, b, c, log.d=FALSE) {
  n <- 0:10000
  if(b==1) {
    sump <- sapply(x, function(x) {
      sum((-1)^n*(2*n + 1)*exp(-(2*n + 1)^2/(8*x)))})
    
    const <- cosh(c/2)*exp(-0.5*c^2*x)/sqrt(2*pi*x^3)
  } else {
    sump <- sapply(x, function(x) {
      sum((-1)^n*gamma(n + b)*(2*n + b)*exp(-(2*n + b)^2/(8*x))/gamma(n + 1))})
    if(all(c==0)) {
      const <- 1/sqrt(2*pi*x^3)
    } else {
      const <- cosh(c/2)^b*2^(b - 1)*exp(-0.5*c^2*x)/(gamma(b)*sqrt(2*pi*x^3))
    }  
  }
  dens <- ifelse(log.d, log(const) + log(sump), const*sump)
  return(dens)
}

# elastic net prior density
denbeta <- function(x, lambda1, lambda2, log.d=FALSE) {
  if(log.d) {
    C <- 0.5*log(lambda2) + dnorm(lambda1/(2*sqrt(lambda2)), log=TRUE) - log(2) - 
      pnorm(-lambda1/(2*sqrt(lambda2)), log.p=TRUE)
    dens <- C - 0.5*(lambda1*abs(x) + lambda2*x^2)
  } else {
    C <- sqrt(lambda2)*dnorm(lambda1/(2*sqrt(lambda2)))/(2*pnorm(-lambda1/(2*sqrt(lambda2))))
    dens <- C*exp(-0.5*(lambda1*abs(x) + lambda2*x^2))
  }
  return(dens)
}

### simulation
set.seed(123)
n <- 50
p <- 12
x <- matrix(rnorm(n*p), ncol=p, nrow=n)
lambda1 <- 1
lambda2 <- 1
beta <- renbeta(p, lambda1=lambda1, lambda2=lambda2)
m <- rep(1, n)
y <- rbinom(n, m, exp(x %*% beta)/(1 + exp(x %*% beta)))

seq.lambda2 <- c(seq(0.05, 0.45, 0.05), seq(0.5, 16, 0.5))[1]
out.mll <- matrix(NA, ncol=2, nrow=length(seq.lambda2))
for(l in 1:length(seq.lambda2)) {
  fit.VB <- ENVB(y, x, lambda1=1, lambda2=seq.lambda2[l], b0=beta, S0=diag(p), eps=1e-06, 
                 maxiter=500)
  fit.Gibbs <- ENgibbs(x, y, m, lambda1=1, lambda2=seq.lambda2[l], b0=c(0, beta), intercept=TRUE, 
                       K=50000, mll=TRUE)
  fit.pen <- penalized(y, x, unpenalized=~0, lambda1=1, lambda2=seq.lambda2[l], model="logistic")
  out.mll[l, ] <- c(fit.VB$mll, fit.Gibbs$mll)
}


plot(seq.lambda2, (out.mll[, 2] - min(out.mll[, 2]))/(max(out.mll[, 2]) - min(out.mll[, 2])), 
     type="l")
lines(seq.lambda2, (out.mll[, 1] - min(out.mll[, 1]))/(max(out.mll[, 1]) - min(out.mll[, 1])))


bgibbs <- apply(fit.Gibbs$beta[-c(1:5000), ], 2, function(b) {
  d <- density(b); return(d$x[which.max(d$y)])})
windows()
par(mfrow=c(3, 4))
for(j in 1:p) {
  h <- hist(fit.Gibbs$beta[-c(1:5000), j], plot=FALSE, breaks=80)
  hist(fit.Gibbs$beta[-c(1:5000), j], freq=FALSE, breaks=80, col=128,
       ylim=c(0, max(c(h$density, dnorm(fit.VB$mu[j], mean=fit.VB$mu[j], 
                                        sd=sqrt(fit.VB$sigma[j, j]))))), 
       xlab=expression(beta), main=paste(letters[j], ")", sep=""))  
  curve(dnorm(x, mean=fit.VB$mu[j], sd=sqrt(fit.VB$sigma[j, j])), add=TRUE, from=h$breaks[1] , 
        to=h$breaks[length(h$breaks)], col=2)
  abline(v=fit.pen@penalized[j], col=4, lty=2)
  abline(v=fit.VB$mu[j], col=2, lty=2)
  abline(v=bgibbs[j], col=1, lty=2)
}
par(mfrow=c(1, 1))

windows()
par(mfrow=c(3, 4))
for(j in 1:p) {
  chi <- seq.lambda2[l]*(fit.VB$sigma[j, j] + fit.VB$mu[j]^2)
  psi <- lambda1^2/(4*seq.lambda2[l])
  h <- hist(fit.Gibbs$tau[-c(1:5000), j] - 1, plot=FALSE, breaks=80)  
  hist(fit.Gibbs$tau[-c(1:5000), j] - 1, freq=FALSE, breaks=80, col=128,
       ylim=c(0, max(c(h$density, (-0.5 + sqrt(0.25 + chi*psi))/chi))), xlab=expression(psi), 
       main=paste(letters[j], ")", sep=""))  
  curve(dgig(x, chi=chi, psi=psi, lambda=0.5), add=TRUE, from=h$breaks[1] , 
        to=h$breaks[length(h$breaks)], col=2)
}
par(mfrow=c(1, 1))








### data

