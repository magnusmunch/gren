##############################  preamble  #############################
# Gibbs sampler implemented in R                                      #
# version: 02                                                         #
# author: Magnus Münch                                                #
# created: 06-12-2016                                                 #
# last edited: 06-12-2016                                             #
#######################################################################

###############################  notes  ###############################
# 06-12-2016: Can include an intercept                                #
#######################################################################

### paths
path.res <- "C:/Users/Magnus/Documents/phd/review/results/"
path.code <- "C:/Users/Magnus/Documents/phd/ENVB/code/"
path.graph <- "C:/Users/Magnus/Documents/phd/ENVB/graphs/"

### libraries
library(mvtnorm)
library(BayesLogit)
library(penalized)
library(Rcpp)
library(pROC)

### functions
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

# source function to estimate VB parameters in C++
sourceCpp(paste(path.code, "ENVB2.cpp", sep=""))

# sampler for omega from PG
romega <- function(beta, x, m, intercept) {
  n <- nrow(x)
  if(intercept) {
    x <- cbind(1, x)
  }
  z <- x %*% as.matrix(beta)
  omega <- rpg(n, m, z)
  return(omega)
}

# sample tau from GIG
rtau <- function(beta, lambda1, lambda2, intercept) {
  if(intercept) {
    beta <- beta[-1]
  }
  psi <- lambda1^2/(4*lambda2)
  chi <- lambda2*beta^2
  tau <- Inf
  while(!all(is.finite(tau))) {
    u <- runif(p)
    y <- rnorm(p)
    z <- sqrt(psi/chi) + 0.5*y^2/chi - sqrt(sqrt(psi)*y^2/(chi^1.5) + 0.25*y^4/(chi^2))
    tau <- ifelse(u <= 1/(1 + z*sqrt(chi/psi)), 1/z + 1, chi*z/psi + 1)
  }
  return(tau)
}

# sample beta from MVN
rbeta <- function(x, kappa, tau, omega, lambda2, intercept) {
  n <- nrow(x)
  p <- ncol(x)
  xtr <- t(x)
  hinv <- 1/lambda2 - 1/(lambda2*tau)
  xtrhinv <- xtr*hinv
  om <- omega
  ominv <- 1/om
  sigma <- matrix(0, nrow=p + intercept, ncol=p + intercept)
  if(intercept) {
    invtrom <- 1/sum(om)
    Omadj <- diag(om) - invtrom*(as.matrix(om) %*% t(as.matrix(om)))
    if(p >= n) {
      xtrxom <- x %*% xtrhinv %*% Omadj + diag(n)
      Ainv <- diag(hinv) - xtrhinv %*% Omadj %*% solve(xtrxom) %*% t(xtrhinv)
    } else {
      Ainv <- solve(xtr %*% Omadj %*% x + diag(lambda2 + lambda2*tau/(tau - 1)))
    }
    Ainvxom <- Ainv %*% xtr %*% as.matrix(om)
    sigma[1, 1] <- invtrom + invtrom^2*(t(as.matrix(om)) %*% x %*% Ainvxom)
    sigma[2:(p + 1), 1] <- sigma[1, 2:(p + 1)] <- -invtrom*Ainvxom
    sigma[2:(p + 1), 2:(p + 1)] <- Ainv
    mu <- sigma %*% rbind(1, xtr) %*% as.matrix(kappa)
  } else {
    if(p >= n) {
      xtrxom <- x %*% xtrhinv + diag(ominv)
      sigma <- diag(hinv) - xtrhinv %*% solve(xtrxom) %*% t(xtrhinv)
    } else {
      sigma <- solve(t(x*om) %*% x + diag(lambda2 + lambda2*tau/(tau - 1)))
    }
    mu <- sigma %*% xtr %*% as.matrix(kappa)
  }
  beta0 <- rnorm(p + intercept)
  beta <- chol(sigma) %*% beta0 + mu
  #beta <- as.numeric(rmvnorm(1, mean=mu, sigma=sigma))
  return(beta)
}

# the Gibbs sampler
ENgibbs <- function(x, y, m, lambda1, lambda2, b0, intercept, K, trace) {
  
  n <- nrow(x)
  p <- ncol(x)
  if(length(lambda2)==1) {
    lambda2 <- rep(lambda2, p)
  }
  beta <- b0
  kappa <- y - m/2
  
  seq.beta <- matrix(NA, ncol=K, nrow=p + intercept)
  seq.tau <- matrix(NA, ncol=K, nrow=p)
  seq.omega <- matrix(NA, ncol=K, nrow=n)
  k <- 0
  while(k < K) {
    k <- k + 1
    if(trace) {cat("\r", k)}
    omega <- seq.omega[, k] <- romega(beta, x, m, intercept)
    tau <- seq.tau[, k] <- rtau(beta, lambda1, lambda2, intercept)
    beta <- seq.beta[, k] <- rbeta(x, kappa, tau, omega, lambda2, intercept)
  }
  
  out <- list(beta=seq.beta, tau=seq.tau, omega=seq.omega)
  return(out)
  
}

# ENVB 
ENVB <- function(x, y, m, lambda1, lambda2, mu0, sigma0, intercept, maxiter, epsilon=1e-07, trace) {
  
  p <- ncol(x)
  n <- nrow(x)
  if(length(lambda2)==1) {
    lambda2 <- rep(lambda2, p)
  }
  phi <- lambda1^2/(4*lambda2)
  kappa <- y - m/2
  
  conv <- FALSE
  niter <- 0
  
  mu <- muold <- mu0
  sigma <- sigmaold <- sigma0
  if(intercept) {
    xaug <- cbind(1, x)
    ciold <- as.numeric(sqrt(colSums(t(xaug) * (sigma %*% t(xaug))) + (colSums(t(xaug)*mu))^2))
    chiold <- as.numeric(lambda2*(diag(sigma)[2:(p + 1)] + mu[2:(p + 1)]^2))
  } else {
    ciold <- as.numeric(sqrt(colSums(t(x) * (sigma %*% t(x))) + (colSums(t(x)*mu))^2))
    chiold <- as.numeric(lambda2*(diag(sigma) + mu^2))
  }
  
  while(!conv & (niter < maxiter)) {
    
    niter <- niter + 1
    if(trace) {cat("\r", niter)}
    
    # estimation of new model parameter (done in Cpp)
    new.param <- est_param(x, kappa, m, n, p, ciold, phi, chiold, lambda2, intercept)
    sigma <- new.param$sigma
    mu <- as.numeric(new.param$mu)
    ci <- as.numeric(new.param$ci)
    chi <- as.numeric(new.param$chi)
    
    # check the convergence of the model parameters
    conv <- max(abs(c(diag(sigma) - diag(sigmaold), mu - muold))) < epsilon
    
    # update old parameters to new ones
    sigmaold <- sigma
    muold <- mu
    ciold <- ci
    chiold <- chi
    
  }
  
  out <- list(niter=niter, conv=conv, sigma=sigma, mu=mu, c=ci, chi=chi)
  return(out)
  
}

### simulation
set.seed(123)
n <- 50
p <- 11
x <- matrix(rnorm(n*p), ncol=p, nrow=n)
lambda1 <- 1
lambda2 <- 1
beta <- renbeta(p, lambda1=lambda1, lambda2=lambda2)
m <- rep(1, n)
y <- rbinom(n, m, exp(x %*% beta)/(1 + exp(x %*% beta)))

fitGibbs <- ENgibbs(x, y, m, lambda1, lambda2, b0=c(0, beta), intercept=TRUE, K=30000, trace=TRUE)
fitVB <- ENVB(x, y, m, lambda1, lambda2, mu0=c(0, beta), sigma0=diag(p + 1), intercept=TRUE, maxiter=100, 
              epsilon=1e-07, trace=TRUE)
fitPen <- penalized(y, x, unpenalized=~1, lambda1=lambda1, lambda2=lambda2, model="logistic")
b.gibbs <- apply(fitGibbs$beta[, -c(1:5000)], 1, function(b) {d <- density(b); return(d$x[which.max(d$y)])})

hist(fitGibbs$beta[1, ], breaks=50)
plot(b.gibbs, fitVB$mu)
plot(b.gibbs, c(fitPen@unpenalized, fitPen@penalized))
plot(c(0, beta), b.gibbs)
plot(c(0, beta), fitVB$mu)
plot(c(0, beta), c(fitPen@unpenalized, fitPen@penalized))

legend.names <- c(expression(beta), expression(hat(beta)^"MCMC"), expression(hat(beta)^"VB"), 
                  expression(hat(beta)^"EN"))
png(paste(path.graph, "VB_MCMC_sim_V01.png", sep=""), width=1920, height=1200)
m <- matrix(c(1:12, rep(13, 4)), nrow=4, ncol=4, byrow=TRUE)
layout(mat=m, heights=c(0.4, 0.4, 0.4, 0.2))
for(j in 1:(p + 1)) {
  h <- hist(fitGibbs$beta[j, -c(1:5000)], plot=FALSE, breaks=80)
  hist(fitGibbs$beta[j, -c(1:5000)], freq=FALSE, breaks=80, col=128,
       ylim=c(0, max(c(h$density, dnorm(fitVB$mu[j], mean=fitVB$mu[j], 
                                        sd=sqrt(fitVB$sigma[j, j]))))), 
       xlab=expression(beta), main=paste(letters[j], ")", sep=""), cex.lab=1.8, cex.main=2, cex.axis=1.3)  
  curve(dnorm(x, mean=fitVB$mu[j], sd=sqrt(fitVB$sigma[j, j])), add=TRUE, from=h$breaks[1] , 
        to=h$breaks[length(h$breaks)], col=2)
  abline(v=c(fitPen@unpenalized, fitPen@penalized)[j], col=4, lty=2)
  abline(v=fitVB$mu[j], col=2, lty=2)
  abline(v=b.gibbs[j], col=1, lty=2)
  abline(v=c(0, beta)[j], col=3, lty=1)
}
par(mar=c(0, 0, 0, 0))
plot(1, type="n", axes=FALSE, xlab="", ylab="")
legend(x="top", inset=0, legend=legend.names, col=c(3, 1, 2, 4), lwd=2, lty=c(1, 2, 2, 2), horiz=TRUE,
       cex=2)
dev.off()




### data
set.seed(123)
data <- read.table("http://www-stat.stanford.edu/~tibs/ElemStatLearn/datasets/SAheart.data", sep=",", head=T,
                   row.names=1)
data$famhist <- ifelse(data$famhist=="Absent", 0, 1)
n <- nrow(data)
ntest <- floor(n*(1/3))
ntrain <- n - ntest
indtest <- sample(1:n, ntest)
ytrain <- data$chd[-indtest]
ytest <- data$chd[indtest]
xtrain <- as.matrix(data[-indtest, names(data)!="chd"])
xtest <- as.matrix(data[indtest, names(data)!="chd"])

p <- ncol(xtrain)
mtrain <- rep(1, n - ntest)
mtest <- rep(1, ntest)
lambda1 <- 1
lambda2 <- 1

fitGibbs2 <- ENgibbs(xtrain, ytrain, mtrain, lambda1, lambda2, b0=c(0, rnorm(p)), intercept=TRUE, K=50000, trace=TRUE)
fitVB2 <- ENVB(xtrain, ytrain, mtrain, lambda1, lambda2, mu0=c(0, rnorm(p)), sigma0=diag(p + 1), intercept=TRUE, 
               maxiter=10000, epsilon=1e-07, trace=TRUE)
fitPen2 <- penalized(ytrain, xtrain, unpenalized=~1, lambda1=lambda1, lambda2=lambda2, model="logistic")
fitGlm2 <- glm(ytrain ~ xtrain, famil="binomial")
b.gibbs2 <- apply(fitGibbs2$beta[, -c(1:5000)], 1, function(b) {d <- density(b); return(d$x[which.max(d$y)])})


hist(fitGibbs2$beta[1, ], breaks=50)
plot(b.gibbs2, fitVB2$mu)
plot(b.gibbs2, c(fitPen2@unpenalized, fitPen2@penalized))
plot(b.gibbs2, coef(fitGlm2))

rocGibbs2 <- roc(ytest, as.numeric(exp(cbind(1, xtest) %*% b.gibbs2)/(1 + exp(cbind(1, xtest) %*% b.gibbs2))))
rocVB2 <- roc(ytest, as.numeric(exp(cbind(1, xtest) %*% fitVB2$mu)/(1 + exp(cbind(1, xtest) %*% fitVB2$mu))))
rocPen2 <- roc(ytest, predict(fitPen2, xtest))
rocGlm2 <- roc(ytest, as.numeric(exp(cbind(1, xtest) %*% coef(fitGlm2))/(1 + exp(cbind(1, xtest) %*% coef(fitGlm2)))))
auc(rocGibbs2)
auc(rocVB2)
auc(rocPen2)
auc(rocGlm2)

titles <- c("intercept", "sbp", "tobacco", "ldl", "adiposity", "famhist", "typea", "obesity", "alcohol",
          "age")
legend.names <- c(expression(hat(beta)^"MCMC"), expression(hat(beta)^"VB"), 
                  expression(hat(beta)^"EN"), expression(hat(beta)^"GLM"))
png(paste(path.graph, "VB_MCMC_data_V02.png", sep=""), width=1920, height=1200)
m <- matrix(c(1:12, rep(13, 4)), nrow=4, ncol=4, byrow=TRUE)
layout(mat=m, heights=c(0.4, 0.4, 0.4, 0.2))
for(j in 1:(p + 1)) {
  h <- hist(fitGibbs2$beta[j, -c(1:5000)], plot=FALSE, breaks=80)
  hist(fitGibbs2$beta[j, -c(1:5000)], freq=FALSE, breaks=80, col=128,
       ylim=c(0, max(c(h$density, dnorm(fitVB2$mu[j], mean=fitVB2$mu[j], 
                                        sd=sqrt(fitVB2$sigma[j, j]))))), 
       xlab=expression(beta), main=titles[j], cex.lab=1.8, cex.main=2, cex.axis=1.3) 
       #xlab=expression(beta), main=paste(letters[j], ")", sep=""), cex.lab=1.8, cex.main=2, cex.axis=1.3)  
  curve(dnorm(x, mean=fitVB2$mu[j], sd=sqrt(fitVB2$sigma[j, j])), add=TRUE, from=h$breaks[1] , 
        to=h$breaks[length(h$breaks)], col=2)
  abline(v=c(fitPen2@unpenalized, fitPen2@penalized)[j], col=4, lty=2)
  abline(v=fitVB2$mu[j], col=2, lty=2)
  abline(v=b.gibbs2[j], col=1, lty=2)
  abline(v=coef(fitGlm2)[j], col=3, lty=2)
}
plot(1, type="n", axes=FALSE, xlab="", ylab="")
plot(1, type="n", axes=FALSE, xlab="", ylab="")
par(mar=c(0, 0, 0, 0))
plot(1, type="n", axes=FALSE, xlab="", ylab="")
legend(x="top", inset=0, legend=legend.names, col=c(1, 2, 4, 3), lwd=2, lty=c(2, 2, 2, 2), horiz=TRUE,
       cex=2)
dev.off()




test <- scan("http://archive.ics.uci.edu/ml/machine-learning-databases/spectrometer/lrs.data", 
             what="character")
test <- gsub(")", "", test)
test <- test[test!="("]
length(test)/103
sapply(1:531, function(i) {})







