##############################  preamble  #############################
# MCEM implemented in R                                               #
# version: 01                                                         #
# author: Magnus Münch                                                #
# created: 03-02-2016                                                 #
# last edited: 07-02-2016                                             #
#######################################################################

###############################  notes  ###############################
# 01-03-2017: This only optimises one lambda2, compares VBEM, MCEM,   #
#             CV and MoM                                              #
#######################################################################

### paths
path.res <- "C:/Users/Magnus/Documents/phd/ENVB/results/"
path.code <- "C:/Users/Magnus/Documents/phd/ENVB/code/"
path.graph <- "C:/Users/Magnus/Documents/phd/ENVB/graphs/"

### libraries
library(Rcpp)
library(penalized)
library(glmnet)
library(GRridge)

### functions
# source function to sample from posterior
sourceCpp(paste(path.code, "ENGibbs.cpp", sep=""))
sourceCpp(paste(path.code, "ENVB2.cpp", sep=""))

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

mll2 <- function(lambda2, sizes, K, lambda1, sum1, sum2) {
  return(-sum(sizes*pnorm(-lambda1/sqrt(4*lambda2), log.p=TRUE)) - 0.5*sum(lambda2*sum1)/K - 
           lambda1^2*sum(sum2/lambda2)/(8*K))
}

lowermll <- function(lambda2, p, lambda1, sum1, sum2) {
  return(-p*pnorm(-lambda1/sqrt(4*lambda2), log.p=TRUE) - 0.5*lambda2*sum1 - lambda1^2*sum2/(8*lambda2))
}

lowermll2 <- function(lambda2, sizes, lambda1, sum1, sum2) {
  return(-sum(sizes*pnorm(-lambda1/sqrt(4*lambda2), log.p=TRUE)) - 0.5*sum(lambda2*sum1) - 
           lambda1^2*sum(sum2/lambda2)/(8))
}

grmagn2 <- function(par, lambda2, sizes, sum1) {
  s <- par[1]
  loglambdag <- par[-1]
  magn <- sum((0.5*lambda2*sum1*exp(loglambdag) - 0.5*sizes + s*sizes)^2) + sum(sizes*loglambdag)^2
  return(magn)
}

grmagn1 <- function(par, lambda2, sizes, sum1) {
  s <- par[1]
  loglambdag <- par[-1]
  magn <- sum((0.5*lambda2*sum1 - 0.5*sizes/lambdag + s*sizes/lambdag)^2) + 
    sum(sizes*log(lambdag))^2
  return(magn)
}

VBEM <- function(x, y, m, n, p, lambda1, lambda2, sigma0, intercept, eps, maxiter) {
  
  # assigning fixed (throughout algorithm) variables
  p <- ncol(x)
  n <- nrow(x)
  kappa <- y - m/2
  lambda <- 0.5
  
  # starting values for lambda2
  lambda2old <- lambda2seq <- lambda2
  if(intercept) {
    xadj <- cbind(1, x)
    lambda2vec <- c(0, rep(lambda2old, p))
  } else {
    xadj <- x
    lambda2vec <- rep(lambda2old, p)
  }
  
  # starting value for model parameters
  sigmaold <- sigma0
  muold <- as.numeric(sigmaold %*% t(xadj) %*% as.matrix(kappa))
  phi <- rep(0.25*lambda1^2/lambda2old, p)
  ci <- as.numeric(sqrt(colSums(t(xadj) * (sigmaold %*% t(xadj))) + (colSums(t(xadj)*muold))^2))
  chi <- as.numeric(lambda2vec*(diag(sigmaold) + muold^2))[(intercept + 1):(intercept + p)]
  
  # sums are needed in optimisation routine
  if(intercept) {
    sum1 <- sum((diag(sigmaold)[-1] + muold[-1]^2)*(1 + sqrt(phi/chi)))
  } else {
    sum1 <- sum((diag(sigmaold) + muold^2)*(1 + sqrt(phi/chi)))
  }
  sum2 <- sum(1/phi + sqrt(chi/phi) + 1)
  
  # keeping track of things:
  lowermllseq <- lowermll(lambda2, p, lambda1, sum1, sum2)
  niter2seq <- vector(mode="list", length=0)
  
  # outer loop of algorithm:
  conv1 <- FALSE
  iter1 <- 0
  while(!conv1 & (iter1 < maxiter)) {
    iter1 <- iter1 + 1
    
    # sums are needed in optimisation routine
    if(intercept) {
      sum1 <- sum((diag(sigmaold)[-1] + muold[-1]^2)*(1 + sqrt(phi/chi)))
    } else {
      sum1 <- sum((diag(sigmaold) + muold^2)*(1 + sqrt(phi/chi)))
    }
    sum2 <- sum(1/phi + sqrt(chi/phi) + 1)
    
    # estimating new lambda2
    opt <- optim(par=lambda2old, fn=lowermll, p=p, lambda1=lambda1, sum1=sum1, sum2=sum2, 
                 control=list(fnscale=-1), method="Brent", lower=0.00001, upper=10000)
    lambda2 <- opt$par
    lowermllseq <- c(lowermllseq, opt$value)
    lambda2seq <- c(lambda2seq, lambda2)
    
    # inner loop of algorithm:
    conv2 <- 0
    iter2 <- 0
    while(!conv2 & (iter2 < maxiter)) {
      iter2 <- iter2 + 1
      
      # estimating new model parameters
      newparam <- est_param(x, kappa, m, n, p, ci, phi, chi, rep(lambda2old, p), intercept)
      sigma <- newparam$sigma
      mu <- newparam$mu
      ci <- newparam$ci
      chi <- newparam$chi
      phi <- rep(0.25*lambda1^2/lambda2, p)
      
      # checking convergence of inner loop
      conv2 <- max(c(abs((mu - muold)/muold), abs(diag((sigma - sigmaold)/sigmaold)))) < eps
      
      muold <- mu
      sigmaold <- sigma
    }
    niter2seq[[iter1]] <- iter2
    
    # checking convergence of outer loop:
    conv1 <- abs((lambda2 - lambda2old)/lambda2old) < eps
    
    lambda2old <- lambda2
  }
  
  out <- list(mu=mu, sigma=sigma, ci=ci, chi=chi, phi=phi, lambda2=lambda2seq, 
              lowermll=lowermllseq, nouteriter=iter1, ninneriter=niter2seq)
  return(out)
  
}

grVBEM1 <- function(x, y, m, groups, lambda1, lambda2, sigma0, intercept, eps, maxiter) {
  
  # assigning fixed (throughout algorithm) variables
  sizes <- rle(groups)$lengths
  G <- length(unique(groups))
  p <- ncol(x)
  n <- nrow(x)
  kappa <- y - m/2
  phi <- 0.25*lambda1^2/lambda2
  
  # starting values for lambdag, s and v
  lambdagold <- lambdagseq <- rep(1, G)
  if(intercept) {
    xadj <- cbind(1, x)
    lambdagvec <- c(0, rep(lambdagold, sizes))
  } else {
    xadj <- x
    lambdagvec <- rep(lambdagold, sizes)
  }
  s <- 0
  #v <- 1
  
  # starting value for model parameters
  sigmaold <- sigma0
  muold <- as.numeric(sigmaold %*% t(xadj) %*% as.matrix(kappa))
  ci <- as.numeric(sqrt(colSums(t(xadj) * (sigmaold %*% t(xadj))) + (colSums(t(xadj)*muold))^2))
  chi <- as.numeric(lambda2*lambdagvec*(diag(sigmaold) + muold^2))[(intercept + 1):(intercept + p)]
  
  # sum is needed in optimisation routine
  sum1 <- sapply(1:G, function(g) {
    ind <- which(groups==g) + intercept; 
    return(sum((diag(sigmaold)[ind] + muold[ind]^2)*(1 + sqrt(phi/chi[ind - intercept]))))})
  
  # keeping track of things:
  lowermllseq <- 0.5*sum(sizes*log(lambdag)) - 0.5*lambda2*sum(lambdag*sum1)
  niter2seq <- vector(mode="list", length=0)
  
  # outer loop of algorithm:
  conv1 <- FALSE
  iter1 <- 0
  while(!conv1 & (iter1 < maxiter)) {
    iter1 <- iter1 + 1
    
    # estimating new lambdag
    # opt <- optim(par=c(s, v, lambdagold), fn=grmagn1, lambda2=lambda2, sum1=sum1, 
    #              p=p, sizes=sizes, method="Nelder-Mead")
    # s <- opt$par[1]
    # v <- opt$par[2]
    # lambdag <- opt$par[-c(1, 2)]
    # lambdagseq <- cbind(lambdagseq, lambdag)
    opt <- optim(par=c(s, log(lambdagold)), fn=grmagn2, lambda2=lambda2, sizes=sizes, sum1=sum1,
                 method="Nelder-Mead", control=list(maxit=1000))
    s <- opt$par[1]
    lambdag <- exp(opt$par[-1])
    lambdagseq <- cbind(lambdagseq, lambdag)
    
    # inner loop of algorithm:
    conv2 <- 0
    iter2 <- 0
    while(!conv2 & (iter2 < maxiter)) {
      iter2 <- iter2 + 1
      
      # estimating new model parameters
      newparam <- est_param(x, kappa, m, n, p, ci, rep(phi, p), chi, rep(lambdag, sizes), intercept)
      sigma <- newparam$sigma
      mu <- newparam$mu
      ci <- newparam$ci
      chi <- newparam$chi
      
      # checking convergence of inner loop
      conv2 <- max(c(abs((mu - muold)/muold), abs(diag((sigma - sigmaold)/sigmaold)))) < eps
      
      muold <- mu
      sigmaold <- sigma
    }
    niter2seq[[iter1]] <- iter2
    
    # sum is needed in optimisation routine
    sum1 <- sapply(1:G, function(g) {
      ind <- which(groups==g) + intercept; 
      return(sum((diag(sigmaold)[ind] + muold[ind]^2)*(1 + sqrt(phi/chi[ind - intercept]))))})
    
    # keeping track of lower bound on marginal log likelihood
    lowermllseq <- c(lowermllseq, 0.5*sum(sizes*log(lambdag)) - 0.5*lambda2*sum(lambdag*sum1))
    
    # checking convergence of outer loop:
    conv1 <- max(abs((lambdag - lambdagold)/lambdagold)) < eps
    
    # updating lambdag for new iteration
    lambdagold <- lambdag
  }
  
  out <- list(mu=mu, sigma=sigma, ci=ci, chi=chi, lambdag=lambdagseq, 
              lowermll=lowermllseq, nouteriter=iter1, ninneriter=niter2seq, conv=conv1)
  return(out)
  
}

MCEM <- function(x, y, m, n, p, lambda1, lambda2, b0, intercept, K, burnin, eps, maxiter) {
  
  sam <- gibbsC(x, y, m, n, p, lambda1, rep(lambda2, p), b0, TRUE, K + burnin)
  beta <- apply(sam$beta[, -c(1:burnin)], 1, function(b) {
    d <- density(b); return(d$x[which.max(d$y)])})
  if(intercept) {
    sum1 <- sum(sam$beta[-1, -c(1:burnin)]^2*sam$tau[, -c(1:burnin)]/(sam$tau[, -c(1:burnin)] - 1))
  } else {
    sum1 <- sum(sam$beta[, -c(1:burnin)]^2*sam$tau[, -c(1:burnin)]/(sam$tau[, -c(1:burnin)] - 1))
  }
  sum2 <- sum(sam$tau[, -c(1:burnin)])
  
  mllseq <- mll(lambda2, p, K, lambda1, sum1, sum2)
  lambda2old <- lambda2seq <- lambda2
  
  conv <- FALSE
  iter <- 0
  while(!conv & (iter < maxiter)) {
    iter <- iter + 1
    print(paste("iteration ", iter, ", current lambda2 is ", round(lambda2, 2), sep=""))
    
    if(intercept) {
      sum1 <- sum(sam$beta[-1, -c(1:3000)]^2*sam$tau[, -c(1:3000)]/(sam$tau[, -c(1:3000)] - 1))
    } else {
      sum1 <- sum(sam$beta[, -c(1:3000)]^2*sam$tau[, -c(1:3000)]/(sam$tau[, -c(1:3000)] - 1))
    }
    sum2 <- sum(sam$tau[, -c(1:3000)])
    
    opt <- optim(par=lambda2old, fn=mll, p=p, K=K, lambda1=lambda1, sum1=sum1, sum2=sum2, control=list(fnscale=-1), 
                 method="Brent", lower=0.00001, upper=10000)
    if(!is.finite(opt$value)) break
    lambda2 <- opt$par
    mllseq <- c(mllseq, opt$value)
    lambda2seq <- c(lambda2seq, lambda2)
    
    conv <- (abs(lambda2 - lambda2old) < eps) #| (abs((mllseq[iter + 1] - mllseq[iter])/mllseq[iter]) < 0.01)
    
    sam <- gibbsC(x, y, m, n, p, lambda1, rep(lambda2, p), beta, TRUE, K + burnin)
    beta <- apply(sam$beta[, -c(1:burnin)], 1, function(b) {
      d <- density(b); return(d$x[which.max(d$y)])})
    
    lambda2old <- lambda2
    
  }
  
  out <- list(beta=sam$beta[, -c(1:burnin)], tau=sam$tau[, -c(1:burnin)], 
              omega=sam$omega[, -c(1:burnin)], lambda2=lambda2seq, mll=mllseq)
  return(out)
  
}

grridge_oneg <- function(x, y, lambda2, intercept, maxiter=10, eps=1e-10) {
  p <- ncol(x)
  lambda2old <- lambda2seq <- lambda2
  if(intercept) {
    xadj <- cbind(1, x) 
    p0 <- p + 1
  } else {
    xadj <- x
    p0 <- p
  }
  conv <- FALSE
  iter <- 0
  while(!conv & (iter < maxiter)) {
    iter <- iter + 1
    if(intercept) {
      fit <- penalized(y, x, unpenalized=~1, lambda1=0, lambda2=lambda2old, model="logistic", trace=FALSE)  
    } else {
      fit <- penalized(y, x, unpenalized=~0, lambda1=0, lambda2=lambda2old, model="logistic", trace=FALSE)  
    }
    best <- coef(fit, which="all")
    linpred <- xadj %*% best
    pred <- as.numeric(exp(linpred)/(1 + linpred))
    # svdxwx <- svd(t(xadj) %*% (xadj*(pred*(1 - pred))))
    # d <- svdxwx$d
    # V <- svdxwx$v
    # ckl <- V %*% (t(V)*(d^2/(d^2 + 2*lambda2old)))
    # uk <- diag(V %*% (t(V)*(d^2/(d^2 + 2*lambda2old)^2)))
    # dkl <- ckl/sqrt(uk)
    # agh <- sum(dkl^2)
    if(intercept) {
      lambda2vec <- c(0, rep(lambda2old, p))
    } else {
      lambda2vec <- rep(lambda2old, p)
    }
    matinv <- solve(t(xadj) %*% xadj + 2*diag(lambda2vec))
    ckl <- (matinv %*% t(xadj) %*% xadj)[-1, -1]
    uk <- diag(matinv %*% t(xadj) %*% xadj %*% matinv)[-1]
    dkl <- diag(1/sqrt(uk)) %*% ckl
    agh <- sum(dkl^2)
    lambda2 <- agh/sum(best[-1]^2/uk - 1)
    conv <- abs(lambda2 - lambda2old) < eps
    lambda2old <- lambda2
    lambda2seq <- c(lambda2seq, lambda2)
  }
  if(intercept) {
    fit <- penalized(y, x, unpenalized=~1, lambda1=0, lambda2=lambda2, model="logistic", trace=FALSE)
  } else {
    fit <- penalized(y, x, unpenalized=~0, lambda1=0, lambda2=lambda2, model="logistic", trace=FALSE)
  }
  best <- coef(fit, which="all")
  out <- list(lambda2=lambda2seq, b=best, niter=iter, conv=conv)
}

### simulation 1 (lambda1=0.1)
set.seed(123)
n <- 100
p <- 50
x <- matrix(rnorm(n*p), ncol=p, nrow=n)
lambda1 <- 0.1
lambda2 <- c(1, 10, 50, 100)
m <- rep(1, n)
b0 <- rnorm(p + 1)
sigma0 <- diag(rchisq(p + 1, 1))

test1.MCEM <- test1.VBEM <- test1.optL2 <- test1.grridge_oneg <- vector(mode="list", length=4)
for(l2 in 1:length(lambda2)) {
  beta <- renbeta(p, lambda1=lambda1, lambda2=lambda2[l2])
  y <- rbinom(n, m, exp(x %*% beta)/(1 + exp(x %*% beta)))
  test1.MCEM[[l2]] <- MCEM(x, y, m, n, p, lambda1, lambda2[l2], b0, intercept=TRUE, K=5000, 
                          burnin=3000, eps=0.001, maxiter=1000)
  test1.VBEM[[l2]] <- VBEM(x, y, m, n, p, lambda1, lambda2[l2], sigma0, intercept=TRUE, eps=0.001, 
                          maxiter=1000)
  test1.optL2[[l2]] <- optL2(y, x, unpenalized=~1, lambda1=lambda1, model="logistic", fold=n)
  test1.grridge_oneg[[l2]] <- grridge_oneg(x, y, 1, intercept=TRUE, maxiter=100, eps=0.001)
}

#save(test1.MCEM, test1.VBEM, test1.optL2, test1.grridge_oneg, 
#     file=paste(path.res, "res_MCEM_lambda1=0.1_V01.Rdata", sep=""))

# lambda2.MCEM <- sapply(1:4, function(l2) {test1.MCEM[[l2]]$lambda2[length(test1.MCEM[[l2]]$lambda2)]})
# lambda2.VBEM <- sapply(1:4, function(l2) {test1.VBEM[[l2]]$lambda2[length(test1.VBEM[[l2]]$lambda2)]})
# lambda2.optL2 <- sapply(1:4, function(l2) {test1.optL2[[l2]]$lambda})
# lambda2.grridge_oneg <- sapply(1:4, function(l2) {
#   test1.grridge_oneg[[l2]]$lambda2[length(test1.grridge_oneg[[l2]]$lambda2)]})
# 
# opar <- par()
# par(mar=par('mar') + c(0, 1, 0, 0))
# plot(lambda2, lambda2.MCEM, xlab=expression(lambda[2]), ylab=expression(hat(lambda)[2]), 
#      ylim=range(c(lambda2.MCEM, lambda2.VBEM, lambda2.optL2, lambda2.grridge_oneg)), col=2, pch=19)
# points(lambda2, lambda2.VBEM, col=3, pch=19)
# points(lambda2, lambda2.optL2, col=4, pch=19)
# points(lambda2, lambda2.grridge_oneg, col=5, pch=19)
# abline(a=0, b=1, lty=2)
# legend("topright", legend=c("MCEM", "VBEM", "CV", "MoM"), pch=19, col=c(2:5))
# par(opar)

### simulation 2 (lambda1=1)
set.seed(123)
n <- 100
p <- 50
x <- matrix(rnorm(n*p), ncol=p, nrow=n)
lambda1 <- 1
lambda2 <- c(1, 10, 50, 100)
m <- rep(1, n)
b0 <- rnorm(p + 1)
sigma0 <- diag(rchisq(p + 1, 1))

test2.MCEM <- test2.VBEM <- test2.optL2 <- test2.grridge_oneg <- vector(mode="list", length=4)
for(l2 in 1:length(lambda2)) {
  beta <- renbeta(p, lambda1=lambda1, lambda2=lambda2[l2])
  y <- rbinom(n, m, exp(x %*% beta)/(1 + exp(x %*% beta)))
  test2.MCEM[[l2]] <- MCEM(x, y, m, n, p, lambda1, lambda2[l2], b0, intercept=TRUE, K=5000, 
                          burnin=3000, eps=0.001, maxiter=1000)
  test2.VBEM[[l2]] <- VBEM(x, y, m, n, p, lambda1, lambda2[l2], sigma0, intercept=TRUE, eps=0.001, 
                          maxiter=1000)
  test2.optL2[[l2]] <- optL2(y, x, unpenalized=~1, lambda1=lambda1, model="logistic", fold=n)
  test2.grridge_oneg[[l2]] <- grridge_oneg(x, y, 1, intercept=TRUE, maxiter=100, eps=0.001)
}

# save(test2.MCEM, test2.VBEM, test2.optL2, test2.grridge_oneg, 
#      file=paste(path.res, "res_MCEM_lambda1=1_V01.Rdata", sep=""))
# 
# lambda2.MCEM <- sapply(1:4, function(l2) {test2.MCEM[[l2]]$lambda2[length(test2.MCEM[[l2]]$lambda2)]})
# lambda2.VBEM <- sapply(1:4, function(l2) {test2.VBEM[[l2]]$lambda2[length(test2.VBEM[[l2]]$lambda2)]})
# lambda2.optL2 <- sapply(1:4, function(l2) {test2.optL2[[l2]]$lambda})
# lambda2.grridge_oneg <- sapply(1:4, function(l2) {
#   test2.grridge_oneg[[l2]]$lambda2[length(test2.grridge_oneg[[l2]]$lambda2)]})
# 
# opar <- par()
# par(mar=par('mar') + c(0, 1, 0, 0))
# plot(lambda2, lambda2.MCEM, xlab=expression(lambda[2]), ylab=expression(hat(lambda)[2]),
#      ylim=range(c(lambda2.MCEM, lambda2.VBEM, lambda2.optL2, lambda2.grridge_oneg)), col=2, pch=19)
# points(lambda2, lambda2.VBEM, col=3, pch=19)
# points(lambda2, lambda2.optL2, col=4, pch=19)
# points(lambda2, lambda2.grridge_oneg, col=5, pch=19)
# abline(a=0, b=1, lty=2)
# legend("topright", legend=c("MCEM", "VBEM", "CV", "MoM"), pch=19, col=c(2:5))
# par(opar)

### simulation 3 (lambda1=10)
set.seed(123)
n <- 100
p <- 50
x <- matrix(rnorm(n*p), ncol=p, nrow=n)
lambda1 <- 10
lambda2 <- c(1, 10, 50, 100)
m <- rep(1, n)
b0 <- rnorm(p + 1)
sigma0 <- diag(rchisq(p + 1, 1))

test3.MCEM <- test3.VBEM <- test3.optL2 <- test3.grridge_oneg <- vector(mode="list", length=4)
for(l2 in 1:length(lambda2)) {
  beta <- renbeta(p, lambda1=lambda1, lambda2=lambda2[l2])
  y <- rbinom(n, m, exp(x %*% beta)/(1 + exp(x %*% beta)))
  test3.MCEM[[l2]] <- MCEM(x, y, m, n, p, lambda1, lambda2[l2], b0, intercept=TRUE, K=5000, 
                           burnin=3000, eps=0.001, maxiter=1000)
  test3.VBEM[[l2]] <- VBEM(x, y, m, n, p, lambda1, lambda2[l2], sigma0, intercept=TRUE, eps=0.001, 
                           maxiter=1000)
  test3.optL2[[l2]] <- optL2(y, x, unpenalized=~1, lambda1=lambda1, model="logistic", fold=n)
  test3.grridge_oneg[[l2]] <- grridge_oneg(x, y, 1, intercept=TRUE, maxiter=100, eps=0.001)
}

# save(test3.MCEM, test3.VBEM, test3.optL2, test3.grridge_oneg,
#      file=paste(path.res, "res_MCEM_lambda1=10_V01.Rdata", sep=""))

# lambda2.MCEM <- sapply(1:4, function(l2) {test3.MCEM[[l2]]$lambda2[length(test3.MCEM[[l2]]$lambda2)]})
# lambda2.VBEM <- sapply(1:4, function(l2) {test3.VBEM[[l2]]$lambda2[length(test3.VBEM[[l2]]$lambda2)]})
# lambda2.optL2 <- sapply(1:4, function(l2) {test3.optL2[[l2]]$lambda})
# lambda2.optL2[2] <- NA
# lambda2.grridge_oneg <- sapply(1:4, function(l2) {
#   test3.grridge_oneg[[l2]]$lambda2[length(test3.grridge_oneg[[l2]]$lambda2)]})
# 
# opar <- par()
# par(mar=par('mar') + c(0, 1, 0, 0))
# plot(lambda2, lambda2.MCEM, xlab=expression(lambda[2]), ylab=expression(hat(lambda)[2]),
#      ylim=range(c(lambda2.MCEM, lambda2.optL2, lambda2.grridge_oneg), na.rm=TRUE), col=2, pch=19)
# #     ylim=range(c(lambda2.MCEM, lambda2.VBEM, lambda2.optL2, lambda2.grridge_oneg), na.rm=TRUE), col=2, pch=19)
# #points(lambda2, lambda2.VBEM, col=3, pch=19)
# points(lambda2, lambda2.optL2, col=4, pch=19)
# points(lambda2, lambda2.grridge_oneg, col=5, pch=19)
# abline(a=0, b=1, lty=2)
# legend("topright", legend=c("MCEM", "VBEM", "CV", "MoM"), pch=19, col=c(2:5))
# par(opar)

### simulation 4 (lambda1=5)
set.seed(123)
n <- 100
p <- 50
x <- matrix(rnorm(n*p), ncol=p, nrow=n)
lambda1 <- 5
lambda2 <- c(1, 10, 50, 100)
m <- rep(1, n)
b0 <- rnorm(p + 1)
sigma0 <- diag(rchisq(p + 1, 1))

test4.MCEM <- test4.VBEM <- test4.optL2 <- test4.grridge_oneg <- vector(mode="list", length=4)
for(l2 in 1:length(lambda2)) {
  beta <- renbeta(p, lambda1=lambda1, lambda2=lambda2[l2])
  y <- rbinom(n, m, exp(x %*% beta)/(1 + exp(x %*% beta)))
  test4.MCEM[[l2]] <- MCEM(x, y, m, n, p, lambda1, lambda2[l2], b0, intercept=TRUE, K=5000, 
                           burnin=3000, eps=0.001, maxiter=1000)
  test4.VBEM[[l2]] <- VBEM(x, y, m, n, p, lambda1, lambda2[l2], sigma0, intercept=TRUE, eps=0.001, 
                           maxiter=1000)
  test4.optL2[[l2]] <- optL2(y, x, unpenalized=~1, lambda1=lambda1, model="logistic", fold=n)
  test4.grridge_oneg[[l2]] <- grridge_oneg(x, y, 1, intercept=TRUE, maxiter=100, eps=0.001)
}

# save(test4.MCEM, test4.VBEM, test4.optL2, test4.grridge_oneg,
#      file=paste(path.res, "res_MCEM_lambda1=10_V01.Rdata", sep=""))
# 
# lambda2.MCEM <- sapply(1:4, function(l2) {test4.MCEM[[l2]]$lambda2[length(test4.MCEM[[l2]]$lambda2)]})
# lambda2.VBEM <- sapply(1:4, function(l2) {test4.VBEM[[l2]]$lambda2[length(test4.VBEM[[l2]]$lambda2)]})
# lambda2.optL2 <- sapply(1:4, function(l2) {test4.optL2[[l2]]$lambda})
# lambda2.optL2[2] <- NA
# lambda2.grridge_oneg <- sapply(1:4, function(l2) {
#   test4.grridge_oneg[[l2]]$lambda2[length(test4.grridge_oneg[[l2]]$lambda2)]})
# 
# opar <- par()
# par(mar=par('mar') + c(0, 1, 0, 0))
# plot(lambda2, lambda2.MCEM, xlab=expression(lambda[2]), ylab=expression(hat(lambda)[2]),
#      ylim=range(c(lambda2.MCEM, lambda2.optL2, lambda2.grridge_oneg), na.rm=TRUE), col=2, pch=19)
# #     ylim=range(c(lambda2.MCEM, lambda2.VBEM, lambda2.optL2, lambda2.grridge_oneg), na.rm=TRUE), col=2, pch=19)
# points(lambda2, lambda2.VBEM, col=3, pch=19)
# points(lambda2, lambda2.optL2, col=4, pch=19)
# points(lambda2, lambda2.grridge_oneg, col=5, pch=19)
# abline(a=0, b=1, lty=2)
# legend("topright", legend=c("MCEM", "VBEM", "CV", "MoM"), pch=19, col=c(2:5))
# par(opar)

### simulation 5 (group wise l1 and l2 penalization)
set.seed(123)
n <- 200
p <- 80
G <- 2
groups <- rep(1:G, each=p/G)
x <- matrix(rnorm(n*p), ncol=p, nrow=n)
lambda1 <- 1
lambda2 <- 1
lambdag <- c(0.9, 1/0.9)
m <- rep(1, n)
b0 <- rnorm(p + 1)
sigma0 <- diag(rchisq(p + 1, 1))
beta <- c(renbeta(p/G, lambda1*sqrt(lambdag[1]), lambda2*lambdag[1]), 
          renbeta(p/G, lambda1*sqrt(lambdag[2]), lambda2*lambdag[2]))
y <- rbinom(n, m, exp(x %*% beta)/(1 + exp(x %*% beta)))

test5.grVBEM1 <- grVBEM1(x, y, m, groups, lambda1, lambda2, sigma0, intercept=TRUE, eps=0.001, 
                         maxiter=50)
plot(test5.grVBEM1$lambdag[1, -c(1:2)], type="l")
plot(test5.grVBEM1$lambdag[2, -c(1:2)], type="l")
str(test5.grVBEM1)
plot(test5.grVBEM1$lowermll, type="l")

a <- 100
b <- 1
x <- seq(0.01, 2, 0.01)
y <- a*x - b*log(x)
plot(x, y, type="l")
