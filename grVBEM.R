##############################  preamble  #############################
# grMCEM implemented in R and C++                                     #
# version: 01                                                         #
# author: Magnus Münch                                                #
# created: 14-03-2017                                                 #
# last edited: 15-03-2017                                             #
#######################################################################

###############################  notes  ###############################
# 15-03-2017: First implementation of group-regularized elastic net   #
#             with lambda1 penalty multiplier, the square root of     #
#             lambda2 penalty multiplier                              #
#######################################################################

# paths
path.ccode <- "C:/Users/Magnus/Documents/phd/ENVB/code/"

### libraries
library(Rcpp)
library(penalized)
library(glmnet)

### functions
# source function for variational Bayes
sourceCpp(paste(path.ccode, "ENVB2.cpp", sep=""))

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

grmagn <- function(par, lambda2, sizes, sum1) {
  s <- par[1]
  loglambdag <- par[-1]
  magn <- sum((0.5*lambda2*sum1*exp(loglambdag) - 0.5*sizes + s*sizes)^2) + sum(sizes*loglambdag)^2
  return(magn)
}

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
  lambda1 <- 2*n*seq.alpha*seq.lam
  #lambda1 <- n*seq.alpha*seq.lam
  lambda2 <- n*(1 - seq.alpha)*seq.lam
  
  out <- list(lambda1=lambda1, lambda2=lambda2, cvll=seq.cvll)
  return(out)
}

grVBEM <- function(x, y, m, groups, lambda1=NULL, lambda2=NULL, intercept, eps, maxiter) {
  
  # assigning fixed (throughout algorithm) variables
  sizes <- rle(groups)$lengths
  G <- length(unique(groups))
  p <- ncol(x)
  n <- nrow(x)
  kappa <- y - m/2
  
  # if no penalty parameters are given we estimate them by cross-validation
  if(is.null(lambda1) | is.null(lambda2)) {
    cat("\r", "Estimating global lambda1 and lambda2 by cross-validation", sep="")
    opt.glob <- cv.pen(x, y, intercept)
    lambda1 <- opt.glob$lambda1[which.min(opt.glob$cvll)]
    lambda2 <- opt.glob$lambda2[which.min(opt.glob$cvll)]
    cat("\n", "Global lambda1 and lambda2 estimated at ", round(lambda1, 2), " and ", 
        round(lambda2, 2), sep="")
  }
  
  # in the multiplier setting phi does not change
  phi <- 0.25*lambda1^2/lambda2
  
  # starting values for lambdag and lagrange multiplier s
  lambdagold <- lambdagseq <- rep(1, G)
  if(intercept) {
    xadj <- cbind(1, x)
    lambdagvec <- c(0, rep(lambdagold, sizes))
  } else {
    xadj <- x
    lambdagvec <- rep(lambdagold, sizes)
  }
  s <- 0
  
  # starting values for mu and sigma
  fit.pen <- penalized(y, x, unpenalized=formula(ifelse(intercept, "~1", "~0")), 
                       model="logistic", lambda1=0, 
                       lambda2=0.5*(lambda1 + lambda2), trace=FALSE)
  #muold <- coef(fit.pen, which="all")
  b0 <- coef(fit.pen, which="all")
  pred0 <- as.numeric(exp(xadj %*% b0)/(1 + exp(xadj %*% b0)))
  w <- sqrt(pred0*(1 - pred0))
  xw <- xadj*w
  svdxw <- svd(xw)
  d <- svdxw$d
  v <- svdxw$v
  invmat <- d^2/(d^2 + lambda1 + lambda2)^2
  sigmaold <- t(t(v)*invmat) %*% t(v)
  
  # rest of the starting values follow from that
  muold <- as.numeric(sigmaold %*% (t(xadj) %*% as.matrix(kappa)))
  ci <- as.numeric(sqrt(colSums(t(xadj) * (sigmaold %*% t(xadj))) + (colSums(t(xadj)*muold))^2))
  chi <- as.numeric(0.5*(lambda1 + lambda2)*(diag(sigmaold) + muold^2))[(intercept + 1):(intercept + p)]
  
  # sum is needed in optimisation routine
  sum1 <- sapply(1:G, function(g) {
    ind <- which(groups==g) + intercept; 
    return(sum((diag(sigmaold)[ind] + muold[ind]^2)*(1 + sqrt(phi/chi[ind - intercept]))))})
  
  # keeping track of things:
  lowermllseq <- 0.5*sum(sizes*log(lambdagold)) - 0.5*lambda2*sum(lambdagold*sum1)
  niter2seq <- numeric(0)
  
  # outer loop of algorithm:
  conv1 <- FALSE
  iter1 <- 0
  cat("\n", "Estimating penalty multipliers by empirical Bayes", "\n", sep="")
  while(!conv1 & (iter1 < maxiter)) {
    iter1 <- iter1 + 1
    
    # estimating new lambdag
    # local_opts <- list(algorithm="NLOPT_LD_MMA", xtol_rel= 1.0e-7)
    # opts <- list(algorithm="NLOPT_LD_AUGLAG", xtol_rel=1.0e-7, maxeval=1000,
    #              local_opts=local_opts)
    # opt <- nloptr(x0=log(lambdag), eval_f=grlowermll, eval_g_eq=grconstr,
    #               opts=opts, lambda2=lambda2, sizes=as.numeric(sizes), sum1=sum1)
    opt <- optim(par=c(s, log(lambdagold)), fn=grmagn, lambda2=lambda2, sizes=sizes, sum1=sum1,
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
      newparam <- est_param(x, kappa, m, n, p, ci, rep(phi, p), chi, rep(lambdag*lambda2, sizes), intercept)
      sigma <- newparam$sigma
      mu <- newparam$mu
      ci <- newparam$ci
      chi <- newparam$chi
      
      # checking convergence of inner loop
      conv2 <- max(c(abs((mu - muold)/muold), abs(diag((sigma - sigmaold)/sigmaold)))) < eps
      
      muold <- mu
      sigmaold <- sigma
    }
    niter2seq <- c(niter2seq, iter2)
    
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
    
    cat("\r", "Penalty multipliers estimated at ", paste(round(lambdag[-G], 2), collapse=", "), 
        " and ", round(lambdag[G], 2), "      ", sep="")
    
  }
  
  out <- list(mu=mu, sigma=sigma, ci=ci, chi=chi, lambda1=lambda1, lambda2=lambda2, 
              lambdag=lambdagseq, lowermll=lowermllseq, nouteriter=iter1, ninneriter=niter2seq, 
              conv=conv1)
  return(out)
  
}







