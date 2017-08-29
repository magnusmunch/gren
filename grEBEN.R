##############################  preamble  #############################
# grMCEM implemented in R and C++                                     #
# version: 01                                                         #
# author: Magnus M?nch                                                #
# created: 14-03-2017                                                 #
# last edited: 28-03-2017                                             #
#######################################################################

###############################  notes  ###############################
# 28-03-2017: Added multiple partitions                               #
# 15-03-2017: First implementation of group-regularized elastic net   #
#             with lambda1 penalty multiplier, the square root of     #
#             lambda2 penalty multiplier                              #
#######################################################################

# paths
path.code <- as.character(ifelse(Sys.info()[1]=="Darwin","/Users/magnusmunch/Documents/PhD/EBEN/code/" ,"~/EBEN/code/"))

### libraries
library(Rcpp)
library(glmnet)
library(penalized)

### functions
# source function for variational Bayes
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

grmagn1 <- function(par, lambda2, sizes, sum1) {
  s <- par[1]
  loglambdag <- par[-1]
  magn <- sum((0.5*lambda2*sum1*exp(loglambdag) - 0.5*sizes + s*sizes)^2) + sum(sizes*loglambdag)^2
  return(magn)
}

grmagn2 <- function(par, lambda2, nparts, sizes, G, uintsec, sum1) {
  
  s <- par[1:nparts]
  loglambdag <- par[-c(1:nparts)]
  partind <- rep(1:nparts, times=unlist(G))
  
  ulambdagprod <- sapply(1:length(uintsec), function(cursec) {exp(sum(sapply(1:nparts, function(part) {
    loglambdag[partind==part][as.numeric(strsplit(uintsec, split=" ")[[cursec]][part])]})))})
  partsum <- sum(sapply(1:nparts, function(part) {sum((sapply(1:G[[part]], function(g) {
    0.5*lambda2*sum(ulambdagprod*sum1*sapply(1:length(uintsec), function(cursec) {
      as.numeric(strsplit(uintsec, split=" ")[[cursec]][part])==g}))}) - 0.5*sizes[[part]] + s[part]*sizes[[part]])^2)}))
  constrsum <- sum(sapply(1:nparts, function(part) {sum(sizes[[part]]*loglambdag[partind==part])^2}))
  magn <- partsum + constrsum
  return(magn)
  
}

# lambdag=unlist(lambdag, use.names=FALSE)
# lambdagold=unlist(lambdagold, use.names=FALSE)

QNupdate <- function(lambdag, lambdagold, Bold, sum1, sizes, eps) {
  
  Qgrad <- as.matrix(0.5*sizes[[1]] - 0.5*lambda2*lambdag*sum1)
  Qhess <- diag(- 0.5*lambda2*lambdag*sum1)
  
  g <- as.matrix(-0.5*sum1*(lambdagold - lambdag))
  s <- as.matrix(log(lambdagold) - log(lambdag))
  inprod <- g - Bold %*% s
  nom <- as.numeric(t(inprod) %*% s)
  if(is.nan(nom/sqrt(sum(inprod^2)*sum(s^2))) | abs(nom/sqrt(sum(inprod^2)*sum(s^2))) < eps) {
    B <- Bold
  } else {
    B <- Bold + inprod %*% t(inprod)/nom
  }

  t <- 0
  hess <- Qhess - 0.5^t*B
  while(any(eigen(hess)$values >= 0)) {
    t <- t + 1
    hess <- Qhess - 0.5^t*B
  }
  
  invmat <- solve(rbind(cbind(hess, sizes[[1]]), c(sizes[[1]], 0)))
  parnew <- c(log(lambdag), 0) - as.numeric(invmat %*% rbind(Qgrad, 0))
  out <- list(lambdagnew=exp(parnew[-length(parnew)]), lambdag=lambdag, B=B, 
              s=parnew[length(parnew)])
  return(out)
  
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
  
  out <- list(lambda1=lambda1, lambda2=lambda2, alpha=seq.alpha, lambda=seq.lam, cvll=seq.cvll)
  return(out)
}

# nparts <- 2
# partsmat <- matrix(c(1, 1, 1, 2, 2, 1, 2, 3, 1, 2), ncol=2)
# partsind <- c(1, 1, 2, 2, 2)
# loglambdag <- c(1, 2, 3, 4, 5)
# sum1 <- c(10, 11, 12, 13, 14)
# sizes <- c(6, 7, 8, 9, 10)
# s <- 1
# 
# loglambdamult <- rowSums(sapply(1:2, function(part) {loglambdag[partsind==part][partsmat[, part]]}))
# 
# c(sum((exp(loglambdamult)*sum1)[c(1:3)]), sum((exp(loglambdamult)*sum1)[c(4:5)]))
# tapply(exp(loglambdamult)*sum1, partsmat[, 1], sum)
# 
# c(sum((exp(loglambdamult)*sum1)[c(1, 4)]), sum((exp(loglambdamult)*sum1)[c(2, 5)]), 
#   sum((exp(loglambdamult)*sum1)[3]))
# tapply(exp(loglambdamult)*sum1, partsmat[, 2], sum)
# 
# unlist(sapply(1:nparts, function(part) {tapply(exp(loglambdamult)*sum1, partsmat[, part], sum)}))
# c(tapply(exp(loglambdamult)*sum1, partsmat[, 1], sum), tapply(exp(loglambdamult)*sum1, partsmat[, 2], sum))
# 
# sum((0.5*lambda2*unlist(sapply(1:nparts, function(part) {tapply(exp(loglambdamult)*sum1, partsmat[, part], sum)})) + 
#   (s - 0.5)*sizes*loglambdag)^2)

# grmagn3 <- function(par, lambda2, nparts, partsind, partsmat, sizes, intsizes, sum1) {
#   
#   s <- par[1]
#   loglambdag <- par[-1]
#   
#   loglambdamult <- rowSums(sapply(1:nparts, function(part) {
#     loglambdag[partsind==part][partsmat[, part]]}))
#   
#   partsum <- sum((0.5*lambda2*unlist(sapply(1:nparts, function(part) {
#     tapply(exp(loglambdamult)*sum1, partsmat[, part], sum)})) + (s - 0.5)*sizes*exp(loglambdag))^2)
#   constrsum <- sum(sizes*loglambdag)^2
#   magn <- partsum + constrsum
#   return(magn)
#   
# }


grmagn3 <- function(par, lambda2, nparts, partsind, partsmat, sizes, G, sum1) {
  
  # s <- rep(par[1:nparts], times=unlist(G))
  # loglambdag <- par[-c(1:nparts)]
  s <- par[1]
  loglambdag <- par[-1]
  
  loglambdamult <- rowSums(sapply(1:nparts, function(part) {
    loglambdag[partsind==part][partsmat[, part]]}))
  
  partsum <- sum((0.5*lambda2*unlist(sapply(1:nparts, function(part) {
    tapply(exp(loglambdamult)*sum1, partsmat[, part], sum)})) + (s - 0.5)*sizes)^2)
  constrsum <- sum(sapply(1:nparts, function(part) {
    sum(loglambdag[partsind==part]*sizes[partsind==part])^2}))
  magn <- partsum + constrsum
  return(magn)
  
}

grEBEN <- function(x, y, m, partitions, lambda1=NULL, lambda2=NULL, intercept, eps, maxiter,
                   trace=TRUE) {
  
  if(!is.list(partitions) & is.null(names(partitions))) {
    partitions <- list(partition1=partitions)
  } else if(is.null(names(partitions))) {
    names(partitions) <- paste("partition", 1:length(partitions), sep="")
  }
  
  partnames <- names(partitions)
  
  # assigning fixed (throughout algorithm) variables
  nparts <- length(partitions)
  sizes <- lapply(partitions, function(part) {rle(sort(part))$lengths})
  G <- lapply(partitions, function(part) {length(unique(part))})
  p <- ncol(x)
  n <- nrow(x)
  kappa <- y - m/2
  
  # if no penalty parameters are given we estimate them by cross-validation
  if(is.null(lambda1) | is.null(lambda2)) {
    if(trace) {cat("\r", "Estimating global lambda1 and lambda2 by cross-validation", sep="")}
    srt <- proc.time()[3]
    opt.glob <- cv.pen(x, y, intercept)
    cv.time <- proc.time()[3] - srt
    lambda1 <- opt.glob$lambda1[which.min(opt.glob$cvll)]
    lambda2 <- opt.glob$lambda2[which.min(opt.glob$cvll)]
    if(trace) {cat("\n", "Global lambda1 and lambda2 estimated at ", round(lambda1, 2), " and ", 
        round(lambda2, 2), " in ", round(cv.time, 2), " seconds", sep="")}
  }
  
  # in the multiplier setting phi does not change
  phi <- 0.25*lambda1^2/lambda2
  
  # starting values for lambdag, lagrange multiplier s and possible approximate hessian B
  lambdagnew <- lambdag <- lambdagold <- lambdagseq <- lapply(G, function(gpart) {rep(1, gpart)})
  # s <- rep(0, times=nparts)
  s <- 0
  
  # starting values for mu and sigma
  fit.pen <- penalized(y, x, unpenalized=formula(ifelse(intercept, "~1", "~0")), 
                       model="logistic", lambda1=0, 
                       lambda2=2*(lambda1 + lambda2), trace=FALSE)
  #muold <- coef(fit.pen, which="all")
  if(intercept) {
    xadj <- cbind(1, x)
  } else {
    xadj <- x
  }
  b0 <- coef(fit.pen, which="all")
  pred0 <- as.numeric(exp(xadj %*% b0)/(1 + exp(xadj %*% b0)))
  w <- sqrt(pred0*(1 - pred0))
  xw <- xadj*w
  svdxw <- svd(xw)
  d <- svdxw$d
  v <- svdxw$v
  invmat <- d^2/(d^2 + 4*(lambda1 + lambda2))^2
  sigmaold <- t(t(v)*invmat) %*% t(v)
  dsigmaold <- diag(sigmaold)
  
  # rest of the starting values follow from that
  muold <- as.numeric(sigmaold %*% (t(xadj) %*% as.matrix(kappa)))
  ci <- as.numeric(sqrt(colSums(t(xadj) * (sigmaold %*% t(xadj))) + (colSums(t(xadj)*muold))^2))
  chi <- as.numeric(0.5*(lambda1 + lambda2)*(diag(sigmaold) + muold^2))[(intercept + 1):(intercept + p)]
  
  # sum is needed in optimisation routine
  intsec <- do.call("paste", c(partitions, sep=" "))
  uintsec <- unique(intsec)
  sum1 <- sapply(1:length(uintsec), function(cursec) {
    ind <- which(intsec==uintsec[[cursec]]) + intercept;
    sum((diag(sigmaold)[ind] + muold[ind]^2)*(1 + sqrt(phi/chi[ind - intercept])))})
  intsizes <- as.numeric(t(table(as.data.frame(matrix(unlist(partitions), ncol=nparts)))))
  partsmat <- unique(matrix(unlist(partitions), ncol=nparts))
  partsind <- rep(1:nparts, times=unlist(G))
  
  # sum1 <- sapply(1:G, function(g) {
  #   ind <- which(groups==g) + intercept; 
  #   return(sum((diag(sigmaold)[ind] + muold[ind]^2)*(1 + sqrt(phi/chi[ind - intercept]))))})
  
  # keeping track of things:
  lambdamult <- exp(rowSums(sapply(1:nparts, function(part) {
    log(unlist(lambdag))[partsind==part][partsmat[, part]]})))
  # lambdagprod <- sapply(1:length(uintsec), function(cursec) {prod(sapply(1:nparts, function(part) {
  #   lambdag[[part]][as.numeric(strsplit(uintsec[cursec], split=" ")[[1]][part])]}))})
  lowermllseq <- 0.5*sum(sapply(1:nparts, function(part) {sum(sizes[[part]]*log(lambdag[[part]]))})) - 
    0.5*lambda2*sum(sum1*lambdamult)
  # lowermllseq <- sum(sizes*log(lambdagold)) - 0.5*lambda2*sum(lambdagold*sum1)
  niter2seq <- numeric(0)
  
  # outer loop of algorithm:
  conv1 <- FALSE
  iter1 <- 0
  if(trace) {cat("\n", "Estimating penalty multipliers by empirical Bayes", "\n", sep="")}
  srt <- proc.time()[3]
  while(!conv1 & (iter1 < maxiter)) {
    
    iter1 <- iter1 + 1
    
    # estimating new lambdag
    if(length(partitions)==1) {
      opt <- optim(par=c(s, log(unlist(lambdag, use.names=FALSE))), fn=grmagn1, lambda2=lambda2, 
                   sizes=sizes[[1]], sum1=sum1, method="BFGS", control=list(maxit=5000))
      s <- opt$par[1]
      lambdagnew <- list(exp(opt$par[-1]))
      lambdagseq <- sapply(1:nparts, function(part) {cbind(lambdagseq[[part]], lambdagnew[[part]])}, simplify=FALSE)
      names(lambdagseq) <- partnames  
    } else {
      opt <- optim(par=c(s, log(unlist(lambdag, use.names=FALSE))), fn=grmagn3, lambda2=lambda2, nparts=nparts, 
                   partsind=partsind, partsmat=partsmat, sizes=unlist(sizes), G=G, sum1=sum1, 
                   method="Nelder-Mead", control=list(maxit=5000))
      s <- opt$par[1]
      lambdagnew <- split(exp(opt$par[-1]), factor(rep(partnames, unlist(G)), levels=partnames))
      # s <- opt$par[1:nparts]
      # lambdagnew <- split(exp(opt$par[-c(1:nparts)]), factor(rep(partnames, unlist(G)), levels=partnames))
      lambdagseq <- sapply(1:nparts, function(part) {cbind(lambdagseq[[part]], lambdagnew[[part]])}, simplify=FALSE)
      names(lambdagseq) <- partnames
    }
    
    # inner loop of algorithm:
    conv2 <- 0
    iter2 <- 0
    while(!conv2 & (iter2 < maxiter)) {
      iter2 <- iter2 + 1
      
      # estimating new model parameters
      lambdamultvec <- apply(sapply(1:nparts, function(part) {lambdagnew[[part]][partitions[[part]]]}), 1, prod)
      newparam <- est_param(x, kappa, m, n, p, ciold, phi, chiold, lambdamultvec*lambda2, intercept)
      dsigma <- newparam$dsigma
      mu <- newparam$mu
      ci <- newparam$ci
      chi <- newparam$chi
      
      # checking convergence of inner loop
      conv2 <- max(c(abs((mu - muold)/muold), abs((dsigma - dsigmaold)/dsigmaold))) < eps
      
      muold <- mu
      dsigmaold <- dsigma
    }
    niter2seq <- c(niter2seq, iter2)
    
    # sum is needed in optimisation routine
    sum1 <- sapply(1:length(uintsec), function(cursec) {
      ind <- which(intsec==uintsec[[cursec]]);
      sum((dsigmaold + muold^2)*(1 + sqrt(phi/chi)))})
    
    # keeping track of lower bound on marginal log likelihood
    lambdamult <- exp(rowSums(sapply(1:nparts, function(part) {
      log(unlist(lambdagnew))[partsind==part][partsmat[, part]]})))
    lowermllseq <- c(lowermllseq, 0.5*sum(sapply(1:nparts, function(part) {sum(sizes[[part]]*log(lambdag[[part]]))})) - 
                       0.5*lambda2*sum(sum1*lambdamult))
    # lowermllseq <- c(lowermllseq, sum(sizes*log(lambdag)) - 0.5*lambda2*sum(lambdag*sum1))
    
    # checking convergence of outer loop:
    conv1 <- max(abs((unlist(lambdagnew) - unlist(lambdag))/unlist(lambdag))) < eps
    
    # updating lambdag for new iteration
    lambdagold <- lambdag
    lambdag <- lambdagnew
    
    # printing progress
    if(trace) {cat("\r", "Penalty multipliers estimated at ", paste(partnames, lapply(lambdag, function(part) {
      paste(round(part, 2), collapse=", ")}), sep=": ", collapse=" and "), "      ", sep="")}
    # cat("\r", "Penalty multipliers estimated at ", paste(round(lambdag[-G], 2), collapse=", "), 
    #     " and ", round(lambdag[G], 2), "      ", sep="")
    
  }
  eb.time <- proc.time()[3] - srt
  if(trace) {cat("\r", "Penalty multipliers estimated at ", paste(partnames, lapply(lambdag, function(part) {
    paste(round(part, 2), collapse=", ")}), sep=": ", collapse=" and "), " in ", 
    round(eb.time, 2), " seconds ", sep="")}
  
  out <- list(mu=mu, dsigma=dsigma, ci=ci, chi=chi, lambda1=lambda1, lambda2=lambda2, 
              lambdag=lambdagseq, lowermll=lowermllseq, nouteriter=iter1, ninneriter=niter2seq, 
              conv=conv1)
  
  return(out)
  
}

# prediction function
predict.grVBEM <- function(object, newx, intercept=TRUE) {
  if(intercept) {
    x <- cbind(1, newx)
  } else {
    x <- newx
  }
  prob <- as.numeric(exp(x %*% object$mu)/(1 + exp(x %*% object$mu)))
  return(prob)
}

# # Variational Bayes estimation of posterior
# VBen <- function(x, y, m, lambda1=NULL, lambda2=NULL, intercept, eps, maxiter, trace=TRUE) {
#   
#   # data characteristics
#   p <- ncol(x)
#   n <- nrow(x)
#   kappa <- y - m/2
#   
#   # if no penalty parameters are given we estimate them by cross-validation
#   cv.time <- NULL
#   if(is.null(lambda1) | is.null(lambda2)) {
#     if(trace) {cat("\r", "Estimating lambda1 and lambda2 by cross-validation", sep="")}
#     srt <- proc.time()[3]
#     opt.glob <- cv.pen(x, y, intercept)
#     cv.time <- proc.time()[3] - srt
#     lambda1 <- opt.glob$lambda1[which.min(opt.glob$cvll)]
#     lambda2 <- opt.glob$lambda2[which.min(opt.glob$cvll)]
#     if(trace) {cat("\n", "Lambda1 and lambda2 estimated at ", round(lambda1, 2), " and ", 
#                    round(lambda2, 2), " in ", round(cv.time, 2), " seconds", sep="")}
#   }
#   
#   # fixed parameters
#   phi <- 0.25*lambda1^2/lambda2
#   
#   # starting values for mu and sigma
#   fit.pen <- penalized(y, x, unpenalized=formula(ifelse(intercept, "~1", "~0")), 
#                        model="logistic", lambda1=0, 
#                        lambda2=2*(lambda1 + lambda2), trace=FALSE)
#   if(intercept) {
#     xadj <- cbind(1, x)
#   } else {
#     xadj <- x
#   }
#   b0 <- coef(fit.pen, which="all")
#   pred0 <- as.numeric(exp(xadj %*% b0)/(1 + exp(xadj %*% b0)))
#   w <- sqrt(pred0*(1 - pred0))
#   xw <- xadj*w
#   svdxw <- svd(xw)
#   d <- svdxw$d
#   v <- svdxw$v
#   invmat <- d^2/(d^2 + 4*(lambda1 + lambda2))^2
#   sigmaold <- t(t(v)*invmat) %*% t(v)
#   
#   # starting values
#   muold <- as.numeric(sigmaold %*% (t(xadj) %*% as.matrix(kappa)))
#   ci <- as.numeric(sqrt(colSums(t(xadj) * (sigmaold %*% t(xadj))) + (colSums(t(xadj)*muold))^2))
#   chi <- as.numeric(0.5*(lambda1 + lambda2)*(diag(sigmaold) + muold^2))[(intercept + 1):(intercept + p)]
#   
#   conv <- FALSE
#   iter <- 0
#   if(trace) {cat("\n", "Estimating model parameters by variational Bayes", "\n", sep="")}
#   srt <- proc.time()[3]
#   while(!conv & (iter < maxiter)) {
#     
#     iter <- iter + 1
#     
#     # estimating new model parameters
#     newparam <- est_param(x, kappa, m, n, p, ci, rep(phi, p), chi, rep(lambda2, p), intercept)
#     sigma <- newparam$sigma
#     mu <- newparam$mu
#     ci <- newparam$ci
#     chi <- newparam$chi
#       
#     # checking convergence
#     conv <- max(c(abs((mu - muold)/muold), abs(diag((sigma - sigmaold)/sigmaold)))) < eps
#       
#     muold <- mu
#     sigmaold <- sigma
#     
#     if(trace) {cat("\r", "Iteration ", iter, "  ", sep="")}
#   }
#   vb.time <- proc.time()[3] - srt
#   
#   # output
#   out <- list(mu=mu, sigma=sigma, ci=ci, chi=chi, lambda1=lambda1, lambda2=lambda2, 
#               niter=iter, conv=conv, cv.time=cv.time, vb.time=vb.time)
#   
#   return(out)
#    
# }
# 
# sourceCpp(paste(path.code, "ENgibbs.cpp", sep=""))
# MCMCen <- function(x, y, m, lambda1=NULL, lambda2=NULL, intercept, K=5000, trace=TRUE) {
#   
#   # data characteristics
#   p <- ncol(x)
#   n <- nrow(x)
#   
#   # if no penalty parameters are given we estimate them by cross-validation
#   cv.time <- NULL
#   if(is.null(lambda1) | is.null(lambda2)) {
#     if(trace) {cat("\r", "Estimating lambda1 and lambda2 by cross-validation", sep="")}
#     srt <- proc.time()[3]
#     opt.glob <- cv.pen(x, y, intercept)
#     cv.time <- proc.time()[3] - srt
#     lambda1 <- opt.glob$lambda1[which.min(opt.glob$cvll)]
#     lambda2 <- opt.glob$lambda2[which.min(opt.glob$cvll)]
#     if(trace) {cat("\n", "Lambda1 and lambda2 estimated at ", round(lambda1, 2), " and ", 
#                    round(lambda2, 2), " in ", round(cv.time, 2), " seconds", sep="")}
#   }
#   
#   # starting values for beta
#   fit.pen <- penalized(y, x, unpenalized=formula(ifelse(intercept, "~1", "~0")), 
#                        model="logistic", lambda1=0, 
#                        lambda2=2*(lambda1 + lambda2), trace=FALSE)
#   if(intercept) {
#     xadj <- cbind(1, x)
#   } else {
#     xadj <- x
#   }
#   b0 <- coef(fit.pen, which="all")
#   
#   # sampling from posterior
#   if(trace) {cat("\n", "Sampling from posterior", sep="")}
#   srt <- proc.time()[3]
#   gsamp <- gibbsC(x, y, m, n, p, lambda1, rep(lambda2, p), b0, intercept, K)
#   mcmc.time <- proc.time()[3] - srt
#   
#   # output
#   out <- list(beta=gsamp$beta, tau=gsamp$tau, omega=gsamp$omega, lambda1=lambda1, 
#               lambda2=lambda2, cv.time=cv.time, mcmc.time=mcmc.time)
#   
#   return(out)
#   
# }
# 
# 
#
# set.seed(9001)
# n <- 200
# p <- 200
# G <- 2
# x <- matrix(rnorm(n*p), ncol=p, nrow=n)
# lambda1 <- 1
# lambda2 <- 1
# lambdag <- c(0.2, 5)
# m <- rep(1, n)
# b0 <- rnorm(p + 1)
# sigma0 <- diag(rchisq(p + 1, 1))
# beta <- c(renbeta(p/G, sqrt(lambdag[1])*lambda1, lambdag[1]*lambda2), 
#           renbeta(p/G, sqrt(lambdag[2])*lambda1, lambdag[2]*lambda2))
# y <- rbinom(n, m, exp(x %*% beta)/(1 + exp(x %*% beta)))
# 
# test1.vb <- VBen(x, y, m, lambda1=1, lambda2=1, intercept=FALSE, eps=0.001, 
#                  maxiter=200, trace=TRUE)
# test1.mcmc <- MCMCen(x, y, m, lambda1=1, lambda2=1, intercept=FALSE, K=10000, trace=TRUE)
# 
# VBmu <- test1.vb$mu
# VBmcmc <- apply(test1.mcmc$beta, 1, mean)
# plot(VBmu, VBmcmc)
# 
# VBsigma <- diag(test1.vb$sigma)
# MCMCsigma <- apply(test1.mcmc$beta, 1, var)
# plot(VBsigma, MCMCsigma, col=rep(1:G, each=p/G))
# 
# VBmult <- VBsigma/exp(mean(log(VBsigma)))
# MCMCmult <- MCMCsigma/exp(mean(log(MCMCsigma)))
# plot(VBmult, MCMCmult)
# 
# VBgmsigmag <- sapply(1:G, function(g) {exp(mean(log(VBsigma[c(((g - 1)*p/G + 1):(g*(p/G)))])))})
# MCMCgmsigmag <- sapply(1:G, function(g) {exp(mean(log(MCMCsigma[c(((g - 1)*p/G + 1):(g*(p/G)))])))})
# 
# VBamsigmag <- sapply(1:G, function(g) {mean(VBsigma[c(((g - 1)*p/G + 1):(g*(p/G)))])})
# MCMCamsigmag <- sapply(1:G, function(g) {mean(MCMCsigma[c(((g - 1)*p/G + 1):(g*(p/G)))])})
# 
# (1/VBgmsigmag)/exp(mean(log(1/VBgmsigmag))) 
# (1/MCMCgmsigmag)/exp(mean(log(1/MCMCgmsigmag))) 
# 
# 1/VBamsigmag
# 1/MCMCamsigmag
# 
# dim(test1.mcmc$beta)
# cov(test1.mcmc$beta[1:100, ])
# cov(test1.mcmc$beta[1:100, ])
# 
# boxplot(cbind(VBsigma[1:100], VBsigma[101:200]))
# boxplot(cbind(MCMCsigma[1:100], MCMCsigma[101:200]))
# boxplot(cbind(beta[1:100], beta[101:200]))
# 
# exp(sum(log(VBmsigma))/3)/VBmsigma
# exp(sum(log(MCMCmsigma))/3)/MCMCmsigma
# 
# cbind(apply(test1.mcmc$beta, 1, var), diag(test1.vb$sigma))
# cbind(VBmult, MCMCmult)



