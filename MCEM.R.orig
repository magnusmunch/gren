<<<<<<< HEAD
##############################  preamble  #############################
# MCEM implemented in R                                               #
# version: 02                                                         #
# author: Magnus Münch                                                #
# created: 03-02-2016                                                 #
# last edited: 15-02-2016                                             #
#######################################################################

###############################  notes  ###############################
# 07-02-2017: added quasi-Newton acceleration                         #
# 15-02-2017: Gives negative lambda2???                               #
#######################################################################

### paths
path.res <- "C:/Users/Magnus/Documents/phd/review/results/"
path.code <- "C:/Users/Magnus/Documents/phd/ENVB/code/"
path.graph <- "C:/Users/Magnus/Documents/phd/ENVB/graphs/"

### libraries
library(Rcpp)
library(penalized)
library(glmnet)
library(pracma)

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

cv.glmnet2 <- function(x, y, family, lambda1, lambda2, nfolds, standardize=FALSE, intercept=TRUE) {
  
  n <- nrow(x)
  alpha <- lambda1/(2*lambda2 + lambda1)
  lambda <- (2*lambda2 + lambda1)/n
  seq.cvll <- numeric(0)
  for(a in 1:length(alpha)) {
    cv.fit <- cv.glmnet(x, y, family=family, alpha=alpha[a], lambda=c(lambda[a], lambda[a] + 1), 
                        nfolds=nfolds, standardize=standardize, intercept=intercept,
                        grouped=FALSE)
    cat("lambda2=", round(lambda2[a], 7), "\t", "cvl=", round(cv.fit$cvm[1], 5), "\n")
    seq.cvll <- c(seq.cvll, cv.fit$cvm[1])
  }
  
  alpha.cv <- alpha[which.min(seq.cvll)]
  lambda.cv <- lambda[which.min(seq.cvll)]
  
  fit <- glmnet(x, y, family=family, alpha=alpha.cv, lambda=lambda.cv, standardize=standardize, 
                intercept=intercept)
  out <- list(fit=fit, cvll=seq.cvll, lambda1=lambda1, lambda2=lambda2)
  return(out)
  
}

mll <- function(lambda2, p, lambda1, c1, c2) {
  return(-p*pnorm(-lambda1/sqrt(4*lambda2), log.p=TRUE) - 
           0.5*lambda2*c1 - lambda1^2*c2/(8*lambda2))
}

# functions for quasi-Newton EM
new_lambda2 <- function(lambda2, lambda2old, Q1old, Bold, lambda1, p, c1, c2, eps, iter) {

  drat <- dnorm(0.5*lambda1/sqrt(lambda2))
  prat <- pnorm(-0.5*lambda1/sqrt(lambda2))
  Q1 <- 0.125*lambda1^2*c2/lambda2^2 - 0.5*c1 - 0.25*p*lambda1*drat/(lambda2^(1.5)*prat)
  Q2 <- 0.375*p*lambda1*drat/(sqrt(lambda2)*prat) - 0.0625*lambda1^2*c2/lambda2^3.5 -
    0.015625*p*lambda1^3*drat*prat/lambda2^3.5 - 0.0625*p*lambda1^2*drat^2/lambda2^3
  Q1b <- 0.125*lambda1^2*c2/lambda2old^2 - 0.5*c1 - 0.25*p*lambda1*drat/(lambda2old^(1.5)*prat)

  if(iter==1) {
    ahess <- Q2
    B <- Bold
  } else {
    gc <- Q1b - Q1old
    sc <- lambda2old - lambda2
    ip <- sc*(gc - Bold*sc)

    if(abs(ip) < eps) {
      B <- Bold
    } else {
      cc <- 1/ip
      vc <- gc - Bold*sc
      B <- Bold + cc*vc^2
    }

    if((Q2 - B) >=0 ) {
      t <- 1
      while((Q2 - 0.5^t*B) >= 0) {
        t <- t + 1
      }
      ahess <- Q2 - 0.5^t*B
    } else {
      ahess <- Q2 - B
    }
  }

  lambda2new <- lambda2 - Q1/ahess

  out <- list(lambda2new=lambda2new, lambda2=lambda2, Q1old=Q1, Bold=B)
  return(out)

}

# Qlambda2 <- function(lambda2, lambda1, p, c1, c2) {
#   -p*pnorm(-0.5*lambda1/sqrt(lambda2), log.p=TRUE) - 0.5*lambda2*c1 - 
#     0.125*lambda1^2*c2/lambda2
# }
# 
# derQlambda2 <- function(lambda2, lambda1, p, c1, c2) {
#   0.125*lambda1^2*c2/lambda2^2 - 0.5*c1 - 
#     p*lambda1*dnorm(0.5*lambda1/sqrt(lambda2))/
#     (4*lambda2^1.5*pnorm(-0.5*lambda1/sqrt(lambda2)))
# }
# 
# EMstep <- function(lambda2, lambda1, p, c1, c2) {
#   opt.rout <- optim(par=lambda2, fn=Qlambda2, gr=derQlambda2, method="Brent", lower=0.001,
#                     upper=1000, control=list(fnscale=-1), lambda1=lambda1, p=p, c1=c1,
#                     c2=c2)
#   gtilde <- opt.rout$par - lambda2
#   return(gtilde)
# }
# 
# line_search <- function(lambda1, lambda2, dir, p, c1, c2, maxa, eps) {
#   
#   au <- ifelse(dir < 0, lambda2/abs(dir), maxa)
#   al <- 0
#   conv <- FALSE
#   while(!conv) {
#     aprop <- (au + al)/2
#     derprop <- derlla(aprop, lambda1, lambda2, dir, p, c1, c2)
#     if(derprop > 0) {
#       au <- aprop
#     } else {
#       al <- aprop
#     }
#     conv <- abs(derprop) < eps
#   }
#   
#   return(aprop)
#   
# }
# 
# derlla <- function(alpha, lambda1, lambda2, dir, p, c1, c2) {
#   -p*lambda1*dir*exp(-lambda1^2/(8*(lambda2 + alpha*dir)))/
#     (sqrt(8*pi)*(lambda2 + alpha*dir)^(3/2)*erfc(lambda1/sqrt(8*(lambda2 + alpha*dir)))) -
#     0.5*dir*c1 + 0.125*lambda1*c2*dir/(lambda2 + alpha*dir)^2
# }
# 
# new_lambda2 <- function(lambda2, gtilde, g, S, c1, c2) {
#   
#   dir <- -gtilde + S*g
#   alpha <- line_search(lambda1, lambda2, dir, p, c1, c2, 1000, eps=eps)
#   step <- alpha*d
#   stepg <- derQlambda2(lambda2 + step, lambda1, p, c1, c2) - g
#   stepgtilde <- EMstep(lambda2 + step, lambda1, p, c1, c2) - gtilde
#   stepstar <- -stepgtilde + S*stepg
#   stepS <- (1 + stepstar/step)*step^2/(stepg*step) - 
#     2*stepstar*step/(stepg*step)
#   
#   out <- list(lambda2=lambda2 + step, g=g + stepg, gtilde=gtilde + stepgtilde, 
#               S=S + stepS)
#   return(out)
#   
# }

MCEM <- function(x, y, m, n, p, lambda1, lambda2, b0, intercept, K, eps, maxiter) {
  
  sam <- gibbsC(x, y, m, n, p, lambda1, rep(lambda2, p), b0, TRUE, 5000)
  beta <- apply(sam$beta, 1, function(b) {d <- density(b); return(d$x[which.max(d$y)])})
  if(intercept) {
    c1 <- sum(sam$beta[-1, ]^2*sam$tau/(sam$tau - 1))/K
  } else {
    c1 <- sum(sam$beta^2*sam$tau/(sam$tau - 1))/K
  }
  c2 <- sum(sam$tau)/K
  
  mllseq <- mll(lambda2, p, lambda1, c1, c2)
  lambda2seq <- lambda2old <- lambda2
  Q1old <- Bold <- 0
  # g <- derQlambda2(lambda2, lambda1, p, c1, c2)
  # gtilde <- EMstep(lambda2, lambda1, p, c1, c2)
  # S <- 0
  
  conv <- FALSE
  iter <- 0
  while(!conv & (iter < maxiter)) {
    iter <- iter + 1
    print(paste("iteration ", iter, ", current lambda2 is ", round(lambda2, 2), sep=""))
    
    qnem_step <- new_lambda2(lambda2, lambda2old, Q1old, Bold, lambda1, p, c1, c2, eps, iter)
    lambda2 <- qnem_step$lambda2new
    lambda2old <- qnem_step$lambda2
    Q1old <- qnem_step$Q1old
    Bold <- qnem_step$Bold
    # qnem_step <- new_lambda2(lambda2, gtilde, g, S, c1, c2)
    # lambda2 <- qnem_step$lambda2
    # gtilde <- qnem_step$gtilde
    # g <- qnem_step$g
    # S <- qnem_step$g
     
    sam <- gibbsC(x, y, m, n, p, lambda1, rep(lambda2, p), beta, TRUE, 5000)
    beta <- apply(sam$beta, 1, function(b) {d <- density(b); return(d$x[which.max(d$y)])})
    
    if(intercept) {
      c1 <- sum(sam$beta[-1, ]^2*sam$tau/(sam$tau - 1))/K
    } else {
      c1 <- sum(sam$beta^2*sam$tau/(sam$tau - 1))/K
    }
    c2 <- sum(sam$tau)/K
    
    lambda2seq <- c(lambda2seq, lambda2)
    mllseq <- c(mllseq, mll(lambda2, p, lambda1, c1, c2))
    conv <- abs(lambda2 - lambda2old) < eps
    # lambda2old <- lambda2
    
  }
  
  out <- list(beta=sam$beta, tau=sam$tau, omega=sam$omega, lambda2=lambda2seq, mll=mllseq)
  return(out)
  
}

### simulation 1
set.seed(123)
n <- 100
p <- 50
m <- rep(1, n)
x <- matrix(rnorm(n*p), ncol=p, nrow=n)
lambda1 <- 1
lambda2 <- 1
#beta <- renbeta(p, lambda1=lambda1, lambda2=lambda2)
beta <- rep(0.1, p)
y <- rbinom(n, m, exp(x %*% beta)/(1 + exp(x %*% beta)))
b0 <- c(0, beta)

test1 <- MCEM(x, y, m, n, p, lambda1, 1, b0, intercept=TRUE, K=5000, eps=0.00001, maxiter=30)
test2 <- optL2(y, x, lambda1=1, model="logistic", fold=n)
test3 <- cv.glmnet2(x, y, family="binomial", lambda1, lambda2=seq(20, 60, length.out=50), 
                    nfolds=n, standardize=FALSE, intercept=TRUE)
plot(test1$lambda2, test1$mll)
plot(test3$lambda2, test3$cvll, type="l")






# simulation 2
lambda2 <- c(1, 10, 20, 30, 50)

mat.lambda2 <- matrix(NA, nrow=length(lambda2), ncol=3)
for(l2 in 1:length(lambda2)) {
  beta <- renbeta(p, lambda1=lambda1, lambda2=lambda2[l2])
  y <- rbinom(n, m, exp(x %*% beta)/(1 + exp(x %*% beta)))
  b0 <- rnorm(p + 1)
  
  test1 <- MCEM(x, y, m, n, p, lambda1, 1, b0, intercept=TRUE, K=5000, eps=0.00001, maxiter=30)
  test2 <- optL2(y, x, lambda1=1, model="logistic", fold=n)
  test3 <- cv.glmnet2(x, y, family="binomial", lambda1, lambda2=seq(1, 100, length.out=50), 
                      nfolds=n, standardize=FALSE, intercept=TRUE)
  mat.lambda2[l2, ] <- c(test1$lambda2[length(test1$lambda2)], test2$lambda, 
                         test3$lambda2[which.min(test3$cvll)])
}






# testing Gibbs sampler
set.seed(123)
n <- 100
p <- 50
m <- rep(1, n)
x <- matrix(rnorm(n*p), ncol=p, nrow=n)
lambda1 <- 1
lambda2 <- 1
#beta <- renbeta(p, lambda1=lambda1, lambda2=lambda2)
beta <- renbeta(p, lambda1, lambda2)
y <- rbinom(n, m, exp(x %*% beta)/(1 + exp(x %*% beta)))
b0 <- c(0, beta)

test1 <- gibbsC(x, y, m, n, p, 1, rep(1, p), b0, TRUE, K)
test2 <- penalized(y, x, lambda1=1, lambda2=2, model="logistic")
hist(test1$beta[2, ], breaks=50)

bmean <- apply(test1$beta, 1, mean)
bmode <- apply(test1$beta, 1, function(b) {d <- density(b); return(d$x[which.max(d$y)])})
plot(c(0, beta), bmean)
plot(c(0, beta), bmode)
plot(c(0, beta), coef(test2, which="all"))
plot(bmean, coef(test2, which="all"))
plot(bmode, coef(test2, which="all"))

str(test2)




=======
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
path.data <- "C:/Users/Magnus/Documents/phd/data/"

### libraries
library(Rcpp)
library(penalized)
library(glmnet)
library(GRridge)
library(mvtnorm)
#library(nloptr)

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

grlowermll <- function(loglambdag, lambda2, sizes, sum1) {
  mll <- 0.5*lambda2*sum(sum1*exp(loglambdag)) - 0.5*sum(sizes*loglambdag)
  grad <- 0.5*lambda2*sum1*exp(loglambdag) - 0.5*sizes
  out <- list(objective=mll, gradient=grad)
  return(out)
}

grconstr <- function(loglambdag, lambda2, sizes, sum1) {
  constr <- sum(sizes*loglambdag)
  grad <- sizes
  out <- list(constraints=constr, jacobian=grad)
  return(out)
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

cv.grVBEM <- function(x, y, m, groups, lambda1=NULL, lambda2=NULL, sigma0, 
                      intercept=TRUE, eps=1e-5, maxiter=100, nfolds, foldid=NULL) {
  
  n <- ncol(x)
  if(is.null(foldid)) {
    rest <- n %% nfolds
    foldsize <- c(rep(n %/% nfolds + as.numeric(rest!=0), times=rest),
                  rep(n %/% nfolds, times=nfolds - rest))
    foldid <- sample(rep(1:nfolds, times=foldsize))
  }
  
  preds <- numeric(n)
  for(k in 1:nfolds) {
    cat("\r", "Fold ", k, sep="")
    xtrain <- x[foldid!=k, ]
    xtest <- x[foldid==k, ]
    ytrain <- y[foldid!=k]
    mtrain <- rep(1, times=sum(foldid!=k))

    fit.grVBEM <- grVBEM(xtrain, ytrain, m=mtrain, groups=groups, lambda1=lambda1, 
                         lambda2=lambda2, sigma0=sigma0, intercept=intercept, 
                         eps=eps, maxiter=maxiter)
    best <- fit.grVBEM$mu
    if(intercept) {
      preds[foldid==k] <- as.numeric(exp(cbind(1, xtest) %*% best)/
                                       (1 + exp(cbind(1, xtest) %*% best)))  
    } else {
      preds[foldid==k] <- as.numeric(exp(xtest %*% best)/
                                       (1 + exp(xtest %*% best)))  
    }
    
  }
  
  return(preds)

}

grVBEM <- function(x, y, m, groups, lambda1=NULL, lambda2=NULL, sigma0, intercept, eps, maxiter) {
  
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
  
  sigmaold <- sigma0
  muold <- as.numeric(sigmaold %*% t(xadj) %*% as.matrix(kappa))
  ci <- as.numeric(sqrt(colSums(t(xadj) * (sigmaold %*% t(xadj))) + (colSums(t(xadj)*muold))^2))
  chi <- as.numeric(lambda2*lambdagvec*(diag(sigmaold) + muold^2))[(intercept + 1):(intercept + p)]
  
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
  cat("\n", "Estimating penalty multipliers by empirical Bayes", sep="")
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
  }
  
  cat("\n", "Penalty multipliers estimated at ", paste(round(lambdag[-G], 2), collapse=", "), 
      " and ", round(lambdag[G], 2), sep="")
  out <- list(mu=mu, sigma=sigma, ci=ci, chi=chi, lambda1=lambda1, lambda2=lambda2, 
              lambdag=lambdagseq, lowermll=lowermllseq, nouteriter=iter1, ninneriter=niter2seq, 
              conv=conv1)
  return(out)
  
}

### simulation 5 (group wise l1 and l2 penalization)
set.seed(123)
n <- 200
p <- 150
G <- 3
groups <- rep(1:G, each=p/G)
x <- matrix(rnorm(n*p), ncol=p, nrow=n)
lambda1 <- 1
lambda2 <- 1
lambdag <- c(0.2, 3, 1/(0.2*3))
m <- rep(1, n)
b0 <- rnorm(p + 1)
sigma0 <- diag(rchisq(p + 1, 1))
beta <- c(renbeta(p/G, lambda1*sqrt(lambdag[1]), lambda2*lambdag[1]), 
          renbeta(p/G, lambda1*sqrt(lambdag[2]), lambda2*lambdag[2]),
          renbeta(p/G, lambda1*sqrt(lambdag[3]), lambda2*lambdag[3]))
y <- rbinom(n, m, exp(x %*% beta)/(1 + exp(x %*% beta)))

test5.grVBEM <- grVBEM(x, y, m, groups, lambda1=NULL, lambda2=NULL, sigma0, intercept=TRUE, 
                       eps=0.001, maxiter=100)
test5.grridge <- grridge(t(x), y, list(group=CreatePartition(as.factor(groups))), unpenal=~1)

ntest <- 1000
xtest <- matrix(rnorm(ntest*p), ncol=p, nrow=ntest)
ytest <- rbinom(ntest, m, exp(xtest %*% beta)/(1 + exp(xtest %*% beta)))

pred5.grVBEM <- as.numeric(exp(cbind(1, xtest) %*% test5.grVBEM$mu)/
                             (1 + exp(cbind(1, xtest) %*% test5.grVBEM$mu)))
pred5.grridge <- predict.grridge(test5.grridge, t(xtest))[, 2]
pred5.truth <- as.numeric(exp(xtest %*% beta)/(1 + exp(xtest %*% beta)))

auc5.grVBEM <- pROC::roc(ytest, pred5.grVBEM)$auc
auc5.grridge <- pROC::roc(ytest, pred5.grridge)$auc
auc5.truth <- pROC::roc(ytest, pred5.truth)$auc

barplot(rbind(test5.grridge$lambdamults$group, lambdag, 
              test5.grVBEM$lambdag[, test5.grVBEM$nouteriter + 1]), beside=TRUE,
        names.arg=c(expression(lambda[1]), expression(lambda[2]), expression(lambda[3])),
        legend.text=paste(c("GRridge", "truth", "VBEM"), 
                          paste(", AUC=", round(c(auc5.grridge, auc5.truth, auc5.grVBEM), 2)), sep=""))
plot(beta, test5.grVBEM$mu[-1])
plot(beta, test5.grridge$betas)

### simulation 6 (group wise l1 and l2 penalization)
set.seed(123)
n <- 200
p <- 150
G <- 3
groups <- rep(1:G, each=p/G)
x <- matrix(rnorm(n*p), ncol=p, nrow=n)
m <- rep(1, n)
b0 <- rnorm(p + 1)
sigma0 <- diag(rchisq(p + 1, 1))
beta <- c(rep(0.5, p/G), rep(0.8, p/G), rep(1.5, p/G))
y <- rbinom(n, m, exp(x %*% beta)/(1 + exp(x %*% beta)))

test6.grVBEM <- grVBEM(x, y, m, groups, lambda1=NULL, lambda2=NULL, sigma0, intercept=TRUE, 
                       eps=0.001, maxiter=100)
test6.grridge <- grridge(t(x), y, list(group=CreatePartition(as.factor(groups))), unpenal=~1)

ntest <- 1000
xtest <- matrix(rnorm(ntest*p), ncol=p, nrow=ntest)
ytest <- rbinom(ntest, m, exp(xtest %*% beta)/(1 + exp(xtest %*% beta)))

pred6.grVBEM <- as.numeric(exp(cbind(1, xtest) %*% test6.grVBEM$mu)/
                             (1 + exp(cbind(1, xtest) %*% test6.grVBEM$mu)))
pred6.grridge <- predict.grridge(test6.grridge, t(xtest))[, 2]
pred6.truth <- as.numeric(exp(xtest %*% beta)/(1 + xtest %*% beta))

auc6.grVBEM <- pROC::roc(ytest, pred6.grVBEM)$auc
auc6.grridge <- pROC::roc(ytest, pred6.grridge)$auc
auc6.truth <- pROC::roc(ytest, pred6.truth)$auc

barplot(rbind(test6.grridge$lambdamults$group,
              test6.grVBEM$lambdag[, test6.grVBEM$nouteriter + 1]), 
        beside=TRUE, names.arg=c(expression(lambda[1]), expression(lambda[2]), expression(lambda[3])),
        legend.text=paste(c("GRridge", "VBEM", "truth"), 
                          paste(", AUC=", round(c(auc6.grridge, auc6.grVBEM, auc6.truth), 2)), sep=""))
plot(beta, test6.grVBEM$mu[-1])
plot(beta, test6.grridge$betas)

### simulation 7 (group wise l1 and l2 penalization)
set.seed(123)
n <- 200
p <- 150
G <- 3
groups <- rep(1:G, each=p/G)
x <- matrix(rnorm(n*p), ncol=p, nrow=n)
m <- rep(1, n)
b0 <- rnorm(p + 1)
sigma0 <- diag(rchisq(p + 1, 1))
beta <- c(rep(0, p/G), rep(0.5, p/G), rep(1.5, p/G))
y <- rbinom(n, m, exp(x %*% beta)/(1 + exp(x %*% beta)))

test7.grVBEM <- grVBEM(x, y, m, groups, lambda1=NULL, lambda2=NULL, sigma0, intercept=TRUE, 
                       eps=0.001, maxiter=100)
test7.grridge <- grridge(t(x), y, list(group=CreatePartition(as.factor(groups))), unpenal=~1)

ntest <- 1000
xtest <- matrix(rnorm(ntest*p), ncol=p, nrow=ntest)
ytest <- rbinom(ntest, m, exp(xtest %*% beta)/(1 + exp(xtest %*% beta)))

pred7.grVBEM <- as.numeric(exp(cbind(1, xtest) %*% test7.grVBEM$mu)/
                             (1 + exp(cbind(1, xtest) %*% test7.grVBEM$mu)))
pred7.grridge <- predict.grridge(test7.grridge, t(xtest))[, 2]
pred7.truth <- as.numeric(exp(xtest %*% beta)/(1 + exp(xtest %*% beta)))

auc7.grVBEM <- pROC::roc(ytest, pred7.grVBEM)$auc
auc7.grridge <- pROC::roc(ytest, pred7.grridge)$auc
auc7.truth <- pROC::roc(ytest, pred7.truth)$auc

barplot(rbind(test7.grridge$lambdamults$group, 
              test7.grVBEM$lambdag[, test7.grVBEM$nouteriter + 1]), 
        beside=TRUE, names.arg=c(expression(lambda[1]), expression(lambda[2]), expression(lambda[3])),
        legend.text=paste(c("GRridge", "VBEM", "truth"), 
                          paste(", AUC=", round(c(auc7.grridge, auc7.grVBEM, auc7.truth), 2)), sep=""))
plot(beta, test7.grVBEM$mu[-1])
plot(beta, test7.grridge$betas)

### simulation 9 (group wise l1 and l2 penalization)
set.seed(123)
set.seed(2345)
n <- 200
p <- 150
G <- 3
groups <- rep(1:G, each=p/G)
rho <- 0.5
sigma <- matrix(rho, ncol=p, nrow=p)
diag(sigma) <- 1
x <- rmvnorm(n, mean=rep(0, p), sigma=sigma)
m <- rep(1, n)
b0 <- rnorm(p + 1)
sigma0 <- diag(rchisq(p + 1, 1))
beta <- c(rep(0, p/G), rep(0.5, p/G), rep(1.5, p/G))
y <- rbinom(n, m, exp(x %*% beta)/(1 + exp(x %*% beta)))

test9.grVBEM <- grVBEM(x, y, m, groups, lambda1=NULL, lambda2=NULL, sigma0, intercept=TRUE, 
                       eps=0.001, maxiter=100)
test9.grridge <- grridge(t(x), y, list(group=CreatePartition(as.factor(groups))), unpenal=~1)

ntest <- 1000
xtest <- matrix(rnorm(ntest*p), ncol=p, nrow=ntest)
ytest <- rbinom(ntest, m, exp(xtest %*% beta)/(1 + exp(xtest %*% beta)))

pred9.grVBEM <- as.numeric(exp(cbind(1, xtest) %*% test9.grVBEM$mu)/
                             (1 + exp(cbind(1, xtest) %*% test9.grVBEM$mu)))
pred9.grridge <- predict.grridge(test9.grridge, t(xtest))[, 2]
pred9.truth <- as.numeric(exp(xtest %*% beta)/(1 + exp(xtest %*% beta)))

auc9.grVBEM <- pROC::roc(ytest, pred9.grVBEM)$auc
auc9.grridge <- pROC::roc(ytest, pred9.grridge)$auc
auc9.truth <- pROC::roc(ytest, pred9.truth)$auc

barplot(rbind(test9.grridge$lambdamults$group, 
              test9.grVBEM$lambdag[, test9.grVBEM$nouteriter + 1]), 
        beside=TRUE, names.arg=c(expression(lambda[1]), expression(lambda[2]), expression(lambda[3])),
        legend.text=paste(c("GRridge", "VBEM", "truth"), 
                          paste(", AUC=", round(c(auc9.grridge, auc9.grVBEM, auc9.truth), 2)), sep=""))
plot(beta, test9.grVBEM$mu[-1])
plot(beta, test9.grridge$betas)

### simulation 10 (group wise l1 and l2 penalization)
set.seed(123)
n <- 100
p <- 150
G <- 3
groups <- rep(1:G, each=p/G)
rho <- 0.3
sigma <- matrix(rho, ncol=p, nrow=p)
diag(sigma) <- 1
x <- rmvnorm(n, mean=rep(0, p), sigma=sigma)
m <- rep(1, n)
b0 <- rnorm(p + 1)
sigma0 <- diag(rchisq(p + 1, 1))
beta <- c(rep(0, p/G), rep(0.5, p/G), rep(1.5, p/G))
y <- rbinom(n, m, exp(x %*% beta)/(1 + exp(x %*% beta)))

test10.grVBEM <- grVBEM(x, y, m, groups, lambda1=NULL, lambda2=NULL, sigma0, intercept=TRUE, 
                        eps=0.001, maxiter=100)
test10.grridge <- grridge(t(x), y, list(group=CreatePartition(as.factor(groups))), unpenal=~1)

ntest <- 1000
xtest <- matrix(rnorm(ntest*p), ncol=p, nrow=ntest)
ytest <- rbinom(ntest, m, exp(xtest %*% beta)/(1 + exp(xtest %*% beta)))

pred10.grVBEM <- as.numeric(exp(cbind(1, xtest) %*% test10.grVBEM$mu)/
                             (1 + exp(cbind(1, xtest) %*% test10.grVBEM$mu)))
pred10.grridge <- predict.grridge(test10.grridge, t(xtest))[, 2]
pred10.truth <- as.numeric(exp(xtest %*% beta)/(1 + exp(xtest %*% beta)))

auc10.grVBEM <- pROC::roc(ytest, pred10.grVBEM)$auc
auc10.grridge <- pROC::roc(ytest, pred10.grridge)$auc
auc10.truth <- pROC::roc(ytest, pred10.truth)$auc

barplot(rbind(test10.grridge$lambdamults$group, 
              test10.grVBEM$lambdag[, test10.grVBEM$nouteriter + 1]), 
        beside=TRUE, names.arg=c(expression(lambda[1]), expression(lambda[2]), expression(lambda[3])),
        legend.text=paste(c("GRridge", "VBEM", "truth"), 
                          paste(", AUC=", round(c(auc10.grridge, auc10.grVBEM, auc10.truth), 2)), sep=""))
plot(beta, test10.grVBEM$mu[-1])
plot(beta, test10.grridge$betas)

### simulation 11 (group wise l1 and l2 penalization)
set.seed(123)
n <- 100
p <- 150
G <- 3
groups <- rep(1:G, each=p/G)
rho <- 0.3
sigma <- matrix(rho, ncol=p, nrow=p)
diag(sigma) <- 1
x <- rmvnorm(n, mean=rep(0, p), sigma=sigma)
m <- rep(1, n)
b0 <- rnorm(p + 1)
sigma0 <- diag(rchisq(p + 1, 1))
beta <- rep(0, p)
y <- rbinom(n, m, exp(x %*% beta)/(1 + exp(x %*% beta)))

test11.grVBEM <- grVBEM(x, y, m, groups, lambda1=NULL, lambda2=NULL, sigma0, intercept=TRUE, 
                        eps=0.001, maxiter=100)
test11.grridge <- grridge(t(x), y, list(group=CreatePartition(as.factor(groups))), unpenal=~1)

ntest <- 1000
xtest <- matrix(rnorm(ntest*p), ncol=p, nrow=ntest)
ytest <- rbinom(ntest, m, exp(xtest %*% beta)/(1 + exp(xtest %*% beta)))

pred11.grVBEM <- as.numeric(exp(cbind(1, xtest) %*% test11.grVBEM$mu)/
                              (1 + exp(cbind(1, xtest) %*% test11.grVBEM$mu)))
pred11.grridge <- predict.grridge(test11.grridge, t(xtest))[, 2]
pred11.truth <- as.numeric(exp(xtest %*% beta)/(1 + exp(xtest %*% beta)))

auc11.grVBEM <- pROC::roc(ytest, pred11.grVBEM)$auc
auc11.grridge <- pROC::roc(ytest, pred11.grridge)$auc
auc11.truth <- pROC::roc(ytest, pred11.truth)$auc

barplot(rbind(test11.grridge$lambdamults$group, 
              test11.grVBEM$lambdag[, test11.grVBEM$nouteriter + 1]), 
        beside=TRUE, names.arg=c(expression(lambda[1]), expression(lambda[2]), expression(lambda[3])),
        legend.text=paste(c("GRridge", "VBEM", "truth"), 
                          paste(", AUC=", round(c(auc11.grridge, auc11.grVBEM, auc11.truth), 2)), sep=""))
plot(beta, test11.grVBEM$mu[-1])
plot(beta, test11.grridge$betas)

### simulation 12 (group wise l1 and l2 penalization)
set.seed(123)
n <- 100
p <- 150
G <- 3
groups <- rep(1:G, each=p/G)
rho <- 0.3
sigma <- matrix(rho, ncol=p, nrow=p)
diag(sigma) <- 1
x <- rmvnorm(n, mean=rep(0, p), sigma=sigma)
m <- rep(1, n)
b0 <- rnorm(p + 1)
sigma0 <- diag(rchisq(p + 1, 1))
beta <- rep(1, p)
y <- rbinom(n, m, exp(x %*% beta)/(1 + exp(x %*% beta)))

test12.grVBEM <- grVBEM(x, y, m, groups, lambda1=NULL, lambda2=NULL, sigma0, intercept=TRUE, 
                        eps=0.001, maxiter=100)
test12.grVBEM2 <- penalized(y, x, unpenalized=~1, lambda1=0.5*test12.grVBEM$lambda1*sqrt(rep(test12.grVBEM$lambdag[, test12.grVBEM$nouteriter + 1], times=rle(groups)$lengths)),
                            lambda2=0.5*test12.grVBEM$lambda2*rep(test12.grVBEM$lambdag[, test12.grVBEM$nouteriter + 1], times=rle(groups)$lengths),
                            model="logistic")
test12.grridge <- grridge(t(x), y, list(group=CreatePartition(as.factor(groups))), unpenal=~1)

ntest <- 1000
xtest <- matrix(rnorm(ntest*p), ncol=p, nrow=ntest)
ytest <- rbinom(ntest, m, exp(xtest %*% beta)/(1 + exp(xtest %*% beta)))

pred12.grVBEM <- as.numeric(exp(cbind(1, xtest) %*% test12.grVBEM$mu)/
                              (1 + exp(cbind(1, xtest) %*% test12.grVBEM$mu)))
pred12.grVBEM2 <- predict(test12.grVBEM2, xtest)
pred12.grridge <- predict.grridge(test12.grridge, t(xtest))[, 2]
pred12.truth <- as.numeric(exp(xtest %*% beta)/(1 + exp(xtest %*% beta)))

auc12.grVBEM <- pROC::roc(ytest, pred12.grVBEM)$auc
auc12.grVBEM2 <- pROC::roc(ytest, pred12.grVBEM2)$auc
auc12.grridge <- pROC::roc(ytest, pred12.grridge)$auc
auc12.truth <- pROC::roc(ytest, pred12.truth)$auc

barplot(rbind(test12.grridge$lambdamults$group, 
              test12.grVBEM$lambdag[, test12.grVBEM$nouteriter + 1]), 
        beside=TRUE, names.arg=c(expression(lambda[1]), expression(lambda[2]), expression(lambda[3])),
        legend.text=paste(c("GRridge", "VBEM", "VBEM + EN", "truth"), 
                          paste(", AUC=", round(c(auc12.grridge, auc12.grVBEM, auc12.grVBEM2, auc12.truth), 2)), sep=""))
plot(beta, test12.grVBEM$mu[-1])
plot(beta, test12.grridge$betas)

# simulation 13
set.seed(456)
n <- 100
p <- 900
G <- 3
groups <- rep(1:G, each=p/G)
rho <- 0.3
sigma <- matrix(rho, ncol=p, nrow=p)
diag(sigma) <- 1
x <- rmvnorm(n, mean=rep(0, p), sigma=sigma)
m <- rep(1, n)
b0 <- rnorm(p + 1)
sigma0 <- diag(rchisq(p + 1, 1))
beta <- rep(0, p)
y <- rbinom(n, m, exp(x %*% beta)/(1 + exp(x %*% beta)))

test13.grVBEM <- grVBEM(x, y, m, groups, lambda1=NULL, lambda2=NULL, sigma0, intercept=TRUE, 
                        eps=0.001, maxiter=100)
test13.grridge <- grridge(t(x), y, unpenal=~1, innfold=10,
                          list(group=CreatePartition(as.factor(groups))))

ntest <- 1000
xtest <- matrix(rnorm(ntest*p), ncol=p, nrow=ntest)
ytest <- rbinom(ntest, m, exp(xtest %*% beta)/(1 + exp(xtest %*% beta)))

pred13.grVBEM <- as.numeric(exp(cbind(1, xtest) %*% test13.grVBEM$mu)/
                              (1 + exp(cbind(1, xtest) %*% test13.grVBEM$mu)))
pred13.grridge <- predict.grridge(test13.grridge, t(xtest))[, 2]
pred13.truth <- as.numeric(exp(xtest %*% beta)/(1 + exp(xtest %*% beta)))

auc13.grVBEM <- pROC::roc(ytest, pred13.grVBEM)$auc
auc13.grridge <- pROC::roc(ytest, pred13.grridge)$auc
auc13.truth <- pROC::roc(ytest, pred13.truth)$auc

barplot(rbind(test13.grridge$lambdamults$group, 
              test13.grVBEM$lambdag[, test13.grVBEM$nouteriter + 1]), 
        beside=TRUE, names.arg=c(expression(lambda[1]), expression(lambda[2]), expression(lambda[3])),
        legend.text=paste(c("GRridge", "VBEM", "truth"), 
                          paste(", AUC=", round(c(auc13.grridge, auc13.grVBEM, auc13.truth), 2)), sep=""))
plot(beta, test13.grVBEM$mu[-1])
plot(beta, test13.grridge$betas)

# simulation 14
set.seed(456)
n <- 100
p <- 900
G <- 3
groups <- rep(1:G, each=p/G)
rho <- 0.3
sigma <- matrix(rho, ncol=p, nrow=p)
diag(sigma) <- 1
x <- rmvnorm(n, mean=rep(0, p), sigma=sigma)
m <- rep(1, n)
b0 <- rnorm(p + 1)
sigma0 <- diag(rchisq(p + 1, 1))
beta <- rep(0.005, p)
y <- rbinom(n, m, exp(x %*% beta)/(1 + exp(x %*% beta)))

test14.grVBEM <- grVBEM(x, y, m, groups, lambda1=NULL, lambda2=NULL, sigma0, intercept=TRUE, 
                        eps=0.001, maxiter=100)
test14.grridge <- grridge(t(x), y, unpenal=~1, #innfold=10,
                          list(group=CreatePartition(as.factor(groups))))

ntest <- 1000
xtest <- matrix(rnorm(ntest*p), ncol=p, nrow=ntest)
ytest <- rbinom(ntest, m, exp(xtest %*% beta)/(1 + exp(xtest %*% beta)))

pred14.grVBEM <- as.numeric(exp(cbind(1, xtest) %*% test14.grVBEM$mu)/
                              (1 + exp(cbind(1, xtest) %*% test14.grVBEM$mu)))
pred14.grridge <- predict.grridge(test14.grridge, t(xtest))[, 2]
pred14.truth <- as.numeric(exp(xtest %*% beta)/(1 + exp(xtest %*% beta)))

auc14.grVBEM <- pROC::roc(ytest, pred14.grVBEM)$auc
auc14.grridge <- pROC::roc(ytest, pred14.grridge)$auc
auc14.truth <- pROC::roc(ytest, pred14.truth)$auc

barplot(rbind(test14.grridge$lambdamults$group, 
              test14.grVBEM$lambdag[, test14.grVBEM$nouteriter + 1]), 
        beside=TRUE, names.arg=c(expression(lambda[1]), expression(lambda[2]), expression(lambda[3])),
        legend.text=paste(c("GRridge", "VBEM", "truth"), 
                          paste(", AUC=", round(c(auc14.grridge, auc14.grVBEM, auc14.truth), 2)), sep=""))
plot(beta, test14.grVBEM$mu[-1])
plot(beta, test14.grridge$betas)

# simulation 15
set.seed(456)
n <- 100
p <- 900
G <- 3
groups <- rep(1:G, each=p/G)
rho <- 0.3
sigma <- matrix(rho, ncol=p, nrow=p)
diag(sigma) <- 1
x <- rmvnorm(n, mean=rep(0, p), sigma=sigma)
m <- rep(1, n)
b0 <- rnorm(p + 1)
sigma0 <- diag(rchisq(p + 1, 1))
beta <- c(rep(0, p/G), rep(0.004, p/G), rep(0.008, p/G))
y <- rbinom(n, m, exp(x %*% beta)/(1 + exp(x %*% beta)))

test15.grVBEM <- grVBEM(x, y, m, groups, lambda1=NULL, lambda2=NULL, sigma0, intercept=TRUE, 
                        eps=0.001, maxiter=100)
test15.grridge <- grridge(t(x), y, unpenal=~1, #innfold=10,
                          list(group=CreatePartition(as.factor(groups))))

ntest <- 1000
xtest <- matrix(rnorm(ntest*p), ncol=p, nrow=ntest)
ytest <- rbinom(ntest, m, exp(xtest %*% beta)/(1 + exp(xtest %*% beta)))

pred15.grVBEM <- as.numeric(exp(cbind(1, xtest) %*% test15.grVBEM$mu)/
                              (1 + exp(cbind(1, xtest) %*% test15.grVBEM$mu)))
pred15.grridge <- predict.grridge(test15.grridge, t(xtest))[, 2]
pred15.truth <- as.numeric(exp(xtest %*% beta)/(1 + exp(xtest %*% beta)))

auc15.grVBEM <- pROC::roc(ytest, pred15.grVBEM)$auc
auc15.grridge <- pROC::roc(ytest, pred15.grridge)$auc
auc15.truth <- pROC::roc(ytest, pred15.truth)$auc

barplot(rbind(test15.grridge$lambdamults$group, 
              test15.grVBEM$lambdag[, test15.grVBEM$nouteriter + 1]), 
        beside=TRUE, names.arg=c(expression(lambda[1]), expression(lambda[2]), expression(lambda[3])),
        legend.text=paste(c("GRridge", "VBEM", "truth"), 
                          paste(", AUC=", round(c(auc15.grridge, auc15.grVBEM, auc15.truth), 2)), sep=""))
plot(beta, test15.grVBEM$mu[-1])
plot(beta, test15.grridge$betas)

### data 1 (group wise l1 and l2 penalization)
load(paste(path.data, "mirsData.RData", sep=""))
parCons <- CreatePartition(mirsData$conservation) # using conservation status as grouping
set.seed(123)
groups <- rep(1:length(parCons), unlist(lapply(parCons, length)))
x <- apply(t(as.matrix(mirsData$transformedData))[, unlist(parCons)], 2, 
           function(x) {(x - mean(x))/sd(x)})       
y <- as.numeric(mirsData$response) - 1
n <- nrow(x)
p <- ncol(x)
m <- rep(1, n)
fit.optL2 <- optL2(y, x, unpenalized=~1, lambda1=0, model="logistic")
b0 <- coef(fit.optL2$fullfit)
pred.b0 <- as.numeric(exp(cbind(1, x) %*% b0)/(1 + exp(cbind(1, x) %*% b0)))
W <- diag(sqrt(pred.b0*(1 - pred.b0)))
Xw <- W %*% cbind(1, x)
invmat <- solve(t(Xw) %*% Xw + diag(c(0, rep(2*fit.optL2$lambda, p))))
sigma0 <- invmat %*% t(Xw) %*% Xw %*% invmat

test16.grVBEM <- grVBEM(x, y, m, groups, lambda1=NULL, lambda2=NULL, sigma0, 
                        intercept=TRUE, eps=0.001, maxiter=500)
test16.grridge <- grridge(t(x), y, list(group=CreatePartition(as.factor(groups))), unpenal=~1)

barplot(rbind(test16.grridge$lambdamults$group, 
              test16.grVBEM$lambdag[, test16.grVBEM$nouteriter + 1]), 
        beside=TRUE, names.arg=c(expression(lambda[1]), expression(lambda[2]), expression(lambda[3])),
        legend.text=c("GRridge", "VBEM"))

pen.fit <- penalized(y, x, unpenalized=~1, model="logistic",
                     lambda1=rep(sqrt(test16.grridge$lambdamults$group)*test16.grVBEM$lambda1, rle(groups)$lengths),
                     lambda2=rep(test16.grridge$lambdamults$group*test16.grVBEM$lambda2, rle(groups)$lengths))
                     

### data 2 (group wise l1 and l2 penalization) permuting the groups (no info)
load(paste(path.data, "mirsData.RData", sep=""))
parCons <- CreatePartition(mirsData$conservation) # using conservation status as grouping
set.seed(123)
n <- ncol(mirsData$transformedData)
p <- nrow(mirsData$transformedData)
groups <- rep(1:length(parCons), unlist(lapply(parCons, length)))
x <- apply(t(as.matrix(mirsData$transformedData))[, unlist(parCons)], 2, 
           function(x) {(x - mean(x))/sd(x)})[, sample(1:p)]       
y <- as.numeric(mirsData$response) - 1
m <- rep(1, n)
fit.optL2 <- optL2(y, x, unpenalized=~1, lambda1=0, model="logistic")
b0 <- coef(fit.optL2$fullfit)
pred.b0 <- as.numeric(exp(cbind(1, x) %*% b0)/(1 + exp(cbind(1, x) %*% b0)))
W <- diag(sqrt(pred.b0*(1 - pred.b0)))
Xw <- W %*% cbind(1, x)
invmat <- solve(t(Xw) %*% Xw + diag(c(0, rep(2*fit.optL2$lambda, p))))
sigma0 <- invmat %*% t(Xw) %*% Xw %*% invmat

test17.grVBEM <- grVBEM(x, y, m, groups, lambda1=NULL, lambda2=NULL, sigma0, 
                        intercept=TRUE, eps=0.001, maxiter=500)
test17.grridge <- grridge(t(x), y, list(group=CreatePartition(as.factor(groups))), unpenal=~1)

barplot(rbind(test17.grridge$lambdamults$group, 
              test17.grVBEM$lambdag[, test17.grVBEM$nouteriter + 1]), 
        beside=TRUE, names.arg=c(expression(lambda[1]), expression(lambda[2]), expression(lambda[3])),
        legend.text=c("GRridge", "VBEM"))

### simulation 18 (more groups)
set.seed(789)
n <- 200
p <- 400
G <- 4
groups <- rep(1:G, each=p/G)
rho <- 0.3
sigma <- matrix(rho, ncol=p, nrow=p)
diag(sigma) <- 1
x <- rmvnorm(n, mean=rep(0, p), sigma=sigma)
m <- rep(1, n)
b0 <- rnorm(p + 1)
sigma0 <- diag(rchisq(p + 1, 1))
beta <- c(rep(0, p/G), rep(0, p/G), rep(0.06, p/G), rep(0.06, p/G))
y <- rbinom(n, m, exp(x %*% beta)/(1 + exp(x %*% beta)))

test18.grVBEM <- grVBEM(x, y, m, groups, lambda1=NULL, lambda2=NULL, sigma0, 
                        intercept=TRUE, eps=0.001, maxiter=500)
test18.grridge <- grridge(t(x), y, list(group=CreatePartition(as.factor(groups))), unpenal=~1)

ntest <- 1000
xtest <- matrix(rnorm(ntest*p), ncol=p, nrow=ntest)
ytest <- rbinom(ntest, m, exp(xtest %*% beta)/(1 + exp(xtest %*% beta)))

pred18.grVBEM <- as.numeric(exp(cbind(1, xtest) %*% test18.grVBEM$mu)/
                              (1 + exp(cbind(1, xtest) %*% test18.grVBEM$mu)))
pred18.grridge <- predict.grridge(test18.grridge, t(xtest))[, 2]
pred18.truth <- as.numeric(exp(xtest %*% beta)/(1 + exp(xtest %*% beta)))

auc18.grVBEM <- pROC::roc(ytest, pred18.grVBEM)$auc
auc18.grridge <- pROC::roc(ytest, pred18.grridge)$auc
auc18.truth <- pROC::roc(ytest, pred18.truth)$auc

barplot(rbind(test18.grridge$lambdamults$group, 
              test18.grVBEM$lambdag[, test18.grVBEM$nouteriter + 1]), 
        beside=TRUE, names.arg=c(expression(lambda[1]), expression(lambda[2]), 
                                 expression(lambda[3]), expression(lambda[4])),
        legend.text=paste(c("GRridge", "VBEM", "truth"), 
                          paste(", AUC=", round(c(auc18.grridge, auc18.grVBEM, auc18.truth), 2)), sep=""))
plot(beta, test18.grVBEM$mu[-1])
plot(beta, test18.grridge$betas)

### data 3 (group wise l1 and l2 penalization) permuting the groups (no info)
load(paste(path.data, "mirsData.RData", sep=""))
parAbund <- CreatePartition(rowSums(mirsData$countData), mingr=25, ngroup=10, decreasing=TRUE) # using abundance as grouping
set.seed(123)
n <- ncol(mirsData$transformedData)
p <- nrow(mirsData$transformedData)
groups <- rep(1:length(parAbund), unlist(lapply(parAbund, length)))
G <- length(unique(groups))
x <- apply(t(as.matrix(mirsData$transformedData))[, unlist(parAbund)], 2, 
           function(x) {(x - mean(x))/sd(x)})[, sample(1:p)]       
y <- as.numeric(mirsData$response) - 1
m <- rep(1, n)
fit.optL2 <- optL2(y, x, unpenalized=~1, lambda1=0, model="logistic")
b0 <- coef(fit.optL2$fullfit)
pred.b0 <- as.numeric(exp(cbind(1, x) %*% b0)/(1 + exp(cbind(1, x) %*% b0)))
W <- diag(sqrt(pred.b0*(1 - pred.b0)))
Xw <- W %*% cbind(1, x)
invmat <- solve(t(Xw) %*% Xw + diag(c(0, rep(2*fit.optL2$lambda, p))))
sigma0 <- invmat %*% t(Xw) %*% Xw %*% invmat

test19.grVBEM <- grVBEM(x, y, m, groups, lambda1=NULL, lambda2=NULL, sigma0, 
                        intercept=TRUE, eps=0.001, maxiter=500)
test19.grridge <- grridge(t(x), y, list(group=CreatePartition(as.factor(groups))), unpenal=~1)

names <- vector(mode="expression", G)
for(g in 1:G) {
  names[g] <- substitute(expression(lambda[i]), list(i=g))[2]
}
barplot(rbind(test19.grridge$lambdamults$group, 
              test19.grVBEM$lambdag[, test19.grVBEM$nouteriter + 1]), 
        beside=TRUE, names.arg=names,
        legend.text=c("GRridge", "VBEM"))

# data 4 (permuted data (no info))
load(paste(path.data, "mirsData.RData", sep=""))
parCons <- CreatePartition(mirsData$conservation) # using conservation status as grouping
set.seed(123)
n <- ncol(mirsData$transformedData)
p <- nrow(mirsData$transformedData)
groups <- rep(1:length(parCons), unlist(lapply(parCons, length)))
x <- apply(t(as.matrix(mirsData$transformedData))[, unlist(parCons)], 2, 
           function(x) {(x - mean(x))/sd(x)})[, sample(1:p)]       
y <- as.numeric(mirsData$response) - 1
m <- rep(1, n)
fit.optL2 <- optL2(y, x, unpenalized=~1, lambda1=0, model="logistic")
b0 <- coef(fit.optL2$fullfit)
pred.b0 <- as.numeric(exp(cbind(1, x) %*% b0)/(1 + exp(cbind(1, x) %*% b0)))
W <- diag(sqrt(pred.b0*(1 - pred.b0)))
Xw <- W %*% cbind(1, x)
invmat <- solve(t(Xw) %*% Xw + diag(c(0, rep(2*fit.optL2$lambda, p))))
sigma0 <- invmat %*% t(Xw) %*% Xw %*% invmat
nfolds <- 10
rest <- n %% nfolds
foldsize <- c(rep(n %/% nfolds + as.numeric(rest!=0), times=rest),
              rep(n %/% nfolds, times=nfolds - rest))
foldid <- sample(rep(1:nfolds, times=foldsize))

pred20.grVBEM <- numeric(n)
for(k in 1:nfolds) {
  print(paste("Fold", k, sep=" "))
  xtrain <- x[foldid!=k, ]
  xtest <- x[foldid==k, ]
  ytrain <- y[foldid!=k]
  mtrain <- rep(1, times=sum(foldid!=k))
  
  fit.grVBEM <- grVBEM(xtrain, ytrain, m=mtrain, groups=groups, lambda1=NULL, 
                       lambda2=NULL, sigma0=sigma0, intercept=TRUE, 
                       eps=0.001, maxiter=500)
  best <- fit.grVBEM$mu
  pred20.grVBEM[foldid==k] <- as.numeric(exp(cbind(1, xtest) %*% best)/
                                           (1 + exp(cbind(1, xtest) %*% best)))
  
}

# test20.grVBEM <- cv.grVBEM(x, y, m, groups, lambda1=NULL, lambda2=NULL, sigma0, 
#                            intercept=TRUE, eps=0.001, maxiter=500, nfolds=nfolds,
#                            foldid=foldid)
test20a.grridge <- grridge(t(x), y, list(group=CreatePartition(as.factor(groups))), unpenal=~1)
test20b.grridge <- grridgeCV(test20a.grridge, t(x), y, outerfold=foldid, fixedfolds=TRUE)

auc20.grVBEM <- pROC::roc(y, pred20.grVBEM)$auc
auc20.grridge <- pROC::roc(y, test20b.grridge[, 3])$auc

# data 5 (checking predictive performance by CV)
load(paste(path.data, "mirsData.RData", sep=""))
parCons <- CreatePartition(mirsData$conservation) # using conservation status as grouping
set.seed(123)
groups <- rep(1:length(parCons), unlist(lapply(parCons, length)))
x <- apply(t(as.matrix(mirsData$transformedData))[, unlist(parCons)], 2, 
           function(x) {(x - mean(x))/sd(x)})       
y <- as.numeric(mirsData$response) - 1
n <- nrow(x)
p <- ncol(x)
m <- rep(1, n)
fit.optL2 <- optL2(y, x, unpenalized=~1, lambda1=0, model="logistic")
b0 <- coef(fit.optL2$fullfit)
pred.b0 <- as.numeric(exp(cbind(1, x) %*% b0)/(1 + exp(cbind(1, x) %*% b0)))
W <- diag(sqrt(pred.b0*(1 - pred.b0)))
Xw <- W %*% cbind(1, x)
invmat <- solve(t(Xw) %*% Xw + diag(c(0, rep(2*fit.optL2$lambda, p))))
sigma0 <- invmat %*% t(Xw) %*% Xw %*% invmat
nfolds <- 10
rest <- n %% nfolds
foldsize <- c(rep(n %/% nfolds + as.numeric(rest!=0), times=rest),
              rep(n %/% nfolds, times=nfolds - rest))
foldid <- sample(rep(1:nfolds, times=foldsize))

pred21.grVBEM <- numeric(n)
for(k in 1:nfolds) {
  print(paste("Fold", k, sep=" "))
  xtrain <- x[foldid!=k, ]
  xtest <- x[foldid==k, ]
  ytrain <- y[foldid!=k]
  mtrain <- rep(1, times=sum(foldid!=k))
  
  fit.grVBEM <- grVBEM(xtrain, ytrain, m=mtrain, groups=groups, lambda1=NULL, 
                       lambda2=NULL, sigma0=sigma0, intercept=TRUE, 
                       eps=0.001, maxiter=500)
  best <- fit.grVBEM$mu
  pred21.grVBEM[foldid==k] <- as.numeric(exp(cbind(1, xtest) %*% best)/
                                           (1 + exp(cbind(1, xtest) %*% best)))
  
}

test21a.grridge <- grridge(t(x), y, list(group=CreatePartition(as.factor(groups))), unpenal=~1)
test21b.grridge <- grridgeCV(test21a.grridge, t(x), y, outerfold=foldid, fixedfolds=TRUE)

pred21.en <- numeric(n)
for(k in 1:nfolds) {
  print(paste("Fold", k, sep=" "))
  xtrain <- x[foldid!=k, ]
  xtest <- x[foldid==k, ]
  ytrain <- y[foldid!=k]
  
  fit.cvpen <- cv.pen(xtrain, ytrain, intercept=TRUE)
  fit.en <- 
  best <- fit.grVBEM$mu
  pred21.grVBEM[foldid==k] <- as.numeric(exp(cbind(1, xtest) %*% best)/
                                           (1 + exp(cbind(1, xtest) %*% best)))
  
}




auc21.grVBEM <- pROC::roc(y, pred21.grVBEM)$auc
auc21.grridge <- pROC::roc(y, test21b.grridge[, 3])$auc
cbind(grVBEM=round(pred21.grVBEM, 2), grridge=round(test21b.grridge[, 3], 2))


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


testdens <- function(x, lambda1, lambda2) {
  c <- 0.5*sqrt(2*lambda2)*dnorm(lambda1/sqrt(2*lambda2))/pnorm(-lambda1/sqrt(2*lambda2))
  dens <- c*exp(-lambda1*abs(x) - lambda2*x^2)
  return(dens)
}
integrate(testdens, lower=-10, upper=10, lambda1=4, lambda2=1)








>>>>>>> grEBEN
