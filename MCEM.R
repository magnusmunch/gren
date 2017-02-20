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




