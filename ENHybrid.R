##############################  preamble  #############################
# code belonging to ENVB2.pdf                                         #
# version: 01                                                         #
# author: Magnus Münch                                                #
# created: 15-12-2016                                                 #
# last edited: 11-01-2016                                             #
#######################################################################



###############################  notes  ###############################
# 11-01-2017: run VBEM until convergence, then optimise               #
#             hyperparameters once                                    #
#######################################################################

### paths
path.data <- "C:/Users/Magnus/Documents/phd/data/"
path.code <- "C:/Users/Magnus/Documents/phd/ENVB/code/"

### libraries
library(Rcpp)
library(glmnet)
library(GRridge)
library(pROC)

### functions
# source function to estimate parameters in C++
sourceCpp(paste(path.code, "ENVB2.cpp", sep=""))

# these three functions are used in sampling from the elastic net prior
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

magnitude <- function(param, lambda1, lambda2, sum1, sum2, sizes) {
  
  # optimises lambdags
  s <- param[1]
  lambdag <- param[-1]
  p <- sum(sizes)
  part1 <- 0.5*lambda2*sum2 - lambda1^2*sum1/(8*lambda2*lambdag^2) + 
    lambda1*sizes*dnorm(lambda1/sqrt(4*lambda2*lambdag))/
    (4*sqrt(lambda2)*lambdag^1.5*pnorm(-lambda1/sqrt(4*lambda2*lambdag))) - s*sizes/lambdag
  part2 <- sum(sizes*log(lambdag))
  
  magn <- sqrt(sum(part1^2) + part2^2)
  return(magn)
  
}

# magnitude2 <- function(param, lambda1, lambda2, sum1, sum2, sizes) {
#   
#   # optimises lambda2gs
#   s <- param[1]
#   lambda2g <- param[-1]
#   G <- length(lambda2g)
#   part1 <- 0.5*sum2 - lambda1^2*sum1/(8*lambda2g^2) + 
#     lambda1*sizes*dnorm(lambda1/sqrt(4*lambda2g))/
#     (4*lambda2g^1.5*pnorm(-lambda1/sqrt(4*lambda2*lambda2g))) - s/(G*lambda2g)
#   part2 <- lambda2 - sum(log(lambda2g))/G
#   
#   magn <- sqrt(sum(part1^2) + part2^2)
#   return(magn)
#   
# }

magnitude3 <- function(param, lambda1, lambda2, sum1, sum2, sizes) {
  
  # constraint on arithmetic mean of penalties
  s <- param[1]
  lambdag <- param[-1]
  p <- sum(sizes)
  part1 <- 0.5*lambda2*sum2 - lambda1^2*sum1/(8*lambda2*lambdag^2) + 
    lambda1*sizes*dnorm(lambda1/sqrt(4*lambda2*lambdag))/
    (sqrt(16*lambda2)*lambdag^1.5*pnorm(-lambda1/sqrt(4*lambda2*lambdag))) + s*sizes
  part2 <- sum(sizes*lambdag) - p
  
  magn <- sqrt(sum(part1^2) + part2^2)
  return(magn)
  
}

# f.nmll <- function(lambdag, lambda1, lambda2, sum1, sum2, sizes) {
#   
#   nmll <- sum(sizes*log(lambda1)) - 0.5*lambda2*sum(lambdag*sum2) - 
#     sum(sizes*pnorm(-lambda1/sqrt(4*lambda2*lambdag), log.p=TRUE)) - 
#     0.125*lambda1^2*sum(sum1/(lambda2*lambdag))
#   return(nmll)
#   
# }
# 
# f.eqcon <- function(lambdag, lambda1, lambda2, sum1, sum2, sizes) {
#   
#   G <- length(lambdag)
#   return(sum(log(lambdag))/G)
#   
# }
# 
# f.ineqcon <- function(param, lambda1, lambda2, sum1, sum2, sizes) {
#   
#   return(param)
#   
# }
# 
# f.uniroot <- function(s, lambda1, lambda2, sum1, sum2) {
#   
#   return(sum(log(-s + sqrt(s^2 + lambda1^2*sum1*sum2/4))) - length(sum1)*log(lambda2) - sum(log(sum2)))
#   
# }

marg.ll <- function(lambda2g, lambda1, sum1, sum2, sizes) {
  # this functions calculates the marginal likelihood 
  # (the part proportional to lambda)
  
  G <- length(lambda2g)
  p <- sum(sizes)
  
  part1 <- p*log(lambda1)
  part2 <- 0.5*sum(sum2*lambda2g)
  part3 <- sum(sizes*pnorm(-lambda1/sqrt(4*lambda2g), log.p=TRUE))
  part4 <- 0.125*lambda1^2*sum(sum1/lambda2g)
  
  ll <- part1 - part2 - part3 - part4
  return(ll)
  
}

enfit <- function(x, y, intercept=TRUE, standardize=TRUE, alpha=NULL, psel, nlambda=1000,
                  nfolds=10, foldid=NULL) {
  
  n <- length(y)
  
  if(is.null(alpha)) {
    seq.alpha <- seq(0.01, 0.99, length.out=50)  
  } else {
    seq.alpha <- alpha
  }
  seq.lam <- numeric(length(seq.alpha))
  seq.df <- numeric(length(seq.alpha))
  seq.cvll <- numeric(length(seq.alpha))
  if(is.null(foldid)) {
    foldid <- sample(rep(1:nfolds, length.out=n))  
  }
  # for(a in 1:length(seq.alpha)) {
  #   fit <- glmnet(x, y, family="binomial", alpha=seq.alpha[a], dfmax=psel, standardize=standardize, 
  #                 intercept=intercept, nlambda=nlambda)
  #   cv.fit <- cv.glmnet(x, y, family="binomial", alpha=seq.alpha[a], lambda=fit$lambda, standardize=standardize, 
  #                       intercept=intercept, foldid=foldid)
  #   seq.lam[a] <- cv.fit$lambda[which.min(cv.fit$cvm)]
  #   seq.df[a] <- fit$df[which.min(cv.fit$cvm)]
  #   seq.cvll[a] <- min(cv.fit$cvm)
  # }
  for(a in 1:length(seq.alpha)) {
    fit <- glmnet(x, y, family="binomial", alpha=seq.alpha[a], dfmax=psel, standardize=standardize, 
                  intercept=intercept, nlambda=nlambda)
    # cv.fit <- cv.glmnet(x, y, family="binomial", alpha=seq.alpha[a], lambda=fit$lambda, standardize=standardize, 
    #                     intercept=intercept, foldid=foldid)
    # seq.lam[a] <- cv.fit$lambda[which.min(cv.fit$cvm)]
    # seq.df[a] <- fit$df[which.min(cv.fit$cvm)]
    # seq.cvll[a] <- min(cv.fit$cvm)
    
    cv.fit <- cv.glmnet(x, y, family="binomial", alpha=seq.alpha[a], 
                        lambda=tail(fit$lambda, n=2L), standardize=standardize, 
                        intercept=intercept, foldid=foldid)
    seq.lam[a] <- tail(fit$lambda, n=1L)
    seq.df[a] <- tail(fit$df, n=1L)
    seq.cvll[a] <- tail(cv.fit$cvm, n=1L)
  }
  
  alpha <- seq.alpha[which.min(seq.cvll)]
  lambda <- seq.lam[which.min(seq.cvll)]
  
  final.fit <- glmnet(x, y, family="binomial", alpha=alpha, lambda=lambda, standardize=standardize, 
                      intercept=intercept)
  final.fit$alpha <- alpha
  return(final.fit)
  
}


intercept=TRUE
psel=50
alpha=NULL
standardize=FALSE 
nlambda=1000
nfolds=10
foldid=NULL
maxiter=1000
epsilon=1e-07
trace=TRUE

# general fitting function
enhybrid <- function(x, y, groupid, intercept=TRUE, psel, alpha=NULL, standardize=FALSE, nlambda=1000, nfolds=10, 
                     foldid=NULL, maxiter=1000, epsilon=1e-07, trace=TRUE) {
  
  # setting general parameters and variables
  p <- ncol(x)
  n <- nrow(x)
  G <- length(unique(groupid))
  sizes <- rle(groupid)$lengths
  m <- rep(1, n)
  kappa <- y - m/2
  if(intercept) {
    xaug <- cbind(1, x)
  } else {
    xaug <- x
  }
  
  # cross-validating to find global lambda values
  if(is.null(alpha)) {
    seq.alpha <- seq(0.01, 0.99, length.out=50)  
  } else {
    seq.alpha <- alpha
  }
  seq.lam <- numeric(length(seq.alpha))
  seq.df <- numeric(length(seq.alpha))
  seq.cvll <- numeric(length(seq.alpha))
  if(is.null(foldid)) {
    foldid <- sample(rep(1:nfolds, length.out=n))  
  }
  for(a in 1:length(seq.alpha)) {
    fit <- glmnet(x, y, family="binomial", alpha=seq.alpha[a], dfmax=psel, standardize=standardize, 
                  intercept=intercept, nlambda=nlambda)
    cv.fit <- cv.glmnet(x, y, family="binomial", alpha=seq.alpha[a], lambda=fit$lambda, standardize=standardize,
                        intercept=intercept, foldid=foldid)
    seq.lam[a] <- cv.fit$lambda[which.min(cv.fit$cvm)]
    seq.df[a] <- fit$df[which.min(cv.fit$cvm)]
    seq.cvll[a] <- min(cv.fit$cvm)
    
    # cv.fit <- cv.glmnet(x, y, family="binomial", alpha=seq.alpha[a], 
    #                     lambda=tail(fit$lambda, n=2L), standardize=standardize, 
    #                     intercept=intercept, foldid=foldid)
    # seq.lam[a] <- tail(fit$lambda, n=1L)
    # seq.df[a] <- tail(fit$df, n=1L)
    # seq.cvll[a] <- tail(cv.fit$cvm, n=1L)
  }
  
  # cross-validated penalty parameters
  alpha <- seq.alpha[which.min(seq.cvll)]
  lambda <- seq.lam[which.min(seq.cvll)]
  lambda1 <- n*alpha*lambda # the n term is due to a 1/n term in the objective function in glmnet
  lambda2 <- n*(1 - alpha)*lambda # the n term is due to a 1/n term in the objective function in glmnet

  if(trace) {print(paste("Global lambda1 and lambda2: [", round(lambda1, 3), " ", round(lambda2, 3), "]"))}
  
  # starting values
  fit <- glmnet(x, y, family="binomial", alpha=alpha, lambda=lambda, standardize=standardize, intercept=intercept)
  mu <- muold <- as.numeric(coef(fit))
  
  # calculation of starting value for sigma
  xmu <- xaug %*% mu
  phat <- as.numeric(exp(xmu)/(1 + exp(xmu)))
  w <- phat*(1 - phat)

  if(intercept) {
    invtrw <- 1/sum(w)
    Wadj <- diag(w) - invtrw*as.matrix(w) %*% t(as.matrix(w))
    Ainv <- 0.5/lambda2*diag(p) - 0.25/lambda2^2*t(x) %*% Wadj %*% solve(diag(n) + (0.5/lambda2)*x %*% t(x) %*% Wadj) %*% x
    xainvxw <- x %*% Ainv %*% t(x) %*% as.matrix(w)
    xainv <- x %*% Ainv
    sigma <- matrix(0, nrow=p + 1, ncol=p + 1)
    sigma[1, 1] <- invtrw + invtrw^2*t(xainvxw) %*% Wadj %*% xainvxw
    sigma[1, 2:(p + 1)] <- sigma[2:(p + 1), 1] <- -invtrw*t(xainv) %*% Wadj %*% xainvxw
    sigma[2:(p + 1), 2:(p + 1)] <- t(xainv) %*% Wadj %*% xainv
    sigmaold <- sigma
  } else {
    Xw <- t(t(x)*sqrt(w))
    svdxw <- svd(Xw)
    U <- svdxw$u
    V <- svdxw$v
    d <- svdxw$d
    invmat <- 1/(d^2 + 2*(lambda1 + lambda2))
    part1 <- invmat^2*d^2
    sigma <- sigmaold <- t(t(V)*part1) %*% t(V)
  }
  
  # starting values ci and chi
  ci <- ciold <- as.numeric(sqrt(colSums(t(xaug) * (sigma %*% t(xaug))) + (colSums(t(xaug)*mu))^2))
  chi <- chiold <- as.numeric(lambda2*(diag(sigma)[(intercept + 1):(p + intercept)] + 
                                         mu[(intercept + 1):(p + intercept)]^2))
  
  lambdag <- lambdagold <- rep(1, G)
  lambda2g <- lambda2gold <- rep(lambda2, G)
  lambda2vec <- rep(lambda2*lambdag, times=sizes)
  phi <- lambda1^2/(4*lambda2vec)
  
  seq.mll <- mllold <- s <- seq.s <- niter1 <- seq.val <- 0
  conv1 <- FALSE
  
  while(!conv1 & (niter1 < maxiter)) {
    
    niter1 <- niter1 + 1
    
    if(trace) {
      cat("\r", "Iteration: ", niter1, ", lambda multipliers: [", paste(round(lambdag, 6), collapse=", "), 
          "]", sep="")
    }
    
    # VBEM updates (until convergence)
    conv2 <- FALSE
    niter2 <- 0
    while(!conv2 & (niter2 < maxiter)) {
      
      niter2 <- niter2 + 1
      
      new.param <- est_param(x, kappa, m, n, p, ciold, phi, chiold, lambda2vec, intercept)
      sigma <- new.param$sigma
      mu <- as.numeric(new.param$mu)
      ci <- as.numeric(new.param$ci)
      chi <- as.numeric(new.param$chi)
      
      conv2 <- max(abs(c(sigma - sigmaold, mu - muold))) < epsilon
      
      sigmaold <- sigma
      muold <- mu
      ciold <- ci
      chiold <- chi
      
    }
    
    # fixed parameters needed in mll calculation
    e.beta <- mu[(intercept + 1):(p + intercept)]
    v.beta <- diag(sigma)[(intercept + 1):(p + intercept)]
    e.psi.inv <- sqrt(phi/chi)
    e.psi <- 1/phi + sqrt(chi/phi)
    vec1 <- e.psi + 1
    vec2 <- (v.beta + e.beta^2)*(1 + e.psi.inv)
    sum1 <- sapply(1:G, function(g) {sum(vec1[groupid==g])})
    sum2 <- sapply(1:G, function(g) {sum(vec2[groupid==g])})
    
    # estimate new lambda multipliers
    # opt <- optim(lambda2gold, marg.ll, method="Nelder-Mead", control=list(maxit=10000),
    #              lambda1=lambda1, sum1=sum1, sum2=sum2, sizes=sizes)
    # lambda2g <- opt$par
    # opt <- optim(c(s, lambdagold), magnitude3, method="Nelder-Mead",
    #              control=list(maxit=10000),
    #              lambda1=lambda1, lambda2=lambda2, sum1=sum1, sum2=sum2, sizes=sizes)
    # lambdag <- opt$par[-1]
    # s <- opt$par[1]
    # s <- uniroot(f.uniroot, interval=c(-1000, 1000), lambda1=lambda1, lambda2=lambda2, sum1=sum1, sum2=sum2)$root
    # lambdag <- (-s + sqrt(s^2 + lambda1^2*sum1*sum2/4))/(lambda2*sum2)
    # lambda2gold <- lambdagold*lambda2
    # opt <- optim(c(s, lambda2gold), magnitude2, method="Nelder-Mead", 
    #              control=list(maxit=10000),
    #              lambda1=lambda1, lambda2=lambda2, sum1=sum1, sum2=sum2, sizes=sizes)
    # lambda2g <- opt$par[-1]
    # lambdag <- lambda2g/lambda2
    # s <- opt$par[1]
    opt <- optim(c(s, lambdagold), magnitude, method="Nelder-Mead",
                 control=list(maxit=10000),
                 lambda1=lambda1, lambda2=lambda2, sum1=sum1, sum2=sum2, sizes=sizes)
    lambdag <- opt$par[-1]
    s <- opt$par[1]
    seq.s <- c(seq.s, s)
    seq.val <- c(seq.val, opt$value)
    # opt <- optim(c(s, lambdagold), magnitude, method="L-BFGS-B", lower=c(-Inf, rep(0.001, times=G)),
    #              upper=rep(Inf, times=G + 1), control=list(maxit=10000),
    #              lambda1=lambda1, lambda2=lambda2, sum1=sum1, sum2=sum2, sizes=sizes)
    # lambdag <- opt$par[-1]
    # s <- opt$par[1]
    # opt <- solnp(c(s, lambdagold), magnitude, ineqfun=f.ineqcon, ineqLB=c(-Inf, rep(0, G)),
    #              ineqUB=rep(Inf, G + 1),
    #              lambda1=lambda1, lambda2=lambda2, sum1=sum1, sum2=sum2, sizes=sizes)
    # lambdag <- opt$pars[-1]
    # opt <- optim(lambdagold, f.nmll, 
    #              lambda1=lambda1, lambda2=lambda2, sum1=sum1, sum2=sum2, 
    #              sizes=sizes)
    # lambdag <- opt$par
    
    # calculate the marginal likelihood of the model
    mll <- p*log(lambda1) - 0.5*lambda2*sum(lambdag*sum2) - lambda1^2/(8*lambda2)*sum(sum1/lambdag) -
      sum(sizes*pnorm(-lambda1/sqrt(4*lambda2*lambdag), log.p=TRUE))
    seq.mll <- c(seq.mll, mll)
    
    # check the convergence of the marginal likelihood
    #conv1 <- abs(seq.mll[niter1 + 1] - seq.mll[niter1]) < epsilon
    conv1 <- max(abs(c(lambdag - lambdagold))) < epsilon
    
    # update old penalty parameters and mll to new ones
    lambdagold <- lambdag
    lambda2vec <- rep(lambdagold*lambda2, times=sizes)
    mll <- mllold
  
  }
  
  # lambda2g <- lambda2*lambdag
  unpenal <- ifelse(intercept, "~1", "~0")
  final.fit <- penalized(y, x, unpenalized=as.formula(unpenal), lambda1=lambda1, lambda2=rep(lambda2g, sizes),
                         model="logistic", standardize=standardize)
  
  out <- list(niter=niter1, conv=conv1, lambda1=lambda1, lambda2=lambda2g, final.fit=final.fit, mu=mu, 
              sigma=sigma, ci=ci, chi=chi, mll=seq.mll)
  return(out)
  
}

### test simulation
set.seed(123)
n <- 100
ntest <- 1000
p <- 1000
G <- 2
groupid <- rep(c(1, 2), each=p/G)
groupid2 <- CreatePartition(as.factor(groupid))
x <- matrix(rnorm(n*p), ncol=p, nrow=n)
beta <- rep(c(0, 2), each=p/G)
prob <- exp(x %*% beta)/(1 + exp(x %*% beta))
y <- as.numeric(runif(n) < prob)

test1.10 <- enhybrid(x, y, groupid, intercept=TRUE, psel=10, alpha=NULL, standardize=FALSE, 
                     nlambda=1000, nfolds=10, foldid=NULL, maxiter=1000, epsilon=1e-07, trace=TRUE)
test1.20 <- enhybrid(x, y, groupid, intercept=TRUE, psel=20, alpha=NULL, standardize=FALSE, 
                     nlambda=1000, nfolds=10, foldid=NULL, maxiter=1000, epsilon=1e-07, trace=TRUE)
test1.30 <- enhybrid(x, y, groupid, intercept=TRUE, psel=30, alpha=NULL, standardize=FALSE, 
                     nlambda=1000, nfolds=10, foldid=NULL, maxiter=1000, epsilon=1e-07, trace=TRUE)
test1.40 <- enhybrid(x, y, groupid, intercept=TRUE, psel=40, alpha=NULL, standardize=FALSE, 
                     nlambda=1000, nfolds=10, foldid=NULL, maxiter=1000, epsilon=1e-07, trace=TRUE)
test1.50 <- enhybrid(x, y, groupid, intercept=TRUE, psel=50, alpha=NULL, standardize=FALSE, 
                     nlambda=1000, nfolds=10, foldid=NULL, maxiter=1000, epsilon=1e-07, trace=TRUE)
test2.10 <- grridge(t(x), y, partitions=list(groupid=groupid2), unpenal=~1, selectionEN=TRUE, 
                    maxsel=10)
test2.20 <- grridge(t(x), y, partitions=list(groupid=groupid2), unpenal=~1, selectionEN=TRUE, 
                    maxsel=20)
test2.30 <- grridge(t(x), y, partitions=list(groupid=groupid2), unpenal=~1, selectionEN=TRUE, 
                    maxsel=30)
test2.40 <- grridge(t(x), y, partitions=list(groupid=groupid2), unpenal=~1, selectionEN=TRUE, 
                    maxsel=40)
test2.50 <- grridge(t(x), y, partitions=list(groupid=groupid2), unpenal=~1, selectionEN=TRUE, 
                    maxsel=50)
test3.10 <- enfit(x, y, psel=10, intercept=TRUE, standardize=FALSE, alpha=NULL, nlambda=1000,
                  nfolds=10, foldid=NULL)
test3.20 <- enfit(x, y, psel=20, intercept=TRUE, standardize=FALSE, alpha=NULL, nlambda=1000,
                  nfolds=10, foldid=NULL)
test3.30 <- enfit(x, y, psel=30, intercept=TRUE, standardize=FALSE, alpha=NULL, nlambda=1000,
                  nfolds=10, foldid=NULL)
test3.40 <- enfit(x, y, psel=40, intercept=TRUE, standardize=FALSE, alpha=NULL, nlambda=1000,
                  nfolds=10, foldid=NULL)
test3.50 <- enfit(x, y, psel=50, intercept=TRUE, standardize=FALSE, alpha=NULL, nlambda=1000,
                  nfolds=10, foldid=NULL)

# test4 <- enfit(x, y, psel=50, intercept=TRUE, standardize=FALSE, alpha=NULL, nlambda=1000,
#                nfolds=10, foldid=NULL)
# test5 <- penalized(y, x, lambda1=test4$lambda*n*test4$alpha, 
#                    lambda2=0.5*test4$lambda*n*(1 - test4$alpha), model="logistic", standardize=FALSE)

# sel1 <- which(test1$final.fit@penalized!=0)
# sel2 <- test2$whichsel
# sel3 <- which(as.numeric(coef(test3, s=test3$lambda[which.min(test3$cvm)]))[-1]!=0)
# sel4 <- which(as.numeric(coef(test4))[-1]!=0)
# sel5 <- which(as.numeric(test5@penalized)!=0)
# 
# sel1pg <- c(sum(sel1 %in% c(1:(p/G))), sum(sel1 %in% c((p/G + 1):(2*p/G))))
# sel2pg <- c(sum(sel2 %in% c(1:(p/G))), sum(sel2 %in% c((p/G + 1):(2*p/G))))
# sel3pg <- c(sum(sel3 %in% c(1:(p/G))), sum(sel3 %in% c((p/G + 1):(2*p/G))))
# sel4pg <- c(sum(sel4 %in% c(1:(p/G))), sum(sel4 %in% c((p/G + 1):(2*p/G))))
# sel5pg <- c(sum(sel5 %in% c(1:(p/G))), sum(sel5 %in% c((p/G + 1):(2*p/G))))

nsel1 <- nsel2 <- nsel3 <- matrix(NA, ncol=2, nrow=5)
wsel1 <- wsel2 <- wsel3 <- vector(mode="list", length=5)
for(m in 1:5) {
  test1 <- get(paste0("test1.", m*10))
  test2 <- get(paste0("test2.", m*10))
  test3 <- get(paste0("test3.", m*10))
  wsel1[[m]] <- which(test1$final.fit@penalized!=0)
  wsel2[[m]] <- test2$resEN$whichEN
  wsel3[[m]] <- which(as.numeric(coef(test3, s=test3$lambda[which.min(test3$cvm)]))[-1]!=0)
  nsel1[m, ] <- c(sum(wsel1[[m]] %in% c(1:(p/G))), sum(wsel1[[m]] %in% c((p/G + 1):(2*p/G))))
  nsel2[m, ] <- c(sum(wsel2[[m]] %in% c(1:(p/G))), sum(wsel2[[m]] %in% c((p/G + 1):(2*p/G))))
  nsel3[m, ] <- c(sum(wsel3[[m]] %in% c(1:(p/G))), sum(wsel3[[m]] %in% c((p/G + 1):(2*p/G))))
}

xtest <- matrix(rnorm(ntest*p), ncol=p, nrow=ntest)
probtest <- exp(xtest %*% beta)/(1 + exp(xtest %*% beta))
ytest <- as.numeric(runif(ntest) < probtest)

# fit1 <- cv.glmnet(x[, sel1], y, family="binomial", alpha=0, standardize=FALSE, intercept=TRUE)
# fit2 <- cv.glmnet(x[, sel2], y, family="binomial", alpha=0, standardize=FALSE, intercept=TRUE)
# fit3 <- cv.glmnet(x[, sel3], y, family="binomial", alpha=0, standardize=FALSE, intercept=TRUE)
# fit4 <- cv.glmnet(x[, sel4], y, family="binomial", alpha=0, standardize=FALSE, intercept=TRUE)
# fit5 <- cv.glmnet(x[, sel5], y, family="binomial", alpha=0, standardize=FALSE, intercept=TRUE)
# 
# pred1 <- as.numeric(predict(fit1, xtest[, sel1], s="lambda.min"))
# pred2 <- as.numeric(predict(fit2, xtest[, sel2], s="lambda.min"))
# pred3 <- as.numeric(predict(fit3, xtest[, sel3], s="lambda.min"))
# pred4 <- as.numeric(predict(fit4, xtest[, sel4], s="lambda.min"))
# pred5 <- as.numeric(predict(fit5, xtest[, sel5], s="lambda.min"))

auc1 <- auc2 <- auc3 <- brier1 <- brier2 <- brier3 <- numeric(5)
for(m in 1:5) {
  lpred1 <- cbind(1, xtest[, wsel1[[m]]]) %*% coef(glm(y ~ 1 + x[, wsel1[[m]]])) 
  lpred2 <- cbind(1, xtest[, wsel2[[m]]]) %*% coef(glm(y ~ 1 + x[, wsel2[[m]]])) 
  lpred3 <- cbind(1, xtest[, wsel3[[m]]]) %*% coef(glm(y ~ 1 + x[, wsel3[[m]]])) 
  pred1 <- as.numeric(exp(lpred1)/(1 + exp(lpred1)))
  pred2 <- as.numeric(exp(lpred1)/(1 + exp(lpred2)))
  pred3 <- as.numeric(exp(lpred1)/(1 + exp(lpred3)))
  auc1[m] <- pROC::auc(ytest, pred1)
  auc2[m] <- pROC::auc(ytest, pred2)
  auc3[m] <- pROC::auc(ytest, pred3)
  brier1[m] <- sum((pred1 - probtest)^2)
  brier2[m] <- sum((pred2 - probtest)^2)
  brier3[m] <- sum((pred3 - probtest)^2)
}

# pred1 <- predict(glm(ytest ~ 1 + xtest[, sel1]), as.data.frame(xtest))
# pred2 <- predict(glm(ytest ~ 1 + xtest[, sel2]), as.data.frame(xtest))
# pred3 <- predict(glm(ytest ~ 1 + xtest[, sel3]), as.data.frame(xtest))
# pred4 <- predict(glm(ytest ~ 1 + xtest[, sel4]), as.data.frame(xtest))
# pred5 <- predict(glm(ytest ~ 1 + xtest[, sel5]), as.data.frame(xtest))
# 
# auc1 <- pROC::auc(ytest, pred1)
# auc2 <- pROC::auc(ytest, pred2)
# auc3 <- pROC::auc(ytest, pred3)
# auc4 <- pROC::auc(ytest, pred4)
# auc5 <- pROC::auc(ytest, pred5)
# 
# brier1 <- sum((pred1 - probtest)^2)
# brier2 <- sum((pred2 - probtest)^2)
# brier3 <- sum((pred3 - probtest)^2)
# brier4 <- sum((pred4 - probtest)^2)
# brier5 <- sum((pred5 - probtest)^2)

# tab2 <- data.frame(method=c("ENhybrid", "GRridge", "lasso", "glmnet", "penalized"), 
#                    total=rowSums(rbind(sel1pg, sel2pg, sel3pg, sel4pg, sel5pg)), 
#                    group1=rbind(sel1pg, sel2pg, sel3pg, sel4pg, sel5pg)[, 1],
#                    group2=rbind(sel1pg, sel2pg, sel3pg, sel4pg, sel5pg)[, 2],
#                    AUC=c(auc1, auc2, auc3, auc4, auc5),
#                    brier=c(brier1, brier2, brier3, brier4, brier5)) 
# rownames(tab2) <- NULL
# tab1 # n <- 100, p <- 1000, beta <- c(0, 0.1)
# tab2 # n <- 100, p <- 1000, beta <- c(0, 2)

plot(seq(10, 50, 10), auc1, type="l", ylim=range(c(auc1, auc2, auc3)))
lines(seq(10, 50, 10), auc2, type="l", col=2)
lines(seq(10, 50, 10), auc3, type="l", col=3)

### test data
data(dataVerlaat)
groupid <- as.numeric(CpGann)[which(CpGann %in% c("Distant", "Island"))]
x <- t(datcenVerlaat)[, which(CpGann %in% c("Distant", "Island"))]
x <- x[, order(groupid)]
x <- apply(x, 2, function(j) {(j - mean(j)/sd(j))})
groupid <- sort(groupid)
groupid2 <- CreatePartition(as.factor(groupid))
y <- respVerlaat
# groupid=1 is distant, groupid=2 is island

test1 <- enhybrid(x, y, groupid, intercept=TRUE, psel=50, alpha=NULL, standardize=FALSE, 
                  nlambda=1000, nfolds=10, foldid=NULL, maxiter=1000, epsilon=1e-07, trace=TRUE)
test2 <- grridge(t(x), y, partitions=list(groupid=groupid2), unpenal=~1, selectionForward=TRUE, 
                 maxsel=52)
test3 <- cv.glmnet(x, y, family="binomial", alpha=1, dfmax=52, nlambda=1000, intercept=TRUE, 
                   standardize=FALSE)
test4 <- enfit(x, y, psel=52, intercept=TRUE, standardize=FALSE, alpha=NULL, nlambda=1000,
               nfolds=10, foldid=NULL)
test5 <- penalized(y, x, lambda1=test4$lambda*n*test4$alpha, 
                   lambda2=0.5*test4$lambda*n*(1 - test4$alpha), model="logistic", standardize=FALSE)

sel1 <- which(as.numeric(c(test1$final.fit@unpenalized, test1$final.fit@penalized))!=0)
sel2 <- test2$whichsel
sel3 <- which(as.numeric(coef(test3, s=test3$lambda[which.min(test4$cvm)]))!=0)
sel4 <- which(as.numeric(coef(test4))!=0)

sel2 %in% sel1
sel4 %in% sel1









library(lattice)
set.seed(123)
n <- 100
p <- 500
x <- matrix(rnorm(n*p), ncol=p, nrow=n)
beta <- rep(0.05, p)
prob <- exp(x %*% beta)/(1 + exp(x %*% beta))
y <- as.numeric(runif(n) < prob)
foldid <- sample(rep(1:30, length.out=n)) 

fit1 <- glmnet(x, y, family="binomial", alpha=0.01, nlambda=1000, dfmax=50, standardize=TRUE, intercept=TRUE)
fit2 <- glmnet(x, y, family="binomial", alpha=0.99, nlambda=1000, dfmax=50, standardize=TRUE, intercept=TRUE)
# seq.lam1 <- seq(min(c(fit1$lambda*0.01, fit2$lambda*0.99)),
#                 max(c(fit1$lambda*0.01, fit2$lambda*0.99)), length.out=30)
# seq.lam2 <- seq(min(c(0.5*fit1$lambda*0.99, 0.5*fit2$lambda*0.01)),
#                 max(c(0.5*fit1$lambda*0.99, 0.5*fit2$lambda*0.01)), length.out=30)
seq.lam1 <- seq(min(c(fit1$lambda*0.01, fit2$lambda*0.99)),
                0.102, length.out=30)
seq.lam2 <- seq(min(c(0.5*fit1$lambda*0.99, 0.5*fit2$lambda*0.01)),
                4.46, length.out=30)
mat.cvll <- mat.df <- matrix(NA, nrow=30, ncol=30)


for(l1 in 1:30) { # number of columns in the matrix
  for(l2 in 1:30) { # rows in the matrix
    alpha <- seq.lam1[l1]/(2*seq.lam2[l2] + seq.lam1[l1])
    lambda <- 2*seq.lam2[l2] + seq.lam1[l1]
    cv.fit <- cv.glmnet(x, y, family="binomial", alpha=alpha, lambda=c(lambda, lambda - 0.001), 
                     standardize=TRUE, intercept=TRUE, foldid=foldid)
    mat.cvll[l2, l1] <- cv.fit$cvm[2]
    mat.df[l2, l1] <- cv.fit$nzero[2]
  }
}

rownames(mat.cvll) <- round(seq.lam2, 3)
colnames(mat.cvll) <- round(seq.lam1, 3)
mat.df[is.na(mat.df)] <- 0
pts <- which(mat.df==48 | mat.df==49 | mat.df==50 | mat.df==51 | mat.df==52 , arr.ind=TRUE)
windows()
theseCol <- heat.colors(150)
# at.levels <- exp(seq(0.5*log(min(mat.cvll[-c(1, 17:30), -c(17:30)])), 
#                      0.5*log(as.numeric(quantile(mat.cvll[-c(1, 17:30), -c(17:30)], probs=0.8))), 
#                      length.out=50))^2
# heatp <- levelplot(mat.cvll[-c(1, 17:30), -c(1, 17:30)], col.regions=theseCol, region=TRUE, 
#                    at=at.levels, ylab=list(label=expression(lambda[1]), rot=-360), 
#                    xlab=expression(lambda[2]), scales=list(x=list(rot=45)))
at.levels <- exp(seq(0.5*log(min(mat.cvll)), 
                     0.5*log(as.numeric(quantile(mat.cvll, probs=0.85))), 
                     length.out=50))^2
heatp <- levelplot(mat.cvll, col.regions=theseCol, region=TRUE, 
                   at=at.levels, ylab=list(label=expression(lambda[1]), rot=-360), 
                   xlab=expression(lambda[2]), scales=list(x=list(rot=45)))
print(heatp)
trellis.focus("panel", 1, 1, highlight=FALSE)
llines(pts[, 2], pts[, 1], col=1)
trellis.unfocus()

wirep <- wireframe(mat.cvll, drape=TRUE, colorkey=FALSE, col.regions=theseCol, at=at.levels)
windows()                
print(wirep)

