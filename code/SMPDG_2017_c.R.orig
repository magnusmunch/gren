<<<<<<< HEAD
##############################  preamble  #############################
# code belonging to abstract_SMPGD_2017_V01.pdf cloud version         #
# version: 02 (cloud)                                                 #
# author: Magnus Münch                                                #
# created: 15-11-2016                                                 #
# last edited: 30-11-2016                                             #
#######################################################################

###############################  notes  ###############################
# 30-11-2016: Now runs parallel, not tested                           #
# 15-11-2016: Building up from basic                                  #
#######################################################################

### paths
path.results <- "/home/ubuntu/ENVB/SMPGD_2017/results/"
path.code <- "/home/ubuntu/ENVB/code/"

### libraries
library(penalized)
library(GRridge)
library(pROC)
library(Rcpp)
library(snowfall)

### functions
# source function to estimate parameters in C++
sourceCpp(paste(path.code, "ENVB2.cpp", sep=""))

# below three functions for marginal likelihood calculations
marg.ll1.2g <- function(lambda, e.beta, v.beta, e.psi.inv, e.psi, p, G, sizes, modmat) {
  # marginal likelihood in lambda1 and lambda2g
  
  lambda1 <- lambda[1]
  lambda2 <- lambda[2:(G + 1)]
  lambda2vec <- rep(lambda2, times=sizes)
  part1 <- p*log(lambda1)
  part2 <- 0.5*sum(lambda2vec*(v.beta + e.beta^2)*(1 + e.psi.inv))
  part3 <- lambda1^2*sum((e.psi + 1)/lambda2vec)/8
  part4 <- sum(sizes*pnorm(-lambda1/(sqrt(4*lambda2)), log.p=TRUE))
  ll <- part1 - part2 - part3 - part4
  
  return(ll)
  
}

gr.marg.ll1.2g <- function(lambda, e.beta, v.beta, e.psi.inv, e.psi, p, G, sizes, modmat) {
  # gradient of marginal likelihood in lambda1 and lambda2g
  
  lambda1 <- lambda[1]
  lambda2 <- lambda[2:(G + 1)]
  lambda2vec <- rep(lambda2, times=sizes)
  part1.1 <- p/lambda1
  part1.2 <- lambda1*sum((e.psi + 1)/lambda2vec)/4
  part1.3 <- sum(sizes*dnorm(lambda1/sqrt(4*lambda2))/
                   (sqrt(lambda2)*pnorm(-lambda1/sqrt(4*lambda2))))/sqrt(4)
  comp1 <- part1.1 - part1.2 + part1.3
  part2.1 <- lambda1^2*(t(as.matrix(e.psi + 1)) %*% modmat)/(8*lambda2^2)
  part2.2 <- 0.5*as.numeric((t(as.matrix((v.beta + e.beta^2)*(1 + e.psi.inv))) %*% modmat))
  part2.3 <- lambda1*sizes*dnorm(lambda1/sqrt(4*lambda2))/
    (4*lambda2^(1.5)*pnorm(-lambda1/sqrt(4*lambda2)))
  comp2 <- part2.1 - part2.2 - part2.3
  return(c(comp1, comp2))
  
}

# the fitting function
envb2 <- function(x, y, groups, lambda1, lambda2, intercept=TRUE, maxiter=1000, epsilon=1e-07, trace=TRUE) {
  
  p <- ncol(x)
  n <- nrow(x)
  G <- length(unique(groups))
  sizes <- rle(groups)$lengths
  modmat <- matrix(0, ncol=G, nrow=p)
  modmat <- sapply(1:G, function(g) {as.numeric(groups==g)})
  m <- rep(1, n)
  kappa <- y - 0.5*m
  xaug <- x
  
  # starting values
  if(intercept) {
    fit <- penalized(y, x, unpenalized=~1, lambda1=0, lambda2=lambda1 + lambda2, model="logistic", trace=FALSE)
    xaug <- cbind(1, x)
  } else {
    fit <- penalized(y, x, unpenalized=~0, lambda1=0, lambda2=lambda1 + lambda2, model="logistic", trace=FALSE)
  }
  
  mu <- muold <- c(fit@unpenalized, fit@penalized)
  
  # calculation of starting value for sigma
  xmu <- xaug %*% mu
  phat <- as.numeric(exp(xmu)/(1 + exp(xmu)))
  
  if(intercept) {
    w <- phat*(1 - phat)
    invtrw <- 1/sum(w)
    Wadj <- diag(w) - invtrw*as.matrix(w) %*% t(as.matrix(w))
    Ainv <- 0.5/lambda2*diag(p) - 0.25/lambda2^2*t(x) %*% Wadj %*% solve(diag(n) + (0.5/lambda2)*x %*% t(x) %*% Wadj) %*% x
    xainvxw <- x %*% Ainv %*% t(x) %*% as.matrix(w)
    xainv <- x %*% Ainv
    sigma <- sigmaold <- matrix(0, nrow=p + 1, ncol=p + 1)
    sigma[1, 1] <- invtrw + invtrw^2*t(xainvxw) %*% Wadj %*% xainvxw
    sigma[1, 2:(p + 1)] <- sigma[2:(p + 1), 1] <- -invtrw*t(xainv) %*% Wadj %*% xainvxw
    sigma[2:(p + 1), 2:(p + 1)] <- t(xainv) %*% Wadj %*% xainv
  } else {
    W <- diag(sqrt(phat*(1 - phat)))
    Xw <- W %*% x
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
  
  lambda2old <- rep(lambda2, G)
  lambda1old <- lambda1
  lambda2vec <- rep(lambda2old, time=sizes)
  phi <- lambda1^2/(4*lambda2vec)
  
  conv <- FALSE
  niter <- 0
  
  while(!conv & (niter < maxiter)) {
    
    niter <- niter + 1
    
    if(trace) {
      cat("\r", "Iteration: ", niter, ", ", 
          paste("lambda1=[", paste(round(lambda1old, 2), collapse=", "), 
                sep=""), "], ", 
          paste("lambda2=[", paste(round(lambda2old, 2), collapse=", "), 
                sep=""), "]", sep="")
    }
    
    # estimation of new model parameter (done in Cpp)
    new.param <- est_param(x, kappa, m, n, p, ciold, phi, chiold, lambda2vec, intercept)
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
    
    # recalculate phi, since it depends on lambda1 and lambda2
    phi <- lambda1old^2/(4*lambda2vec)
    
    # fixed parameters needed in function optimisation
    e.beta <- mu[(intercept + 1):(p + intercept)]
    v.beta <- diag(sigma)[(intercept + 1):(p + intercept)]
    e.psi.inv <- sqrt(phi/chi)
    e.psi <- 1/phi + sqrt(chi/phi)
    
    lambdaold <- c(lambda1old, lambda2old)
    opt.rout <- tryCatch({
      optim(par=lambdaold, fn=marg.ll1.2g, gr=gr.marg.ll1.2g,
            method="L-BFGS-B", lower=rep(0.001, G + 1), 
            upper=rep(Inf, G + 1), control=list(fnscale=-1), 
            e.beta=e.beta, v.beta=v.beta, e.psi.inv=e.psi.inv, 
            e.psi=e.psi, p=p, G=G, sizes=sizes, modmat=modmat)},
      error=function(war) {
        optim(par=lambdaold, fn=marg.ll1.2g,
              method="L-BFGS-B", lower=rep(0.001, G + 1), 
              upper=rep(Inf, G + 1), control=list(fnscale=-1), 
              e.beta=e.beta, v.beta=v.beta, e.psi.inv=e.psi.inv, 
              e.psi=e.psi, p=p, G=G, sizes=sizes, modmat=modmat)
      })
    lambda1 <- opt.rout$par[1]
    lambda2 <- opt.rout$par[2:(G + 1)]
    
    # update old hyperparameters to new ones
    lambda1old <- lambda1
    lambda2old <- lambda2
    lambda2vec <- rep(lambda2old, times=sizes)
    
  }
  
  out <- list(niter=niter, conv=conv, sigma=sigma, mu=mu, c=ci, chi=chi, lambda1=lambda1, lambda2=lambda2)
  return(out)
  
}

# function to sparsify coefficient vector in simulations
sparsify <- function(vec, frac){
  if(frac==0){
    return(vec)
  } else {
    N0 <- round(frac*length(vec))
    vecnew <- vec
    vecnew[1:N0] <- 0
    return(vecnew)
  }
}

### simulations
n <- 100           # Nb of observations    
ntest <- 1000          ## Nb of test set observations     
p <- 200            # Nb of variables per group
G <- 10            # Nb of groups
meanBeta <- 0.01   # Beta variances per group are VarBeta*(1:G); use this for CorX=0.5
CorX <- 0.5       # correlation within variable block
Nblock <- 10*G    # number of correlation blocks
settings <- c(ntrain=n, p=p, G=G, meanBeta=meanBeta, CorX=CorX, Nblock=Nblock)

nrep <- 50  #number of repeats per simulation setting
facvec <- c(1.3, 1.6, 2)   #tunes how much weaker each next group is. The '2' means that the second group is twice as weak as the first, etc
fractvec <- c(0, 0.7, 0.9) #tunes the sparsity per group. E.g. 0.9 means that 9/10 betas in a group are set to 0

aucmat <- matrix(NA, ncol=5, nrow=nrep)
briermat <- matrix(NA, ncol=5, nrow=nrep)
msemat <- matrix(NA, ncol=5, nrow=nrep)

vsettings <- expand.grid(fac=facvec, fract=fractvec)
vsettings <- vsettings[c(1, 5, 9, 2:4, 6:8), ]
parfun <- function(r) {
  sourceCpp(paste(path.code, "ENVB2.cpp", sep=""))
  fac <- vsettings$fac[r]
  fract <- vsettings$fract[r]
  print(paste("fac=", fac))
  print(paste("fract=", fract))
  for(reptit in 1:nrep){
    print(paste("repeat=", reptit))
    set.seed(reptit)
    reps <- rev(sapply(0:(G - 1), function(i) {fac^(-i)}))
    meansB <- rep(reps, each=p)*meanBeta/mean(rep(reps, each=p))
    Beta <- meansB
    Beta <- rev(sparsify(Beta, frac=fract))
    
    ### FITTING THE MODELS
    pblock <- G*p/Nblock
    grs <- rep(1:G, each=p)
    P <- G*p #Complete number of variables
    X <- Reduce(cbind, lapply(1:Nblock, function(z) {
      matrix(rep(rnorm(n, sd=sqrt(CorX/(1 - CorX))), times=pblock), n, pblock)})) + matrix(rnorm(n*G*p), n, G*p)
    X <- t((t(X) - apply(t(X), 1, mean))/apply(t(X), 1, sd))
    
    lpred <- X %*% Beta 
    logisticintercept <- 0
    prob <- 1/(1 + exp(-(lpred + logisticintercept)))
    Y <- rbinom(length(prob), 1, prob)
      
    # ENVB
    vbSim <- tryCatch({
      envb2(x=X, y=Y, groups=grs, lambda1=1, lambda2=500, maxiter=1000, epsilon=1e-06, intercept=TRUE, trace=FALSE)},
      error=function(war) {return(NULL)})
      
    # GRridge
    groups <- CreatePartition(grs, grsize=p, uniform=T, decreasing=F)
    partsim <- list(grouping=groups)
    grSim <- tryCatch({
      grridge(t(X), Y, unpenal=~1, partsim, savepredobj="all", innfold=10, method="stable")},
      error=function(war) {return(NULL)})
      
    # ridge
    rrSim <- tryCatch({optL2(Y, X, unpenalized=~1, lambda1=0, model="logistic", fold=10, trace=FALSE)},
                      error=function(war) {return(NULL)})
      
    # lasso
    lrSim <- tryCatch({optL1(Y, X, unpenalized=~1, lambda2=0, model="logistic", fold=10, trace=FALSE)},
                      error=function(war) {return(NULL)})
    
    # elastic net
    enSim <- tryCatch({optL2(Y, X, unpenalized=~1, lambda1=vbSim$lambda1, model="logistic", fold=10, trace=FALSE)},
                      error=function(war) {return(NULL)})
    
    # calculating mse
    if(is.null(vbSim)) {
      vbMse <- NA
    } else {
      vbMse <- var((as.numeric(vbSim$mu[-1]) - Beta)^2)
    }
    if(is.null(grSim)) {
      grMse <- NA
    } else {
      grMse <- var((grSim$betas - Beta)^2)
    }
    if(is.null(rrSim)) {
      rrMse <- NA
    } else {
      rrMse <- var((rrSim$fullfit@penalized - Beta)^2)
    }
    if(is.null(lrSim)) {
      lrMse <- NA
    } else {
      lrMse <- var((lrSim$fullfit@penalized - Beta)^2)
    }
    if(is.null(enSim)) {
      enMse <- NA
    } else {
      enMse <- var((enSim$fullfit@penalized - Beta)^2)
    }
    msemat[reptit, ] <- c(vbMse, grMse, rrMse, lrMse, enMse)
      
    ### TESTING THE MODELS
    # making the test data
    Xtest <- Reduce(cbind, lapply(1:Nblock, function(z) {
      matrix(rep(rnorm(ntest, sd=sqrt(CorX/(1 - CorX))), times=pblock), ntest, pblock)})) + 
      matrix(rnorm(ntest*G*p), ntest, G*p)
    lpredtest <- Xtest %*% Beta
    logisticintercept <- 0
    
    probtest <- as.numeric(1/(1 + exp(-(lpredtest + logisticintercept))))
    Ytest <- rbinom(length(probtest), 1, probtest)
    
    # making predictions, calculating AUC and mean of brier residuals
    cutoffs <- rev(seq(0, 1, by=0.005))
    if(is.null(vbSim)) {
      vbBrier <- NA
      vbAuc <- NA
    } else {
      vbPred <- as.numeric(1/(1 + exp(-(cbind(1, Xtest) %*% vbSim$mu))))
      vbBrier <- mean((vbPred - probtest)^2)
      vbRoc <- GRridge::roc(probs=as.numeric(vbPred), true=Ytest, cutoffs)
      vbAuc <- GRridge::auc(vbRoc)
    }
    if(is.null(grSim)) {
      grBrier <- NA
      grAuc <- NA
    } else {
      grPred <- predict.grridge(grSim, t(Xtest))
      grBrier <- mean((grPred[, 2] - probtest)^2)
      grRoc <- GRridge::roc(probs=as.numeric(grPred[, 2]), true=Ytest, cutoffs)
      grAuc <- GRridge::auc(grRoc)
    }
    if(is.null(rrSim)) {
      rrBrier <- NA
      rrAuc <- NA
    } else {
      rrPred <- predict(rrSim$fullfit, Xtest)
      rrBrier <- mean((rrPred - probtest)^2)
      rrRoc <- GRridge::roc(probs=as.numeric(rrPred), true=Ytest, cutoffs)
      rrAuc <- GRridge::auc(rrRoc)
    }
    if(is.null(lrSim)) {
      lrBrier <- NA
      lrAuc <- NA
    } else {
      lrPred <- predict(lrSim$fullfit, Xtest)
      lrBrier <- mean((lrPred - probtest)^2)
      lrRoc <- GRridge::roc(probs=as.numeric(lrPred), true=Ytest, cutoffs)
      lrAuc <- GRridge::auc(lrRoc)
    }
    if(is.null(enSim)) {
      enBrier <- NA
      enAuc <- NA
    } else {
      enPred <- predict(enSim$fullfit, Xtest)
      enBrier <- mean((enPred - probtest)^2)
      enRoc <- GRridge::roc(probs=as.numeric(enPred), true=Ytest, cutoffs)
      enAuc <- GRridge::auc(enRoc)
    }
      
    briermat[reptit, ] <- c(vbBrier, grBrier, rrBrier, lrBrier, enBrier)
    aucmat[reptit, ] <- c(vbAuc, grAuc, rrAuc, lrAuc, enAuc)
     
  }
    
  colnames(msemat) <- colnames(briermat) <- colnames(aucmat) <- c("ENVB", "GRridge", "RR", "LR", "EN")
  out <- list(msemat=msemat, briermat=briermat, aucmat=aucmat)
  return(out)
    
}   

#sfInit(parallel=TRUE, cpus=9)
#sfLibrary(penalized)
#sfLibrary(GRridge)
#sfLibrary(pROC)
#sfLibrary(Rcpp)
#sfLibrary(snowfall)
#sfExportAll()
#out <- sfSapply(1:nrow(vsettings), parfun)
out <- parfun(9)
#names(out) <- apply(vsettings, 1, paste, c("fac", "fract"), collapse=", ")
save(out, file=paste(path.results, "SMPDG_2017_c_v02_run10.Rdata", sep=""))
#sfStop()



#path.res <- "C:/Users/Magnus/Documents/phd/ENVB/abstract_SMPGD_2017/results/"
#load(paste(path.res, "SMPDG_2017_c_v02_run01.Rdata", sep=""))
=======
##############################  preamble  #############################
# code belonging to abstract_SMPGD_2017_V01.pdf cloud version         #
# version: 02 (cloud)                                                 #
# author: Magnus Münch                                                #
# created: 15-11-2016                                                 #
# last edited: 30-11-2016                                             #
#######################################################################

###############################  notes  ###############################
# 30-11-2016: Now runs parallel, not tested                           #
# 15-11-2016: Building up from basic                                  #
#######################################################################

### paths
path.results <- "/home/ubuntu/ENVB/SMPGD_2017/results/"
path.code <- "/home/ubuntu/ENVB/code/"

### libraries
library(penalized)
library(GRridge)
library(pROC)
library(Rcpp)
library(snowfall)

### functions
# source function to estimate parameters in C++
sourceCpp(paste(path.code, "ENVB2.cpp", sep=""))

# below three functions for marginal likelihood calculations
marg.ll1.2g <- function(lambda, e.beta, v.beta, e.psi.inv, e.psi, p, G, sizes, modmat) {
  # marginal likelihood in lambda1 and lambda2g
  
  lambda1 <- lambda[1]
  lambda2 <- lambda[2:(G + 1)]
  lambda2vec <- rep(lambda2, times=sizes)
  part1 <- p*log(lambda1)
  part2 <- 0.5*sum(lambda2vec*(v.beta + e.beta^2)*(1 + e.psi.inv))
  part3 <- lambda1^2*sum((e.psi + 1)/lambda2vec)/8
  part4 <- sum(sizes*pnorm(-lambda1/(sqrt(4*lambda2)), log.p=TRUE))
  ll <- part1 - part2 - part3 - part4
  
  return(ll)
  
}

gr.marg.ll1.2g <- function(lambda, e.beta, v.beta, e.psi.inv, e.psi, p, G, sizes, modmat) {
  # gradient of marginal likelihood in lambda1 and lambda2g
  
  lambda1 <- lambda[1]
  lambda2 <- lambda[2:(G + 1)]
  lambda2vec <- rep(lambda2, times=sizes)
  part1.1 <- p/lambda1
  part1.2 <- lambda1*sum((e.psi + 1)/lambda2vec)/4
  part1.3 <- sum(sizes*dnorm(lambda1/sqrt(4*lambda2))/
                   (sqrt(lambda2)*pnorm(-lambda1/sqrt(4*lambda2))))/sqrt(4)
  comp1 <- part1.1 - part1.2 + part1.3
  part2.1 <- lambda1^2*(t(as.matrix(e.psi + 1)) %*% modmat)/(8*lambda2^2)
  part2.2 <- 0.5*as.numeric((t(as.matrix((v.beta + e.beta^2)*(1 + e.psi.inv))) %*% modmat))
  part2.3 <- lambda1*sizes*dnorm(lambda1/sqrt(4*lambda2))/
    (4*lambda2^(1.5)*pnorm(-lambda1/sqrt(4*lambda2)))
  comp2 <- part2.1 - part2.2 - part2.3
  return(c(comp1, comp2))
  
}

# the fitting function
envb2 <- function(x, y, groups, lambda1, lambda2, intercept=TRUE, maxiter=1000, epsilon=1e-07, trace=TRUE) {
  
  p <- ncol(x)
  n <- nrow(x)
  G <- length(unique(groups))
  sizes <- rle(groups)$lengths
  modmat <- matrix(0, ncol=G, nrow=p)
  modmat <- sapply(1:G, function(g) {as.numeric(groups==g)})
  m <- rep(1, n)
  kappa <- y - 0.5*m
  xaug <- x
  
  # starting values
  if(intercept) {
    fit <- penalized(y, x, unpenalized=~1, lambda1=0, lambda2=lambda1 + lambda2, model="logistic", trace=FALSE)
    xaug <- cbind(1, x)
  } else {
    fit <- penalized(y, x, unpenalized=~0, lambda1=0, lambda2=lambda1 + lambda2, model="logistic", trace=FALSE)
  }
  
  mu <- muold <- c(fit@unpenalized, fit@penalized)
  
  # calculation of starting value for sigma
  xmu <- xaug %*% mu
  phat <- as.numeric(exp(xmu)/(1 + exp(xmu)))
  
  if(intercept) {
    w <- phat*(1 - phat)
    invtrw <- 1/sum(w)
    Wadj <- diag(w) - invtrw*as.matrix(w) %*% t(as.matrix(w))
    Ainv <- 0.5/lambda2*diag(p) - 0.25/lambda2^2*t(x) %*% Wadj %*% solve(diag(n) + (0.5/lambda2)*x %*% t(x) %*% Wadj) %*% x
    xainvxw <- x %*% Ainv %*% t(x) %*% as.matrix(w)
    xainv <- x %*% Ainv
    sigma <- sigmaold <- matrix(0, nrow=p + 1, ncol=p + 1)
    sigma[1, 1] <- invtrw + invtrw^2*t(xainvxw) %*% Wadj %*% xainvxw
    sigma[1, 2:(p + 1)] <- sigma[2:(p + 1), 1] <- -invtrw*t(xainv) %*% Wadj %*% xainvxw
    sigma[2:(p + 1), 2:(p + 1)] <- t(xainv) %*% Wadj %*% xainv
  } else {
    W <- diag(sqrt(phat*(1 - phat)))
    Xw <- W %*% x
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
  
  lambda2old <- rep(lambda2, G)
  lambda1old <- lambda1
  lambda2vec <- rep(lambda2old, time=sizes)
  phi <- lambda1^2/(4*lambda2vec)
  
  conv <- FALSE
  niter <- 0
  
  while(!conv & (niter < maxiter)) {
    
    niter <- niter + 1
    
    if(trace) {
      cat("\r", "Iteration: ", niter, ", ", 
          paste("lambda1=[", paste(round(lambda1old, 2), collapse=", "), 
                sep=""), "], ", 
          paste("lambda2=[", paste(round(lambda2old, 2), collapse=", "), 
                sep=""), "]", sep="")
    }
    
    # estimation of new model parameter (done in Cpp)
    new.param <- est_param(x, kappa, m, n, p, ciold, phi, chiold, lambda2vec, intercept)
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
    
    # recalculate phi, since it depends on lambda1 and lambda2
    phi <- lambda1old^2/(4*lambda2vec)
    
    # fixed parameters needed in function optimisation
    e.beta <- mu[(intercept + 1):(p + intercept)]
    v.beta <- diag(sigma)[(intercept + 1):(p + intercept)]
    e.psi.inv <- sqrt(phi/chi)
    e.psi <- 1/phi + sqrt(chi/phi)
    
    lambdaold <- c(lambda1old, lambda2old)
    opt.rout <- tryCatch({
      optim(par=lambdaold, fn=marg.ll1.2g, gr=gr.marg.ll1.2g,
            method="L-BFGS-B", lower=rep(0.001, G + 1), 
            upper=rep(Inf, G + 1), control=list(fnscale=-1), 
            e.beta=e.beta, v.beta=v.beta, e.psi.inv=e.psi.inv, 
            e.psi=e.psi, p=p, G=G, sizes=sizes, modmat=modmat)},
      error=function(war) {
        optim(par=lambdaold, fn=marg.ll1.2g,
              method="L-BFGS-B", lower=rep(0.001, G + 1), 
              upper=rep(Inf, G + 1), control=list(fnscale=-1), 
              e.beta=e.beta, v.beta=v.beta, e.psi.inv=e.psi.inv, 
              e.psi=e.psi, p=p, G=G, sizes=sizes, modmat=modmat)
      })
    lambda1 <- opt.rout$par[1]
    lambda2 <- opt.rout$par[2:(G + 1)]
    
    # update old hyperparameters to new ones
    lambda1old <- lambda1
    lambda2old <- lambda2
    lambda2vec <- rep(lambda2old, times=sizes)
    
  }
  
  out <- list(niter=niter, conv=conv, sigma=sigma, mu=mu, c=ci, chi=chi, lambda1=lambda1, lambda2=lambda2)
  return(out)
  
}

# function to sparsify coefficient vector in simulations
sparsify <- function(vec, frac){
  if(frac==0){
    return(vec)
  } else {
    N0 <- round(frac*length(vec))
    vecnew <- vec
    vecnew[1:N0] <- 0
    return(vecnew)
  }
}

### simulations
n <- 100           # Nb of observations    
ntest <- 1000          ## Nb of test set observations     
p <- 200            # Nb of variables per group
G <- 10            # Nb of groups
meanBeta <- 0.01   # Beta variances per group are VarBeta*(1:G); use this for CorX=0.5
CorX <- 0.5       # correlation within variable block
Nblock <- 10*G    # number of correlation blocks
settings <- c(ntrain=n, p=p, G=G, meanBeta=meanBeta, CorX=CorX, Nblock=Nblock)

nrep <- 50  #number of repeats per simulation setting
facvec <- c(1.3, 1.6, 2)   #tunes how much weaker each next group is. The '2' means that the second group is twice as weak as the first, etc
fractvec <- c(0, 0.7, 0.9) #tunes the sparsity per group. E.g. 0.9 means that 9/10 betas in a group are set to 0

aucmat <- matrix(NA, ncol=5, nrow=nrep)
briermat <- matrix(NA, ncol=5, nrow=nrep)
msemat <- matrix(NA, ncol=5, nrow=nrep)

vsettings <- expand.grid(fac=facvec, fract=fractvec)
vsettings <- vsettings[c(1, 5, 9, 2:4, 6:8), ]
parfun <- function(r) {
  sourceCpp(paste(path.code, "ENVB2.cpp", sep=""))
  fac <- vsettings$fac[r]
  fract <- vsettings$fract[r]
  print(paste("fac=", fac))
  print(paste("fract=", fract))
  for(reptit in 1:nrep){
    print(paste("repeat=", reptit))
    set.seed(reptit)
    reps <- rev(sapply(0:(G - 1), function(i) {fac^(-i)}))
    meansB <- rep(reps, each=p)*meanBeta/mean(rep(reps, each=p))
    Beta <- meansB
    Beta <- rev(sparsify(Beta, frac=fract))
    
    ### FITTING THE MODELS
    pblock <- G*p/Nblock
    grs <- rep(1:G, each=p)
    P <- G*p #Complete number of variables
    X <- Reduce(cbind, lapply(1:Nblock, function(z) {
      matrix(rep(rnorm(n, sd=sqrt(CorX/(1 - CorX))), times=pblock), n, pblock)})) + matrix(rnorm(n*G*p), n, G*p)
    X <- t((t(X) - apply(t(X), 1, mean))/apply(t(X), 1, sd))
    
    lpred <- X %*% Beta 
    logisticintercept <- 0
    prob <- 1/(1 + exp(-(lpred + logisticintercept)))
    Y <- rbinom(length(prob), 1, prob)
      
    # ENVB
    vbSim <- tryCatch({
      envb2(x=X, y=Y, groups=grs, lambda1=1, lambda2=500, maxiter=1000, epsilon=1e-06, intercept=TRUE, trace=FALSE)},
      error=function(war) {return(NULL)})
      
    # GRridge
    groups <- CreatePartition(grs, grsize=p, uniform=T, decreasing=F)
    partsim <- list(grouping=groups)
    grSim <- tryCatch({
      grridge(t(X), Y, unpenal=~1, partsim, savepredobj="all", innfold=10, method="stable")},
      error=function(war) {return(NULL)})
      
    # ridge
    rrSim <- tryCatch({optL2(Y, X, unpenalized=~1, lambda1=0, model="logistic", fold=10, trace=FALSE)},
                      error=function(war) {return(NULL)})
      
    # lasso
    lrSim <- tryCatch({optL1(Y, X, unpenalized=~1, lambda2=0, model="logistic", fold=10, trace=FALSE)},
                      error=function(war) {return(NULL)})
    
    # elastic net
    enSim <- tryCatch({optL2(Y, X, unpenalized=~1, lambda1=vbSim$lambda1, model="logistic", fold=10, trace=FALSE)},
                      error=function(war) {return(NULL)})
    
    # calculating mse
    if(is.null(vbSim)) {
      vbMse <- NA
    } else {
      vbMse <- var((as.numeric(vbSim$mu[-1]) - Beta)^2)
    }
    if(is.null(grSim)) {
      grMse <- NA
    } else {
      grMse <- var((grSim$betas - Beta)^2)
    }
    if(is.null(rrSim)) {
      rrMse <- NA
    } else {
      rrMse <- var((rrSim$fullfit@penalized - Beta)^2)
    }
    if(is.null(lrSim)) {
      lrMse <- NA
    } else {
      lrMse <- var((lrSim$fullfit@penalized - Beta)^2)
    }
    if(is.null(enSim)) {
      enMse <- NA
    } else {
      enMse <- var((enSim$fullfit@penalized - Beta)^2)
    }
    msemat[reptit, ] <- c(vbMse, grMse, rrMse, lrMse, enMse)
      
    ### TESTING THE MODELS
    # making the test data
    Xtest <- Reduce(cbind, lapply(1:Nblock, function(z) {
      matrix(rep(rnorm(ntest, sd=sqrt(CorX/(1 - CorX))), times=pblock), ntest, pblock)})) + 
      matrix(rnorm(ntest*G*p), ntest, G*p)
    lpredtest <- Xtest %*% Beta
    logisticintercept <- 0
    
    probtest <- as.numeric(1/(1 + exp(-(lpredtest + logisticintercept))))
    Ytest <- rbinom(length(probtest), 1, probtest)
    
    # making predictions, calculating AUC and mean of brier residuals
    cutoffs <- rev(seq(0, 1, by=0.005))
    if(is.null(vbSim)) {
      vbBrier <- NA
      vbAuc <- NA
    } else {
      vbPred <- as.numeric(1/(1 + exp(-(cbind(1, Xtest) %*% vbSim$mu))))
      vbBrier <- mean((vbPred - probtest)^2)
      vbRoc <- GRridge::roc(probs=as.numeric(vbPred), true=Ytest, cutoffs)
      vbAuc <- GRridge::auc(vbRoc)
    }
    if(is.null(grSim)) {
      grBrier <- NA
      grAuc <- NA
    } else {
      grPred <- predict.grridge(grSim, t(Xtest))
      grBrier <- mean((grPred[, 2] - probtest)^2)
      grRoc <- GRridge::roc(probs=as.numeric(grPred[, 2]), true=Ytest, cutoffs)
      grAuc <- GRridge::auc(grRoc)
    }
    if(is.null(rrSim)) {
      rrBrier <- NA
      rrAuc <- NA
    } else {
      rrPred <- predict(rrSim$fullfit, Xtest)
      rrBrier <- mean((rrPred - probtest)^2)
      rrRoc <- GRridge::roc(probs=as.numeric(rrPred), true=Ytest, cutoffs)
      rrAuc <- GRridge::auc(rrRoc)
    }
    if(is.null(lrSim)) {
      lrBrier <- NA
      lrAuc <- NA
    } else {
      lrPred <- predict(lrSim$fullfit, Xtest)
      lrBrier <- mean((lrPred - probtest)^2)
      lrRoc <- GRridge::roc(probs=as.numeric(lrPred), true=Ytest, cutoffs)
      lrAuc <- GRridge::auc(lrRoc)
    }
    if(is.null(enSim)) {
      enBrier <- NA
      enAuc <- NA
    } else {
      enPred <- predict(enSim$fullfit, Xtest)
      enBrier <- mean((enPred - probtest)^2)
      enRoc <- GRridge::roc(probs=as.numeric(enPred), true=Ytest, cutoffs)
      enAuc <- GRridge::auc(enRoc)
    }
      
    briermat[reptit, ] <- c(vbBrier, grBrier, rrBrier, lrBrier, enBrier)
    aucmat[reptit, ] <- c(vbAuc, grAuc, rrAuc, lrAuc, enAuc)
     
  }
    
  colnames(msemat) <- colnames(briermat) <- colnames(aucmat) <- c("ENVB", "GRridge", "RR", "LR", "EN")
  out <- list(msemat=msemat, briermat=briermat, aucmat=aucmat)
  return(out)
    
}   

#sfInit(parallel=TRUE, cpus=9)
#sfLibrary(penalized)
#sfLibrary(GRridge)
#sfLibrary(pROC)
#sfLibrary(Rcpp)
#sfLibrary(snowfall)
#sfExportAll()
#out <- sfSapply(1:nrow(vsettings), parfun)
out <- parfun(8)
#names(out) <- apply(vsettings, 1, paste, c("fac", "fract"), collapse=", ")
save(out, file=paste(path.results, "SMPDG_2017_c_v02_run09.Rdata", sep=""))
#sfStop()



#path.res <- "C:/Users/Magnus/Documents/phd/ENVB/abstract_SMPGD_2017/results/"
#load(paste(path.res, "SMPDG_2017_c_v02_run01.Rdata", sep=""))
>>>>>>> grEBEN
