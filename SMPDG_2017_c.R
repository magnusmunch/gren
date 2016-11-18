##############################  preamble  #############################
# code belonging to abstract_SMPGD_2017_V01.pdf cloud version         #
# version: 02 (cloud)                                                 #
# author: Magnus Münch                                                #
# created: 15-11-2016                                                 #
# last edited: 18-11-2016                                             #
#######################################################################

###############################  notes  ###############################
# 15-11-2016: Building up from basic                                  #
#######################################################################

### paths
path.results <- "/home/magnusmunch/SMPDG_2017/"

### libraries
library(penalized)
library(GRridge)
library(pROC)

### functions
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
envb2 <- function(x, y, groups, lambda1, lambda2, maxiter=1000, epsilon=1e-07, trace=TRUE) {
  
  p <- ncol(x)
  n <- nrow(x)
  G <- length(unique(groups))
  sizes <- rle(groups)$lengths
  modmat <- matrix(0, ncol=G, nrow=p)
  modmat <- sapply(1:G, function(g) {as.numeric(groups==g)})
  m <- 1
  kappa <- y - 0.5*m
  
  # starting values
  fit <- penalized(y, x, unpenalized=~0, lambda1=0, lambda2=lambda1 + lambda2, model="logistic")
  
  mu <- muold <- fit@penalized
  
  xmu <- x %*% mu
  phat <- as.numeric(exp(xmu)/(1 + exp(xmu)))
  W <- diag(sqrt(phat*(1 - phat)))
  Xw <- W %*% x
  svdxw <- svd(Xw)
  U <- svdxw$u
  V <- svdxw$v
  d <- svdxw$d
  td <- min(n, p)
  invmat <- 1/(d^2 + 2*(lambda1 + lambda2))
  part1 <- invmat^2*d^2
  sigma <- sigmaold <- t(t(V)*part1) %*% t(V)
  
  ci <- ciold <- as.numeric(sqrt(colSums(t(x) * (sigma %*% t(x))) + (colSums(t(x)*mu))^2))
  chi <- chiold <- as.numeric(lambda2*(diag(sigma) + mu^2))
  
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
          paste("lambda1=[", paste(round(lambda1, 2), collapse=", "), 
                sep=""), "], ", 
          paste("lambda2=[", paste(round(lambda2, 2), collapse=", "), 
                sep=""), "]", sep="")
    }
    
    # estimation of new model parameters
    h <- (1 + sqrt(phi/chiold))*lambda2vec
    om <- (0.5*m/ciold)*tanh(ciold/2)
    omsq <- sqrt(om)
    hinvsq <- 1/sqrt(h)
    
    # using the Woodbury identity 
    V <- x*omsq
    U <- t(V)
    hinv <- 1/h
    vhinv <- t(t(V)*hinv)
    sigma <- diag(hinv) - (hinv*U) %*% solve(diag(n) + vhinv %*% U) %*% vhinv
    
    mu <- as.numeric(sigma %*% (t(x) %*% kappa))
    ci <- as.numeric(sqrt(colSums(t(x) * (sigma %*% t(x))) + (colSums(t(x)*mu))^2))
    chi <- as.numeric(lambda2vec*(diag(sigma) + mu^2))
    
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
    e.beta <- mu
    v.beta <- diag(sigma)
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
  
  out <- list(niter=niter, conv=conv, sigma=sigma, mu=mu, c=ci, chi=chi, lambda1=lambda1)
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
set.seed(123)
n <- 100           # Nb of observations    
ntest <- 1000          ## Nb of test set observations     
p <- 200            # Nb of variables per group
G <- 10            # Nb of groups
meanBeta <- 0.01   # Beta variances per group are VarBeta*(1:G); use this for CorX=0.5
CorX <- 0.5       # correlation within variable block
Nblock <- 10*G    # number of correlation blocks
settings <- c(ntrain=n, p=p, G=G, meanBeta=meanBeta, CorX=CorX, Nblock=Nblock)

nrep <- 100  #number of repeats per simulation setting
facvec <- c(1.3, 1.6, 2)   #tunes how much weaker each next group is. The '2' means that the second group is twice as weak as the first, etc
fractvec <- c(0, 0.7, 0.9) #tunes the sparsity per group. E.g. 0.9 means that 9/10 betas in a group are set to 0

out <- vector(mode="list", length=length(facvec)*length(fractvec))
names(out) <- as.vector(t(outer(paste("factor ", facvec, ",", sep=""), paste("fraction", fractvec), paste)))

aucmat <- matrix(NA, ncol=5, nrow=nrep)
briermat <- matrix(NA, ncol=5, nrow=nrep)
msemat <- matrix(NA, ncol=5, nrow=nrep)
for(fac in facvec){
  for(fract in fractvec){
    for(reptit in 1:nrep){
      print(paste("fac=", fac))
      print(paste("fract=", fract))
      print(paste("repeat=", reptit))
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
      vbSim <- envb2(x=X, y=Y, groups=grs, lambda1=1, lambda2=1, maxiter=1000, epsilon=1e-06, trace=TRUE)

      # GRridge
      groups <- CreatePartition(grs, grsize=p, uniform=T, decreasing=F)
      partsim <- list(grouping=groups)
      grSim <- grridge(t(X), Y, unpenal=~0, partsim, savepredobj="all", innfold=10, method="stable")
      
      # ridge
      rrSim <- optL2(Y, X, unpenalized=~0, lambda1=0, model="logistic", fold=10)
      
      # lasso
      lrSim <- optL1(Y, X, unpenalized=~0, lambda2=0, model="logistic", fold=10)
      
      # elastic net
      enSim <- optL2(Y, X, unpenalized=~0, lambda1=vbSim$lambda1, model="logistic", fold=10)
      
      # calculating mse
      grMse <- var((grSim$betas - Beta)^2)
      vbMse <- var((vbSim$mu - Beta)^2)
      rrMse <- var((rrSim$fullfit@penalized - Beta)^2)
      lrMse <- var((lrSim$fullfit@penalized - Beta)^2)
      enMse <- var((enSim$fullfit@penalized - Beta)^2)
      msemat[reptit, ] <- c(vbMse, grMse, rrMse, lrMse, enMse)
      
      ### TESTING THE MODELS
      # making the test data
      Xtest <- Reduce(cbind, lapply(1:Nblock, function(z) {
        matrix(rep(rnorm(ntest, sd=sqrt(CorX/(1 - CorX))), times=pblock), ntest, pblock)})) + 
        matrix(rnorm(ntest*G*p), ntest, G*p)
      lpredtest <- Xtest %*% Beta 
      logisticintercept <- 0
      
      probtest <- 1/(1 + exp(-(lpredtest + logisticintercept)))
      Ytest <- rbinom(length(probtest), 1, probtest)
      
      # making predictions
      vbPred <- 1/(1 + exp(-(Xtest %*% vbSim$mu)))
      grPred <- predict.grridge(grSim, t(Xtest))
      rrPred <- predict(rrSim$fullfit, Xtest)
      lrPred <- predict(lrSim$fullfit, Xtest)
      enPred <- predict(enSim$fullfit, Xtest)
      
      # calculating brier residuals and auc
      grBrier <- mean((grPred[, 2] - probtest)^2)
      vbBrier <- mean((vbPred - probtest)^2)
      rrBrier <- mean((rrPred - probtest)^2)
      lrBrier <- mean((lrPred - probtest)^2)
      enBrier <- mean((enPred - probtest)^2)
      briermat[reptit, ] <- c(vbBrier, grBrier, rrBrier, lrBrier, enBrier)
      
      cutoffs <- rev(seq(0, 1, by=0.005))
      grRoc <- GRridge::roc(probs=as.numeric(grPred[, 2]), true=Ytest, cutoffs) #ridge, sel
      vbRoc <- GRridge::roc(probs=as.numeric(vbPred), true=Ytest, cutoffs)
      rrRoc <- GRridge::roc(probs=as.numeric(rrPred), true=Ytest, cutoffs)
      lrRoc <- GRridge::roc(probs=as.numeric(lrPred), true=Ytest, cutoffs)
      enRoc <- GRridge::roc(probs=as.numeric(enPred), true=Ytest, cutoffs)
      aucmat[reptit, ] <- c(GRridge::auc(vbRoc), GRridge::auc(grRoc), GRridge::auc(rrRoc), GRridge::auc(lrRoc),
                            GRridge::auc(enRoc))
      
    }
    
    colnames(msemat) <- colnames(briermat) <- colnames(aucmat) <- c("ENVB", "GRridge", "RR", "LR", "EN")
    out[[3*fac + fract - 3]] <- list(msemat=msemat, briermat=briermat, aucmat=aucmat)
    
  }   
}

save(out, file=paste(path.results, "SMPDG_2017_c_v01.Rdata", sep=""))

