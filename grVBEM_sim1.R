##############################  preamble  #############################
# simulations for grVBEM                                              #
# version: 01                                                         #
# author: Magnus Münch                                                #
# created: 15-03-2017                                                 #
# last edited: 15-03-2017                                             #
#######################################################################

###############################  notes  ###############################
#                                                                     #
#######################################################################

### paths
path.rcode <- "C:/Users/Magnus/Documents/phd/ENVB/code/"
path.res <- "C:/Users/Magnus/Documents/phd/ENVB/results/"

### libraries
library(GRridge)
library(mvtnorm)
library(grpreg)

### source grENVB functions
source(paste(path.rcode, "grVBEM.R", sep=""))

### the simulation
settings <- cbind(expand.grid(f=c(1.3, 1.6, 2), q=c(0, 0.7, 0.9)), G=10, 
                  pg=200, meanbeta=0.01, nblock=20, rho=0.3, sigma2=1, m=1)
n <- 100
ntest <- 1000
reps <- 5
briermat <- aucmat <- matrix(NA, ncol=5, nrow=reps)
for(s in 1:nrow(settings)) {
  
  # specifying parameters for current setting
  p <- settings$G[s]*settings$pg[s]
  groups1 <- rep(1:settings$G[s], each=settings$pg[s])
  groups2 <- CreatePartition(as.factor(groups1))
  pblock <- p/settings$nblock[s]
  ubeta <- rep(rev(sapply(0:(settings$G[s] - 1), function(g) {settings$f[s]^(-g)})), 
               each=settings$pg[s])
  beta <- ubeta*settings$meanbeta[s]/mean(ubeta)
  bsigma <- matrix(settings$rho[s], ncol=pblock, nrow=pblock)
  diag(bsigma) <- settings$sigma2[s]
  
  for(r in 1:reps) {
    
    # creating the data
    x <- do.call(cbind, replicate(settings$nblock[s], 
                                  rmvnorm(n, mean=rep(0, pblock), 
                                          sigma=bsigma), simplify=FALSE))
    prob <- as.numeric(exp(x %*% beta)/(1 + exp(x %*% beta)))
    y <- rbinom(n, settings$m[s], prob)
    
    # fitting the models
    fit.grridge <- grridge(t(x), y, groups2, unpenal=~1)
    fit.grVBEM <- grVBEM(x, y, m=rep(settings$m[s], n), groups1, lambda1=NULL, 
                         lambda2=NULL, intercept=TRUE, eps=0.001, maxiter=500)
    fit.grpreg <- cv.grpreg(x, y, group=groups1, penalty="grLasso", family="binomial")
    
    # creating test data
    xtest <- do.call(cbind, replicate(settings$nblock[s], 
                                      rmvnorm(ntest, mean=rep(0, pblock), 
                                              sigma=bsigma), simplify=FALSE))
    probtest <- as.numeric(exp(xtest %*% beta)/(1 + exp(xtest %*% beta)))
    ytest <- rbinom(ntest, rep(settings$m[s], ntest), probtest)
    
    # computing predictions
    pred.truth <- as.numeric(exp(xtest %*% beta)/(1 + exp(xtest %*% beta)))
    pred.grridge <- predict.grridge(fit.grridge, t(xtest))[, 2]
    pred.grVBEM <- as.numeric(exp(cbind(1, xtest) %*% fit.grVBEM$mu)/
                                (1 + exp(cbind(1, xtest) %*% fit.grVBEM$mu)))
    pred.grpreg <- predict(fit.grpreg, xtest, type="response")
    pred.ridge <- predict.grridge(fit.grridge, t(xtest))[, 1]
    
    # calculating AUCs on the predictions
    auc.truth <- pROC::roc(ytest, pred.truth)$auc
    auc.grridge <- pROC::roc(ytest, pred.grridge)$auc
    auc.grVBEM <- pROC::roc(ytest, pred.grVBEM)$auc
    auc.grpreg <- pROC::roc(ytest, pred.grpreg)$auc
    auc.ridge <- pROC::roc(ytest, pred.ridge)$auc
    
    # calculating Brier residuals on the predictions
    brier.truth <- mean((pred.truth - probtest)^2)
    brier.grridge <- mean((pred.grridge - probtest)^2)
    brier.grVBEM <- mean((pred.grVBEM - probtest)^2)
    brier.grpreg <- mean((pred.grpreg - probtest)^2)
    brier.ridge <- mean((pred.ridge - probtest)^2)
    
    # filling in the matrices with results
    aucmat[r, ] <- c(auc.truth, auc.grridge, auc.grVBEM, auc.grpreg, auc.ridge)
    briermat[r, ] <- c(brier.truth, brier.grridge, brier.grVBEM, brier.grpreg, brier.ridge)
  
  }
  
  # naming the columns and saving the objects
  colnames(aucmat) <- colnames(briermat) <- c("truth", "grridge", "grVBEM", "grpred", "ridge")
  assign(paste("setting", s, sep=""), list(auc=aucmat, brier=briermat))
  save(paste("setting", s, sep=""), file=paste(path.res, "grVBEM_res1_setting", s, ".Rdata", sep=""))
  
}

