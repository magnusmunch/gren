##############################  preamble  #############################
# code belonging to abstract_SMPGD_2017_V01.pdf                       #
# version: 02                                                         #
# author: Magnus Münch                                                #
# created: 15-11-2016                                                 #
# last edited: 18-11-2016                                             #
#######################################################################

###############################  notes  ###############################
# 15-11-2016: Building up from basic                                  #
#######################################################################

### paths
path.results <- "C:/Users/Magnus/Documents/phd/ENVB/abstract_SMPGD_2017/results/"
path.code <- "C:/Users/Magnus/Documents/phd/ENVB/code/"
path.graph <- "C:/Users/Magnus/Documents/phd/ENVB/graphs/"

### libraries
library(penalized)
library(mvtnorm)
library(GRridge)
library(pROC)
library(Rcpp)

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
    fit <- penalized(y, x, unpenalized=~1, lambda1=0, lambda2=lambda1 + lambda2, model="logistic")
    xaug <- cbind(1, x)
  } else {
    fit <- penalized(y, x, unpenalized=~0, lambda1=0, lambda2=lambda1 + lambda2, model="logistic")
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

nrep <- 3  #number of repeats per simulation setting
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
      set.seed(as.numeric(gsub("\\.", "", paste(fac, fract, reptit, sep=""))))
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
      vbSim <- envb2(x=X, y=Y, groups=grs, lambda1=1, lambda2=500, maxiter=1000, epsilon=1e-06, trace=TRUE, 
                     intercept=TRUE)
      
      # GRridge
      groups <- CreatePartition(grs, grsize=p, uniform=T, decreasing=F)
      partsim <- list(grouping=groups)
      grSim <- tryCatch({
        grridge(t(X), Y, unpenal=~1, partsim, savepredobj="all", innfold=10, method="stable")},
        error=function(war) {return(NULL)})
      
      # ridge
      rrSim <- optL2(Y, X, unpenalized=~1, lambda1=0, model="logistic", fold=10)
      
      # lasso
      lrSim <- optL1(Y, X, unpenalized=~1, lambda2=0, model="logistic", fold=10)
      
      # elastic net
      enSim <- optL2(Y, X, unpenalized=~1, lambda1=vbSim$lambda1, model="logistic", fold=10)
      
      # calculating mse
      if(is.null(grSim)) {
        grMse <- NA
      } else {
        grMse <- var((grSim$betas - Beta)^2)
      }
      vbMse <- var((vbSim$mu[-1] - Beta)^2)
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
      vbPred <- 1/(1 + exp(-(cbind(1, Xtest) %*% vbSim$mu)))
      rrPred <- predict(rrSim$fullfit, Xtest)
      lrPred <- predict(lrSim$fullfit, Xtest)
      enPred <- predict(enSim$fullfit, Xtest)
      
      # GRridge stuff
      cutoffs <- rev(seq(0, 1, by=0.005))
      if(is.null(grSim)) {
        grBrier <- NA
        grAuc <- NA
      } else {
        grPred <- predict.grridge(grSim, t(Xtest))
        grBrier <- mean((grPred[, 2] - probtest)^2)
        grRoc <- GRridge::roc(probs=as.numeric(grPred[, 2]), true=Ytest, cutoffs)
        grAuc <- GRridge::auc(grRoc)
      }
      
      # calculating brier residuals and auc
      
      vbBrier <- mean((vbPred - probtest)^2)
      rrBrier <- mean((rrPred - probtest)^2)
      lrBrier <- mean((lrPred - probtest)^2)
      enBrier <- mean((enPred - probtest)^2)
      briermat[reptit, ] <- c(vbBrier, grBrier, rrBrier, lrBrier, enBrier)
      
      vbRoc <- GRridge::roc(probs=as.numeric(vbPred), true=Ytest, cutoffs)
      rrRoc <- GRridge::roc(probs=as.numeric(rrPred), true=Ytest, cutoffs)
      lrRoc <- GRridge::roc(probs=as.numeric(lrPred), true=Ytest, cutoffs)
      enRoc <- GRridge::roc(probs=as.numeric(enPred), true=Ytest, cutoffs)
      aucmat[reptit, ] <- c(GRridge::auc(vbRoc), grAuc, GRridge::auc(rrRoc), GRridge::auc(lrRoc),
                            GRridge::auc(enRoc))
      
    }
    
    colnames(msemat) <- colnames(briermat) <- colnames(aucmat) <- c("ENVB", "GRridge", "RR", "LR", "EN")
    out[[3*which(fac==facvec) + which(fract==fractvec) - 3]] <- list(msemat=msemat, briermat=briermat, aucmat=aucmat)
    
  }   
}

save(out, file=paste(path.results, "SMPDG_2017_v01_run02.Rdata", sep=""))

### reading in results from surfsara
load(paste(path.results, "SMPDG_2017_c_v01.Rdata", sep=""))
boxplot(out[[1]]$briermat[, -4])
boxplot(out[[2]]$msemat[, -4])
boxplot(out[[3]]$msemat[, -4])
boxplot(out[[4]]$msemat[, -4])
boxplot(out[[5]]$msemat[, -4])
boxplot(out[[5]]$msemat[, 1])

### comparing to Gibbs sampling 
library(metabolomics)
data(treated)


str(treated)

### data example
data(dataVerlaat)
cpganngroup <- CreatePartition(CpGann)
grs <- rep(1:6, times=sapply(1:length(cpganngroup), function(g) {length(cpganngroup[[g]])}))
grs2 <- CreatePartition(as.factor(grs))
X <- t(datcenVerlaat)[, order(unlist(cpganngroup))]
X <- apply(X, 2, function(j) {(j - mean(j)/sd(j))})
Y <- respVerlaat

test1 <- envb2(x=X, y=Y, groups=grs, lambda1=1, lambda2=100, intercept=TRUE, maxiter=500, epsilon=1e-06, trace=TRUE)
test2 <- grridge(t(X), Y, partitions=list(cpg=grs2), unpenal=~1, innfold=10)
test3 <- grridge(datcenVerlaat, respVerlaat, partition=list(gpg=CreatePartition(CpGann)),
                 unpenal=~1, innfold=10)


### creating graphs with performance in simulations
library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)
library(grid)

# results from surfsara setting 1
load(paste(path.results, "SMPDG_2017_c_v02_run02.Rdata", sep=""))
data1.1 <- melt(as.data.frame(out[[1]]))
data1.2 <- melt(as.data.frame(out[[2]]))
data1.3 <- melt(as.data.frame(out[[3]]))

th1 <- theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
            panel.background=element_blank(), axis.line=element_line(colour="black"),
            plot.title=element_text(hjust=0.5))
th2 <- theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             panel.background=element_blank(), axis.line=element_line(colour="black"),
             axis.title.y=element_blank(), axis.title.x=element_blank())

plot.big1.1 <- ggplot() + geom_boxplot(aes(x=variable, y=value), data1.1) + th1
plot.zoom1.1 <- ggplot() + geom_boxplot(aes(x=variable, y=value),
                                        data1.1[!(data1.1$variable %in% c("GRridge", "LR")), ]) + th2
ylim1.1 <- c(min(boxplot(out[[1]][, -c(2, 4)], plot=FALSE)$stats),
             max(boxplot(out[[1]][, -c(2, 4)], plot=FALSE)$stats))
ylim1.2 <- ggplot_build(plot.big1.1)$layout$panel_ranges[[1]]$y.range[2] - 
  c(diff(ggplot_build(plot.big1.1)$layout$panel_ranges[[1]]$y.range)/1.5,
    diff(ggplot_build(plot.big1.1)$layout$panel_ranges[[1]]$y.range)*0.08)
plot1.1 <- ggplot() + geom_boxplot(aes(x=variable, y=value), data1.1) + th1 + 
  annotation_custom(grob=ggplotGrob(plot.zoom1.1 + coord_cartesian(ylim=ylim1.1)), xmax=3, 
                    ymin=ylim1.2[1], ymax=ylim1.2[2]) + labs(title="a)", x="", y="Average MSE")
plot1.2 <- ggplot() + geom_boxplot(aes(x=variable, y=value), data1.2) + th1 + 
  labs(title="b)", x="", y="Mean Brier residuals")
plot1.3 <- ggplot() + geom_boxplot(aes(x=variable, y=value), data1.3) + th1 + 
  labs(title="c)", x="f=1.3, q=0", y="AUC")

# reading in results from surfsara setting 2
load(paste(path.results, "SMPDG_2017_c_v02_run03.Rdata", sep=""))
data2.1 <- melt(as.data.frame(out[[1]]))
data2.2 <- melt(as.data.frame(out[[2]]))
data2.3 <- melt(as.data.frame(out[[3]]))

plot.big2.1 <- ggplot() + geom_boxplot(aes(x=variable, y=value), data2.1) + th1
plot.zoom2.1 <- ggplot() + geom_boxplot(aes(x=variable, y=value),
                                        data2.1[!(data2.1$variable %in% c("GRridge", "LR")), ]) + th2
ylim2.1 <- c(min(boxplot(out[[1]][, -c(2, 4)], plot=FALSE)$stats),
             max(boxplot(out[[1]][, -c(2, 4)], plot=FALSE)$stats))
ylim2.2 <- ggplot_build(plot.big2.1)$layout$panel_ranges[[1]]$y.range[2] - 
  c(diff(ggplot_build(plot.big2.1)$layout$panel_ranges[[1]]$y.range)/1.5,
    diff(ggplot_build(plot.big2.1)$layout$panel_ranges[[1]]$y.range)*0.08)
plot2.1 <- ggplot() + geom_boxplot(aes(x=variable, y=value), data2.1) + th1 + 
  annotation_custom(grob=ggplotGrob(plot.zoom2.1 + coord_cartesian(ylim=ylim2.1)), xmax=3, ymin=ylim2.2[1],
                    ymax=ylim2.2[2]) + labs(title="d)", x="", y="")
plot2.2 <- ggplot() + geom_boxplot(aes(x=variable, y=value), data2.2) + th1 + 
  labs(title="e)", x="", y="")
plot2.3 <- ggplot() + geom_boxplot(aes(x=variable, y=value), data2.3) + th1 + 
  labs(title="f)", x="f=1.6, q=0.7", y="")

# reading in results from surfsara setting 3
load(paste(path.results, "SMPDG_2017_c_v02_run04.Rdata", sep=""))
data3.1 <- melt(as.data.frame(out[[1]]))
data3.2 <- melt(as.data.frame(out[[2]]))
data3.3 <- melt(as.data.frame(out[[3]]))

plot.big3.1 <- ggplot() + geom_boxplot(aes(x=variable, y=value), data3.1) + th1
plot.zoom3.1 <- ggplot() + geom_boxplot(aes(x=variable, y=value),
                                        data3.1[!(data3.1$variable %in% c("GRridge", "LR")), ]) + th2
ylim3.1 <- c(min(boxplot(out[[1]][, -c(2, 4)], plot=FALSE)$stats),
             max(boxplot(out[[1]][, -c(2, 4)], plot=FALSE)$stats))
ylim3.2 <- ggplot_build(plot.big3.1)$layout$panel_ranges[[1]]$y.range[2] - 
  c(diff(ggplot_build(plot.big3.1)$layout$panel_ranges[[1]]$y.range)/1.5,
    diff(ggplot_build(plot.big3.1)$layout$panel_ranges[[1]]$y.range)*0.08)
plot3.1 <- ggplot() + geom_boxplot(aes(x=variable, y=value), data3.1) + th1 + 
  annotation_custom(grob=ggplotGrob(plot.zoom3.1 + coord_cartesian(ylim=ylim3.1)), xmax=3, ymin=ylim3.2[1],
                    ymax=ylim3.2[2]) + labs(title="g)", x="", y="")
plot3.2 <- ggplot() + geom_boxplot(aes(x=variable, y=value), data3.2) + th1 + 
  labs(title="h)", x="", y="")
plot3.3 <- ggplot() + geom_boxplot(aes(x=variable, y=value), data3.3) + th1 + 
  labs(title="i)", x="f=2, q=0.9", y="")

# reading in results from surfsara setting 4
load(paste(path.results, "SMPDG_2017_c_v02_run05.Rdata", sep=""))
data4.1 <- melt(as.data.frame(out[[1]]))
data4.2 <- melt(as.data.frame(out[[2]]))
data4.3 <- melt(as.data.frame(out[[3]]))

plot.big4.1 <- ggplot() + geom_boxplot(aes(x=variable, y=value), data4.1) + th1
plot.zoom4.1 <- ggplot() + geom_boxplot(aes(x=variable, y=value),
                                        data4.1[!(data4.1$variable %in% c("GRridge", "LR")), ]) + th2
ylim4.1 <- c(min(boxplot(out[[1]][, -c(2, 4)], plot=FALSE)$stats),
             max(boxplot(out[[1]][, -c(2, 4)], plot=FALSE)$stats))
ylim4.2 <- ggplot_build(plot.big4.1)$layout$panel_ranges[[1]]$y.range[2] - 
  c(diff(ggplot_build(plot.big4.1)$layout$panel_ranges[[1]]$y.range)/1.5,
    diff(ggplot_build(plot.big4.1)$layout$panel_ranges[[1]]$y.range)*0.08)
plot4.1 <- ggplot() + geom_boxplot(aes(x=variable, y=value), data4.1) + th1 + 
  annotation_custom(grob=ggplotGrob(plot.zoom4.1 + coord_cartesian(ylim=ylim4.1)), xmax=3, ymin=ylim4.2[1],
                    ymax=ylim4.2[2]) + labs(title="a)", x="", y="")
plot4.2 <- ggplot() + geom_boxplot(aes(x=variable, y=value), data4.2) + th1 + 
  labs(title="b)", x="", y="")
plot4.3 <- ggplot() + geom_boxplot(aes(x=variable, y=value), data4.3) + th1 + 
  labs(title="c)", x="f=1.3, q=0.7", y="")

# reading in results from surfsara setting 5
load(paste(path.results, "SMPDG_2017_c_v02_run06.Rdata", sep=""))
data5.1 <- melt(as.data.frame(out[[1]]))
data5.2 <- melt(as.data.frame(out[[2]]))
data5.3 <- melt(as.data.frame(out[[3]]))

plot.big5.1 <- ggplot() + geom_boxplot(aes(x=variable, y=value), data5.1) + th1
plot.zoom5.1 <- ggplot() + geom_boxplot(aes(x=variable, y=value),
                                        data5.1[!(data5.1$variable %in% c("GRridge", "LR", "EN")), ]) + th2
ylim5.1 <- c(min(boxplot(out[[1]][, -c(2, 4, 5)], plot=FALSE)$stats),
             max(boxplot(out[[1]][, -c(2, 4, 5)], plot=FALSE)$stats))
ylim5.2 <- ggplot_build(plot.big5.1)$layout$panel_ranges[[1]]$y.range[2] - 
  c(diff(ggplot_build(plot.big5.1)$layout$panel_ranges[[1]]$y.range)/1.5,
    diff(ggplot_build(plot.big5.1)$layout$panel_ranges[[1]]$y.range)*0.08)
plot5.1 <- ggplot() + geom_boxplot(aes(x=variable, y=value), data5.1) + th1 + 
  annotation_custom(grob=ggplotGrob(plot.zoom5.1 + coord_cartesian(ylim=ylim5.1)), xmax=3, ymin=ylim5.2[1],
                    ymax=ylim5.2[2]) + labs(title="d)", x="", y="")
plot5.2 <- ggplot() + geom_boxplot(aes(x=variable, y=value), data5.2) + th1 + 
  labs(title="e)", x="", y="")
plot5.3 <- ggplot() + geom_boxplot(aes(x=variable, y=value), data5.3) + th1 + 
  labs(title="f)", x="f=1.6, q=0", y="")

boxplot(out[[1]][, -c(2, 4, 5)], outline=FALSE)

# visualize
windows()
grid.arrange(plot1.1, plot2.1, plot3.1,
             plot1.2, plot2.2, plot3.2,
             plot1.3, plot2.3, plot3.3, ncol=3, nrow=3)

windows()
grid.arrange(plot4.1, plot5.1,
             plot4.2, plot5.2,
             plot4.3, plot5.3, ncol=2, nrow=3)

# save to png
g1 <- arrangeGrob(plot1.1, plot2.1, plot3.1,
                  plot1.2, plot2.2, plot3.2,
                  plot1.3, plot2.3, plot3.3, ncol=3, nrow=3)
g2 <- arrangeGrob(plot4.1, plot5.1,
                  plot4.2, plot5.2,
                  plot4.3, plot5.3, ncol=2, nrow=3)
ggsave(paste(path.graph, "ENVB_sim01_V01.png", sep=""), g1, device="png", width=297, height=210,
       units="mm", dpi=300)
ggsave(paste(path.graph, "ENVB_sim02_V01.png", sep=""), g2, device="png", width=297, height=210,
       units="mm", dpi=300)



















################################## TESTING SPARSE GROUP LASSO ##################################
library(SGL)

set.seed(123)
n <- 100           # Nb of observations    
ntest <- 1000          ## Nb of test set observations     
p <- 200            # Nb of variables per group
G <- 10            # Nb of groups
meanBeta <- 0.01   # Beta variances per group are VarBeta*(1:G); use this for CorX=0.5
CorX <- 0.5       # correlation within variable block
Nblock <- 10*G    # number of correlation blocks

nrep <- 3  #number of repeats per simulation setting
fac <- 1.3   #tunes how much weaker each next group is. The '2' means that the second group is twice as weak as the first, etc
fract <- 0 #tunes the sparsity per group. E.g. 0.9 means that 9/10 betas in a group are set to 0

reps <- rev(sapply(0:(G - 1), function(i) {fac^(-i)}))
meansB <- rep(reps, each=p)*meanBeta/mean(rep(reps, each=p))
Beta <- meansB
Beta <- rev(sparsify(Beta, frac=fract))
      
### FITTING THE MODELS
pblock <- G*p/Nblock
grs <- rep(1:G, each=p)
P <- G*p #Complete number of variables
X <- Reduce(cbind, lapply(1:Nblock, function(z) {
  matrix(rep(rnorm(n, sd=sqrt(CorX/(1 - CorX))), times=pblock), n, pblock)})) + 
  matrix(rnorm(n*G*p), n, G*p)
X <- t((t(X) - apply(t(X), 1, mean))/apply(t(X), 1, sd))
      
lpred <- X %*% Beta 
logisticintercept <- 0
prob <- 1/(1 + exp(-(lpred + logisticintercept)))
Y <- rbinom(length(prob), 1, prob)
      
# ENVB
vbSim <- envb2(x=X, y=Y, groups=grs, lambda1=1, lambda2=500, maxiter=1000, epsilon=1e-06, trace=TRUE, 
               intercept=TRUE)
vb2Sim <- penalized(Y, X, unpenalized=~1, lambda1=vbSim$lambda1, lambda2=rep(vbSim$lambda2, each=p),
                    model="logistic")
      
# GRridge
groups <- CreatePartition(grs, grsize=p, uniform=T, decreasing=F)
partsim <- list(grouping=groups)
grSim <- tryCatch({
  grridge(t(X), Y, unpenal=~1, partsim, savepredobj="all", innfold=10, method="stable")},
  error=function(war) {return(NULL)})
      
# ridge
rrSim <- optL2(Y, X, unpenalized=~1, lambda1=0, model="logistic", fold=10)
      
# lasso
lrSim <- optL1(Y, X, unpenalized=~1, lambda2=0, model="logistic", fold=10)
      
# elastic net
enSim <- optL2(Y, X, unpenalized=~1, lambda1=vbSim$lambda1, model="logistic", fold=10)

# sparse group lasso
sglSim <- cvSGL(list(x=X, y=Y), grs, type="logit", standardize=FALSE, verbose=TRUE)

# calculating mse
grMse <- var((grSim$betas - Beta)^2)
vbMse <- var((vbSim$mu[-1] - Beta)^2)
vb2Mse <- var((vb2Sim@penalized - Beta)^2)
rrMse <- var((rrSim$fullfit@penalized - Beta)^2)
lrMse <- var((lrSim$fullfit@penalized - Beta)^2)
enMse <- var((enSim$fullfit@penalized - Beta)^2)
sglMse <- var((sglSim$fit$beta[, which.min(sglSim$lldiff)] - Beta)^2)
msemat <- c(vbMse, vb2Mse, grMse, rrMse, lrMse, enMse, sglMse)

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
grPred <- predict.grridge(grSim, t(Xtest))
vbPred <- 1/(1 + exp(-(cbind(1, Xtest) %*% vbSim$mu)))
vb2Pred <- predict(vb2Sim, Xtest)
rrPred <- predict(rrSim$fullfit, Xtest)
lrPred <- predict(lrSim$fullfit, Xtest)
enPred <- predict(enSim$fullfit, Xtest)
sglPred <- 1/(1 + exp(-(cbind(1, Xtest) %*% c(sglSim$fit$intercepts[which.min(sglSim$lldiff)], 
                                              sglSim$fit$beta[, which.min(sglSim$lldiff)]))))

# brier residuals
grBrier <- mean((grPred[, 2] - probtest)^2)
vbBrier <- mean((vbPred - probtest)^2)
vb2Brier <- mean((vb2Pred - probtest)^2)
rrBrier <- mean((rrPred - probtest)^2)
lrBrier <- mean((lrPred - probtest)^2)
enBrier <- mean((enPred - probtest)^2)
sglBrier <- mean((sglPred - probtest)^2)
briermat <- c(vbBrier, vb2Brier, grBrier, rrBrier, lrBrier, enBrier, sglBrier)

# AUC
cutoffs <- rev(seq(0, 1, by=0.005))
grRoc <- GRridge::roc(probs=as.numeric(grPred), true=Ytest, cutoffs)
vbRoc <- GRridge::roc(probs=as.numeric(vbPred), true=Ytest, cutoffs)
vb2Roc <- GRridge::roc(probs=as.numeric(vb2Pred), true=Ytest, cutoffs)
rrRoc <- GRridge::roc(probs=as.numeric(rrPred), true=Ytest, cutoffs)
lrRoc <- GRridge::roc(probs=as.numeric(lrPred), true=Ytest, cutoffs)
enRoc <- GRridge::roc(probs=as.numeric(enPred), true=Ytest, cutoffs)
sglRoc <- GRridge::roc(probs=as.numeric(sglPred), true=Ytest, cutoffs)
aucmat <- c(GRridge::auc(vbRoc), GRridge::auc(vb2Roc), GRridge::auc(grRoc), GRridge::auc(rrRoc), 
            GRridge::auc(lrRoc), GRridge::auc(enRoc), GRridge::auc(sglRoc))


