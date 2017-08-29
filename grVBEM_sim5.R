##############################  preamble  #############################
# simulations for grVBEM to assess variable selection                 #
# version: 01                                                         #
# author: Magnus M?nch                                                #
# created: 11-07-2017                                                 #
# last edited: 11-07-2017                                             #
#######################################################################

###############################  notes  ###############################
#                                                                     #
#######################################################################

### paths
path.code <- as.character(ifelse(Sys.info()[1]=="Darwin","/Users/magnusmunch/Documents/PhD/EBEN/code/" ,"~/EBEN/code/"))
path.res <- as.character(ifelse(Sys.info()[1]=="Darwin","/Users/magnusmunch/Documents/PhD/EBEN/results/" ,"~/EBEN/results/"))
path.data <- as.character(ifelse(Sys.info()[1]=="Darwin","/Users/magnusmunch/Documents/PhD/EBEN/data/" ,"~/EBEN/data/"))
path.graph <- "/Users/magnusmunch/Documents/PhD/EBEN/graphs/"

### libraries
library(Rcpp)
library(glmnet)
library(penalized)
library(GRridge)
library(mvtnorm)
library(pROC)
library(psych)

### functions
# source grENVB functions
source(paste(path.code, "grVBEM.R", sep=""))

# function to cross-validate elastic net
cv.en <- function(x, y, intercept=intercept) {
  fit.pen <- cv.pen(x, y, intercept=intercept)
  fit.en <- glmnet(x, y, family="binomial", alpha=fit.pen$alpha[which.min(fit.pen$cvll)],
                   lambda=fit.pen$lambda[which.min(fit.pen$cvll)])
  return(fit.en)
}

# function to cross-validate lambda1 and lambda2 in enet with fixed number of selected vars
# cv.pensel <- function(x, y, intercept, nsel, type=c("ridge", "lasso", "enet")) {
#   p <- ncol(x)
#   n <- nrow(x)
#   if(type=="enet") {
#     seq.alpha <- seq(0.01, 0.99, length.out=50)
#   } else if(type=="lasso"){
#     seq.alpha <- 1
#   } else if(type=="ridge") {
#     seq.alpha <- 0
#     df.fit <- list(df=NA)
#   }
#   nlambda <- 100
#   seq.lam <- seq.df <- seq.cvll <- numeric(length(seq.alpha))
#   for(a in 1:length(seq.alpha)) {
#     # first find the correct lambda for every alpha
#     suppressWarnings(df.fit <- glmnet(x, y, family="binomial", alpha=seq.alpha[a], 
#                                       standardize=FALSE, dfmax=nsel, intercept=intercept))
#     if(type=="ridge") {
#       lambda <- NULL
#     } else {
#       lambda <- tail(df.fit$lambda[df.fit$df <= nsel], n=2L)
#     }
#     found <- any(df.fit$df==nsel)
#     lambdamax <- max(abs(t(x) %*% (y - mean(y)*(1 - mean(y)))))/(n*seq.alpha[a])
#     while(!found & type!="ridge") {
#       found <- any(df.fit$df==nsel)
#       if(!found) {
#         lambdamax <- ifelse(length(tail(df.fit$lambda[df.fit$df < nsel], n=1L))==0,
#                             lambdamax, 
#                             ifelse(!all(df.fit$df < nsel),
#                                    ifelse(df.fit$df[1] > nsel,
#                                           tail(df.fit$lambda[which(!(df.fit$df < nsel))], n=1L),
#                                           df.fit$lambda[which(!(df.fit$df < nsel))[1] - 1]),
#                                    tail(df.fit$lambda, n=1L)))
#                             # tail(df.fit$lambda[df.fit$df < nsel], n=1L))
#         lambdamin <- ifelse(length(df.fit$lambda[df.fit$df > nsel])==0, 
#                             lambdamax*ifelse(n < p, 0.01, 0.0001),
#                             ifelse(df.fit$df[1] > nsel,
#                                    df.fit$lambda[df.fit$df < nsel][1],
#                                    df.fit$lambda[df.fit$df > nsel][1]))
#         lambda <- exp(seq(log(lambdamax), log(lambdamin), length.out=nlambda))
#         if(length(unique(lambda))==1) {
#           found <- TRUE
#           lambda <- c(max(abs(t(x) %*% (y - mean(y)*(1 - mean(y)))))/(n*seq.alpha[a]),
#                       unique(lambda))
#         }
#         suppressWarnings(df.fit <- glmnet(x, y, family="binomial", alpha=seq.alpha[a], 
#                                           standardize=FALSE, lambda=lambda, 
#                                           intercept=intercept))
#       } else {
#         lambda <- tail(df.fit$lambda[1:tail(which(df.fit$df==nsel), n=1L)], n=2L)
#         if(length(unique(lambda))==1) {
#           lambda <- c(max(abs(t(x) %*% (y - mean(y)*(1 - mean(y)))))/(n*seq.alpha[a]), 
#                       unique(lambda))
#         }
#       }
#     }
#     
#     # cross-validate the lambda and alpha combination
#     cv.fit <- cv.glmnet(x, y, family="binomial", lambda=lambda, alpha=seq.alpha[a], 
#                         standardize=FALSE, intercept=intercept)
#     seq.cvll[a] <- ifelse(type=="ridge", min(cv.fit$cvm), tail(cv.fit$cvm, n=1L))
#     seq.lam[a] <- ifelse(type=="ridge", cv.fit$lambda.min, tail(cv.fit$lambda, n=1L))
#     seq.df[a] <- ifelse(type=="ridge", p, ifelse(!any(df.fit$df==nsel),
#                                                  max(df.fit$df[which.min(abs(df.fit$df - nsel))]),
#                                                  df.fit$df[tail(which(df.fit$df==nsel), n=1L)]))
#   }
#   cvll <- min(seq.cvll)
#   df <- seq.df[which.min(seq.cvll)]
#   
#   lambdaglmnet <- seq.lam[which.min(seq.cvll)]
#   alphaglmnet <- seq.alpha[which.min(seq.cvll)]
#   lambda1pen <- lambdaglmnet*alphaglmnet*n
#   lambda2pen <- lambdaglmnet*(1 - alphaglmnet)*n*0.5
#   lambda1bayes <- lambdaglmnet*alphaglmnet*n*2
#   lambda2bayes <- lambdaglmnet*(1 - alphaglmnet)*n
#   
#   out <- list(lambdaglmnet=lambdaglmnet, alphaglmnet=alphaglmnet, lambda1pen=lambda1pen, 
#               lambda2pen=lambda2pen, lambda1bayes=lambda1bayes, lambda2bayes=lambda2bayes, 
#               cvll=cvll, df=df)
#   return(out)
# }

cv.pensel <- function(x, y, intercept, nsel, type=c("ridge", "lasso", "enet")) {
  p <- ncol(x)
  n <- nrow(x)
  lambda.min.ratio <- ifelse(n < p, 0.01, 0.0001)
  nlambda <- 100
  if(type=="enet") {
    seq.alpha <- seq(0.01, 0.99, length.out=50)
  } else if(type=="lasso"){
    seq.alpha <- 1
  } else if(type=="ridge") {
    seq.alpha <- 0
    lambdamax <- 1000*max(abs(t(x) %*% (y - mean(y)*(1 - mean(y)))))/n
    lambdamin <- lambdamax*lambda.min.ratio
    lambda <- exp(seq(log(lambdamax), log(lambdamin), length.out=nlambda))
  }
  seq.lam <- seq.df <- seq.cvll <- numeric(length(seq.alpha))
  for(a in 1:length(seq.alpha)) {
    df <- p
    lambdaseq <- NULL
    found <- ifelse(type=="ridge", TRUE, FALSE)
    inbet <- FALSE
    count <- 0
    while(!found) {
      suppressWarnings(fit.df <- glmnet(x, y, family="binomial", alpha=seq.alpha[a], 
                                        nlambda=nlambda, lambda=lambdaseq, standardize=FALSE, 
                                        intercept=intercept, 
                                        dfmax=ifelse(is.null(lambdaseq), nsel, p + 1)))
      if(any(fit.df$df==nsel)) {
        if(tail(fit.df$df, n=1L)==nsel & !inbet) {
          lambdamax <- tail(fit.df$lambda[fit.df$df==nsel], n=1L)
          lambdamin <- lambdamax*lambda.min.ratio
        } else {
          found <- TRUE
          lambda <- c(max(abs(t(x) %*% (y - mean(y)*(1 - mean(y)))))/(n*seq.alpha[a]), 
                     tail(fit.df$lambda[fit.df$df==nsel], n=1L))
        }
      } else if(all(fit.df$df < nsel)) {
        lambdamax <- tail(fit.df$lambda, n=1L)
        lambdamin <- lambdamax*lambda.min.ratio
      } else if(all(fit.df$df > nsel)) {
        lambdamax <- max(abs(t(x) %*% (y - mean(y)*(1 - mean(y)))))/(n*seq.alpha[a])
        lambdamin <- max(fit.df$lambda)
        inbet <- TRUE
        count <- count + 1
        if(count==5) {
          found <- TRUE
          lambda <- c(max(abs(t(x) %*% (y - mean(y)*(1 - mean(y)))))/(n*seq.alpha[a]),
                      fit.df$lambda[which.max(fit.df$df[which.min(abs(fit.df$df - nsel))])])
        }
      } else if(length(unique(lambdaseq))==1) {
        found <- TRUE
        lambda <- c(max(abs(t(x) %*% (y - mean(y)*(1 - mean(y)))))/(n*seq.alpha[a]),
                    fit.df$lambda[which.max(fit.df$df[which.min(abs(fit.df$df - nsel))])])
      } else {
        ind <- c(rle(fit.df$df < nsel)$lengths[1], rle(fit.df$df < nsel)$lengths[1] + 1)
        lambdamax <- max(fit.df$lambda[ind])
        lambdamin <- min(fit.df$lambda[ind])
        inbet <- TRUE
      }
      if(!found) {
        lambdaseq <- exp(seq(log(lambdamax), log(lambdamin), length.out=nlambda))
      }    
    }
    fit.cv <- cv.glmnet(x, y, family="binomial", lambda=lambda, alpha=seq.alpha[a], 
                        standardize=FALSE, intercept=intercept)
    seq.cvll[a] <- ifelse(type=="ridge", min(fit.cv$cvm), tail(fit.cv$cvm, n=1L))
    seq.lam[a] <- ifelse(type=="ridge", fit.cv$lambda.min, tail(fit.cv$lambda, n=1L))
    if(type=="ridge") {
      seq.df[a] <- p
    } else {
      seq.df[a] <- ifelse(any(fit.df$df==nsel),
                          fit.df$df[fit.df$lambda==lambda[2]],
                          max(fit.df$df[which.min(abs(fit.df$df - nsel))]))
    }
  }
  indmin <- which.min(seq.cvll)
  
  cvll <- seq.cvll[indmin]
  df <- seq.df[indmin]
  lambdaglmnet <- seq.lam[indmin]
  alphaglmnet <- seq.alpha[indmin]
  
  lambda1pen <- lambdaglmnet*alphaglmnet*n
  lambda2pen <- lambdaglmnet*(1 - alphaglmnet)*n*0.5
  lambda1bayes <- lambdaglmnet*alphaglmnet*n*2
  lambda2bayes <- lambdaglmnet*(1 - alphaglmnet)*n
  
  out <- list(lambdaglmnet=lambdaglmnet, alphaglmnet=alphaglmnet, lambda1pen=lambda1pen, 
              lambda2pen=lambda2pen, lambda1bayes=lambda1bayes, lambda2bayes=lambda2bayes, 
              cvll=cvll, df=df)
  
}

# function to run grVBEM with post hoc variable selection
# grVBEM2sel <- function(x, y, m, partitions, lambda1=NULL, lambda2=NULL, nsel, intercept, 
#                        monotone, posterior, eps, maxiter, trace=TRUE, alphastart) {
#   
#   nparts <- length(partitions)
#   sizes <- lapply(partitions, function(part) {rle(sort(part))$lengths})
#   fit.grVBEM <- grVBEM2(x, y, m, partitions, lambda1, lambda2, intercept, monotone, 
#                         posterior, eps, maxiter, trace, alphastart)
#   lambdagvec <- exp(colSums(log(sapply(1:nparts, function(part) {
#     rep(fit.grVBEM$lambdag[[part]][, fit.grVBEM$nouteriter + 1], times=sizes[[part]])}))))
#   
#   lambda1 <- lambda1/2
#   lambda2 <- lambda2/2
#   fit.grensel <- penalized(y, x, ~1, lambda1*sqrt(lambdagvec), lambda2*lambdagvec, 
#                            model="logistic", trace=FALSE)
#   lambda1min <- 0
#   lambda1max <- max(abs(as.numeric(t(x) %*% (y - 0.5)))/sqrt(lambdagvec))
#   nselcur <- sum(fit.grensel@penalized!=0)
#   
#   while(nselcur!=nsel) {
#     if(nselcur < nsel) {
#       lambda1max <- lambda1
#     } else if(nselcur > nsel) {
#       lambda1min <- lambda1
#     }
#     lambda1 <- (lambda1max + lambda1min)/2
#     fit.grensel <- penalized(y, x, ~1, lambda1*sqrt(lambdagvec), lambda2*lambdagvec, 
#                              model="logistic", trace=FALSE)
#     nselcur <- sum(fit.grensel@penalized!=0)
#   }
#   
#   return(fit.grensel)
#   
# }

grVBEM2sel <- function(x, y, m, partitions, lambda1=NULL, lambda2=NULL, nsel, intercept, 
                       monotone, posterior, eps, maxiter, trace=TRUE, alphastart) {
  
  nparts <- length(partitions)
  sizes <- lapply(partitions, function(part) {rle(sort(part))$lengths})
  fit.grVBEM <- grVBEM2(x, y, m, partitions, lambda1, lambda2, intercept, monotone, 
                        posterior, eps, maxiter, trace, alphastart)
  lambdagvec <- exp(colSums(log(sapply(1:nparts, function(part) {
    rep(fit.grVBEM$lambdag[[part]][, fit.grVBEM$nouteriter + 1], times=sizes[[part]])}))))
  
  alpha <- lambda1/(2*lambda2 + lambda1)
  
  # selecting variables
  p <- ncol(x)
  n <- nrow(x)
  lambda.min.ratio <- ifelse(n < p, 0.01, 0.0001)
  nlambda <- 100
  

  lambdaseq <- NULL
  found <- inbet <- FALSE
  count <- 0
  while(!found) {
    suppressWarnings(fit.df <- glmnet(x, y, family="binomial", alpha=alpha, 
                                      nlambda=nlambda, lambda=lambdaseq, standardize=FALSE, 
                                      intercept=intercept, 
                                      dfmax=ifelse(is.null(lambdaseq), nsel, p + 1)))
    if(any(fit.df$df==nsel)) {
      if(tail(fit.df$df, n=1L)==nsel & !inbet) {
      
        lambdamax <- tail(fit.df$lambda[fit.df$df==nsel], n=1L)
        lambdamin <- lambdamax*lambda.min.ratio
      } else {
        found <- TRUE
      }
    } else if(all(fit.df$df < nsel)) {
      lambdamax <- tail(fit.df$lambda, n=1L)
      lambdamin <- lambdamax*lambda.min.ratio
    } else if(all(fit.df$df > nsel)) {
      lambdamax <- max(abs(t(x) %*% (y - mean(y)*(1 - mean(y)))))/(n*alpha)
      lambdamin <- max(fit.df$lambda)
      inbet <- TRUE
      count <- count + 1
      if(count==5) {
        found <- TRUE
      }
    } else if(length(unique(lambdaseq))==1) {
      found <- TRUE
    } else {
      ind <- c(rle(fit.df$df < nsel)$lengths[1], rle(fit.df$df < nsel)$lengths[1] + 1)
      lambdamax <- max(fit.df$lambda[ind])
      lambdamin <- min(fit.df$lambda[ind])
      inbet <- TRUE
    }
    if(!found) {
      lambdaseq <- exp(seq(log(lambdamax), log(lambdamin), length.out=nlambda))
    }    
  }
  ind <- ifelse(any(fit.df$df==nsel), tail(which(fit.df$df==nsel), n=1L),
                which.max(fit.df$df[which.min(abs(fit.df$df - nsel))]))
  indsel <- which(as.numeric(fit.df$beta[, ind])!=0)
  return(indsel)
}

### the simulation
set.seed(2017)
p <- 500
G <- 5
pg <- p/G
meanbeta <- 0.03
nblock <- 20
n <- 50
rho <- 0.3
sigma2 <- 1
m <- rep(1, n)
q <- 0.9
f <- 1.7
groups <- rep(1:G, each=pg)
partition1 <- list(groups=groups)
partition2 <- list(groups=CreatePartition(as.factor(groups)))
pblock <- p/nblock
ubeta <- rep(rev(sapply(0:(G - 1), function(g) {c(f^(-g), 0)})), 
             times=rep(c(pg*q, pg*(1 - q) + 1), times=G))
beta <- ubeta*meanbeta/mean(ubeta)
sigma <- matrix(rho, ncol=pblock, nrow=pblock)
diag(sigma) <- sigma2

ntest <- 1000
nreps <- 50
nselseq <- c(1, 2, 3, 5, 10, 15, 20, 30)

aucmat <- kappamat <- precmat <- recmat <- matrix(NA, ncol=4, nrow=nreps*length(nselseq))
nselmat <- matrix(NA, ncol=4*G, nrow=nreps*length(nselseq))

for(r in 1:nreps) {
  
  # creating the data
  x <- do.call(cbind, replicate(nblock, rmvnorm(n, mean=rep(0, pblock), sigma=sigma), 
                                simplify=FALSE))
  prob <- as.numeric(exp(x %*% beta)/(1 + exp(x %*% beta)))
  y <- rbinom(n, m, prob)
  
  # creating test data
  xtest <- do.call(cbind, replicate(nblock, rmvnorm(ntest, mean=rep(0, pblock), 
                                                    sigma=sigma), simplify=FALSE))
  probtest <- as.numeric(exp(xtest %*% beta)/(1 + exp(xtest %*% beta)))
  ytest <- rbinom(ntest, rep(m, ntest), probtest)
    
  # different number of selected variables
  for(cursel in 1:length(nselseq)) {
    
    nsel <- nselseq[cursel]
    
    print(paste("################### rep ", r, ", nsel ", nsel, " ###################", sep=""))
    
    # cross-validating the various penalty parameters
    penparsen <- cv.pensel(x, y, intercept=TRUE, nsel=nsel, type="enet")
    penparsridge <- cv.pensel(x, y, intercept=TRUE, nsel=nsel, type="ridge")
    penparslasso <- cv.pensel(x, y, intercept=TRUE, nsel=nsel, type="lasso")
    
    # fitting the different models
    fit.grridge <- grridge(t(x), y, partition2, selectionEN=TRUE, maxsel=nsel, stepsel=nsel,
                           optl=penparsridge$lambda2pen, trace=FALSE)
    fit.lasso <- glmnet(x, y, family="binomial", alpha=1, lambda=penparslasso$lambdaglmnet,
                        standardize=FALSE)
    fit.en <- glmnet(x, y, family="binomial", alpha=penparsen$alphaglmnet, 
                     lambda=penparsen$lambdaglmnet, standardize=FALSE)
    fit.gren <- grVBEM2sel(x, y, m, partition1, penparsen$lambda1bayes, penparsen$lambda2bayes, 
                           nsel, TRUE, FALSE, FALSE, 0.001, 1000, FALSE, NULL)
    
    # which variables are selected per method
    sel.grridge <- as.numeric(fit.grridge$resEN$whichEN)
    sel.lasso <- which(as.numeric(coef(fit.lasso)[-1])!=0)
    sel.en <- which(as.numeric(coef(fit.en)[-1])!=0)
    sel.gren <- as.numeric(fit.gren)
    
    # refit for predictions
    fit.grridge2 <- cv.glmnet(cbind(1, x[, sel.grridge]), y, family="binomial", alpha=0, standardize=FALSE)
    fit.lasso2 <- cv.glmnet(cbind(1, x[, sel.lasso]), y, family="binomial", alpha=0, standardize=FALSE)
    fit.en2 <- cv.glmnet(cbind(1, x[, sel.en]), y, family="binomial", alpha=0, standardize=FALSE)
    fit.gren2 <- cv.glmnet(cbind(1, x[, sel.gren]), y, family="binomial", alpha=0, standardize=FALSE)
    
    # predict using new data
    pred.grridge <- as.numeric(predict(fit.grridge2, newx=cbind(1, xtest[, sel.grridge]), s="lambda.min"))
    pred.lasso <- as.numeric(predict(fit.lasso2, newx=cbind(1, xtest[, sel.lasso]), s="lambda.min"))
    pred.en <- as.numeric(predict(fit.en2, newx=cbind(1, xtest[, sel.en]), s="lambda.min"))
    pred.gren <- as.numeric(predict(fit.gren2, newx=cbind(1, xtest[, sel.gren]), s="lambda.min"))
    
    # calculating the auc's on the test data
    auc.grridge <- as.numeric(pROC::roc(ytest, pred.grridge)$auc)
    auc.lasso <- as.numeric(pROC::roc(ytest, pred.lasso)$auc)
    auc.en <- as.numeric(pROC::roc(ytest, pred.en)$auc)
    auc.gren <- as.numeric(pROC::roc(ytest, pred.gren)$auc)
    
    # confusion tables
    tab.grridge <- table(as.numeric(c(1:p) %in% sel.grridge), as.numeric(beta!=0))
    tab.lasso <- table(as.numeric(c(1:p) %in% sel.lasso), as.numeric(beta!=0))
    tab.en <- table(as.numeric(c(1:p) %in% sel.en), as.numeric(beta!=0))
    tab.gren <- table(as.numeric(c(1:p) %in% sel.gren), as.numeric(beta!=0))
    
    # calculating the kappa of variable selection
    kappa.grridge <- cohen.kappa(tab.grridge)$kappa
    kappa.lasso <- cohen.kappa(tab.lasso)$kappa
    kappa.en <- cohen.kappa(tab.en)$kappa
    kappa.gren <- cohen.kappa(tab.gren)$kappa
    
    # precisions
    prec.grridge <- tab.grridge[2, 2]/(tab.grridge[2, 2] + tab.grridge[2, 1])
    prec.lasso <- tab.lasso[2, 2]/(tab.lasso[2, 2] + tab.lasso[2, 1])
    prec.en <- tab.en[2, 2]/(tab.en[2, 2] + tab.en[2, 1])
    prec.gren <- tab.gren[2, 2]/(tab.gren[2, 2] + tab.gren[2, 1])
    
    # recalls
    rec.grridge <- tab.grridge[2, 2]/(tab.grridge[2, 2] + tab.grridge[1, 2])
    rec.lasso <- tab.lasso[2, 2]/(tab.lasso[2, 2] + tab.lasso[1, 2])
    rec.en <- tab.en[2, 2]/(tab.en[2, 2] + tab.en[1, 2])
    rec.gren <- tab.gren[2, 2]/(tab.gren[2, 2] + tab.gren[1, 2])
    
    # determining the number of selected variables per group
    nsel.grridge <- sapply(1:G, function(g) {sum(sel.grridge %in% c(((g - 1)*pg + 1):(g*pg)))})
    nsel.lasso <- sapply(1:G, function(g) {sum(sel.lasso %in% c(((g - 1)*pg + 1):(g*pg)))})
    nsel.en <- sapply(1:G, function(g) {sum(sel.en %in% c(((g - 1)*pg + 1):(g*pg)))})
    nsel.gren <- sapply(1:G, function(g) {sum(sel.gren %in% c(((g - 1)*pg + 1):(g*pg)))})
    
    # assigning the metrics to a matrix row
    aucmat[(cursel - 1)*nreps + r, ] <- c(auc.grridge, auc.lasso, auc.en, auc.gren)
    kappamat[(cursel - 1)*nreps + r, ]  <- c(kappa.grridge, kappa.lasso, kappa.en, kappa.gren) 
    precmat[(cursel - 1)*nreps + r, ]  <- c(prec.grridge, prec.lasso, prec.en, prec.gren) 
    recmat[(cursel - 1)*nreps + r, ]  <- c(rec.grridge, rec.lasso, rec.en, rec.gren) 
    nselmat[(cursel - 1)*nreps + r, ] <- c(nsel.grridge, nsel.lasso, nsel.en, nsel.gren)
    
    # res1 is with q=0.7
    # res2 is with q=0.9 and meanbeta=0.03
    res3 <- list(auc=aucmat, kappa=kappamat, precision=precmat, recall=recmat, nsel=nselmat)
    save(res3, file=paste(path.res, "grVBEM_sim5_res3.Rdata", sep=""))
  }
  
}

# load(paste(path.res, "grVBEM_sim5_res1.Rdata", sep=""))
# boxplot(res1$kappa[c(1:13), ], names=c("grridge", "lasso", "enet", "gren"),
#         main=paste(nselseq[1], "variable selected"))
# boxplot(res1$kappa[c((1*nreps + 1):(1*nreps + 13)), ], names=c("grridge", "lasso", "enet", "gren"),
#         main=paste(nselseq[2], "variable selected"))
# boxplot(res1$kappa[c((2*nreps + 1):(2*nreps + 13)), ], names=c("grridge", "lasso", "enet", "gren"),
#         main=paste(nselseq[3], "variable selected"))
# boxplot(res1$kappa[c((3*nreps + 1):(3*nreps + 13)), ], names=c("grridge", "lasso", "enet", "gren"),
#         main=paste(nselseq[4], "variable selected"))
# boxplot(res1$kappa[c((4*nreps + 1):(4*nreps + 13)), ], names=c("grridge", "lasso", "enet", "gren"),
#         main=paste(nselseq[5], "variable selected"))
# boxplot(res1$kappa[c((5*nreps + 1):(5*nreps + 13)), ], names=c("grridge", "lasso", "enet", "gren"),
#         main=paste(nselseq[6], "variable selected"))
# boxplot(res1$kappa[c((6*nreps + 1):(6*nreps + 13)), ], names=c("grridge", "lasso", "enet", "gren"),
#         main=paste(nselseq[7], "variable selected"))
# boxplot(res1$kappa[c((7*nreps + 1):(7*nreps + 13)), ], names=c("grridge", "lasso", "enet", "gren"),
#         main=paste(nselseq[8], "variable selected"))
# 
# tabkap <- cbind(nselseq, t(sapply(1:length(nselseq), function(cursel) {
#   apply(res1$kappa[c(((cursel - 1)*nreps + 1):((cursel - 1)*nreps + nreps)), ], 2,
#         median, na.rm=TRUE)})))
# colnames(tabkap) <- c("nsel", "grridge", "lasso", "enet", "gren")
# tabkap

load(paste(path.res, "grVBEM_sim5_res2.Rdata", sep=""))
f1score <- 2*res2$precision*res2$recall/(res2$precision + res2$recall)

tabkap <- cbind(nselseq, t(sapply(1:length(nselseq), function(cursel) {
  apply(res2$kappa[c(((cursel - 1)*nreps + 1):((cursel - 1)*nreps + nreps)), ], 2,
        median, na.rm=TRUE)})))
colnames(tabkap) <- c("nsel", "grridge", "lasso", "enet", "gren")

tabprec <- cbind(nselseq, t(sapply(1:length(nselseq), function(cursel) {
  apply(res2$precision[c(((cursel - 1)*nreps + 1):((cursel - 1)*nreps + nreps)), ], 2,
        median, na.rm=TRUE)})))
colnames(tabprec) <- c("nsel", "grridge", "lasso", "enet", "gren")

tabrec <- cbind(nselseq, t(sapply(1:length(nselseq), function(cursel) {
  apply(res2$recall[c(((cursel - 1)*nreps + 1):((cursel - 1)*nreps + nreps)), ], 2,
        median, na.rm=TRUE)})))
colnames(tabrec) <- c("nsel", "grridge", "lasso", "enet", "gren")

tabf1 <- cbind(nselseq, t(sapply(1:length(nselseq), function(cursel) {
  apply(f1score[c(((cursel - 1)*nreps + 1):((cursel - 1)*nreps + nreps)), ], 2,
        median, na.rm=TRUE)})))
colnames(tabf1) <- c("nsel", "grridge", "lasso", "enet", "gren")

tabauc <- cbind(nselseq, t(sapply(1:length(nselseq), function(cursel) {
  apply(res2$auc[c(((cursel - 1)*nreps + 1):((cursel - 1)*nreps + nreps)), ], 2,
        median, na.rm=TRUE)})))
colnames(tabauc) <- c("nsel", "grridge", "lasso", "enet", "gren")

png(paste(path.graph, "grVBEM_sim5_res2_kappa.png", sep=""), width=6, height=6, units='in', res=300)
plot(grridge ~ nsel, data=tabkap, type="l", ylim=range(tabkap[, -1]), 
     main="Cohen's kappa against number of selected variables",
     xlab=expression(hat(p)), ylab=expression(kappa))
lines(lasso ~ nsel, data=tabkap, type="l", col=2)
lines(enet ~ nsel, data=tabkap, type="l", col=3)
lines(gren ~ nsel, data=tabkap, type="l", col=4)
legend("bottomright", legend=c("GRridge", "lasso", "enet", "GRenet"), lty=1, col=c(1:4))
dev.off()

png(paste(path.graph, "grVBEM_sim5_res2_f1.png", sep=""), width=6, height=6, units='in', res=300)
plot(grridge ~ nsel, data=tabf1, type="l", ylim=range(tabf1[, -1]), 
     main=expression(bold(F[1]~"score against number of selected variables")),
     xlab=expression(hat(p)), ylab="F1")
lines(lasso ~ nsel, data=tabf1, type="l", col=2)
lines(enet ~ nsel, data=tabf1, type="l", col=3)
lines(gren ~ nsel, data=tabf1, type="l", col=4)
legend("bottomright", legend=c("GRridge", "lasso", "enet", "GRenet"), lty=1, col=c(1:4))
dev.off()

png(paste(path.graph, "grVBEM_sim5_res2_auc.png", sep=""), width=6, height=6, units='in', res=300)
plot(grridge ~ nsel, data=tabauc, type="l", ylim=range(tabauc[, -1]), 
     main="AUC after selection against number of selected variables",
     xlab=expression(hat(p)), ylab=expression(F[1]))
lines(lasso ~ nsel, data=tabauc, type="l", col=2)
lines(enet ~ nsel, data=tabauc, type="l", col=3)
lines(gren ~ nsel, data=tabauc, type="l", col=4)
legend("bottomright", legend=c("GRridge", "lasso", "enet", "GRenet"), lty=1, col=c(1:4))
dev.off()

plot(grridge ~ nsel, data=tabprec, type="l", ylim=range(tabprec[, -1]), 
     main="Precision against number of selected variables",
     xlab=expression(hat(p)), ylab="precision")
lines(lasso ~ nsel, data=tabprec, type="l", col=2)
lines(enet ~ nsel, data=tabprec, type="l", col=3)
lines(gren ~ nsel, data=tabprec, type="l", col=4)
legend("topright", legend=c("GRridge", "lasso", "enet", "GRenet"), lty=1, col=c(1:4))

plot(grridge ~ nsel, data=tabrec, type="l", ylim=range(tabrec[, -1]), 
     main="Recall against number of selected variables",
     xlab=expression(hat(p)), ylab="recall")
lines(lasso ~ nsel, data=tabrec, type="l", col=2)
lines(enet ~ nsel, data=tabrec, type="l", col=3)
lines(gren ~ nsel, data=tabrec, type="l", col=4)
legend("bottomright", legend=c("GRridge", "lasso", "enet", "GRenet"), lty=1, col=c(1:4))

plot(grridge ~ nsel, data=tabf1, type="l", ylim=range(tabf1[, -1]), 
     main=expression(bold(F[1]~"score against number of selected variables")),
     xlab=expression(hat(p)), ylab="F1")
lines(lasso ~ nsel, data=tabf1, type="l", col=2)
lines(enet ~ nsel, data=tabf1, type="l", col=3)
lines(gren ~ nsel, data=tabf1, type="l", col=4)
legend("bottomright", legend=c("GRridge", "lasso", "enet", "GRenet"), lty=1, col=c(1:4))

plot(grridge ~ nsel, data=tabauc, type="l", ylim=range(tabauc[, -1]), 
     main="AUC after selection against number of selected variables",
     xlab=expression(hat(p)), ylab=expression(F[1]))
lines(lasso ~ nsel, data=tabauc, type="l", col=2)
lines(enet ~ nsel, data=tabauc, type="l", col=3)
lines(gren ~ nsel, data=tabauc, type="l", col=4)
legend("bottomright", legend=c("GRridge", "lasso", "enet", "GRenet"), lty=1, col=c(1:4))





load(paste(path.res, "grVBEM_sim5_res3.Rdata", sep=""))
f1score <- 2*res3$precision*res3$recall/(res3$precision + res3$recall)

tabkap <- cbind(nselseq, t(sapply(1:length(nselseq), function(cursel) {
  apply(res3$kappa[c(((cursel - 1)*nreps + 1):((cursel - 1)*nreps + nreps)), ], 2,
        median, na.rm=TRUE)})))
colnames(tabkap) <- c("nsel", "grridge", "lasso", "enet", "gren")

tabauc <- cbind(nselseq, t(sapply(1:length(nselseq), function(cursel) {
  apply(res3$auc[c(((cursel - 1)*nreps + 1):((cursel - 1)*nreps + nreps)), ], 2,
        median, na.rm=TRUE)})))
colnames(tabauc) <- c("nsel", "grridge", "lasso", "enet", "gren")

tabf1 <- cbind(nselseq, t(sapply(1:length(nselseq), function(cursel) {
  apply(f1score[c(((cursel - 1)*nreps + 1):((cursel - 1)*nreps + nreps)), ], 2,
        median, na.rm=TRUE)})))
colnames(tabf1) <- c("nsel", "grridge", "lasso", "enet", "gren")

png(paste(path.graph, "grVBEM_sim5_res3_kappa.png", sep=""), width=6, height=6, units='in', res=300)
plot(grridge ~ nsel, data=tabkap, type="l", ylim=range(tabkap[, -1]), 
     main="Cohen's kappa against number of selected variables",
     xlab=expression(hat(p)), ylab=expression(kappa))
lines(lasso ~ nsel, data=tabkap, type="l", col=2)
lines(enet ~ nsel, data=tabkap, type="l", col=3)
lines(gren ~ nsel, data=tabkap, type="l", col=4)
legend("bottomright", legend=c("GRridge", "lasso", "enet", "GRenet"), lty=1, col=c(1:4))
dev.off()

png(paste(path.graph, "grVBEM_sim5_res3_f1.png", sep=""), width=6, height=6, units='in', res=300)
plot(grridge ~ nsel, data=tabf1, type="l", ylim=range(tabf1[, -1]), 
     main=expression(bold(F[1]~"score against number of selected variables")),
     xlab=expression(hat(p)), ylab="F1")
lines(lasso ~ nsel, data=tabf1, type="l", col=2)
lines(enet ~ nsel, data=tabf1, type="l", col=3)
lines(gren ~ nsel, data=tabf1, type="l", col=4)
legend("bottomright", legend=c("GRridge", "lasso", "enet", "GRenet"), lty=1, col=c(1:4))
dev.off()

png(paste(path.graph, "grVBEM_sim5_res3_auc.png", sep=""), width=6, height=6, units='in', res=300)
plot(grridge ~ nsel, data=tabauc, type="l", ylim=range(tabauc[, -1]), 
     main="AUC after selection against number of selected variables",
     xlab=expression(hat(p)), ylab=expression(F[1]))
lines(lasso ~ nsel, data=tabauc, type="l", col=2)
lines(enet ~ nsel, data=tabauc, type="l", col=3)
lines(gren ~ nsel, data=tabauc, type="l", col=4)
legend("bottomright", legend=c("GRridge", "lasso", "enet", "GRenet"), lty=1, col=c(1:4))
dev.off()




fint <- function(s, gamma) {
  return(exp(-0.5*log(s/(1 - s))^2/gamma^2)/(1 - s))
}

vary <- function(gamma) {
  # gamma <- sqrt(t(beta) %*% sigma %*% as.matrix(beta))
  const <- 1/(gamma*sqrt(2*pi))
  int <- sapply(gamma, function(g) {integrate(fint, lower=0, upper=1, gamma=g)$value})
  return(const*int - const^2*int^2)
}

seq.gam <- seq(8.5, 9, 0.01)
plot(seq.gam, vary(seq.gam), type="l")



rho <- 0.7
sigma2 <- 0.1
sigma <- matrix(rho, ncol=pblock, nrow=pblock)
diag(sigma) <- sigma2
p <- 25
n <- 10000
beta <- c(1:p)/10000
x <- rmvnorm(n, rep(0, p), sigma=sigma)
y <- rbinom(n, rep(1, n), prob=exp(x %*% beta)/(1 + exp(x %*% beta)))
var(y)



vary(beta, sigma)

beta <- c(1:p)
seq.m <- exp(seq(log(1000), log(10000), length.out=100))
plot(seq.m, sapply(seq.m, function(m) {vary(beta/m, sigma)}), type="l")

sigma2 <- 1
sigma <- matrix(rho, ncol=pblock, nrow=pblock)
diag(sigma) <- sigma2
sqrt(t(beta) %*% sigma %*% as.matrix(beta))
beta <- 1:25
seq.m <- 1:10
plot(seq.m, sapply(seq.m, function(m) {sqrt(t(m*beta) %*% sigma %*% as.matrix(m*beta))}), type="l")



