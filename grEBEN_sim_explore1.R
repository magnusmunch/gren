##############################  preamble  #############################
# Simulations                                                         #
# version: 01                                                         #
# author: Magnus Munch                                                #
# created: 26-09-2017                                                 #
# last edited: 26-09-2017                                             #
#######################################################################

###############################  notes  ###############################
#                                                                     #
#######################################################################

# paths
path.code <- as.character(ifelse(Sys.info()[1]=="Darwin","/Users/magnusmunch/Documents/PhD/EBEN/code/" ,"~/EBEN/code/"))
path.graph <- "/Users/magnusmunch/Documents/PhD/EBEN/graphs/"
path.res <- as.character(ifelse(Sys.info()[1]=="Darwin","/Users/magnusmunch/Documents/PhD/EBEN/results/" ,"~/EBEN/results/"))

### libraries
library(psych)
library(mvtnorm)
library(pROC)

### functions
# source function for variational Bayes
source(paste(path.code, "grVBEM.R", sep=""))

# source version of grridge that works
source(paste(path.code, "mygrridge.R", sep=""))

# various functions of the elastic net prior
dtau <- function(x, lambda1, lambda2, log=FALSE) {
  scale <- 8*lambda2/lambda1^2
  if(log) {
    dens <- -(log(0.5) - 0.5*log(x*(x > 1)) - 0.5*log(scale) - 0.5*log(pi) - x/scale - 
                pnorm(-sqrt(2/scale), log.p=TRUE))
  } else {
    dens <- 0.5*sqrt((x > 1)/(ifelse(x==0, 1, x)*scale*pi))*exp(-x/scale)/pnorm(-sqrt(2/scale))
  }
  return(dens)
}

ptau <- function(q, lambda1, lambda2, log.p=FALSE) {
  scale <- 8*lambda2/lambda1^2
  denom <- pnorm(-sqrt(2/scale))
  if(log.p) {
    prob <- log(1 - pnorm(-sqrt(2*ifelse(q <= 1, 1, q)/scale)/denom))
  } else {
    prob <- 1 - pnorm(-sqrt(2*ifelse(q <= 1, 1, q)/scale)/denom)
  }
  return(prob)
}

qtau <- function(p, lambda1, lambda2, log.p=FALSE) {
  scale <- 8*lambda2/lambda1^2
  if(log.p) {
    aux <- pnorm(-sqrt(2/scale), log.p=log.p)
    tau <- 0.5*scale*qnorm(exp(aux) - exp(aux + p))^2
  } else {
    tau <- 0.5*scale*qnorm((1 - p)*pnorm(-sqrt(2/scale)))^2
  }
  return(tau)
}

rtau <- function(n, lambda1, lambda2, log.p=FALSE) {
  u <- runif(n)
  if(log.p) {
    tau <- qtau(log(u), lambda1, lambda2, log.p=log.p)
  } else {
    tau <- qtau(u, lambda1, lambda2, log.p=log.p)
  }
  return(tau)
}

denet <- function(x, lambda1=1, lambda2=1, log=FALSE) {
  if(log) {
    dens <- 0.5*log(lambda2) - log(2) - 0.5*lambda1*abs(x) - 0.5*lambda2*x^2 +
      dnorm(0.5*lambda1/sqrt(lambda2), log=TRUE) - pnorm(-0.5*lambda1/sqrt(lambda2), log.p=TRUE)
  } else {
    dens <- 0.5*sqrt(lambda2)*exp(-0.5*(lambda1*abs(x) + lambda2*x^2))*
      exp(dnorm(0.5*lambda1/sqrt(lambda2), log=TRUE) - 
            pnorm(-0.5*lambda1/sqrt(lambda2), log.p=TRUE))
  }
  return(dens)
}

renet <- function(n, lambda1=1, lambda2=1, log.p=FALSE) {
  if(lambda1==0) {
    beta <- rnorm(n, 0, sqrt(1/lambda2))
  } else if(lambda2==0) {
    u <- runif(n) - 0.5
    beta <- -2*sign(u)*log(1 - 2*abs(u))/lambda1
  } else {
    tau <- rtau(n, lambda1, lambda2, log.p=log.p)
    sigma <- (tau - 1)/(lambda2*tau)
    beta <- rnorm(n, 0, sqrt(sigma))
  }
  return(beta)
}

varenet <- function(lambda1, lambda2) {
  return(lambda1^2/(4*lambda2^2) + 1/lambda2 - lambda1/(2*lambda2^(3/2))*
           exp(dnorm(lambda1/(2*sqrt(lambda2)), log=TRUE) - 
                 pnorm(-lambda1/(2*sqrt(lambda2)), log.p=TRUE)))
}

# ### setting 1
# n <- 200
# p <- 500
# G <- 5
# 
# lambda1 <- 1
# lambda2 <- 1
# lambdag <- exp(seq(-2, 2, length.out=G))
# groups <- rep(1:G, each=p/G)
# partitions1 <- list(groups=CreatePartition(as.factor(groups)))
# partitions2 <- list(groups=groups)
# 
# methods <- c("enet", "enettrue", "grridge", "greben.beta", "greben.mu", "grebentrue.beta", 
#              "grebentrue.mu", "greben0.05.beta", "greben0.05.mu", "ridge", "lasso")
# ntest <- 1000
# nreps <- 50
# msemat <- cormat <- aucmat <- briermat <- matrix(NA, nrow=length(methods), ncol=nreps)
# varbetamat <- varbenetmat <- varbenettruemat <- matrix(NA, nrow=G, ncol=nreps)
# for(k in 1:nreps) {
#   
#   set.seed(2000 + k)
#   cat(paste("Iteration ", k, "\n"))
#   beta <- as.numeric(sapply(1:length(lambdag), function(g) {
#     renet(p/G, lambda1*sqrt(lambdag[g]), lambda2*lambdag[g], log.p=FALSE)}))
#   x <- matrix(rnorm(n*p), ncol=p, nrow=n)
#   m <- rep(1, n)
#   y <- rbinom(n, m, as.numeric(exp(x %*% beta)/(1 + exp(x %*% beta))))
#   
#   xtest <- matrix(rnorm(ntest*p), ncol=p, nrow=ntest)
#   mtest <- rep(1, ntest)
#   ytest <- rbinom(ntest, mtest, as.numeric(exp(xtest %*% beta)/(1 + exp(xtest %*% beta))))
# 
#   est1.lambda <- cv.pen(x, y, intercept=TRUE)
#   est2.lambda <- cv.glmnet(x, y, family="binomial", alpha=0.05, standardize=FALSE)
#   est3.lambda <- cv.glmnet(x, y, family="binomial", alpha=1, standardize=FALSE)
#   est4.lambda <- cv.glmnet(x, y, family="binomial", alpha=0, standardize=FALSE)
#   
#   fit1.enet <- glmnet(x, y, family="binomial", alpha=est1.lambda$alpha[which.min(est1.lambda$cvll)], 
#                       lambda=est1.lambda$lambda[which.min(est1.lambda$cvll)], standardize=FALSE)
#   fit2.enet <- glmnet(x, y, family="binomial", alpha=lambda1/(lambda1 + 2*lambda2), 
#                       lambda=(0.5*lambda1 + lambda2)/n, standardize=FALSE)
#   fit1.grridge <- grridge(t(x), y, partitions1)
#   fit1.greben <- grEBEN(x, y, m, partitions=partitions2, 
#                         lambda1=est1.lambda$lambda1bayes[which.min(est1.lambda$cvll)],
#                         lambda2=est1.lambda$lambda2bayes[which.min(est1.lambda$cvll)])
#   fit2.greben <- grEBEN(x, y, m, partitions=partitions2, lambda1=lambda2, lambda2=lambda2)
#   fit3.greben <- grEBEN(x, y, m, partitions=partitions2, 
#                         lambda1=2*n*est2.lambda$lambda.min*0.05,
#                         lambda2=n*est2.lambda$lambda.min*0.95)
#   fit1.ridge <- glmnet(x, y, family="binomial", alpha=0, lambda=est4.lambda$lambda.min, 
#                        standardize=FALSE)
#   fit1.lasso <- glmnet(x, y, family="binomial", alpha=1, lambda=est3.lambda$lambda.min, 
#                        standardize=FALSE)
# 
#   msemat[, k] <- c(mean((c(0, beta) - as.numeric(coef(fit1.enet)))^2),
#                    mean((c(0, beta) - as.numeric(coef(fit2.enet)))^2),
#                    mean((c(0, beta) - c(predict.grridge(fit1.grridge, t(matrix(0, ncol=p)))[, 2], 
#                                         fit1.grridge$betas))^2),
#                    mean((c(0, beta) - as.numeric(fit1.greben$beta))^2),
#                    mean((c(0, beta) - fit1.greben$mu)^2),
#                    mean((c(0, beta) - as.numeric(fit2.greben$beta))^2),
#                    mean((c(0, beta) - fit2.greben$mu)^2),
#                    mean((c(0, beta) - as.numeric(fit3.greben$beta))^2),
#                    mean((c(0, beta) - fit3.greben$mu)^2),
#                    mean((c(0, beta) - as.numeric(coef(fit1.ridge)))^2),
#                    mean((c(0, beta) - as.numeric(coef(fit1.lasso)))^2))
#   
#   cormat[, k] <- c(cor(c(0, beta), as.numeric(coef(fit1.enet))),
#                    cor(c(0, beta), as.numeric(coef(fit2.enet))),
#                    cor(c(0, beta), c(predict.grridge(fit1.grridge, t(matrix(0, ncol=p)))[, 2], 
#                                      fit1.grridge$betas)),
#                    cor(c(0, beta), as.numeric(fit1.greben$beta)),
#                    cor(c(0, beta), fit1.greben$mu),
#                    cor(c(0, beta), as.numeric(fit2.greben$beta)),
#                    cor(c(0, beta), fit2.greben$mu),
#                    cor(c(0, beta), as.numeric(fit3.greben$beta)),
#                    cor(c(0, beta), fit3.greben$mu),
#                    cor(c(0, beta), as.numeric(coef(fit1.ridge))),
#                    cor(c(0, beta), as.numeric(coef(fit1.lasso))))
#   
#   predmat <- rbind(as.numeric(predict(fit1.enet, xtest, type="response")),
#                    as.numeric(predict(fit2.enet, xtest, type="response")),
#                    predict.grridge(fit1.grridge, t(xtest))[, 2],
#                    predict.grEBEN(fit1.greben, xtest, type="penalized"),
#                    predict.grEBEN(fit1.greben, xtest, type="VB"),
#                    predict.grEBEN(fit2.greben, xtest, type="penalized"),
#                    predict.grEBEN(fit2.greben, xtest, type="VB"),
#                    predict.grEBEN(fit3.greben, xtest, type="penalized"),
#                    predict.grEBEN(fit3.greben, xtest, type="VB"),
#                    as.numeric(predict(fit1.ridge, xtest, type="response")),
#                    as.numeric(predict(fit1.lasso, xtest, type="response")))
#   
#   aucmat[, k] <- sapply(1:nrow(aucmat), function(method) {
#     pROC::roc(ytest, predmat[method, ])$auc})
#   
#   briermat[, k] <- sapply(1:nrow(aucmat), function(method) {
#     mean((ytest - predmat[method, ])^2)})
#   
#   varbetamat[, k] <- sapply(1:G, function(g) {var(beta[((g - 1)*p/G + 1):(g*p/G)])})
#   
#   varbenetmat[, k] <- sapply(1:G, function(g) {
#     var(as.numeric(coef(fit1.enet))[-1][((g - 1)*p/G + 1):(g*p/G)])})
#  
#   varbenettruemat[, k] <- sapply(1:G, function(g) {
#     var(as.numeric(coef(fit2.enet))[-1][((g - 1)*p/G + 1):(g*p/G)])})
#   
# }
# 
# rownames(msemat) <- rownames(cormat) <- rownames(aucmat) <- rownames(briermat) <- methods
# rownames(varbetamat) <- rownames(varbenetmat) <- rownames(varbenettruemat) <- 
#   paste("group", c(1:G), sep="")
# results1 <- list(mse=msemat, cor=cormat, auc=aucmat, brier=briermat, varbeta=varbetamat, 
#                  varbenet=varbenetmat, varbenettrue=varbenettruemat)
# # save(results1, file=paste(path.res, "grEBEN_sim_explore1_res1.Rdata", sep=""))




# ### setting 2
# n <- 200
# p <- 500
# G <- 5
# 
# rho <- 0.9
# lambda1 <- 1
# lambda2 <- 1
# lambdag <- exp(seq(-2, 2, length.out=G))
# groups <- rep(1:G, each=p/G)
# partitions1 <- list(groups=CreatePartition(as.factor(groups)))
# partitions2 <- list(groups=groups)
# 
# methods <- c("enet", "enettrue", "grridge", "greben.beta", "greben.mu", "grebentrue.beta", 
#              "grebentrue.mu", "greben0.05.beta", "greben0.05.mu", "ridge", "lasso")
# ntest <- 1000
# nreps <- 50
# msemat <- cormat <- aucmat <- briermat <- matrix(NA, nrow=length(methods), ncol=nreps)
# pencormat  <- penrankmat <- penmsemat <- matrix(NA, nrow=4, ncol=nreps)
# varbetamat <- varbenetmat <- varbenettruemat <- matrix(NA, nrow=G, ncol=nreps)
# for(k in 1:nreps) {
#   
#   set.seed(3000 + k)
#   cat(paste("Iteration ", k, "\n"))
#   beta <- as.numeric(sapply(1:length(lambdag), function(g) {
#     renet(p/G, lambda1*sqrt(lambdag[g]), lambda2*lambdag[g], log.p=FALSE)}))
#   sigma <- matrix(rho, ncol=100, nrow=100)
#   diag(sigma) <- 1
#   x <- matrix(NA, ncol=p, nrow=n)
#   x[, c(1:50, 451:500)] <- rmvnorm(n, mean=rep(0, 100), sigma=sigma)
#   x[, c(101:150, 351:400)] <- rmvnorm(n, mean=rep(0, 100), sigma=sigma)
#   x[, c(201:300)] <- rmvnorm(n, mean=rep(0, 100), sigma=sigma)
#   x[, c(51:100, 151:200, 301:350, 401:450)] <- matrix(rnorm(200*n), ncol=200, nrow=n)
#   m <- rep(1, n)
#   y <- rbinom(n, m, as.numeric(exp(x %*% beta)/(1 + exp(x %*% beta))))
#   
#   xtest <- matrix(NA, ncol=p, nrow=ntest)
#   xtest[, c(1:50, 451:500)] <- rmvnorm(ntest, mean=rep(0, 100), sigma=sigma)
#   xtest[, c(101:150, 351:400)] <- rmvnorm(ntest, mean=rep(0, 100), sigma=sigma)
#   xtest[, c(201:300)] <- rmvnorm(ntest, mean=rep(0, 100), sigma=sigma)
#   xtest[, c(51:100, 151:200, 301:350, 401:450)] <- matrix(rnorm(200*ntest), ncol=200, nrow=ntest)
#   mtest <- rep(1, ntest)
#   ytest <- rbinom(ntest, mtest, as.numeric(exp(xtest %*% beta)/(1 + exp(xtest %*% beta))))
#   
#   est1.lambda <- cv.pen(x, y, intercept=TRUE)
#   est2.lambda <- cv.glmnet(x, y, family="binomial", alpha=0.05, standardize=FALSE)
#   est3.lambda <- cv.glmnet(x, y, family="binomial", alpha=1, standardize=FALSE)
#   est4.lambda <- cv.glmnet(x, y, family="binomial", alpha=0, standardize=FALSE)
#   
#   fit1.enet <- glmnet(x, y, family="binomial", alpha=est1.lambda$alpha[which.min(est1.lambda$cvll)], 
#                       lambda=est1.lambda$lambda[which.min(est1.lambda$cvll)], standardize=FALSE)
#   fit2.enet <- glmnet(x, y, family="binomial", alpha=lambda1/(lambda1 + 2*lambda2), 
#                       lambda=(0.5*lambda1 + lambda2)/n, standardize=FALSE)
#   fit1.grridge <- grridge(t(x), y, partitions1)
#   fit1.greben <- grEBEN(x, y, m, partitions=partitions2, 
#                         lambda1=est1.lambda$lambda1bayes[which.min(est1.lambda$cvll)],
#                         lambda2=est1.lambda$lambda2bayes[which.min(est1.lambda$cvll)])
#   fit2.greben <- grEBEN(x, y, m, partitions=partitions2, lambda1=lambda2, lambda2=lambda2)
#   fit3.greben <- grEBEN(x, y, m, partitions=partitions2, 
#                         lambda1=2*n*est2.lambda$lambda.min*0.05,
#                         lambda2=n*est2.lambda$lambda.min*0.95)
#   fit1.ridge <- glmnet(x, y, family="binomial", alpha=0, lambda=est4.lambda$lambda.min, 
#                        standardize=FALSE)
#   fit1.lasso <- glmnet(x, y, family="binomial", alpha=1, lambda=est3.lambda$lambda.min, 
#                        standardize=FALSE)
#   
#   msemat[, k] <- c(mean((c(0, beta) - as.numeric(coef(fit1.enet)))^2),
#                    mean((c(0, beta) - as.numeric(coef(fit2.enet)))^2),
#                    mean((c(0, beta) - c(predict.grridge(fit1.grridge, t(matrix(0, ncol=p)))[, 2], 
#                                         fit1.grridge$betas))^2),
#                    mean((c(0, beta) - as.numeric(fit1.greben$beta))^2),
#                    mean((c(0, beta) - fit1.greben$mu)^2),
#                    mean((c(0, beta) - as.numeric(fit2.greben$beta))^2),
#                    mean((c(0, beta) - fit2.greben$mu)^2),
#                    mean((c(0, beta) - as.numeric(fit3.greben$beta))^2),
#                    mean((c(0, beta) - fit3.greben$mu)^2),
#                    mean((c(0, beta) - as.numeric(coef(fit1.ridge)))^2),
#                    mean((c(0, beta) - as.numeric(coef(fit1.lasso)))^2))
#   
#   cormat[, k] <- c(cor(c(0, beta), as.numeric(coef(fit1.enet))),
#                    cor(c(0, beta), as.numeric(coef(fit2.enet))),
#                    cor(c(0, beta), c(predict.grridge(fit1.grridge, t(matrix(0, ncol=p)))[, 2], 
#                                      fit1.grridge$betas)),
#                    cor(c(0, beta), as.numeric(fit1.greben$beta)),
#                    cor(c(0, beta), fit1.greben$mu),
#                    cor(c(0, beta), as.numeric(fit2.greben$beta)),
#                    cor(c(0, beta), fit2.greben$mu),
#                    cor(c(0, beta), as.numeric(fit3.greben$beta)),
#                    cor(c(0, beta), fit3.greben$mu),
#                    cor(c(0, beta), as.numeric(coef(fit1.ridge))),
#                    cor(c(0, beta), as.numeric(coef(fit1.lasso))))
#   
#   predmat <- rbind(as.numeric(predict(fit1.enet, xtest, type="response")),
#                    as.numeric(predict(fit2.enet, xtest, type="response")),
#                    predict.grridge(fit1.grridge, t(xtest))[, 2],
#                    predict.grEBEN(fit1.greben, xtest, type="penalized"),
#                    predict.grEBEN(fit1.greben, xtest, type="VB"),
#                    predict.grEBEN(fit2.greben, xtest, type="penalized"),
#                    predict.grEBEN(fit2.greben, xtest, type="VB"),
#                    predict.grEBEN(fit3.greben, xtest, type="penalized"),
#                    predict.grEBEN(fit3.greben, xtest, type="VB"),
#                    as.numeric(predict(fit1.ridge, xtest, type="response")),
#                    as.numeric(predict(fit1.lasso, xtest, type="response")))
#   
#   aucmat[, k] <- sapply(1:nrow(aucmat), function(method) {
#     pROC::roc(ytest, predmat[method, ])$auc})
#   
#   briermat[, k] <- sapply(1:nrow(aucmat), function(method) {
#     mean((ytest - predmat[method, ])^2)})
#   
#   penmat <- cbind(lambdag, fit1.grridge$lambdamults$groups, 
#                   fit1.greben$lambdag$groups[, fit1.greben$nouteriter + 1],
#                   fit2.greben$lambdag$groups[, fit2.greben$nouteriter + 1],
#                   fit3.greben$lambdag$groups[, fit3.greben$nouteriter + 1])
#   pencormat[, k] <- cor(penmat)[-1, 1]
#   penrankmat[, k] <- cor(penmat, method="spearman")[-1, 1]
#   penmsemat[, k] <- apply(penmat[, -1], 2, function(x) {mean((x - penmat[, 1])^2)})
#   
#   varbetamat[, k] <- sapply(1:G, function(g) {var(beta[((g - 1)*p/G + 1):(g*p/G)])})
#   
#   varbenetmat[, k] <- sapply(1:G, function(g) {
#     var(as.numeric(coef(fit1.enet))[-1][((g - 1)*p/G + 1):(g*p/G)])})
#   
#   varbenettruemat[, k] <- sapply(1:G, function(g) {
#     var(as.numeric(coef(fit2.enet))[-1][((g - 1)*p/G + 1):(g*p/G)])})
#   
# }
# 
# rownames(msemat) <- rownames(cormat) <- rownames(aucmat) <- rownames(briermat) <- methods
# rownames(pencormat) <- rownames(penrankmat) <- rownames(penmsemat) <- methods[c(3, 4, 6, 8)]
# rownames(varbetamat) <- rownames(varbenetmat) <- rownames(varbenettruemat) <- 
#   paste("group", c(1:G), sep="")
# results2 <- list(mse=msemat, cor=cormat, auc=aucmat, brier=briermat, pencor=pencormat,
#                  penrank=penrankmat, penmse=penmsemat, varbeta=varbetamat, 
#                  varbenet=varbenetmat, varbenettrue=varbenettruemat)
# save(results2, file=paste(path.res, "grEBEN_sim_explore1_res2.Rdata", sep=""))





# ### setting 3
# n <- 100
# p <- 200
# G <- 4
# rho <- 0.8
# lambda1 <- 1
# lambda2 <- 1
# lambdag <- exp(seq(-1.5, 1.5, length.out=G))
# groups <- rep(1:G, p/G)
# partitions1 <- list(groups=CreatePartition(as.factor(groups)))
# partitions2 <- list(groups=groups)
# 
# methods <- c("enet", "enettrue", "grridge", "greben.beta", "greben.mu", "grebentrue.beta", 
#              "grebentrue.mu", "greben0.05.beta", "greben0.05.mu", "ridge", "lasso")
# ntest <- 1000
# nreps <- 50
# msemat <- cormat <- aucmat <- briermat <- matrix(NA, nrow=length(methods), ncol=nreps)
# pencormat  <- penrankmat <- penmsemat <- matrix(NA, nrow=4, ncol=nreps)
# varbetamat <- varbenetmat <- varbenettruemat <- matrix(NA, nrow=G, ncol=nreps)
# for(k in 1:nreps) {
#   
#   set.seed(4000 + k)
#   cat(paste("Iteration ", k, "\n"))
#   beta <- as.numeric(sapply(1:length(lambdag), function(g) {
#     renet(p/G, lambda1*sqrt(lambdag[g]), lambda2*lambdag[g], log.p=FALSE)}))
#   x <- matrix(rnorm(n*p/2), ncol=p/2, nrow=n)
#   x <- cbind(x, rho*x + sqrt(1 - rho^2)*rnorm(n*p/2))
#   m <- rep(1, n)
#   y <- rbinom(n, m, as.numeric(exp(x %*% beta)/(1 + exp(x %*% beta))))
# 
#   xtest <- matrix(rnorm(n*p/2), ncol=p/2, nrow=ntest)
#   xtest <- cbind(xtest, rho*xtest + sqrt(1 - rho^2)*rnorm(ntest*p/2))
#   mtest <- rep(1, ntest)
#   ytest <- rbinom(ntest, mtest, as.numeric(exp(xtest %*% beta)/(1 + exp(xtest %*% beta))))
# 
#   est1.lambda <- cv.pen(x, y, intercept=TRUE)
#   est2.lambda <- cv.glmnet(x, y, family="binomial", alpha=0.05, standardize=FALSE)
#   est3.lambda <- cv.glmnet(x, y, family="binomial", alpha=1, standardize=FALSE)
#   est4.lambda <- cv.glmnet(x, y, family="binomial", alpha=0, standardize=FALSE)
# 
#   fit1.enet <- glmnet(x, y, family="binomial", alpha=est1.lambda$alpha[which.min(est1.lambda$cvll)], 
#                       lambda=est1.lambda$lambda[which.min(est1.lambda$cvll)], standardize=FALSE)
#   fit2.enet <- glmnet(x, y, family="binomial", alpha=lambda1/(lambda1 + 2*lambda2), 
#                       lambda=(0.5*lambda1 + lambda2)/n, standardize=FALSE)
#   fit1.grridge <- grridge(t(x), y, partitions1)
#   fit1.greben <- grEBEN(x, y, m, partitions=partitions2, 
#                         lambda1=est1.lambda$lambda1bayes[which.min(est1.lambda$cvll)],
#                         lambda2=est1.lambda$lambda2bayes[which.min(est1.lambda$cvll)])
#   fit2.greben <- grEBEN(x, y, m, partitions=partitions2, lambda1=lambda2, lambda2=lambda2)
#   fit3.greben <- grEBEN(x, y, m, partitions=partitions2, 
#                         lambda1=2*n*est2.lambda$lambda.min*0.05,
#                         lambda2=n*est2.lambda$lambda.min*0.95)
#   fit1.ridge <- glmnet(x, y, family="binomial", alpha=0, lambda=est4.lambda$lambda.min, 
#                        standardize=FALSE)
#   fit1.lasso <- glmnet(x, y, family="binomial", alpha=1, lambda=est3.lambda$lambda.min, 
#                        standardize=FALSE)
# 
#   msemat[, k] <- c(mean((c(0, beta) - as.numeric(coef(fit1.enet)))^2),
#                    mean((c(0, beta) - as.numeric(coef(fit2.enet)))^2),
#                    mean((c(0, beta) - c(predict.grridge(fit1.grridge, t(matrix(0, ncol=p)))[, 2], 
#                                         fit1.grridge$betas))^2),
#                    mean((c(0, beta) - as.numeric(fit1.greben$beta))^2),
#                    mean((c(0, beta) - fit1.greben$mu)^2),
#                    mean((c(0, beta) - as.numeric(fit2.greben$beta))^2),
#                    mean((c(0, beta) - fit2.greben$mu)^2),
#                    mean((c(0, beta) - as.numeric(fit3.greben$beta))^2),
#                    mean((c(0, beta) - fit3.greben$mu)^2),
#                    mean((c(0, beta) - as.numeric(coef(fit1.ridge)))^2),
#                    mean((c(0, beta) - as.numeric(coef(fit1.lasso)))^2))
#   
#   cormat[, k] <- c(cor(c(0, beta), as.numeric(coef(fit1.enet))),
#                    cor(c(0, beta), as.numeric(coef(fit2.enet))),
#                    cor(c(0, beta), c(predict.grridge(fit1.grridge, t(matrix(0, ncol=p)))[, 2], 
#                                      fit1.grridge$betas)),
#                    cor(c(0, beta), as.numeric(fit1.greben$beta)),
#                    cor(c(0, beta), fit1.greben$mu),
#                    cor(c(0, beta), as.numeric(fit2.greben$beta)),
#                    cor(c(0, beta), fit2.greben$mu),
#                    cor(c(0, beta), as.numeric(fit3.greben$beta)),
#                    cor(c(0, beta), fit3.greben$mu),
#                    cor(c(0, beta), as.numeric(coef(fit1.ridge))),
#                    cor(c(0, beta), as.numeric(coef(fit1.lasso))))
#   
#   predmat <- rbind(as.numeric(predict(fit1.enet, xtest, type="response")),
#                    as.numeric(predict(fit2.enet, xtest, type="response")),
#                    predict.grridge(fit1.grridge, t(xtest))[, 2],
#                    predict.grEBEN(fit1.greben, xtest, type="penalized"),
#                    predict.grEBEN(fit1.greben, xtest, type="VB"),
#                    predict.grEBEN(fit2.greben, xtest, type="penalized"),
#                    predict.grEBEN(fit2.greben, xtest, type="VB"),
#                    predict.grEBEN(fit3.greben, xtest, type="penalized"),
#                    predict.grEBEN(fit3.greben, xtest, type="VB"),
#                    as.numeric(predict(fit1.ridge, xtest, type="response")),
#                    as.numeric(predict(fit1.lasso, xtest, type="response")))
#   
#   aucmat[, k] <- sapply(1:nrow(aucmat), function(method) {
#     pROC::roc(ytest, predmat[method, ])$auc})
#   
#   briermat[, k] <- sapply(1:nrow(aucmat), function(method) {
#     mean((ytest - predmat[method, ])^2)})
#   
#   penmat <- cbind(lambdag, fit1.grridge$lambdamults$groups, 
#                   fit1.greben$lambdag$groups[, fit1.greben$nouteriter + 1],
#                   fit2.greben$lambdag$groups[, fit2.greben$nouteriter + 1],
#                   fit3.greben$lambdag$groups[, fit3.greben$nouteriter + 1])
#   pencormat[, k] <- cor(penmat)[-1, 1]
#   penrankmat[, k] <- cor(penmat, method="spearman")[-1, 1]
#   penmsemat[, k] <- apply(penmat[, -1], 2, function(x) {mean((x - penmat[, 1])^2)})
#   
#   varbetamat[, k] <- sapply(1:G, function(g) {var(beta[((g - 1)*p/G + 1):(g*p/G)])})
#   
#   varbenetmat[, k] <- sapply(1:G, function(g) {
#     var(as.numeric(coef(fit1.enet))[-1][((g - 1)*p/G + 1):(g*p/G)])})
#   
#   varbenettruemat[, k] <- sapply(1:G, function(g) {
#     var(as.numeric(coef(fit2.enet))[-1][((g - 1)*p/G + 1):(g*p/G)])})
# }
# 
# rownames(msemat) <- rownames(cormat) <- rownames(aucmat) <- rownames(briermat) <- methods
# rownames(pencormat) <- rownames(penrankmat) <- rownames(penmsemat) <- methods[c(3, 4, 6, 8)]
# rownames(varbetamat) <- rownames(varbenetmat) <- rownames(varbenettruemat) <- 
#   paste("group", c(1:G), sep="")
# results3 <- list(mse=msemat, cor=cormat, auc=aucmat, brier=briermat, pencor=pencormat,
#                  penrank=penrankmat, penmse=penmsemat, varbeta=varbetamat, 
#                  varbenet=varbenetmat, varbenettrue=varbenettruemat)
# save(results3, file=paste(path.res, "grEBEN_sim_explore1_res3.Rdata", sep=""))


  

# ### setting 4
# n <- 100
# p <- 500
# G <- 10
# sigma2 <- 1
# rho <- 0.8
# lambda1 <- 1
# lambda2 <- 1
# lambdag <- exp(seq(-1.5, 1.5, length.out=G))
# groups <- rep(1:G, p/G)
# partitions1 <- list(groups=CreatePartition(as.factor(groups)))
# partitions2 <- list(groups=groups)
# 
# methods <- c("enet", "enettrue", "grridge", "greben.beta", "greben.mu", "grebentrue.beta", 
#              "grebentrue.mu", "greben0.05.beta", "greben0.05.mu", "ridge", "lasso")
# ntest <- 1000
# nreps <- 50
# msemat <- cormat <- aucmat <- briermat <- matrix(NA, nrow=length(methods), ncol=nreps)
# pencormat  <- penrankmat <- penmsemat <- matrix(NA, nrow=4, ncol=nreps)
# varbetamat <- varbenetmat <- varbenettruemat <- matrix(NA, nrow=G, ncol=nreps)
# pselmat <- matrix(NA, nrow=6, ncol=nreps)
# for(k in 1:nreps) {
#   
#   set.seed(5000 + k)
#   cat(paste("Iteration ", k, "\n"))
#   beta <- as.numeric(sapply(1:length(lambdag), function(g) {
#     renet(p/G, lambda1*sqrt(lambdag[g]), lambda2*lambdag[g], log.p=FALSE)}))
#   x <- matrix(rnorm(n*p/2, rep(0, n*p/2), sigma2), ncol=p/2, nrow=n)
#   x <- cbind(x, rho*x[, c((p/2):1)] + sqrt(1 - rho^2)*rnorm(n*p/2, rep(0, n*p/2), sigma2))
#   m <- rep(1, n)
#   y <- rbinom(n, m, as.numeric(exp(x %*% beta)/(1 + exp(x %*% beta))))
#   
#   xtest <- matrix(rnorm(ntest*p/2, rep(0, n*p/2), sigma2), ncol=p/2, nrow=ntest)
#   xtest <- cbind(xtest, rho*xtest[, c((p/2):1)] + 
#                    sqrt(1 - rho^2)*rnorm(ntest*p/2, rep(0, n*p/2), sigma2))
#   mtest <- rep(1, ntest)
#   ytest <- rbinom(ntest, mtest, as.numeric(exp(xtest %*% beta)/(1 + exp(xtest %*% beta))))
#   
#   est1.lambda <- cv.pen(x, y, intercept=TRUE)
#   est2.lambda <- cv.glmnet(x, y, family="binomial", alpha=0.05, standardize=FALSE)
#   est3.lambda <- cv.glmnet(x, y, family="binomial", alpha=1, standardize=FALSE)
#   est4.lambda <- cv.glmnet(x, y, family="binomial", alpha=0, standardize=FALSE)
#   
#   fit1.enet <- glmnet(x, y, family="binomial", alpha=est1.lambda$alpha[which.min(est1.lambda$cvll)], 
#                       lambda=est1.lambda$lambda[which.min(est1.lambda$cvll)], standardize=FALSE)
#   fit2.enet <- glmnet(x, y, family="binomial", alpha=lambda1/(lambda1 + 2*lambda2), 
#                       lambda=(0.5*lambda1 + lambda2)/n, standardize=FALSE)
#   fit1.grridge <- grridge(t(x), y, partitions1)
#   fit1.greben <- grEBEN(x, y, m, partitions=partitions2, 
#                         lambda1=est1.lambda$lambda1bayes[which.min(est1.lambda$cvll)],
#                         lambda2=est1.lambda$lambda2bayes[which.min(est1.lambda$cvll)])
#   fit2.greben <- grEBEN(x, y, m, partitions=partitions2, lambda1=lambda2, lambda2=lambda2)
#   fit3.greben <- grEBEN(x, y, m, partitions=partitions2, 
#                         lambda1=2*n*est2.lambda$lambda.min*0.05,
#                         lambda2=n*est2.lambda$lambda.min*0.95)
#   fit1.ridge <- glmnet(x, y, family="binomial", alpha=0, lambda=est4.lambda$lambda.min, 
#                        standardize=FALSE)
#   fit1.lasso <- glmnet(x, y, family="binomial", alpha=1, lambda=est3.lambda$lambda.min, 
#                        standardize=FALSE)
#   
#   msemat[, k] <- c(mean((c(0, beta) - as.numeric(coef(fit1.enet)))^2),
#                    mean((c(0, beta) - as.numeric(coef(fit2.enet)))^2),
#                    mean((c(0, beta) - c(predict.grridge(fit1.grridge, t(matrix(0, ncol=p)))[, 2], 
#                                         fit1.grridge$betas))^2),
#                    mean((c(0, beta) - as.numeric(fit1.greben$beta))^2),
#                    mean((c(0, beta) - fit1.greben$mu)^2),
#                    mean((c(0, beta) - as.numeric(fit2.greben$beta))^2),
#                    mean((c(0, beta) - fit2.greben$mu)^2),
#                    mean((c(0, beta) - as.numeric(fit3.greben$beta))^2),
#                    mean((c(0, beta) - fit3.greben$mu)^2),
#                    mean((c(0, beta) - as.numeric(coef(fit1.ridge)))^2),
#                    mean((c(0, beta) - as.numeric(coef(fit1.lasso)))^2))
#   
#   cormat[, k] <- c(cor(c(0, beta), as.numeric(coef(fit1.enet))),
#                    cor(c(0, beta), as.numeric(coef(fit2.enet))),
#                    cor(c(0, beta), c(predict.grridge(fit1.grridge, t(matrix(0, ncol=p)))[, 2], 
#                                      fit1.grridge$betas)),
#                    cor(c(0, beta), as.numeric(fit1.greben$beta)),
#                    cor(c(0, beta), fit1.greben$mu),
#                    cor(c(0, beta), as.numeric(fit2.greben$beta)),
#                    cor(c(0, beta), fit2.greben$mu),
#                    cor(c(0, beta), as.numeric(fit3.greben$beta)),
#                    cor(c(0, beta), fit3.greben$mu),
#                    cor(c(0, beta), as.numeric(coef(fit1.ridge))),
#                    cor(c(0, beta), as.numeric(coef(fit1.lasso))))
#   
#   predmat <- rbind(as.numeric(predict(fit1.enet, xtest, type="response")),
#                    as.numeric(predict(fit2.enet, xtest, type="response")),
#                    predict.grridge(fit1.grridge, t(xtest))[, 2],
#                    predict.grEBEN(fit1.greben, xtest, type="penalized"),
#                    predict.grEBEN(fit1.greben, xtest, type="VB"),
#                    predict.grEBEN(fit2.greben, xtest, type="penalized"),
#                    predict.grEBEN(fit2.greben, xtest, type="VB"),
#                    predict.grEBEN(fit3.greben, xtest, type="penalized"),
#                    predict.grEBEN(fit3.greben, xtest, type="VB"),
#                    as.numeric(predict(fit1.ridge, xtest, type="response")),
#                    as.numeric(predict(fit1.lasso, xtest, type="response")))
#   
#   aucmat[, k] <- sapply(1:nrow(aucmat), function(method) {
#     pROC::roc(ytest, predmat[method, ])$auc})
#   
#   briermat[, k] <- sapply(1:nrow(aucmat), function(method) {
#     mean((ytest - predmat[method, ])^2)})
#   
#   penmat <- cbind(lambdag, fit1.grridge$lambdamults$groups, 
#                   fit1.greben$lambdag$groups[, fit1.greben$nouteriter + 1],
#                   fit2.greben$lambdag$groups[, fit2.greben$nouteriter + 1],
#                   fit3.greben$lambdag$groups[, fit3.greben$nouteriter + 1])
#   pencormat[, k] <- cor(penmat)[-1, 1]
#   penrankmat[, k] <- cor(penmat, method="spearman")[-1, 1]
#   penmsemat[, k] <- apply(penmat[, -1], 2, function(x) {mean((x - penmat[, 1])^2)})
#   
#   pselmat[, k] <- c(fit1.enet$df, fit2.enet$df, sum(fit1.greben$beta[-1]!=0),
#                     sum(fit2.greben$beta[-1]!=0), sum(fit3.greben$beta[-1]!=0), fit1.lasso$df)
#   
#   varbetamat[, k] <- sapply(1:G, function(g) {var(beta[((g - 1)*p/G + 1):(g*p/G)])})
#   
#   varbenetmat[, k] <- sapply(1:G, function(g) {
#     var(as.numeric(coef(fit1.enet))[-1][((g - 1)*p/G + 1):(g*p/G)])})
#   
#   varbenettruemat[, k] <- sapply(1:G, function(g) {
#     var(as.numeric(coef(fit2.enet))[-1][((g - 1)*p/G + 1):(g*p/G)])})
# }
# 
# rownames(msemat) <- rownames(cormat) <- rownames(aucmat) <- rownames(briermat) <- methods
# rownames(pencormat) <- rownames(penrankmat) <- rownames(penmsemat) <- methods[c(3, 4, 6, 8)]
# rownames(pselmat) <- methods[c(1, 2, 4, 6, 8, 11)]
# rownames(varbetamat) <- rownames(varbenetmat) <- rownames(varbenettruemat) <- 
#   paste("group", c(1:G), sep="")
# results4 <- list(mse=msemat, cor=cormat, auc=aucmat, brier=briermat, pencor=pencormat,
#                  penrank=penrankmat, penmse=penmsemat, psel=pselmat, varbeta=varbetamat, 
#                  varbenet=varbenetmat, varbenettrue=varbenettruemat)
# save(results4, file=paste(path.res, "grEBEN_sim_explore1_res4.Rdata", sep=""))




# ### setting 5
# n <- 150
# p <- 500
# G <- 10
# sigma2 <- 1
# rho <- 0.8
# lambda1 <- 1
# lambda2 <- 1
# lambdag <- exp(seq(-2, 2, length.out=G))
# groups <- rep(1:G, p/G)
# partitions1 <- list(groups=CreatePartition(as.factor(groups)))
# partitions2 <- list(groups=groups)
# 
# methods <- c("enet", "enettrue", "grridge", "greben.beta", "greben.mu", "grebentrue.beta", 
#              "grebentrue.mu", "greben0.05.beta", "greben0.05.mu", "ridge", "lasso")
# ntest <- 1000
# nreps <- 50
# msemat <- cormat <- aucmat <- briermat <- matrix(NA, nrow=length(methods), ncol=nreps)
# pencormat  <- penrankmat <- penmsemat <- matrix(NA, nrow=4, ncol=nreps)
# varbetamat <- varbenetmat <- varbenettruemat <- matrix(NA, nrow=G, ncol=nreps)
# pselmat <- matrix(NA, nrow=6, ncol=nreps)
# for(k in 1:nreps) {
#   
#   set.seed(5000 + k)
#   cat(paste("Iteration ", k, "\n"))
#   beta <- as.numeric(sapply(1:length(lambdag), function(g) {
#     renet(p/G, lambda1*sqrt(lambdag[g]), lambda2*lambdag[g], log.p=FALSE)}))
#   x <- matrix(rnorm(n*p/2, rep(0, n*p/2), sigma2), ncol=p/2, nrow=n)
#   x <- cbind(x, rho*x[, c((p/2):1)] + sqrt(1 - rho^2)*rnorm(n*p/2, rep(0, n*p/2), sigma2))
#   m <- rep(1, n)
#   y <- rbinom(n, m, as.numeric(exp(x %*% beta)/(1 + exp(x %*% beta))))
#   
#   xtest <- matrix(rnorm(ntest*p/2, rep(0, n*p/2), sigma2), ncol=p/2, nrow=ntest)
#   xtest <- cbind(xtest, rho*xtest[, c((p/2):1)] + 
#                    sqrt(1 - rho^2)*rnorm(ntest*p/2, rep(0, n*p/2), sigma2))
#   mtest <- rep(1, ntest)
#   ytest <- rbinom(ntest, mtest, as.numeric(exp(xtest %*% beta)/(1 + exp(xtest %*% beta))))
#   
#   est1.lambda <- cv.pen(x, y, intercept=TRUE)
#   est2.lambda <- cv.glmnet(x, y, family="binomial", alpha=0.05, standardize=FALSE)
#   est3.lambda <- cv.glmnet(x, y, family="binomial", alpha=1, standardize=FALSE)
#   est4.lambda <- cv.glmnet(x, y, family="binomial", alpha=0, standardize=FALSE)
#   
#   fit1.enet <- glmnet(x, y, family="binomial", alpha=est1.lambda$alpha[which.min(est1.lambda$cvll)], 
#                       lambda=est1.lambda$lambda[which.min(est1.lambda$cvll)], standardize=FALSE)
#   fit2.enet <- glmnet(x, y, family="binomial", alpha=lambda1/(lambda1 + 2*lambda2), 
#                       lambda=(0.5*lambda1 + lambda2)/n, standardize=FALSE)
#   fit1.grridge <- grridge(t(x), y, partitions1)
#   fit1.greben <- grEBEN(x, y, m, partitions=partitions2, 
#                         lambda1=est1.lambda$lambda1bayes[which.min(est1.lambda$cvll)],
#                         lambda2=est1.lambda$lambda2bayes[which.min(est1.lambda$cvll)])
#   fit2.greben <- grEBEN(x, y, m, partitions=partitions2, lambda1=lambda2, lambda2=lambda2)
#   fit3.greben <- grEBEN(x, y, m, partitions=partitions2, 
#                         lambda1=2*n*est2.lambda$lambda.min*0.05,
#                         lambda2=n*est2.lambda$lambda.min*0.95)
#   fit1.ridge <- glmnet(x, y, family="binomial", alpha=0, lambda=est4.lambda$lambda.min, 
#                        standardize=FALSE)
#   fit1.lasso <- glmnet(x, y, family="binomial", alpha=1, lambda=est3.lambda$lambda.min, 
#                        standardize=FALSE)
#   
#   msemat[, k] <- c(mean((c(0, beta) - as.numeric(coef(fit1.enet)))^2),
#                    mean((c(0, beta) - as.numeric(coef(fit2.enet)))^2),
#                    mean((c(0, beta) - c(predict.grridge(fit1.grridge, t(matrix(0, ncol=p)))[, 2], 
#                                         fit1.grridge$betas))^2),
#                    mean((c(0, beta) - as.numeric(fit1.greben$beta))^2),
#                    mean((c(0, beta) - fit1.greben$mu)^2),
#                    mean((c(0, beta) - as.numeric(fit2.greben$beta))^2),
#                    mean((c(0, beta) - fit2.greben$mu)^2),
#                    mean((c(0, beta) - as.numeric(fit3.greben$beta))^2),
#                    mean((c(0, beta) - fit3.greben$mu)^2),
#                    mean((c(0, beta) - as.numeric(coef(fit1.ridge)))^2),
#                    mean((c(0, beta) - as.numeric(coef(fit1.lasso)))^2))
#   
#   cormat[, k] <- c(cor(c(0, beta), as.numeric(coef(fit1.enet))),
#                    cor(c(0, beta), as.numeric(coef(fit2.enet))),
#                    cor(c(0, beta), c(predict.grridge(fit1.grridge, t(matrix(0, ncol=p)))[, 2], 
#                                      fit1.grridge$betas)),
#                    cor(c(0, beta), as.numeric(fit1.greben$beta)),
#                    cor(c(0, beta), fit1.greben$mu),
#                    cor(c(0, beta), as.numeric(fit2.greben$beta)),
#                    cor(c(0, beta), fit2.greben$mu),
#                    cor(c(0, beta), as.numeric(fit3.greben$beta)),
#                    cor(c(0, beta), fit3.greben$mu),
#                    cor(c(0, beta), as.numeric(coef(fit1.ridge))),
#                    cor(c(0, beta), as.numeric(coef(fit1.lasso))))
#   
#   predmat <- rbind(as.numeric(predict(fit1.enet, xtest, type="response")),
#                    as.numeric(predict(fit2.enet, xtest, type="response")),
#                    predict.grridge(fit1.grridge, t(xtest))[, 2],
#                    predict.grEBEN(fit1.greben, xtest, type="penalized"),
#                    predict.grEBEN(fit1.greben, xtest, type="VB"),
#                    predict.grEBEN(fit2.greben, xtest, type="penalized"),
#                    predict.grEBEN(fit2.greben, xtest, type="VB"),
#                    predict.grEBEN(fit3.greben, xtest, type="penalized"),
#                    predict.grEBEN(fit3.greben, xtest, type="VB"),
#                    as.numeric(predict(fit1.ridge, xtest, type="response")),
#                    as.numeric(predict(fit1.lasso, xtest, type="response")))
#   
#   aucmat[, k] <- sapply(1:nrow(aucmat), function(method) {
#     pROC::roc(ytest, predmat[method, ])$auc})
#   
#   briermat[, k] <- sapply(1:nrow(aucmat), function(method) {
#     mean((ytest - predmat[method, ])^2)})
#   
#   penmat <- cbind(lambdag, fit1.grridge$lambdamults$groups, 
#                   fit1.greben$lambdag$groups[, fit1.greben$nouteriter + 1],
#                   fit2.greben$lambdag$groups[, fit2.greben$nouteriter + 1],
#                   fit3.greben$lambdag$groups[, fit3.greben$nouteriter + 1])
#   pencormat[, k] <- cor(penmat)[-1, 1]
#   penrankmat[, k] <- cor(penmat, method="spearman")[-1, 1]
#   penmsemat[, k] <- apply(penmat[, -1], 2, function(x) {mean((x - penmat[, 1])^2)})
#   
#   pselmat[, k] <- c(fit1.enet$df, fit2.enet$df, sum(fit1.greben$beta[-1]!=0),
#                     sum(fit2.greben$beta[-1]!=0), sum(fit3.greben$beta[-1]!=0), fit1.lasso$df)
#   
#   varbetamat[, k] <- sapply(1:G, function(g) {var(beta[((g - 1)*p/G + 1):(g*p/G)])})
#   
#   varbenetmat[, k] <- sapply(1:G, function(g) {
#     var(as.numeric(coef(fit1.enet))[-1][((g - 1)*p/G + 1):(g*p/G)])})
#   
#   varbenettruemat[, k] <- sapply(1:G, function(g) {
#     var(as.numeric(coef(fit2.enet))[-1][((g - 1)*p/G + 1):(g*p/G)])})
# }
# 
# rownames(msemat) <- rownames(cormat) <- rownames(aucmat) <- rownames(briermat) <- methods
# rownames(pencormat) <- rownames(penrankmat) <- rownames(penmsemat) <- methods[c(3, 4, 6, 8)]
# rownames(pselmat) <- methods[c(1, 2, 4, 6, 8, 11)]
# rownames(varbetamat) <- rownames(varbenetmat) <- rownames(varbenettruemat) <- 
#   paste("group", c(1:G), sep="")
# results5 <- list(mse=msemat, cor=cormat, auc=aucmat, brier=briermat, pencor=pencormat,
#                  penrank=penrankmat, penmse=penmsemat, psel=pselmat, varbeta=varbetamat, 
#                  varbenet=varbenetmat, varbenettrue=varbenettruemat)
# save(results5, file=paste(path.res, "grEBEN_sim_explore1_res5.Rdata", sep=""))



# ### setting 6
# n <- 150
# p <- 600
# G <- 6
# sigma2 <- 1
# rho <- 0.5
# lambda1 <- 1
# lambda2 <- 1
# lambdag <- exp(seq(-2, 2, length.out=G))
# groups <- rep(1:G, p/G)
# partitions1 <- list(groups=CreatePartition(as.factor(groups)))
# partitions2 <- list(groups=groups)
# 
# methods <- c("enet", "enettrue", "grridge", "greben.beta", "greben.mu", "grebentrue.beta", 
#              "grebentrue.mu", "greben0.05.beta", "greben0.05.mu", "ridge", "lasso")
# ntest <- 1000
# nreps <- 50
# msemat <- cormat <- aucmat <- briermat <- matrix(NA, nrow=length(methods), ncol=nreps)
# pencormat  <- penrankmat <- penmsemat <- matrix(NA, nrow=4, ncol=nreps)
# varbetamat <- varbenetmat <- varbenettruemat <- matrix(NA, nrow=G, ncol=nreps)
# pselmat <- matrix(NA, nrow=6, ncol=nreps)
# for(k in 1:nreps) {
#   
#   set.seed(6000 + k)
#   cat(paste("Iteration ", k, "\n"))
#   beta <- as.numeric(sapply(1:length(lambdag), function(g) {
#     renet(p/G, lambda1*sqrt(lambdag[g]), lambda2*lambdag[g], log.p=FALSE)}))
#   x <- matrix(rnorm(n*p/2, rep(0, n*p/2), sigma2), ncol=p/2, nrow=n)
#   x <- cbind(x, rho*x[, c((p/2):1)] + sqrt(1 - rho^2)*rnorm(n*p/2, rep(0, n*p/2), sigma2))
#   m <- rep(1, n)
#   y <- rbinom(n, m, as.numeric(exp(x %*% beta)/(1 + exp(x %*% beta))))
#   
#   xtest <- matrix(rnorm(ntest*p/2, rep(0, n*p/2), sigma2), ncol=p/2, nrow=ntest)
#   xtest <- cbind(xtest, rho*xtest[, c((p/2):1)] + 
#                    sqrt(1 - rho^2)*rnorm(ntest*p/2, rep(0, n*p/2), sigma2))
#   mtest <- rep(1, ntest)
#   ytest <- rbinom(ntest, mtest, as.numeric(exp(xtest %*% beta)/(1 + exp(xtest %*% beta))))
#   
#   est1.lambda <- cv.pen(x, y, intercept=TRUE)
#   est2.lambda <- cv.glmnet(x, y, family="binomial", alpha=0.05, standardize=FALSE)
#   est3.lambda <- cv.glmnet(x, y, family="binomial", alpha=1, standardize=FALSE)
#   est4.lambda <- cv.glmnet(x, y, family="binomial", alpha=0, standardize=FALSE)
#   
#   fit1.enet <- glmnet(x, y, family="binomial", alpha=est1.lambda$alpha[which.min(est1.lambda$cvll)], 
#                       lambda=est1.lambda$lambda[which.min(est1.lambda$cvll)], standardize=FALSE)
#   fit2.enet <- glmnet(x, y, family="binomial", alpha=lambda1/(lambda1 + 2*lambda2), 
#                       lambda=(0.5*lambda1 + lambda2)/n, standardize=FALSE)
#   fit1.grridge <- grridge(t(x), y, partitions1)
#   fit1.greben <- grEBEN(x, y, m, partitions=partitions2, 
#                         lambda1=est1.lambda$lambda1bayes[which.min(est1.lambda$cvll)],
#                         lambda2=est1.lambda$lambda2bayes[which.min(est1.lambda$cvll)])
#   fit2.greben <- grEBEN(x, y, m, partitions=partitions2, lambda1=lambda2, lambda2=lambda2)
#   fit3.greben <- grEBEN(x, y, m, partitions=partitions2, 
#                         lambda1=2*n*est2.lambda$lambda.min*0.05,
#                         lambda2=n*est2.lambda$lambda.min*0.95)
#   fit1.ridge <- glmnet(x, y, family="binomial", alpha=0, lambda=est4.lambda$lambda.min, 
#                        standardize=FALSE)
#   fit1.lasso <- glmnet(x, y, family="binomial", alpha=1, lambda=est3.lambda$lambda.min, 
#                        standardize=FALSE)
#   
#   msemat[, k] <- c(mean((c(0, beta) - as.numeric(coef(fit1.enet)))^2),
#                    mean((c(0, beta) - as.numeric(coef(fit2.enet)))^2),
#                    mean((c(0, beta) - c(predict.grridge(fit1.grridge, t(matrix(0, ncol=p)))[, 2], 
#                                         fit1.grridge$betas))^2),
#                    mean((c(0, beta) - as.numeric(fit1.greben$beta))^2),
#                    mean((c(0, beta) - fit1.greben$mu)^2),
#                    mean((c(0, beta) - as.numeric(fit2.greben$beta))^2),
#                    mean((c(0, beta) - fit2.greben$mu)^2),
#                    mean((c(0, beta) - as.numeric(fit3.greben$beta))^2),
#                    mean((c(0, beta) - fit3.greben$mu)^2),
#                    mean((c(0, beta) - as.numeric(coef(fit1.ridge)))^2),
#                    mean((c(0, beta) - as.numeric(coef(fit1.lasso)))^2))
#   
#   cormat[, k] <- c(cor(c(0, beta), as.numeric(coef(fit1.enet))),
#                    cor(c(0, beta), as.numeric(coef(fit2.enet))),
#                    cor(c(0, beta), c(predict.grridge(fit1.grridge, t(matrix(0, ncol=p)))[, 2], 
#                                      fit1.grridge$betas)),
#                    cor(c(0, beta), as.numeric(fit1.greben$beta)),
#                    cor(c(0, beta), fit1.greben$mu),
#                    cor(c(0, beta), as.numeric(fit2.greben$beta)),
#                    cor(c(0, beta), fit2.greben$mu),
#                    cor(c(0, beta), as.numeric(fit3.greben$beta)),
#                    cor(c(0, beta), fit3.greben$mu),
#                    cor(c(0, beta), as.numeric(coef(fit1.ridge))),
#                    cor(c(0, beta), as.numeric(coef(fit1.lasso))))
#   
#   predmat <- rbind(as.numeric(predict(fit1.enet, xtest, type="response")),
#                    as.numeric(predict(fit2.enet, xtest, type="response")),
#                    predict.grridge(fit1.grridge, t(xtest))[, 2],
#                    predict.grEBEN(fit1.greben, xtest, type="penalized"),
#                    predict.grEBEN(fit1.greben, xtest, type="VB"),
#                    predict.grEBEN(fit2.greben, xtest, type="penalized"),
#                    predict.grEBEN(fit2.greben, xtest, type="VB"),
#                    predict.grEBEN(fit3.greben, xtest, type="penalized"),
#                    predict.grEBEN(fit3.greben, xtest, type="VB"),
#                    as.numeric(predict(fit1.ridge, xtest, type="response")),
#                    as.numeric(predict(fit1.lasso, xtest, type="response")))
#   
#   aucmat[, k] <- sapply(1:nrow(aucmat), function(method) {
#     pROC::roc(ytest, predmat[method, ])$auc})
#   
#   briermat[, k] <- sapply(1:nrow(aucmat), function(method) {
#     mean((ytest - predmat[method, ])^2)})
#   
#   penmat <- cbind(lambdag, fit1.grridge$lambdamults$groups, 
#                   fit1.greben$lambdag$groups[, fit1.greben$nouteriter + 1],
#                   fit2.greben$lambdag$groups[, fit2.greben$nouteriter + 1],
#                   fit3.greben$lambdag$groups[, fit3.greben$nouteriter + 1])
#   pencormat[, k] <- cor(penmat)[-1, 1]
#   penrankmat[, k] <- cor(penmat, method="spearman")[-1, 1]
#   penmsemat[, k] <- apply(penmat[, -1], 2, function(x) {mean((x - penmat[, 1])^2)})
#   
#   pselmat[, k] <- c(fit1.enet$df, fit2.enet$df, sum(fit1.greben$beta[-1]!=0),
#                     sum(fit2.greben$beta[-1]!=0), sum(fit3.greben$beta[-1]!=0), fit1.lasso$df)
#   
#   varbetamat[, k] <- sapply(1:G, function(g) {var(beta[((g - 1)*p/G + 1):(g*p/G)])})
#   
#   varbenetmat[, k] <- sapply(1:G, function(g) {
#     var(as.numeric(coef(fit1.enet))[-1][((g - 1)*p/G + 1):(g*p/G)])})
#   
#   varbenettruemat[, k] <- sapply(1:G, function(g) {
#     var(as.numeric(coef(fit2.enet))[-1][((g - 1)*p/G + 1):(g*p/G)])})
# }
# 
# rownames(msemat) <- rownames(cormat) <- rownames(aucmat) <- rownames(briermat) <- methods
# rownames(pencormat) <- rownames(penrankmat) <- rownames(penmsemat) <- methods[c(3, 4, 6, 8)]
# rownames(pselmat) <- methods[c(1, 2, 4, 6, 8, 11)]
# rownames(varbetamat) <- rownames(varbenetmat) <- rownames(varbenettruemat) <- 
#   paste("group", c(1:G), sep="")
# results6 <- list(mse=msemat, cor=cormat, auc=aucmat, brier=briermat, pencor=pencormat,
#                  penrank=penrankmat, penmse=penmsemat, psel=pselmat, varbeta=varbetamat, 
#                  varbenet=varbenetmat, varbenettrue=varbenettruemat)
# save(results6, file=paste(path.res, "grEBEN_sim_explore1_res6.Rdata", sep=""))



# ### setting 7
# n <- 200
# p <- 1000
# G <- 10
# sigma2 <- 1
# rho <- 0.7
# pblock <- 10
# bSigma <- matrix(c(rep(rep(c(sigma2, rho), times=c(1, pblock)), times=9), sigma2), ncol=pblock,
#                  nrow=pblock, byrow=TRUE)
# lambda1 <- 1
# lambda2 <- 1
# lambdag <- exp(seq(-2, 2, length.out=G))
# groups <- rep(1:G, p/G)
# partitions1 <- list(groups=CreatePartition(as.factor(groups)))
# partitions2 <- list(groups=groups)
# 
# methods <- c("enet", "enettrue", "grridge", "greben.beta", "greben.mu", "grebentrue.beta", 
#              "grebentrue.mu", "greben0.05.beta", "greben0.05.mu", "ridge", "lasso")
# ntest <- 1000
# nreps <- 50
# msemat <- cormat <- aucmat <- briermat <- matrix(NA, nrow=length(methods), ncol=nreps)
# pencormat  <- penrankmat <- penmsemat <- matrix(NA, nrow=4, ncol=nreps)
# varbetamat <- varbenetmat <- varbenettruemat <- matrix(NA, nrow=G, ncol=nreps)
# pselmat <- matrix(NA, nrow=6, ncol=nreps)
# for(k in 1:nreps) {
#   
#   set.seed(7000 + k)
#   cat(paste("Iteration ", k, "\n"))
#   beta <- as.numeric(sapply(1:length(lambdag), function(g) {
#     renet(p/G, lambda1*sqrt(lambdag[g]), lambda2*lambdag[g], log.p=FALSE)}))
#   x <- do.call(cbind, replicate(p/pblock, rmvnorm(n, rep(0, pblock), bSigma), simplify=FALSE))
#   m <- rep(1, n)
#   y <- rbinom(n, m, as.numeric(exp(x %*% beta)/(1 + exp(x %*% beta))))
#   
#   xtest <- matrix(rnorm(ntest*p/2, rep(0, n*p/2), sigma2), ncol=p/2, nrow=ntest)
#   xtest <- cbind(xtest, rho*xtest[, c((p/2):1)] + 
#                    sqrt(1 - rho^2)*rnorm(ntest*p/2, rep(0, n*p/2), sigma2))
#   mtest <- rep(1, ntest)
#   ytest <- rbinom(ntest, mtest, as.numeric(exp(xtest %*% beta)/(1 + exp(xtest %*% beta))))
#   
#   est1.lambda <- cv.pen(x, y, intercept=TRUE)
#   est2.lambda <- cv.glmnet(x, y, family="binomial", alpha=0.05, standardize=FALSE)
#   est3.lambda <- cv.glmnet(x, y, family="binomial", alpha=1, standardize=FALSE)
#   est4.lambda <- cv.glmnet(x, y, family="binomial", alpha=0, standardize=FALSE)
#   
#   fit1.enet <- glmnet(x, y, family="binomial", alpha=est1.lambda$alpha[which.min(est1.lambda$cvll)], 
#                       lambda=est1.lambda$lambda[which.min(est1.lambda$cvll)], standardize=FALSE)
#   fit2.enet <- glmnet(x, y, family="binomial", alpha=lambda1/(lambda1 + 2*lambda2), 
#                       lambda=(0.5*lambda1 + lambda2)/n, standardize=FALSE)
#   fit1.grridge <- grridge(t(x), y, partitions1)
#   fit1.greben <- grEBEN(x, y, m, partitions=partitions2, 
#                         lambda1=est1.lambda$lambda1bayes[which.min(est1.lambda$cvll)],
#                         lambda2=est1.lambda$lambda2bayes[which.min(est1.lambda$cvll)])
#   fit2.greben <- grEBEN(x, y, m, partitions=partitions2, lambda1=lambda2, lambda2=lambda2)
#   fit3.greben <- grEBEN(x, y, m, partitions=partitions2, 
#                         lambda1=2*n*est2.lambda$lambda.min*0.05,
#                         lambda2=n*est2.lambda$lambda.min*0.95)
#   fit1.ridge <- glmnet(x, y, family="binomial", alpha=0, lambda=est4.lambda$lambda.min, 
#                        standardize=FALSE)
#   fit1.lasso <- glmnet(x, y, family="binomial", alpha=1, lambda=est3.lambda$lambda.min, 
#                        standardize=FALSE)
#   
#   msemat[, k] <- c(mean((c(0, beta) - as.numeric(coef(fit1.enet)))^2),
#                    mean((c(0, beta) - as.numeric(coef(fit2.enet)))^2),
#                    mean((c(0, beta) - c(predict.grridge(fit1.grridge, t(matrix(0, ncol=p)))[, 2], 
#                                         fit1.grridge$betas))^2),
#                    mean((c(0, beta) - as.numeric(fit1.greben$beta))^2),
#                    mean((c(0, beta) - fit1.greben$mu)^2),
#                    mean((c(0, beta) - as.numeric(fit2.greben$beta))^2),
#                    mean((c(0, beta) - fit2.greben$mu)^2),
#                    mean((c(0, beta) - as.numeric(fit3.greben$beta))^2),
#                    mean((c(0, beta) - fit3.greben$mu)^2),
#                    mean((c(0, beta) - as.numeric(coef(fit1.ridge)))^2),
#                    mean((c(0, beta) - as.numeric(coef(fit1.lasso)))^2))
#   
#   cormat[, k] <- c(cor(c(0, beta), as.numeric(coef(fit1.enet))),
#                    cor(c(0, beta), as.numeric(coef(fit2.enet))),
#                    cor(c(0, beta), c(predict.grridge(fit1.grridge, t(matrix(0, ncol=p)))[, 2], 
#                                      fit1.grridge$betas)),
#                    cor(c(0, beta), as.numeric(fit1.greben$beta)),
#                    cor(c(0, beta), fit1.greben$mu),
#                    cor(c(0, beta), as.numeric(fit2.greben$beta)),
#                    cor(c(0, beta), fit2.greben$mu),
#                    cor(c(0, beta), as.numeric(fit3.greben$beta)),
#                    cor(c(0, beta), fit3.greben$mu),
#                    cor(c(0, beta), as.numeric(coef(fit1.ridge))),
#                    cor(c(0, beta), as.numeric(coef(fit1.lasso))))
#   
#   predmat <- rbind(as.numeric(predict(fit1.enet, xtest, type="response")),
#                    as.numeric(predict(fit2.enet, xtest, type="response")),
#                    predict.grridge(fit1.grridge, t(xtest))[, 2],
#                    predict.grEBEN(fit1.greben, xtest, type="penalized"),
#                    predict.grEBEN(fit1.greben, xtest, type="VB"),
#                    predict.grEBEN(fit2.greben, xtest, type="penalized"),
#                    predict.grEBEN(fit2.greben, xtest, type="VB"),
#                    predict.grEBEN(fit3.greben, xtest, type="penalized"),
#                    predict.grEBEN(fit3.greben, xtest, type="VB"),
#                    as.numeric(predict(fit1.ridge, xtest, type="response")),
#                    as.numeric(predict(fit1.lasso, xtest, type="response")))
#   
#   aucmat[, k] <- sapply(1:nrow(aucmat), function(method) {
#     pROC::roc(ytest, predmat[method, ])$auc})
#   
#   briermat[, k] <- sapply(1:nrow(aucmat), function(method) {
#     mean((ytest - predmat[method, ])^2)})
#   
#   penmat <- cbind(lambdag, fit1.grridge$lambdamults$groups, 
#                   fit1.greben$lambdag$groups[, fit1.greben$nouteriter + 1],
#                   fit2.greben$lambdag$groups[, fit2.greben$nouteriter + 1],
#                   fit3.greben$lambdag$groups[, fit3.greben$nouteriter + 1])
#   pencormat[, k] <- cor(penmat)[-1, 1]
#   penrankmat[, k] <- cor(penmat, method="spearman")[-1, 1]
#   penmsemat[, k] <- apply(penmat[, -1], 2, function(x) {mean((x - penmat[, 1])^2)})
#   
#   pselmat[, k] <- c(fit1.enet$df, fit2.enet$df, sum(fit1.greben$beta[-1]!=0),
#                     sum(fit2.greben$beta[-1]!=0), sum(fit3.greben$beta[-1]!=0), fit1.lasso$df)
#   
#   varbetamat[, k] <- sapply(1:G, function(g) {var(beta[((g - 1)*p/G + 1):(g*p/G)])})
#   
#   varbenetmat[, k] <- sapply(1:G, function(g) {
#     var(as.numeric(coef(fit1.enet))[-1][((g - 1)*p/G + 1):(g*p/G)])})
#   
#   varbenettruemat[, k] <- sapply(1:G, function(g) {
#     var(as.numeric(coef(fit2.enet))[-1][((g - 1)*p/G + 1):(g*p/G)])})
# }
# 
# rownames(msemat) <- rownames(cormat) <- rownames(aucmat) <- rownames(briermat) <- methods
# rownames(pencormat) <- rownames(penrankmat) <- rownames(penmsemat) <- methods[c(3, 4, 6, 8)]
# rownames(pselmat) <- methods[c(1, 2, 4, 6, 8, 11)]
# rownames(varbetamat) <- rownames(varbenetmat) <- rownames(varbenettruemat) <- 
#   paste("group", c(1:G), sep="")
# results7 <- list(mse=msemat, cor=cormat, auc=aucmat, brier=briermat, pencor=pencormat,
#                  penrank=penrankmat, penmse=penmsemat, psel=pselmat, varbeta=varbetamat, 
#                  varbenet=varbenetmat, varbenettrue=varbenettruemat)
# save(results7, file=paste(path.res, "grEBEN_sim_explore1_res7.Rdata", sep=""))


  
### setting 9
n <- 200
p <- 1000
G <- 4
sigma2 <- 1
rho <- 0.7
pblock <- 10
bSigma <- matrix(c(rep(rep(c(sigma2, rho), times=c(1, pblock)), times=9), sigma2), ncol=pblock,
                 nrow=pblock, byrow=TRUE)
lambda1 <- 1
lambda2 <- 1
lambdag <- exp(seq(-1.5, 1.5, length.out=G))
groups <- rep(1:G, p/G)
partitions1 <- list(groups=CreatePartition(as.factor(groups)))
partitions2 <- list(groups=groups)

methods <- c("enet", "enettrue", "grridge", "greben.beta", "greben.mu", "grebentrue.beta", 
             "grebentrue.mu", "greben0.05.beta", "greben0.05.mu", "ridge", "lasso")
ntest <- 1000
nreps <- 50
msemat <- cormat <- aucmat <- briermat <- matrix(NA, nrow=length(methods), ncol=nreps)
pencormat  <- penrankmat <- penmsemat <- matrix(NA, nrow=4, ncol=nreps)
varbetamat <- varbenetmat <- varbenettruemat <- matrix(NA, nrow=G, ncol=nreps)
pselmat <- matrix(NA, nrow=6, ncol=nreps)
for(k in 1:nreps) {
  
  set.seed(9000 + k)
  cat(paste("Iteration ", k, "\n"))
  beta <- as.numeric(sapply(1:length(lambdag), function(g) {
    renet(p/G, lambda1*sqrt(lambdag[g]), lambda2*lambdag[g], log.p=FALSE)}))
  x <- do.call(cbind, replicate(p/pblock, rmvnorm(n, rep(0, pblock), bSigma), simplify=FALSE))
  m <- rep(1, n)
  y <- rbinom(n, m, as.numeric(exp(x %*% beta)/(1 + exp(x %*% beta))))
  
  xtest <- matrix(rnorm(ntest*p/2, rep(0, n*p/2), sigma2), ncol=p/2, nrow=ntest)
  xtest <- cbind(xtest, rho*xtest[, c((p/2):1)] + 
                   sqrt(1 - rho^2)*rnorm(ntest*p/2, rep(0, n*p/2), sigma2))
  mtest <- rep(1, ntest)
  ytest <- rbinom(ntest, mtest, as.numeric(exp(xtest %*% beta)/(1 + exp(xtest %*% beta))))
  
  est1.lambda <- cv.pen(x, y, intercept=TRUE)
  est2.lambda <- cv.glmnet(x, y, family="binomial", alpha=0.05, standardize=FALSE)
  est3.lambda <- cv.glmnet(x, y, family="binomial", alpha=1, standardize=FALSE)
  est4.lambda <- cv.glmnet(x, y, family="binomial", alpha=0, standardize=FALSE)
  
  fit1.enet <- glmnet(x, y, family="binomial", alpha=est1.lambda$alpha[which.min(est1.lambda$cvll)], 
                      lambda=est1.lambda$lambda[which.min(est1.lambda$cvll)], standardize=FALSE)
  fit2.enet <- glmnet(x, y, family="binomial", alpha=lambda1/(lambda1 + 2*lambda2), 
                      lambda=(0.5*lambda1 + lambda2)/n, standardize=FALSE)
  fit1.grridge <- grridge(t(x), y, partitions1)
  fit1.greben <- grEBEN(x, y, m, partitions=partitions2, 
                        lambda1=est1.lambda$lambda1bayes[which.min(est1.lambda$cvll)],
                        lambda2=est1.lambda$lambda2bayes[which.min(est1.lambda$cvll)])
  fit2.greben <- grEBEN(x, y, m, partitions=partitions2, lambda1=lambda2, lambda2=lambda2)
  fit3.greben <- grEBEN(x, y, m, partitions=partitions2, 
                        lambda1=2*n*est2.lambda$lambda.min*0.05,
                        lambda2=n*est2.lambda$lambda.min*0.95)
  fit1.ridge <- glmnet(x, y, family="binomial", alpha=0, lambda=est4.lambda$lambda.min, 
                       standardize=FALSE)
  fit1.lasso <- glmnet(x, y, family="binomial", alpha=1, lambda=est3.lambda$lambda.min, 
                       standardize=FALSE)
  
  msemat[, k] <- c(mean((c(0, beta) - as.numeric(coef(fit1.enet)))^2),
                   mean((c(0, beta) - as.numeric(coef(fit2.enet)))^2),
                   mean((c(0, beta) - c(predict.grridge(fit1.grridge, t(matrix(0, ncol=p)))[, 2], 
                                        fit1.grridge$betas))^2),
                   mean((c(0, beta) - as.numeric(fit1.greben$beta))^2),
                   mean((c(0, beta) - fit1.greben$mu)^2),
                   mean((c(0, beta) - as.numeric(fit2.greben$beta))^2),
                   mean((c(0, beta) - fit2.greben$mu)^2),
                   mean((c(0, beta) - as.numeric(fit3.greben$beta))^2),
                   mean((c(0, beta) - fit3.greben$mu)^2),
                   mean((c(0, beta) - as.numeric(coef(fit1.ridge)))^2),
                   mean((c(0, beta) - as.numeric(coef(fit1.lasso)))^2))
  
  cormat[, k] <- c(cor(c(0, beta), as.numeric(coef(fit1.enet))),
                   cor(c(0, beta), as.numeric(coef(fit2.enet))),
                   cor(c(0, beta), c(predict.grridge(fit1.grridge, t(matrix(0, ncol=p)))[, 2], 
                                     fit1.grridge$betas)),
                   cor(c(0, beta), as.numeric(fit1.greben$beta)),
                   cor(c(0, beta), fit1.greben$mu),
                   cor(c(0, beta), as.numeric(fit2.greben$beta)),
                   cor(c(0, beta), fit2.greben$mu),
                   cor(c(0, beta), as.numeric(fit3.greben$beta)),
                   cor(c(0, beta), fit3.greben$mu),
                   cor(c(0, beta), as.numeric(coef(fit1.ridge))),
                   cor(c(0, beta), as.numeric(coef(fit1.lasso))))
  
  predmat <- rbind(as.numeric(predict(fit1.enet, xtest, type="response")),
                   as.numeric(predict(fit2.enet, xtest, type="response")),
                   predict.grridge(fit1.grridge, t(xtest))[, 2],
                   predict.grEBEN(fit1.greben, xtest, type="penalized"),
                   predict.grEBEN(fit1.greben, xtest, type="VB"),
                   predict.grEBEN(fit2.greben, xtest, type="penalized"),
                   predict.grEBEN(fit2.greben, xtest, type="VB"),
                   predict.grEBEN(fit3.greben, xtest, type="penalized"),
                   predict.grEBEN(fit3.greben, xtest, type="VB"),
                   as.numeric(predict(fit1.ridge, xtest, type="response")),
                   as.numeric(predict(fit1.lasso, xtest, type="response")))
  
  aucmat[, k] <- sapply(1:nrow(aucmat), function(method) {
    pROC::roc(ytest, predmat[method, ])$auc})
  
  briermat[, k] <- sapply(1:nrow(aucmat), function(method) {
    mean((ytest - predmat[method, ])^2)})
  
  penmat <- cbind(lambdag, fit1.grridge$lambdamults$groups, 
                  fit1.greben$lambdag$groups[, fit1.greben$nouteriter + 1],
                  fit2.greben$lambdag$groups[, fit2.greben$nouteriter + 1],
                  fit3.greben$lambdag$groups[, fit3.greben$nouteriter + 1])
  pencormat[, k] <- cor(penmat)[-1, 1]
  penrankmat[, k] <- cor(penmat, method="spearman")[-1, 1]
  penmsemat[, k] <- apply(penmat[, -1], 2, function(x) {mean((x - penmat[, 1])^2)})
  
  pselmat[, k] <- c(fit1.enet$df, fit2.enet$df, sum(fit1.greben$beta[-1]!=0),
                    sum(fit2.greben$beta[-1]!=0), sum(fit3.greben$beta[-1]!=0), fit1.lasso$df)
  
  varbetamat[, k] <- sapply(1:G, function(g) {var(beta[((g - 1)*p/G + 1):(g*p/G)])})
  
  varbenetmat[, k] <- sapply(1:G, function(g) {
    var(as.numeric(coef(fit1.enet))[-1][((g - 1)*p/G + 1):(g*p/G)])})
  
  varbenettruemat[, k] <- sapply(1:G, function(g) {
    var(as.numeric(coef(fit2.enet))[-1][((g - 1)*p/G + 1):(g*p/G)])})
}

rownames(msemat) <- rownames(cormat) <- rownames(aucmat) <- rownames(briermat) <- methods
rownames(pencormat) <- rownames(penrankmat) <- rownames(penmsemat) <- methods[c(3, 4, 6, 8)]
rownames(pselmat) <- methods[c(1, 2, 4, 6, 8, 11)]
rownames(varbetamat) <- rownames(varbenetmat) <- rownames(varbenettruemat) <- 
  paste("group", c(1:G), sep="")
results9 <- list(mse=msemat, cor=cormat, auc=aucmat, brier=briermat, pencor=pencormat,
                 penrank=penrankmat, penmse=penmsemat, psel=pselmat, varbeta=varbetamat, 
                 varbenet=varbenetmat, varbenettrue=varbenettruemat)
save(results9, file=paste(path.res, "grEBEN_sim_explore1_res9.Rdata", sep=""))



  
# ### graphs
# # setting 1
# load(paste(path.res, "grEBEN_sim_explore1_res1.Rdata", sep=""))
# opar <- par()
# leg.names <- c("enet", "enet with true penalties", "GRridge", "frequentist grEBEN",
#                "Bayesian grEBEN", "frequentist grEBEN with true penalties", 
#                "Bayesian grEBEN with true penalties", expression(paste("frequentist grEBEN with ", alpha==0.05)),
#                expression(paste("Bayesian grEBEN with ", alpha==0.05)), "ridge", "lasso")
# png(paste(path.graph, "grEBEN_sim_explore1_boxplot1.png", sep=""), width=1200, height=720, res=120)
# layout(matrix(c(1, 2, 5, 3, 4, 5), 2, 3, byrow=TRUE))
# boxplot(t(results1$mse), las=2, main="a)", ylab="MSE", col=rainbow(nrow(results1$mse)), names=NA,
#         xaxt="n")
# boxplot(t(results1$cor), las=2, main="b)", ylab=expression(paste("cor(", beta, ",", hat(beta), ")")), 
#         col=rainbow(nrow(results1$cor)), names=NA,
#         xaxt="n")
# boxplot(t(results1$auc), las=2, main="c)", ylab="AUC", col=rainbow(nrow(results1$auc)), names=NA,
#         xaxt="n")
# boxplot(t(results1$brier), las=2, main="d)", ylab="Brier score", col=rainbow(nrow(results1$brier)), 
#         names=NA, xaxt="n")
# plot.new()
# legend("topleft", legend=leg.names, col=rainbow(nrow(results1$brier)), pch=15)
# par(opar)
# dev.off()



# # setting 2
# load(paste(path.res, "grEBEN_sim_explore1_res2.Rdata", sep=""))
# opar <- par()
# leg.names <- c("enet", "enet with true penalties", "GRridge", "frequentist grEBEN",
#                "Bayesian grEBEN", "frequentist grEBEN with true penalties",
#                "Bayesian grEBEN with true penalties", expression(paste("frequentist grEBEN with ", alpha==0.05)),
#                expression(paste("Bayesian grEBEN with ", alpha==0.05)), "ridge", "lasso")
# png(paste(path.graph, "grEBEN_sim_explore1_boxplot2.1.png", sep=""), width=1200, height=720, res=120)
# layout(matrix(c(1, 2, 5, 3, 4, 5), 2, 3, byrow=TRUE))
# boxplot(t(results2$mse), las=2, main="a)", ylab="MSE", col=rainbow(nrow(results2$mse)), names=NA,
#         xaxt="n")
# boxplot(t(results2$cor), las=2, main="b)", ylab=expression(paste("cor(", beta, ",", hat(beta), ")")),
#         col=rainbow(nrow(results2$cor)), names=NA,
#         xaxt="n")
# boxplot(t(results2$auc), las=2, main="c)", ylab="AUC", col=rainbow(nrow(results2$auc)), names=NA,
#         xaxt="n")
# boxplot(t(results2$brier), las=2, main="d)", ylab="Brier score", col=rainbow(nrow(results2$brier)),
#         names=NA, xaxt="n")
# plot.new()
# legend("topleft", legend=leg.names, col=rainbow(nrow(results2$brier)), pch=15)
# par(opar)
# dev.off()
# 
# leg.names <- c("GRridge", "grEBEN", "grEBEN with true penalties", 
#                expression(paste("grEBEN with ", alpha==0.05)))
# png(paste(path.graph, "grEBEN_sim_explore1_boxplot2.2.png", sep=""), width=1200, height=720, res=120)
# layout(matrix(c(1, 2, 3, 4), 2, 2, byrow=TRUE))
# boxplot(t(results2$pencor), las=2, main="d)", ylab=expression(paste("cor(", lambda[g], ",", hat(lambda)[g], ")")), 
#         col=rainbow(nrow(results2$pencor)), names=NA, xaxt="n")
# boxplot(t(results2$penrank), las=2, main="d)", ylab=expression(paste("Spearman cor(", lambda[g], ",", hat(lambda)[g], ")")), 
#         col=rainbow(nrow(results2$penrank)), names=NA, xaxt="n")
# boxplot(ifelse(t(results2$penmse) >= 10, NA, t(results2$penmse)), las=2, main="d)", ylab=expression(paste("MSE(", hat(lambda)[g], ")")), 
#         col=rainbow(nrow(results2$penmse)), names=NA, xaxt="n")
# plot.new()
# legend("topleft", legend=leg.names, col=rainbow(nrow(results2$pencor)), pch=15)
# par(opar)
# dev.off()



# # setting 3
# load(paste(path.res, "grEBEN_sim_explore1_res3.Rdata", sep=""))
# opar <- par()
# leg.names <- c("enet", "enet with true penalties", "GRridge", "frequentist grEBEN",
#                "Bayesian grEBEN", "frequentist grEBEN with true penalties",
#                "Bayesian grEBEN with true penalties", expression(paste("frequentist grEBEN with ", alpha==0.05)),
#                expression(paste("Bayesian grEBEN with ", alpha==0.05)), "ridge", "lasso")
# png(paste(path.graph, "grEBEN_sim_explore1_boxplot3.1.png", sep=""), width=1200, height=720, res=120)
# layout(matrix(c(1, 2, 5, 3, 4, 5), 2, 3, byrow=TRUE))
# boxplot(t(results3$mse), las=2, main="a)", ylab="MSE", col=rainbow(nrow(results3$mse)), names=NA,
#         xaxt="n")
# boxplot(t(results3$cor), las=2, main="b)", ylab=expression(paste("cor(", beta, ",", hat(beta), ")")),
#         col=rainbow(nrow(results3$cor)), names=NA,
#         xaxt="n")
# boxplot(t(results3$auc), las=2, main="c)", ylab="AUC", col=rainbow(nrow(results3$auc)), names=NA,
#         xaxt="n")
# boxplot(t(results3$brier), las=2, main="d)", ylab="Brier score", col=rainbow(nrow(results3$brier)),
#         names=NA, xaxt="n")
# plot.new()
# legend("topleft", legend=leg.names, col=rainbow(nrow(results3$brier)), pch=15)
# par(opar)
# dev.off()
# 
# leg.names <- c("GRridge", "grEBEN", "grEBEN with true penalties", 
#                expression(paste("grEBEN with ", alpha==0.05)))
# png(paste(path.graph, "grEBEN_sim_explore1_boxplot3.2.png", sep=""), width=1200, height=720, res=120)
# layout(matrix(c(1, 2, 3, 4), 2, 2, byrow=TRUE))
# boxplot(t(results3$pencor), las=2, main="d)", ylab=expression(paste("cor(", lambda[g], ",", hat(lambda)[g], ")")), 
#         col=rainbow(nrow(results3$pencor)), names=NA, xaxt="n")
# boxplot(t(results3$penrank), las=2, main="d)", ylab=expression(paste("Spearman cor(", lambda[g], ",", hat(lambda)[g], ")")), 
#         col=rainbow(nrow(results3$penrank)), names=NA, xaxt="n")
# boxplot(t(results3$penrank), las=2, main="d)", ylab=expression(paste("Spearman cor(", lambda[g], ",", hat(lambda)[g], ")")), 
#         col=rainbow(nrow(results3$penrank)), names=NA, xaxt="n")
# plot.new()
# legend("topleft", legend=leg.names, col=rainbow(nrow(results3$pencor)), pch=15)
# par(opar)
# dev.off()




# # setting 4
# load(paste(path.res, "grEBEN_sim_explore1_res4.Rdata", sep=""))
# opar <- par()
# leg.names <- c("enet", "enet with true penalties", "GRridge", "frequentist grEBEN",
#                "Bayesian grEBEN", "frequentist grEBEN with true penalties",
#                "Bayesian grEBEN with true penalties", expression(paste("frequentist grEBEN with ", alpha==0.05)),
#                expression(paste("Bayesian grEBEN with ", alpha==0.05)), "ridge", "lasso")
# png(paste(path.graph, "grEBEN_sim_explore1_boxplot4.1.png", sep=""), width=1200, height=720, res=120)
# layout(matrix(c(1, 2, 5, 3, 4, 5), 2, 3, byrow=TRUE))
# boxplot(t(results4$mse), las=2, main="a)", ylab="MSE", col=rainbow(nrow(results4$mse)), names=NA,
#         xaxt="n")
# boxplot(t(results4$cor), las=2, main="b)", ylab=expression(paste("cor(", beta, ",", hat(beta), ")")),
#         col=rainbow(nrow(results4$cor)), names=NA,
#         xaxt="n")
# boxplot(t(results4$auc), las=2, main="c)", ylab="AUC", col=rainbow(nrow(results4$auc)), names=NA,
#         xaxt="n")
# boxplot(t(results4$brier), las=2, main="d)", ylab="Brier score", col=rainbow(nrow(results4$brier)),
#         names=NA, xaxt="n")
# plot.new()
# legend("topleft", legend=leg.names, col=rainbow(nrow(results4$brier)), pch=15)
# par(opar)
# dev.off()
# 
# leg.names <- c("GRridge", "grEBEN", "grEBEN with true penalties", 
#                expression(paste("grEBEN with ", alpha==0.05)))
# png(paste(path.graph, "grEBEN_sim_explore1_boxplot4.2.png", sep=""), width=1200, height=720, res=120)
# layout(matrix(c(1, 2, 3, 4), 2, 2, byrow=TRUE))
# boxplot(t(results4$pencor), las=2, main="d)", ylab=expression(paste("cor(", lambda[g], ",", hat(lambda)[g], ")")), 
#         col=rainbow(nrow(results4$pencor)), names=NA, xaxt="n")
# boxplot(t(results4$penrank), las=2, main="d)", ylab=expression(paste("Spearman cor(", lambda[g], ",", hat(lambda)[g], ")")), 
#         col=rainbow(nrow(results4$penrank)), names=NA, xaxt="n")
# boxplot(t(results4$penmse), las=2, main="d)", ylab=expression(paste("MSE(", hat(lambda)[g], ")")), 
#         col=rainbow(nrow(results4$penmse)), names=NA, xaxt="n")
# plot.new()
# legend("topleft", legend=leg.names, col=rainbow(nrow(results4$pencor)), pch=15)
# par(opar)
# dev.off()



# # setting 5
# load(paste(path.res, "grEBEN_sim_explore1_res5.Rdata", sep=""))
# opar <- par()
# leg.names <- c("enet", "enet with true penalties", "GRridge", "frequentist grEBEN",
#                "Bayesian grEBEN", "frequentist grEBEN with true penalties",
#                "Bayesian grEBEN with true penalties", expression(paste("frequentist grEBEN with ", alpha==0.05)),
#                expression(paste("Bayesian grEBEN with ", alpha==0.05)), "ridge", "lasso")
# png(paste(path.graph, "grEBEN_sim_explore1_boxplot5.1.png", sep=""), width=1200, height=720, res=120)
# layout(matrix(c(1, 2, 5, 3, 4, 5), 2, 3, byrow=TRUE))
# boxplot(t(results5$mse), las=2, main="a)", ylab="MSE", col=rainbow(nrow(results5$mse)), names=NA,
#         xaxt="n")
# boxplot(t(results5$cor), las=2, main="b)", ylab=expression(paste("cor(", beta, ",", hat(beta), ")")),
#         col=rainbow(nrow(results5$cor)), names=NA,
#         xaxt="n")
# boxplot(t(results5$auc), las=2, main="c)", ylab="AUC", col=rainbow(nrow(results5$auc)), names=NA,
#         xaxt="n")
# boxplot(t(results5$brier), las=2, main="d)", ylab="Brier score", col=rainbow(nrow(results5$brier)),
#         names=NA, xaxt="n")
# plot.new()
# legend("topleft", legend=leg.names, col=rainbow(nrow(results5$brier)), pch=15)
# par(opar)
# dev.off()
# 
# leg.names <- c("GRridge", "grEBEN", "grEBEN with true penalties",
#                expression(paste("grEBEN with ", alpha==0.05)))
# png(paste(path.graph, "grEBEN_sim_explore1_boxplot5.2.png", sep=""), width=1200, height=720, res=120)
# layout(matrix(c(1, 2, 3, 4), 2, 2, byrow=TRUE))
# boxplot(t(results5$pencor), las=2, main="d)", ylab=expression(paste("cor(", lambda[g], ",", hat(lambda)[g], ")")),
#         col=rainbow(nrow(results5$pencor)), names=NA, xaxt="n")
# boxplot(t(results5$penrank), las=2, main="d)", ylab=expression(paste("Spearman cor(", lambda[g], ",", hat(lambda)[g], ")")),
#         col=rainbow(nrow(results5$penrank)), names=NA, xaxt="n")
# boxplot(t(results5$penmse), las=2, main="d)", ylab=expression(paste("MSE(", hat(lambda)[g], ")")),
#         col=rainbow(nrow(results5$penmse)), names=NA, xaxt="n")
# plot.new()
# legend("topleft", legend=leg.names, col=rainbow(nrow(results5$pencor)), pch=15)
# par(opar)
# dev.off()



# # setting 6
# load(paste(path.res, "grEBEN_sim_explore1_res6.Rdata", sep=""))
# opar <- par()
# leg.names <- c("enet", "enet with true penalties", "GRridge", "frequentist grEBEN",
#                "Bayesian grEBEN", "frequentist grEBEN with true penalties",
#                "Bayesian grEBEN with true penalties", expression(paste("frequentist grEBEN with ", alpha==0.05)),
#                expression(paste("Bayesian grEBEN with ", alpha==0.05)), "ridge", "lasso")
# png(paste(path.graph, "grEBEN_sim_explore1_boxplot6.1.png", sep=""), width=1200, height=720, res=120)
# layout(matrix(c(1, 2, 5, 3, 4, 5), 2, 3, byrow=TRUE))
# boxplot(t(results6$mse), las=2, main="a)", ylab="MSE", col=rainbow(nrow(results6$mse)), names=NA,
#         xaxt="n")
# boxplot(t(results6$cor), las=2, main="b)", ylab=expression(paste("cor(", beta, ",", hat(beta), ")")),
#         col=rainbow(nrow(results6$cor)), names=NA,
#         xaxt="n")
# boxplot(t(results6$auc), las=2, main="c)", ylab="AUC", col=rainbow(nrow(results6$auc)), names=NA,
#         xaxt="n")
# boxplot(t(results6$brier), las=2, main="d)", ylab="Brier score", col=rainbow(nrow(results6$brier)),
#         names=NA, xaxt="n")
# plot.new()
# legend("topleft", legend=leg.names, col=rainbow(nrow(results6$brier)), pch=15)
# par(opar)
# dev.off()
# 
# leg.names <- c("GRridge", "grEBEN", "grEBEN with true penalties",
#                expression(paste("grEBEN with ", alpha==0.05)))
# png(paste(path.graph, "grEBEN_sim_explore1_boxplot6.2.png", sep=""), width=1200, height=720, res=120)
# layout(matrix(c(1, 2, 3, 4), 2, 2, byrow=TRUE))
# boxplot(t(results6$pencor), las=2, main="d)", ylab=expression(paste("cor(", lambda[g], ",", hat(lambda)[g], ")")),
#         col=rainbow(nrow(results6$pencor)), names=NA, xaxt="n")
# boxplot(t(results6$penrank), las=2, main="d)", ylab=expression(paste("Spearman cor(", lambda[g], ",", hat(lambda)[g], ")")),
#         col=rainbow(nrow(results6$penrank)), names=NA, xaxt="n")
# boxplot(t(results6$penmse), las=2, main="d)", ylab=expression(paste("MSE(", hat(lambda)[g], ")")),
#         col=rainbow(nrow(results6$penmse)), names=NA, xaxt="n")
# plot.new()
# legend("topleft", legend=leg.names, col=rainbow(nrow(results6$pencor)), pch=15)
# par(opar)
# dev.off()





# # setting 7
# load(paste(path.res, "grEBEN_sim_explore1_res7.Rdata", sep=""))
# opar <- par()
# leg.names <- c("enet", "enet with true penalties", "GRridge", "frequentist grEBEN",
#                "Bayesian grEBEN", "frequentist grEBEN with true penalties",
#                "Bayesian grEBEN with true penalties", expression(paste("frequentist grEBEN with ", alpha==0.05)),
#                expression(paste("Bayesian grEBEN with ", alpha==0.05)), "ridge", "lasso")
# png(paste(path.graph, "grEBEN_sim_explore1_boxplot7.1.png", sep=""), width=1200, height=720, res=120)
# layout(matrix(c(1, 2, 5, 3, 4, 5), 2, 3, byrow=TRUE))
# boxplot(t(results7$mse), las=2, main="a)", ylab="MSE", col=rainbow(nrow(results7$mse)), names=NA,
#         xaxt="n")
# boxplot(t(results7$cor), las=2, main="b)", ylab=expression(paste("cor(", beta, ",", hat(beta), ")")),
#         col=rainbow(nrow(results7$cor)), names=NA,
#         xaxt="n")
# boxplot(t(results7$auc), las=2, main="c)", ylab="AUC", col=rainbow(nrow(results7$auc)), names=NA,
#         xaxt="n")
# boxplot(t(results7$brier), las=2, main="d)", ylab="Brier score", col=rainbow(nrow(results7$brier)),
#         names=NA, xaxt="n")
# plot.new()
# legend("topleft", legend=leg.names, col=rainbow(nrow(results7$brier)), pch=15)
# par(opar)
# dev.off()
# 
# leg.names <- c("GRridge", "grEBEN", "grEBEN with true penalties",
#                expression(paste("grEBEN with ", alpha==0.05)))
# png(paste(path.graph, "grEBEN_sim_explore1_boxplot7.2.png", sep=""), width=1200, height=720, res=120)
# layout(matrix(c(1, 2, 3, 4), 2, 2, byrow=TRUE))
# boxplot(t(results7$pencor), las=2, main="d)", ylab=expression(paste("cor(", lambda[g], ",", hat(lambda)[g], ")")),
#         col=rainbow(nrow(results7$pencor)), names=NA, xaxt="n")
# boxplot(t(results7$penrank), las=2, main="d)", ylab=expression(paste("Spearman cor(", lambda[g], ",", hat(lambda)[g], ")")),
#         col=rainbow(nrow(results7$penrank)), names=NA, xaxt="n")
# boxplot(t(results7$penmse), las=2, main="d)", ylab=expression(paste("MSE(", hat(lambda)[g], ")")),
#         col=rainbow(nrow(results7$penmse)), names=NA, xaxt="n")
# plot.new()
# legend("topleft", legend=leg.names, col=rainbow(nrow(results7$pencor)), pch=15)
# par(opar)
# dev.off()





# # setting 8
# load(paste(path.res, "grEBEN_sim_explore1_res8.Rdata", sep=""))
# opar <- par()
# leg.names <- c("enet", "enet with true penalties", "GRridge", "frequentist grEBEN",
#                "Bayesian grEBEN", "frequentist grEBEN with true penalties",
#                "Bayesian grEBEN with true penalties", expression(paste("frequentist grEBEN with ", alpha==0.05)),
#                expression(paste("Bayesian grEBEN with ", alpha==0.05)), "ridge", "lasso")
# png(paste(path.graph, "grEBEN_sim_explore1_boxplot8.1.png", sep=""), width=1200, height=720, res=120)
# layout(matrix(c(1, 2, 5, 3, 4, 5), 2, 3, byrow=TRUE))
# boxplot(t(results8$mse), las=2, main="a)", ylab="MSE", col=rainbow(nrow(results8$mse)), names=NA,
#         xaxt="n")
# boxplot(t(results8$cor), las=2, main="b)", ylab=expression(paste("cor(", beta, ",", hat(beta), ")")),
#         col=rainbow(nrow(results8$cor)), names=NA,
#         xaxt="n")
# boxplot(t(results8$auc), las=2, main="c)", ylab="AUC", col=rainbow(nrow(results8$auc)), names=NA,
#         xaxt="n")
# boxplot(t(results8$brier), las=2, main="d)", ylab="Brier score", col=rainbow(nrow(results8$brier)),
#         names=NA, xaxt="n")
# plot.new()
# legend("topleft", legend=leg.names, col=rainbow(nrow(results8$brier)), pch=15)
# par(opar)
# dev.off()
# 
# leg.names <- c("GRridge", "grEBEN", "grEBEN with true penalties",
#                expression(paste("grEBEN with ", alpha==0.05)))
# png(paste(path.graph, "grEBEN_sim_explore1_boxplot8.2.png", sep=""), width=1200, height=720, res=120)
# layout(matrix(c(1, 2, 3, 4), 2, 2, byrow=TRUE))
# boxplot(t(results8$pencor), las=2, main="d)", ylab=expression(paste("cor(", lambda[g], ",", hat(lambda)[g], ")")),
#         col=rainbow(nrow(results8$pencor)), names=NA, xaxt="n")
# boxplot(t(results8$penrank), las=2, main="d)", ylab=expression(paste("Spearman cor(", lambda[g], ",", hat(lambda)[g], ")")),
#         col=rainbow(nrow(results8$penrank)), names=NA, xaxt="n")
# boxplot(t(results8$penmse), las=2, main="d)", ylab=expression(paste("MSE(", hat(lambda)[g], ")")),
#         col=rainbow(nrow(results8$penmse)), names=NA, xaxt="n")
# plot.new()
# legend("topleft", legend=leg.names, col=rainbow(nrow(results8$pencor)), pch=15)
# par(opar)
# dev.off()