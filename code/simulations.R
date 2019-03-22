#!/usr/bin/env Rscript

### installation of gren if updated on github
if(substr(system('git log -n 1 --format="%h %aN %s %ad"', intern=TRUE), 1, 7)!=
   substr(packageDescription("gren")$GithubSHA1, 1, 7)) {
  if(!("devtools" %in% installed.packages())) {
    install.packages("devtools")
  }
  library(devtools)
  install_github("magnusmunch/gren/rpackage", local=FALSE,
                 auth_token=Sys.getenv("GITHUB_PAT"))
}

### parallelisation
parallel <- TRUE

### libraries
library(gren)
library(GRridge)
library(grpreg)
library(SGL)
library(foreach)
library(doParallel)
library(irr)
library(mvtnorm)

### loading functions
source("code/functions.R")

################################# simulation 1 #################################
n <- 100
ntest <- 1000
p <- 1000
G <- 5
rho <- 0.5
rg <- 50
beta.mean <- 0.07
f <- 1.6

Sigma <- diag(p)
for(i in 1:p) {
  for(j in 1:p) {
    Sigma[i, j] <- rho^abs(i - j)
  }
}
b <- 2*p*beta.mean/sum(rg*f^(c(1:G) - 1))
beta <- numeric(p)

part <- rep(c(1:G), each=p/G)
csel <- 2^c(1:8)
nreps <- 100

### analysis splits in parallel
ncores <- min(detectCores() - 1, nreps)
cluster <- makeForkCluster(ncores)
if(parallel) {
  registerDoParallel(cluster)
} else {
  registerDoSEQ()
}

res <- foreach(k=c(1:nreps), .errorhandling="pass") %dopar% {
  set.seed(2019 + k)
  
  for(g in 1:G) {
    idz <- c(((g-1)*p/G + 1):(g*p/G - rg))
    ida <- c(((g*p/G - rg) + 1):(g*p/G))
    beta[idz] <- 0
    beta[ida] <- runif(rg, 0, b*f^(g-1))
  }
  xtrain <- rmvnorm(n, rep(0, p), Sigma)
  ytrain <- rbinom(n, 1, as.numeric(1/(1 + exp(-xtrain %*% beta))))
  
  xtest <- rmvnorm(ntest, rep(0, p), Sigma)
  ytest <- rbinom(ntest, 1, as.numeric(1/(1 + exp(-xtest %*% beta))))
  
  # fitting models
  fit.gren1 <- gren(xtrain, ytrain, partitions=list(part=part), alpha=0.05, 
                    standardize=TRUE, trace=FALSE, psel=csel)
  fit.gren2 <- gren(xtrain, ytrain, partitions=list(part=part), alpha=0.5, 
                    standardize=TRUE, trace=FALSE, psel=csel)
  fit.gren3 <- gren(xtrain, ytrain, partitions=list(part=part), alpha=0.95, 
                    standardize=TRUE, trace=FALSE, psel=csel)
  
  fit.grridge <- grridge(t(xtrain), ytrain, list(part=split(1:p, part)))
  
  fit.cmcp1 <- sel.grpreg(xtrain, ytrain, part, penalty="cMCP", 
                          family="binomial", alpha=0.05, psel=csel)
  fit.cmcp2 <- sel.grpreg(xtrain, ytrain, part, penalty="cMCP", 
                          family="binomial", 
                          alpha=0.5, psel=csel)
  fit.cmcp3 <- sel.grpreg(xtrain, ytrain, part, penalty="cMCP", 
                          family="binomial", alpha=0.95, psel=csel)
  
  fit.gel1 <- sel.grpreg(xtrain, ytrain, part, penalty="gel", family="binomial", 
                         alpha=0.05, psel=csel)
  fit.gel2 <- sel.grpreg(xtrain, ytrain, part, penalty="gel", family="binomial", 
                         alpha=0.5, psel=csel)
  fit.gel3 <- sel.grpreg(xtrain, ytrain, part, penalty="gel", family="binomial", 
                         alpha=0.95, psel=csel)
  
  pred <- data.frame(ridge=predict.grridge(fit.grridge, t(xtest))[, 1],
                     grridge=predict.grridge(fit.grridge, t(xtest))[, 2],
                     gren1=predict(fit.gren1, xtest, type="groupreg"),
                     gren2=predict(fit.gren2, xtest, type="groupreg"),
                     gren3=predict(fit.gren3, xtest, type="groupreg"),
                     enet1=predict(fit.gren1, xtest, type="regular"),
                     enet2=predict(fit.gren2, xtest, type="regular"),
                     enet3=predict(fit.gren3, xtest, type="regular"),
                     cmcp1=predict(fit.cmcp1, xtest),
                     cmcp2=predict(fit.cmcp2, xtest),
                     cmcp3=predict(fit.cmcp3, xtest),
                     gel1=predict(fit.gel1, xtest),
                     gel2=predict(fit.gel2, xtest),
                     gel3=predict(fit.gel3, xtest))
  
  auc <- apply(pred, 2, function(s) {pROC::auc(ytest, s)})
  
  const <- sum((ytest - mean(ytest))^2)
  briers <- apply(pred, 2, function(s) {1 - sum((ytest - s)^2)/const})
  
  coef <- data.frame(ridge=coef(fit.grridge$predobj$NoGroups),
                      grridge=coef(fit.grridge$predobj$GroupRegul),
                      gren1=as.matrix(coef(fit.gren1, type="groupreg")),
                      gren2=as.matrix(coef(fit.gren2, type="groupreg")),
                      gren3=as.matrix(coef(fit.gren3, type="groupreg")),
                      enet1=as.matrix(coef(fit.gren1, type="regular")),
                      enet2=as.matrix(coef(fit.gren2, type="regular")),
                      enet3=as.matrix(coef(fit.gren3, type="regular")),
                      cmcp1=coef(fit.cmcp1),
                      cmcp2=coef(fit.cmcp2),
                      cmcp3=coef(fit.cmcp3),
                      gel1=coef(fit.gel1),
                      gel2=coef(fit.gel2),
                      gel3=coef(fit.gel3))
  
  mse <- apply(coef, 2, function(s) {mean((s - c(0, beta))^2)})
  
  kappa <- apply(coef, 2, function(s) {kappa2(cbind(s[-1]!=0, beta!=0))$value})
  
  psel <- c(ridge=p, grridge=p, 
            gren1=fit.gren1$freq.model$groupreg$df,
            gren2=fit.gren2$freq.model$groupreg$df,
            gren3=fit.gren3$freq.model$groupreg$df,
            enet1=fit.gren1$freq.model$regular$df,
            enet2=fit.gren2$freq.model$regular$df,
            enet3=fit.gren3$freq.model$regular$df,
            cmcp1=apply(as.matrix(coef(fit.cmcp1)), 2, function(b) {sum(b!=0)}),
            cmcp2=apply(as.matrix(coef(fit.cmcp2)), 2, function(b) {sum(b!=0)}),
            cmcp3=apply(as.matrix(coef(fit.cmcp3)), 2, function(b) {sum(b!=0)}),
            gel1=apply(as.matrix(coef(fit.gel1)), 2, function(b) {sum(b!=0)}),
            gel2=apply(as.matrix(coef(fit.gel2)), 2, function(b) {sum(b!=0)}),
            gel3=apply(as.matrix(coef(fit.gel3)), 2, function(b) {sum(b!=0)}))
  
  mults <- c(grridge=fit.grridge$lambdamults$part,
             gren1=fit.gren1$lambdag$part,
             gren2=fit.gren2$lambdag$part,
             gren3=fit.gren3$lambdag$part)
  
  list(psel=psel, auc=auc, briers=briers, mse=mse, kappa=kappa, mults=mults)
  
}
if(parallel) {stopCluster(cluster)}
save(res, file="results/simulations_res1.Rdata")



################################# simulation 2 #################################
n <- 100
ntest <- 1000
p <- 900
G <- 300
rho <- 0.5
Gactive <- 50

Sigma <- diag(G)
for(i in 1:G) {
  for(j in 1:G) {
    Sigma[i, j] <- rho^abs(i - j)
  }
}
beta.active.mean <- 0.3

nclass <- 3
q <- c(-Inf, qnorm((c(1:nclass)/nclass)[-nclass]), Inf)
beta <- numeric(p)
beta.active <- seq(beta.active.mean*Gactive/(Gactive/2 + 0.5),
                   beta.active.mean/(Gactive/2 + 0.5), length.out=Gactive)
for(g in 1:Gactive) {
  id <- ((g - 1)*nclass + 1):((g - 1)*nclass + nclass - 1)
  beta[id] <- beta.active[g]
}
beta0 <- as.numeric(rep(1/3, p) %*% beta)

part <- rep(c(1:G), each=p/G)
csel <- 2^c(1:8)
nreps <- 100

### analysis splits in parallel
# ncores <- min(detectCores() - 1, nreps)
ncores <- 50
cluster <- makeForkCluster(ncores)
if(parallel) {
  registerDoParallel(cluster)
} else {
  registerDoSEQ()
}

res <- foreach(k=c(1:nreps), .errorhandling="pass") %dopar% {
  
  print(paste("rep", k))
  set.seed(2019 + k)
  
  ztrain <- rmvnorm(n, rep(0, G), Sigma)
  xtrain <- matrix(0, nrow=n, ncol=p)
  for(g in 1:G) {
    id <- c(((g-1)*p/G + 1):(g*p/G))
    for(class in 1:nclass) {
      xtrain[, id[class]] <- (ztrain[, g] >= q[class] & 
                                ztrain[, g] < q[class + 1])
    }
  }
  ytrain <- rbinom(n, 1, as.numeric(1/(1 + exp(-xtrain %*% beta + beta0))))
  
  ztest <- rmvnorm(n, rep(0, G), Sigma)
  xtest <- matrix(0, nrow=ntest, ncol=p)
  for(g in 1:G) {
    id <- c(((g-1)*p/G + 1):(g*p/G))
    for(class in 1:nclass) {
      xtest[, id[class]] <- (ztest[, g] >= q[class] & ztest[, g] < q[class + 1])
    }
  }
  ytest <- rbinom(ntest, 1, as.numeric(1/(1 + exp(-xtest %*% beta + beta0))))
  
  # fitting models
  fit.gren1 <- gren(xtrain, ytrain, partitions=list(part=part), alpha=0.05, 
                    standardize=TRUE, trace=FALSE, psel=csel,
                    control=list(maxit.vb=10, maxit=200))
  fit.gren2 <- gren(xtrain, ytrain, partitions=list(part=part), alpha=0.5, 
                    standardize=TRUE, trace=FALSE, psel=csel,
                    control=list(maxit.vb=10, maxit=200))
  fit.gren3 <- gren(xtrain, ytrain, partitions=list(part=part), alpha=0.95, 
                    standardize=TRUE, trace=FALSE, psel=csel,
                    control=list(maxit.vb=10, maxit=200))
  
  fit.grridge <- grridge(t(xtrain), ytrain, list(part=split(1:p, part)))
  
  fit.cmcp1 <- sel.grpreg(xtrain, ytrain, part, penalty="cMCP", 
                          family="binomial", alpha=0.05, psel=csel)
  fit.cmcp2 <- sel.grpreg(xtrain, ytrain, part, penalty="cMCP", 
                          family="binomial", 
                          alpha=0.5, psel=csel)
  fit.cmcp3 <- sel.grpreg(xtrain, ytrain, part, penalty="cMCP", 
                          family="binomial", alpha=0.95, psel=csel)
  
  fit.gel1 <- sel.grpreg(xtrain, ytrain, part, penalty="gel", family="binomial", 
                         alpha=0.05, psel=csel)
  fit.gel2 <- sel.grpreg(xtrain, ytrain, part, penalty="gel", family="binomial", 
                         alpha=0.5, psel=csel)
  fit.gel3 <- sel.grpreg(xtrain, ytrain, part, penalty="gel", family="binomial", 
                         alpha=0.95, psel=csel)
  
  pred <- data.frame(ridge=predict.grridge(fit.grridge, t(xtest))[, 1],
                     grridge=predict.grridge(fit.grridge, t(xtest))[, 2],
                     gren1=predict(fit.gren1, xtest, type="groupreg"),
                     gren2=predict(fit.gren2, xtest, type="groupreg"),
                     gren3=predict(fit.gren3, xtest, type="groupreg"),
                     enet1=predict(fit.gren1, xtest, type="regular"),
                     enet2=predict(fit.gren2, xtest, type="regular"),
                     enet3=predict(fit.gren3, xtest, type="regular"),
                     cmcp1=predict(fit.cmcp1, xtest),
                     cmcp2=predict(fit.cmcp2, xtest),
                     cmcp3=predict(fit.cmcp3, xtest),
                     gel1=predict(fit.gel1, xtest),
                     gel2=predict(fit.gel2, xtest),
                     gel3=predict(fit.gel3, xtest))
  
  auc <- apply(pred, 2, function(s) {pROC::auc(ytest, s)})
  
  const <- sum((ytest - mean(ytest))^2)
  briers <- apply(pred, 2, function(s) {1 - sum((ytest - s)^2)/const})
  
  coef <- data.frame(ridge=coef(fit.grridge$predobj$NoGroups),
                     grridge=coef(fit.grridge$predobj$GroupRegul),
                     gren1=as.matrix(coef(fit.gren1, type="groupreg")),
                     gren2=as.matrix(coef(fit.gren2, type="groupreg")),
                     gren3=as.matrix(coef(fit.gren3, type="groupreg")),
                     enet1=as.matrix(coef(fit.gren1, type="regular")),
                     enet2=as.matrix(coef(fit.gren2, type="regular")),
                     enet3=as.matrix(coef(fit.gren3, type="regular")),
                     cmcp1=coef(fit.cmcp1),
                     cmcp2=coef(fit.cmcp2),
                     cmcp3=coef(fit.cmcp3),
                     gel1=coef(fit.gel1),
                     gel2=coef(fit.gel2),
                     gel3=coef(fit.gel3))
  
  mse <- apply(coef, 2, function(s) {mean((s - c(0, beta))^2)})
  
  kappa <- apply(coef, 2, function(s) {kappa2(cbind(s[-1]!=0, beta!=0))$value})
  
  psel <- c(ridge=p, grridge=p, 
            gren1=fit.gren1$freq.model$groupreg$df,
            gren2=fit.gren2$freq.model$groupreg$df,
            gren3=fit.gren3$freq.model$groupreg$df,
            enet1=fit.gren1$freq.model$regular$df,
            enet2=fit.gren2$freq.model$regular$df,
            enet3=fit.gren3$freq.model$regular$df,
            cmcp1=apply(as.matrix(coef(fit.cmcp1)), 2, function(b) {sum(b!=0)}),
            cmcp2=apply(as.matrix(coef(fit.cmcp2)), 2, function(b) {sum(b!=0)}),
            cmcp3=apply(as.matrix(coef(fit.cmcp3)), 2, function(b) {sum(b!=0)}),
            gel1=apply(as.matrix(coef(fit.gel1)), 2, function(b) {sum(b!=0)}),
            gel2=apply(as.matrix(coef(fit.gel2)), 2, function(b) {sum(b!=0)}),
            gel3=apply(as.matrix(coef(fit.gel3)), 2, function(b) {sum(b!=0)}))
  
  mults <- c(grridge=fit.grridge$lambdamults$part,
             gren1=fit.gren1$lambdag$part,
             gren2=fit.gren2$lambdag$part,
             gren3=fit.gren3$lambdag$part)
  
  list(psel=psel, auc=auc, briers=briers, mse=mse, kappa=kappa, mults=mults)
  
}
if(parallel) {stopCluster(cluster)}
save(res, file="results/simulations_res2.Rdata")



################################# simulation 3 #################################
n <- 100
ntest <- 1000
p <- 1000
alpha <- 0.5
lambda <- 100
G <- 5
q <- 0.5

rho <- 0.7
Sigma <- matrix(rho, nrow=p/G, ncol=p/G)
diag(Sigma) <- 1

part <- rep(c(1:G), each=p/G)
csel <- 2^c(1:8)
nreps <- 100

### analysis splits in parallel
ncores <- min(detectCores() - 1, nreps)
cluster <- makeForkCluster(ncores)
if(parallel) {
  registerDoParallel(cluster)
} else {
  registerDoSEQ()
}

res <- foreach(k=c(1:nreps), .errorhandling="pass") %dopar% {
  set.seed(2019 + k)
  
  beta <- as.numeric(sapply(1:G, function(g) {
    b <- renet(p/G, lambda*alpha, 0.5*(1 - alpha)*lambda);
    b[abs(b)<=quantile(abs(b), q)] <- 0
    return(b)}))
  xtrain <- do.call(cbind, replicate(G, rmvnorm(n, mean=rep(0, p/G),
                                                sigma=Sigma), simplify=FALSE))
  ytrain <- rbinom(n, 1, as.numeric(1/(1 + exp(-xtrain %*% beta))))
  
  xtest <- do.call(cbind, replicate(G, rmvnorm(ntest, mean=rep(0, p/G),
                                               sigma=Sigma), simplify=FALSE))
  ytest <- rbinom(ntest, 1, as.numeric(1/(1 + exp(-xtest %*% beta))))
  
  # fitting models
  fit.gren1 <- gren(xtrain, ytrain, partitions=list(part=part), alpha=0.05, 
                    standardize=TRUE, trace=FALSE, psel=csel)
  fit.gren2 <- gren(xtrain, ytrain, partitions=list(part=part), alpha=0.5, 
                    standardize=TRUE, trace=FALSE, psel=csel)
  fit.gren3 <- gren(xtrain, ytrain, partitions=list(part=part), alpha=0.95, 
                    standardize=TRUE, trace=FALSE, psel=csel)
  
  fit.grridge <- grridge(t(xtrain), ytrain, list(part=split(1:p, part)))
  
  fit.cmcp1 <- sel.grpreg(xtrain, ytrain, part, penalty="cMCP", 
                          family="binomial", alpha=0.05, psel=csel)
  fit.cmcp2 <- sel.grpreg(xtrain, ytrain, part, penalty="cMCP", 
                          family="binomial", 
                          alpha=0.5, psel=csel)
  fit.cmcp3 <- sel.grpreg(xtrain, ytrain, part, penalty="cMCP", 
                          family="binomial", alpha=0.95, psel=csel)
  
  fit.gel1 <- sel.grpreg(xtrain, ytrain, part, penalty="gel", family="binomial", 
                         alpha=0.05, psel=csel)
  fit.gel2 <- sel.grpreg(xtrain, ytrain, part, penalty="gel", family="binomial", 
                         alpha=0.5, psel=csel)
  fit.gel3 <- sel.grpreg(xtrain, ytrain, part, penalty="gel", family="binomial", 
                         alpha=0.95, psel=csel)
  
  pred <- data.frame(ridge=predict.grridge(fit.grridge, t(xtest))[, 1],
                     grridge=predict.grridge(fit.grridge, t(xtest))[, 2],
                     gren1=predict(fit.gren1, xtest, type="groupreg"),
                     gren2=predict(fit.gren2, xtest, type="groupreg"),
                     gren3=predict(fit.gren3, xtest, type="groupreg"),
                     enet1=predict(fit.gren1, xtest, type="regular"),
                     enet2=predict(fit.gren2, xtest, type="regular"),
                     enet3=predict(fit.gren3, xtest, type="regular"),
                     cmcp1=predict(fit.cmcp1, xtest),
                     cmcp2=predict(fit.cmcp2, xtest),
                     cmcp3=predict(fit.cmcp3, xtest),
                     gel1=predict(fit.gel1, xtest),
                     gel2=predict(fit.gel2, xtest),
                     gel3=predict(fit.gel3, xtest))
  
  auc <- apply(pred, 2, function(s) {pROC::auc(ytest, s)})
  
  const <- sum((ytest - mean(ytest))^2)
  briers <- apply(pred, 2, function(s) {1 - sum((ytest - s)^2)/const})
  
  coef <- data.frame(ridge=coef(fit.grridge$predobj$NoGroups),
                     grridge=coef(fit.grridge$predobj$GroupRegul),
                     gren1=as.matrix(coef(fit.gren1, type="groupreg")),
                     gren2=as.matrix(coef(fit.gren2, type="groupreg")),
                     gren3=as.matrix(coef(fit.gren3, type="groupreg")),
                     enet1=as.matrix(coef(fit.gren1, type="regular")),
                     enet2=as.matrix(coef(fit.gren2, type="regular")),
                     enet3=as.matrix(coef(fit.gren3, type="regular")),
                     cmcp1=coef(fit.cmcp1),
                     cmcp2=coef(fit.cmcp2),
                     cmcp3=coef(fit.cmcp3),
                     gel1=coef(fit.gel1),
                     gel2=coef(fit.gel2),
                     gel3=coef(fit.gel3))
  
  mse <- apply(coef, 2, function(s) {mean((s - c(0, beta))^2)})
  
  kappa <- apply(coef, 2, function(s) {kappa2(cbind(s[-1]!=0, beta!=0))$value})
  
  psel <- c(ridge=p, grridge=p, 
            gren1=fit.gren1$freq.model$groupreg$df,
            gren2=fit.gren2$freq.model$groupreg$df,
            gren3=fit.gren3$freq.model$groupreg$df,
            enet1=fit.gren1$freq.model$regular$df,
            enet2=fit.gren2$freq.model$regular$df,
            enet3=fit.gren3$freq.model$regular$df,
            cmcp1=apply(as.matrix(coef(fit.cmcp1)), 2, function(b) {sum(b!=0)}),
            cmcp2=apply(as.matrix(coef(fit.cmcp2)), 2, function(b) {sum(b!=0)}),
            cmcp3=apply(as.matrix(coef(fit.cmcp3)), 2, function(b) {sum(b!=0)}),
            gel1=apply(as.matrix(coef(fit.gel1)), 2, function(b) {sum(b!=0)}),
            gel2=apply(as.matrix(coef(fit.gel2)), 2, function(b) {sum(b!=0)}),
            gel3=apply(as.matrix(coef(fit.gel3)), 2, function(b) {sum(b!=0)}))
  
  mults <- c(grridge=fit.grridge$lambdamults$part,
             gren1=fit.gren1$lambdag$part,
             gren2=fit.gren2$lambdag$part,
             gren3=fit.gren3$lambdag$part)
  
  list(psel=psel, auc=auc, briers=briers, mse=mse, kappa=kappa, mults=mults)
  
}
if(parallel) {stopCluster(cluster)}
save(res, file="results/simulations_res3.Rdata")




################################# simulation 4 #################################
n <- 100
p <- 1000
ntest <- 1000

alpha <- 0.5
lambda <- 100
G <- 4
lambdag <- exp(seq(-2, 2, length.out=4))
q <- 0.5
pblock <- 25
rho <- 0.7
Sigma <- matrix(rho, ncol=pblock, nrow=pblock); diag(Sigma) <- 1

part <- rep(c(1:G), each=p/G)
csel <- 2^c(1:8)
nreps <- 100

### analysis splits in parallel
ncores <- min(detectCores() - 1, nreps)
cluster <- makeForkCluster(ncores)
if(parallel) {
  registerDoParallel(cluster)
} else {
  registerDoSEQ()
}

res <- foreach(k=c(1:nreps), .errorhandling="pass") %dopar% {

  print(paste("rep", k))
  set.seed(2019 + k)
  
  beta <- as.numeric(sapply(1:G, function(g) {
    b <- renet(p/G, lambda*lambdag[g]*alpha, 0.5*(1 - alpha)*lambda*lambdag[g]);
    b[abs(b)<=quantile(abs(b), q)] <- 0
    return(b)}))
  xtrain <- do.call(cbind, replicate(p/pblock, rmvnorm(
    n, mean=rep(0, pblock), sigma=Sigma), simplify=FALSE))
  ytrain <- rbinom(n, 1, as.numeric(1/(1 + exp(-xtrain %*% beta))))
  
  xtest <- do.call(cbind, replicate(p/pblock, rmvnorm(
    ntest, mean=rep(0, pblock), sigma=Sigma), simplify=FALSE))
  ytest <- rbinom(ntest, 1, as.numeric(1/(1 + exp(-xtest %*% beta))))
  
  # fitting models
  fit.gren1 <- gren(xtrain, ytrain, partitions=list(part=part), alpha=0.05, 
                    standardize=TRUE, trace=FALSE, psel=csel)
  fit.gren2 <- gren(xtrain, ytrain, partitions=list(part=part), alpha=0.5, 
                    standardize=TRUE, trace=FALSE, psel=csel)
  fit.gren3 <- gren(xtrain, ytrain, partitions=list(part=part), alpha=0.95, 
                    standardize=TRUE, trace=FALSE, psel=csel)
  
  fit.grridge <- grridge(t(xtrain), ytrain, list(part=split(1:p, part)))
  
  fit.cmcp1 <- sel.grpreg(xtrain, ytrain, part, penalty="cMCP", 
                          family="binomial", alpha=0.05, psel=csel)
  fit.cmcp2 <- sel.grpreg(xtrain, ytrain, part, penalty="cMCP", 
                          family="binomial", 
                          alpha=0.5, psel=csel)
  fit.cmcp3 <- sel.grpreg(xtrain, ytrain, part, penalty="cMCP", 
                          family="binomial", alpha=0.95, psel=csel)
  
  fit.gel1 <- sel.grpreg(xtrain, ytrain, part, penalty="gel", family="binomial", 
                         alpha=0.05, psel=csel)
  fit.gel2 <- sel.grpreg(xtrain, ytrain, part, penalty="gel", family="binomial", 
                         alpha=0.5, psel=csel)
  fit.gel3 <- sel.grpreg(xtrain, ytrain, part, penalty="gel", family="binomial", 
                         alpha=0.95, psel=csel)
  
  pred <- data.frame(ridge=predict.grridge(fit.grridge, t(xtest))[, 1],
                     grridge=predict.grridge(fit.grridge, t(xtest))[, 2],
                     gren1=predict(fit.gren1, xtest, type="groupreg"),
                     gren2=predict(fit.gren2, xtest, type="groupreg"),
                     gren3=predict(fit.gren3, xtest, type="groupreg"),
                     enet1=predict(fit.gren1, xtest, type="regular"),
                     enet2=predict(fit.gren2, xtest, type="regular"),
                     enet3=predict(fit.gren3, xtest, type="regular"),
                     cmcp1=predict(fit.cmcp1, xtest),
                     cmcp2=predict(fit.cmcp2, xtest),
                     cmcp3=predict(fit.cmcp3, xtest),
                     gel1=predict(fit.gel1, xtest),
                     gel2=predict(fit.gel2, xtest),
                     gel3=predict(fit.gel3, xtest))
  
  auc <- apply(pred, 2, function(s) {pROC::auc(ytest, s)})
  
  const <- sum((ytest - mean(ytest))^2)
  briers <- apply(pred, 2, function(s) {1 - sum((ytest - s)^2)/const})
  
  coef <- data.frame(ridge=coef(fit.grridge$predobj$NoGroups),
                     grridge=coef(fit.grridge$predobj$GroupRegul),
                     gren1=as.matrix(coef(fit.gren1, type="groupreg")),
                     gren2=as.matrix(coef(fit.gren2, type="groupreg")),
                     gren3=as.matrix(coef(fit.gren3, type="groupreg")),
                     enet1=as.matrix(coef(fit.gren1, type="regular")),
                     enet2=as.matrix(coef(fit.gren2, type="regular")),
                     enet3=as.matrix(coef(fit.gren3, type="regular")),
                     cmcp1=coef(fit.cmcp1),
                     cmcp2=coef(fit.cmcp2),
                     cmcp3=coef(fit.cmcp3),
                     gel1=coef(fit.gel1),
                     gel2=coef(fit.gel2),
                     gel3=coef(fit.gel3))
  
  mse <- apply(coef, 2, function(s) {mean((s - c(0, beta))^2)})
  
  kappa <- apply(coef, 2, function(s) {kappa2(cbind(s[-1]!=0, beta!=0))$value})
  
  psel <- c(ridge=p, grridge=p, 
            gren1=fit.gren1$freq.model$groupreg$df,
            gren2=fit.gren2$freq.model$groupreg$df,
            gren3=fit.gren3$freq.model$groupreg$df,
            enet1=fit.gren1$freq.model$regular$df,
            enet2=fit.gren2$freq.model$regular$df,
            enet3=fit.gren3$freq.model$regular$df,
            cmcp1=apply(as.matrix(coef(fit.cmcp1)), 2, function(b) {sum(b!=0)}),
            cmcp2=apply(as.matrix(coef(fit.cmcp2)), 2, function(b) {sum(b!=0)}),
            cmcp3=apply(as.matrix(coef(fit.cmcp3)), 2, function(b) {sum(b!=0)}),
            gel1=apply(as.matrix(coef(fit.gel1)), 2, function(b) {sum(b!=0)}),
            gel2=apply(as.matrix(coef(fit.gel2)), 2, function(b) {sum(b!=0)}),
            gel3=apply(as.matrix(coef(fit.gel3)), 2, function(b) {sum(b!=0)}))
  
  mults <- c(grridge=fit.grridge$lambdamults$part,
             gren1=fit.gren1$lambdag$part,
             gren2=fit.gren2$lambdag$part,
             gren3=fit.gren3$lambdag$part)
  
  list(psel=psel, auc=auc, briers=briers, mse=mse, kappa=kappa, mults=mults)
  
}
if(parallel) {stopCluster(cluster)}
save(res, file="results/simulations_res4.Rdata")



################################# simulation 5 #################################
n <- 100
p <- 1000
ntest <- 1000

G <- 10
rho <- 0.5
Sigma <- diag(p)
for(i in 1:p) {
  for(j in 1:p) {
    Sigma[i, j] <- rho^abs(i - j)
  }
}

beta.group <- c(rep(0, 8), 0.2, 0.5)
q.zero <- 0.85

part <- rep(c(1:G), each=p/G)
csel <- 2^c(1:8)
nreps <- 100

### analysis splits in parallel
ncores <- min(detectCores() - 1, nreps)
cluster <- makeForkCluster(ncores)
if(parallel) {
  registerDoParallel(cluster)
} else {
  registerDoSEQ()
}

res <- foreach(k=c(1:nreps), .errorhandling="pass") %dopar% {
  
  print(paste("rep", k))
  set.seed(2019 + k)
  
  beta <- numeric(p)
  for(g in 1:G) {
    beta[((g - 1)*p/G + 1):(g*p/G)] <- rep(beta.group[g], p/G)
    if(beta.group[g]!=0) {
      beta[((g - 1)*p/G + 1):((g - 1)*p/G + q.zero*p/G)] <- 0
    }
  }
  
  xtrain <- rmvnorm(n, mean=rep(0, p), sigma=Sigma)
  ytrain <- rbinom(n, 1, as.numeric(1/(1 + exp(-xtrain %*% beta))))
  
  xtest <- rmvnorm(ntest, mean=rep(0, p), sigma=Sigma)
  ytest <- rbinom(ntest, 1, as.numeric(1/(1 + exp(-xtest %*% beta))))
  
  # fitting models
  fit.gren1 <- gren(xtrain, ytrain, partitions=list(part=part), alpha=0.05, 
                    standardize=TRUE, trace=FALSE, psel=csel)
  fit.gren2 <- gren(xtrain, ytrain, partitions=list(part=part), alpha=0.5, 
                    standardize=TRUE, trace=FALSE, psel=csel)
  fit.gren3 <- gren(xtrain, ytrain, partitions=list(part=part), alpha=0.95, 
                    standardize=TRUE, trace=FALSE, psel=csel)
  
  fit.grridge <- grridge(t(xtrain), ytrain, list(part=split(1:p, part)))
  
  fit.cmcp1 <- sel.grpreg(xtrain, ytrain, part, penalty="cMCP", 
                          family="binomial", alpha=0.05, psel=csel)
  fit.cmcp2 <- sel.grpreg(xtrain, ytrain, part, penalty="cMCP", 
                          family="binomial", 
                          alpha=0.5, psel=csel)
  fit.cmcp3 <- sel.grpreg(xtrain, ytrain, part, penalty="cMCP", 
                          family="binomial", alpha=0.95, psel=csel)
  
  fit.gel1 <- sel.grpreg(xtrain, ytrain, part, penalty="gel", family="binomial", 
                         alpha=0.05, psel=csel)
  fit.gel2 <- sel.grpreg(xtrain, ytrain, part, penalty="gel", family="binomial", 
                         alpha=0.5, psel=csel)
  fit.gel3 <- sel.grpreg(xtrain, ytrain, part, penalty="gel", family="binomial", 
                         alpha=0.95, psel=csel)
  
  pred <- data.frame(ridge=predict.grridge(fit.grridge, t(xtest))[, 1],
                     grridge=predict.grridge(fit.grridge, t(xtest))[, 2],
                     gren1=predict(fit.gren1, xtest, type="groupreg"),
                     gren2=predict(fit.gren2, xtest, type="groupreg"),
                     gren3=predict(fit.gren3, xtest, type="groupreg"),
                     enet1=predict(fit.gren1, xtest, type="regular"),
                     enet2=predict(fit.gren2, xtest, type="regular"),
                     enet3=predict(fit.gren3, xtest, type="regular"),
                     cmcp1=predict(fit.cmcp1, xtest),
                     cmcp2=predict(fit.cmcp2, xtest),
                     cmcp3=predict(fit.cmcp3, xtest),
                     gel1=predict(fit.gel1, xtest),
                     gel2=predict(fit.gel2, xtest),
                     gel3=predict(fit.gel3, xtest))
  
  auc <- apply(pred, 2, function(s) {pROC::auc(ytest, s)})
  
  const <- sum((ytest - mean(ytest))^2)
  briers <- apply(pred, 2, function(s) {1 - sum((ytest - s)^2)/const})
  
  coef <- data.frame(ridge=coef(fit.grridge$predobj$NoGroups),
                     grridge=coef(fit.grridge$predobj$GroupRegul),
                     gren1=as.matrix(coef(fit.gren1, type="groupreg")),
                     gren2=as.matrix(coef(fit.gren2, type="groupreg")),
                     gren3=as.matrix(coef(fit.gren3, type="groupreg")),
                     enet1=as.matrix(coef(fit.gren1, type="regular")),
                     enet2=as.matrix(coef(fit.gren2, type="regular")),
                     enet3=as.matrix(coef(fit.gren3, type="regular")),
                     cmcp1=coef(fit.cmcp1),
                     cmcp2=coef(fit.cmcp2),
                     cmcp3=coef(fit.cmcp3),
                     gel1=coef(fit.gel1),
                     gel2=coef(fit.gel2),
                     gel3=coef(fit.gel3))
  
  mse <- apply(coef, 2, function(s) {mean((s - c(0, beta))^2)})
  
  kappa <- apply(coef, 2, function(s) {kappa2(cbind(s[-1]!=0, beta!=0))$value})
  
  psel <- c(ridge=p, grridge=p, 
            gren1=fit.gren1$freq.model$groupreg$df,
            gren2=fit.gren2$freq.model$groupreg$df,
            gren3=fit.gren3$freq.model$groupreg$df,
            enet1=fit.gren1$freq.model$regular$df,
            enet2=fit.gren2$freq.model$regular$df,
            enet3=fit.gren3$freq.model$regular$df,
            cmcp1=apply(as.matrix(coef(fit.cmcp1)), 2, function(b) {sum(b!=0)}),
            cmcp2=apply(as.matrix(coef(fit.cmcp2)), 2, function(b) {sum(b!=0)}),
            cmcp3=apply(as.matrix(coef(fit.cmcp3)), 2, function(b) {sum(b!=0)}),
            gel1=apply(as.matrix(coef(fit.gel1)), 2, function(b) {sum(b!=0)}),
            gel2=apply(as.matrix(coef(fit.gel2)), 2, function(b) {sum(b!=0)}),
            gel3=apply(as.matrix(coef(fit.gel3)), 2, function(b) {sum(b!=0)}))
  
  mults <- c(grridge=fit.grridge$lambdamults$part,
             gren1=fit.gren1$lambdag$part,
             gren2=fit.gren2$lambdag$part,
             gren3=fit.gren3$lambdag$part)
  
  list(psel=psel, auc=auc, briers=briers, mse=mse, kappa=kappa, mults=mults)
  
}
if(parallel) {stopCluster(cluster)}
save(res, file="results/simulations_res5.Rdata")