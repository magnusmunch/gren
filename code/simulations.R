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
nreps <- 2

### analysis splits in parallel
ncores <- min(detectCores() - 1, nreps)
cluster <- makeForkCluster(ncores)
if(parallel) {
  registerDoParallel(cluster)
} else {
  registerDoSEQ()
}

res <- foreach(k=c(1:nreps)) %dopar% {
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
                    standardize=TRUE, trace=FALSE)
  fit.gren2 <- gren(xtrain, ytrain, partitions=list(part=part), alpha=0.5, 
                    standardize=TRUE, trace=FALSE)
  fit.gren3 <- gren(xtrain, ytrain, partitions=list(part=part), alpha=0.95, 
                    standardize=TRUE, trace=FALSE)
  
  fit.grridge <- grridge(t(xtrain), ytrain, list(part=split(1:p, part)))
  
  fit.sgl1 <- SGL(list(x=xtrain, y=ytrain), part, type="logit", alpha=0.05)
  fit.sgl2 <- SGL(list(x=xtrain, y=ytrain), part, type="logit", alpha=0.5)
  fit.sgl3 <- SGL(list(x=xtrain, y=ytrain), part, type="logit", alpha=0.95)
  
  fit.cmcp1 <- grpreg(xtrain, ytrain, part, penalty="cMCP", family="binomial", 
                      alpha=0.05)
  fit.cmcp2 <- grpreg(xtrain, ytrain, part, penalty="cMCP", family="binomial", 
                      alpha=0.5)
  fit.cmcp3 <- grpreg(xtrain, ytrain, part, penalty="cMCP", family="binomial", 
                      alpha=0.95)
  
  fit.gel1 <- grpreg(xtrain, ytrain, part, penalty="gel", family="binomial", 
                     alpha=0.05)
  fit.gel2 <- grpreg(xtrain, ytrain, part, penalty="gel", family="binomial", 
                     alpha=0.5)
  fit.gel3 <- grpreg(xtrain, ytrain, part, penalty="gel", family="binomial", 
                     alpha=0.95)
  
  pred <- data.frame(ridge=predict.grridge(fit.grridge, t(xtest))[, 1],
                     grridge=predict.grridge(fit.grridge, t(xtest))[, 2],
                     gren1=predict(fit.gren1, xtest, type="groupreg"),
                     gren2=predict(fit.gren2, xtest, type="groupreg"),
                     gren3=predict(fit.gren3, xtest, type="groupreg"),
                     enet1=predict(fit.gren1, xtest, type="regular"),
                     enet2=predict(fit.gren2, xtest, type="regular"),
                     enet3=predict(fit.gren3, xtest, type="regular"),
                     sgl1=predictSGL(fit.sgl1, xtest),
                     sgl2=predictSGL(fit.sgl2, xtest),
                     sgl3=predictSGL(fit.sgl3, xtest),
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
                      sgl1=rbind(fit.sgl1$intercept, fit.sgl1$beta),
                      sgl2=rbind(fit.sgl2$intercept, fit.sgl2$beta),
                      sgl3=rbind(fit.sgl3$intercept, fit.sgl3$beta),
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
            sgl1=apply(fit.sgl1$beta, 2, function(b) {sum(b!=0)}),
            sgl2=apply(fit.sgl2$beta, 2, function(b) {sum(b!=0)}),
            sgl3=apply(fit.sgl3$beta, 2, function(b) {sum(b!=0)}),
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
res <- sapply(c("psel", "auc", "briers", "mse", "kappa", "mults"), function(s) {
  sapply(res, function(m) {m[[s]]}, simplify=FALSE)}, simplify=FALSE)
  

sapply(res$psel, function(s) {s[grep("gren1", names(s))]}, simplify=FALSE)




save(results1, file=paste(path.res, "gren_sim1_res1.Rdata", sep=""))