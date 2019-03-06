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

### libraries
library(gren)
library(GRridge)
library(grpreg)
library(SGL)
library(foreach)
library(doParallel)
library(irr)
library(pROC)
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
nreps <- 1
ncores <- min(detectCores() - 1, nreps)

### analysis splits in parallel
cluster <- makeCluster(ncores, type="FORK")
registerDoParallel(cluster)

res1 <- foreach(k=c(1:nreps)) %dopar% {
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
  fit1.gren1 <- gren(xtrain, ytrain, partitions=list(part=part), alpha=0.05, 
                     standardize=TRUE, trace=FALSE)
  fit1.gren2 <- gren(xtrain, ytrain, partitions=list(part=part), alpha=0.5, 
                     standardize=TRUE, trace=FALSE)
  fit1.gren3 <- gren(xtrain, ytrain, partitions=list(part=part), alpha=0.95, 
                     standardize=TRUE, trace=FALSE)
  
  fit1.grridge <- grridge(t(xtrain), ytrain, list(part=split(1:p, part)))
  
  fit1.sgl1 <- SGL(list(x=xtrain, y=ytrain), part, type="logit", alpha=0.05)
  fit1.sgl2 <- SGL(list(x=xtrain, y=ytrain), part, type="logit", alpha=0.5)
  fit1.sgl3 <- SGL(list(x=xtrain, y=ytrain), part, type="logit", alpha=0.95)
  
  fit1.cmcp1 <- grpreg(xtrain, ytrain, part, penalty="cMCP", family="binomial", 
                       alpha=0.05)
  fit1.cmcp2 <- grpreg(xtrain, ytrain, part, penalty="cMCP", family="binomial", 
                       alpha=0.5)
  fit1.cmcp3 <- grpreg(xtrain, ytrain, part, penalty="cMCP", family="binomial", 
                       alpha=0.95)
  
  fit1.gel1 <- grpreg(xtrain, ytrain, part, penalty="gel", family="binomial", 
                      alpha=0.05)
  fit1.gel2 <- grpreg(xtrain, ytrain, part, penalty="gel", family="binomial", 
                      alpha=0.5)
  fit1.gel3 <- grpreg(xtrain, ytrain, part, penalty="gel", family="binomial", 
                      alpha=0.95)
  
  pred1 <- data.frame(ridge=predict.grridge(fit1.grridge, t(xtest))[, 1],
                      grridge=predict.grridge(fit1.grridge, t(xtest))[, 2],
                      gren1=predict(fit1.gren1, xtest, type="groupreg"),
                      gren2=predict(fit1.gren2, xtest, type="groupreg"),
                      gren3=predict(fit1.gren3, xtest, type="groupreg"),
                      enet1=predict(fit1.gren1, xtest, type="regular"),
                      enet2=predict(fit1.gren2, xtest, type="regular"),
                      enet3=predict(fit1.gren3, xtest, type="regular"),
                      sgl1=predictSGL(fit1.sgl1, xtest),
                      sgl2=predictSGL(fit1.sgl2, xtest),
                      sgl3=predictSGL(fit1.sgl3, xtest),
                      cmcp1=predict(fit1.cmcp1, xtest),
                      cmcp2=predict(fit1.cmcp2, xtest),
                      cmcp3=predict(fit1.cmcp3, xtest),
                      gel1=predict(fit1.gel1, xtest),
                      gel2=predict(fit1.gel2, xtest),
                      gel3=predict(fit1.gel3, xtest))
  
  auc1 <- apply(pred1, 2, function(s) {pROC::auc(ytest, s)})
  
  const <- sum((ytest - mean(ytest))^2)
  briers1 <- apply(pred1, 2, function(s) {1 - sum((ytest - s)^2)/const})
  
  coef1 <- data.frame(ridge=coef(fit1.grridge$predobj$NoGroups),
                      grridge=coef(fit1.grridge$predobj$GroupRegul),
                      gren1=as.matrix(coef(fit1.gren1, type="groupreg")),
                      gren2=as.matrix(coef(fit1.gren2, type="groupreg")),
                      gren3=as.matrix(coef(fit1.gren3, type="groupreg")),
                      enet1=as.matrix(coef(fit1.gren1, type="regular")),
                      enet2=as.matrix(coef(fit1.gren2, type="regular")),
                      enet3=as.matrix(coef(fit1.gren3, type="regular")),
                      sgl1=rbind(fit1.sgl1$intercept, fit1.sgl1$beta),
                      sgl2=rbind(fit1.sgl2$intercept, fit1.sgl2$beta),
                      sgl3=rbind(fit1.sgl3$intercept, fit1.sgl3$beta),
                      cmcp1=coef(fit1.cmcp1),
                      cmcp2=coef(fit1.cmcp2),
                      cmcp3=coef(fit1.cmcp3),
                      gel1=coef(fit1.gel1),
                      gel2=coef(fit1.gel2),
                      gel3=coef(fit1.gel3))
  
  mse1 <- apply(coef1, 2, function(s) {mean((s - c(0, beta))^2)})
  
  kappa1 <- apply(coef1, 2, function(s) {
    kappa2(cbind(s[-1]!=0, beta!=0))$value})
  
  psel1 <- c(ridge=p, grridge=p, 
             gren1=fit1.gren1$freq.model$groupreg$df,
             gren2=fit1.gren2$freq.model$groupreg$df,
             gren3=fit1.gren3$freq.model$groupreg$df,
             enet1=fit1.gren1$freq.model$regular$df,
             enet2=fit1.gren2$freq.model$regular$df,
             enet3=fit1.gren3$freq.model$regular$df,
             sgl1=apply(fit1.sgl1$beta, 2, function(b) {sum(b!=0)}),
             sgl2=apply(fit1.sgl2$beta, 2, function(b) {sum(b!=0)}),
             sgl3=apply(fit1.sgl3$beta, 2, function(b) {sum(b!=0)}),
             cmcp1=apply(as.matrix(coef(fit1.cmcp1)), 2, function(b) {
               sum(b!=0)}),
             cmcp2=apply(as.matrix(coef(fit1.cmcp2)), 2, function(b) {
               sum(b!=0)}),
             cmcp3=apply(as.matrix(coef(fit1.cmcp3)), 2, function(b) {
               sum(b!=0)}),
             gel1=apply(as.matrix(coef(fit1.gel1)), 2, function(b) {sum(b!=0)}),
             gel2=apply(as.matrix(coef(fit1.gel2)), 2, function(b) {sum(b!=0)}),
             gel3=apply(as.matrix(coef(fit1.gel3)), 2, function(b) {sum(b!=0)}))
  
  mults1 <- c(grridge=fit1.grridge$lambdamults$part,
              gren1=fit1.gren1$lambdag$part,
              gren2=fit1.gren2$lambdag$part,
              gren3=fit1.gren3$lambdag$part)
  
  list(auc=auc1, briers=briers1, mse=mse1, kappa=kappa1, mults=mults1)
  
}


results1 <- list(auc=auc1, briers=briers1, mse=mse1, kappa=kappa1, psel=psel1,
                 lambdag=lambdagest1)
save(results1, file=paste(path.res, "gren_sim1_res1.Rdata", sep=""))