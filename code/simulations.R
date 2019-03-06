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
nreps <- 100
psel1 <- kappa1 <- data.frame(enet1=matrix(NA, nrow=nreps, ncol=100),
                              enet2=matrix(NA, nrow=nreps, ncol=100),
                              enet3=matrix(NA, nrow=nreps, ncol=100),
                              gren1=matrix(NA, nrow=nreps, ncol=100),
                              gren2=matrix(NA, nrow=nreps, ncol=100),
                              gren3=matrix(NA, nrow=nreps, ncol=100),
                              sgl1=matrix(NA, nrow=nreps, ncol=20),
                              sgl2=matrix(NA, nrow=nreps, ncol=20),
                              sgl3=matrix(NA, nrow=nreps, ncol=20),
                              cmcp1=matrix(NA, nrow=nreps, ncol=100),
                              cmcp2=matrix(NA, nrow=nreps, ncol=100),
                              cmcp3=matrix(NA, nrow=nreps, ncol=100),
                              gel1=matrix(NA, nrow=nreps, ncol=100),
                              gel2=matrix(NA, nrow=nreps, ncol=100),
                              gel3=matrix(NA, nrow=nreps, ncol=100))
auc1 <- mse1 <- briers1 <- data.frame(ridge=numeric(nreps), 
                                      grridge=numeric(nreps),
                                      enet1=matrix(NA, nrow=nreps, ncol=100),
                                      enet2=matrix(NA, nrow=nreps, ncol=100),
                                      enet3=matrix(NA, nrow=nreps, ncol=100),
                                      gren1=matrix(NA, nrow=nreps, ncol=100),
                                      gren2=matrix(NA, nrow=nreps, ncol=100),
                                      gren3=matrix(NA, nrow=nreps, ncol=100),
                                      sgl1=matrix(NA, nrow=nreps, ncol=20),
                                      sgl2=matrix(NA, nrow=nreps, ncol=20),
                                      sgl3=matrix(NA, nrow=nreps, ncol=20),
                                      cmcp1=matrix(NA, nrow=nreps, ncol=100),
                                      cmcp2=matrix(NA, nrow=nreps, ncol=100),
                                      cmcp3=matrix(NA, nrow=nreps, ncol=100),
                                      gel1=matrix(NA, nrow=nreps, ncol=100),
                                      gel2=matrix(NA, nrow=nreps, ncol=100),
                                      gel3=matrix(NA, nrow=nreps, ncol=100))
lambdagest1 <- list(grridge=matrix(NA, nrow=nreps, ncol=G),
                    gren1=matrix(NA, nrow=nreps, ncol=G),
                    gren2=matrix(NA, nrow=nreps, ncol=G),
                    gren3=matrix(NA, nrow=nreps, ncol=G))

### analysis splits in parallel
ncores <- min(detectCores() - 1, nreps)
cluster <- makeCluster(ncores, type="FORK")
registerDoParallel(cluster)

out <- foreach(k=c(1:nreps)) %dopar% {
  set.seed(2019 + k)
  
  for(g in 1:G) {
    idz <- c(((g-1)*p/G + 1):(g*p/G - rg))
    ida <- c(((g*p/G - rg) + 1):(g*p/G))
    beta[idz] <- 0
    beta[ida] <- runif(rg, 0, b1*f^(g-1))
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
  
  const <- sum((ytest - mean(ytest))^2)
  briers1$ridge[r] <- 1 - sum((ytest - pred1.ridge)^2)/const
  
  coef(fit1.gren1, type="groupreg")
  str(coef(fit1.gren1$freq.model$groupreg))
  
  mse1$ridge[r] <- mean((coef(fit1.ridge)[-1] - beta)^2)
  
  mse1$grridge[r, ] <- sapply(fit1.grridge, function(s) {
    mean((replace(rep(0, p), s$resEN$whichEN, s$resEN$betasEN) - beta)^2)})
  
  mse1$gren1[r, ] <- apply(fit1.gren1$freq.model$groupreg$beta, 2, function(b) {
    mean((b - beta)^2)})
  mse1$gren2[r, ] <- apply(fit1.gren2$freq.model$groupreg$beta, 2, function(b) {
    mean((b - beta)^2)})
  mse1$gren3[r, ] <- apply(fit1.gren3$freq.model$groupreg$beta, 2, function(b) {
    mean((b - beta)^2)})
  
  mse1$enet1[r, ] <- apply(fit1.gren1$freq.model$regular$beta, 2, function(b) {
    mean((b - beta)^2)})
  mse1$enet2[r, ] <- apply(fit1.gren2$freq.model$regular$beta, 2, function(b) {
    mean((b - beta)^2)})
  mse1$enet3[r, ] <- apply(fit1.gren3$freq.model$regular$beta, 2, function(b) {
    mean((b - beta)^2)})
  
  mse1$sglasso1[r, ] <- apply(fit1.sglasso1$beta, 2, function(b) {
    mean((b - beta)^2)})
  mse1$sglasso2[r, ] <- apply(fit1.sglasso2$beta, 2, function(b) {
    mean((b - beta)^2)})
  mse1$sglasso3[r, ] <- apply(fit1.sglasso3$beta, 2, function(b) {
    mean((b - beta)^2)})
  
  mse1$cmcp1[r, ] <- apply(fit1.cmcp1$beta[-1, ], 2, function(b) {
    mean((b - beta)^2)})
  mse1$cmcp2[r, ] <- apply(fit1.cmcp2$beta[-1, ], 2, function(b) {
    mean((b - beta)^2)})
  mse1$cmcp3[r, ] <- apply(fit1.cmcp3$beta[-1, ], 2, function(b) {
    mean((b - beta)^2)})
  
  mse1$gelasso1[r, ] <- apply(fit1.gelasso1$beta[-1, ], 2, function(b) {
    mean((b - beta)^2)})
  mse1$gelasso2[r, ] <- apply(fit1.gelasso2$beta[-1, ], 2, function(b) {
    mean((b - beta)^2)})
  mse1$gelasso3[r, ] <- apply(fit1.gelasso3$beta[-1, ], 2, function(b) {
    mean((b - beta)^2)})
  
  kappa1$grridge[r, ] <- sapply(fit1.grridge, function(s) {
    kappa2(cbind(beta!=0, replace(rep(FALSE, p), s$resEN$whichEN,
                                  TRUE)))$value})
  
  kappa1$gren1[r, ] <- apply(fit1.gren1$freq.model$groupreg$beta, 2,
                             function(b) {kappa2(cbind(beta!=0, b!=0))$value})
  kappa1$gren2[r, ] <- apply(fit1.gren2$freq.model$groupreg$beta, 2,
                             function(b) {kappa2(cbind(beta!=0, b!=0))$value})
  kappa1$gren3[r, ] <- apply(fit1.gren3$freq.model$groupreg$beta, 2,
                             function(b) {kappa2(cbind(beta!=0, b!=0))$value})
  
  kappa1$enet1[r, ] <- apply(fit1.gren1$freq.model$regular$beta, 2,
                             function(b) {kappa2(cbind(beta!=0, b!=0))$value})
  kappa1$enet2[r, ] <- apply(fit1.gren2$freq.model$regular$beta, 2,
                             function(b) {kappa2(cbind(beta!=0, b!=0))$value})
  kappa1$enet3[r, ] <- apply(fit1.gren3$freq.model$regular$beta, 2,
                             function(b) {kappa2(cbind(beta!=0, b!=0))$value})
  
  kappa1$sglasso1[r, ] <- apply(fit1.sglasso1$beta, 2, function(b) {
    kappa2(cbind(beta!=0, b!=0))$value})
  kappa1$sglasso2[r, ] <- apply(fit1.sglasso2$beta, 2, function(b) {
    kappa2(cbind(beta!=0, b!=0))$value})
  kappa1$sglasso3[r, ] <- apply(fit1.sglasso3$beta, 2, function(b) {
    kappa2(cbind(beta!=0, b!=0))$value})
  
  kappa1$cmcp1[r, ] <- apply(fit1.cmcp1$beta[-1, ], 2, function(b) {
    kappa2(cbind(beta!=0, b!=0))$value})
  kappa1$cmcp2[r, ] <- apply(fit1.cmcp2$beta[-1, ], 2, function(b) {
    kappa2(cbind(beta!=0, b!=0))$value})
  kappa1$cmcp3[r, ] <- apply(fit1.cmcp3$beta[-1, ], 2, function(b) {
    kappa2(cbind(beta!=0, b!=0))$value})
  
  kappa1$gelasso1[r, ] <- apply(fit1.gelasso1$beta[-1, ], 2, function(b) {
    kappa2(cbind(beta!=0, b!=0))$value})
  kappa1$gelasso2[r, ] <- apply(fit1.gelasso2$beta[-1, ], 2, function(b) {
    kappa2(cbind(beta!=0, b!=0))$value})
  kappa1$gelasso3[r, ] <- apply(fit1.gelasso3$beta[-1, ], 2, function(b) {
    kappa2(cbind(beta!=0, b!=0))$value})
  
  psel1$grridge[r, ] <- sapply(fit1.grridge, function(s) {
    length(s$resEN$whichEN)})
  
  psel1$gren1[r, ] <- fit1.gren1$freq.model$groupreg$df
  psel1$gren2[r, ] <- fit1.gren2$freq.model$groupreg$df
  psel1$gren3[r, ] <- fit1.gren3$freq.model$groupreg$df
  
  psel1$enet1[r, ] <- fit1.gren1$freq.model$regular$df
  psel1$enet2[r, ] <- fit1.gren2$freq.model$regular$df
  psel1$enet3[r, ] <- fit1.gren3$freq.model$regular$df
  
  psel1$sglasso1[r, ] <- apply(fit1.sglasso1$beta, 2, function(b) {sum(b!=0)})
  psel1$sglasso2[r, ] <- apply(fit1.sglasso2$beta, 2, function(b) {sum(b!=0)})
  psel1$sglasso3[r, ] <- apply(fit1.sglasso3$beta, 2, function(b) {sum(b!=0)})
  
  psel1$cmcp1[r, ] <- apply(fit1.cmcp1$beta, 2, function(b) {sum(b!=0)})
  psel1$cmcp2[r, ] <- apply(fit1.cmcp2$beta, 2, function(b) {sum(b!=0)})
  psel1$cmcp3[r, ] <- apply(fit1.cmcp3$beta, 2, function(b) {sum(b!=0)})
  
  psel1$gelasso1[r, ] <- apply(fit1.gelasso1$beta, 2, function(b) {sum(b!=0)})
  psel1$gelasso2[r, ] <- apply(fit1.gelasso2$beta, 2, function(b) {sum(b!=0)})
  psel1$gelasso3[r, ] <- apply(fit1.gelasso3$beta, 2, function(b) {sum(b!=0)})
  
  lambdagest1$grridge[r, ] <- fit1.grridge[[1]]$lambdamults$groups
  lambdagest1$gren1[r, ] <- fit1.gren1$lambdag$groups
  lambdagest1$gren2[r, ] <- fit1.gren2$lambdag$groups
  lambdagest1$gren3[r, ] <- fit1.gren3$lambdag$groups
  
  results1 <- list(auc=auc1, briers=briers1, mse=mse1, kappa=kappa1, psel=psel1,
                   lambdag=lambdagest1)
  save(results1, file=paste(path.res, "gren_sim1_res1.Rdata", sep=""))
  
}