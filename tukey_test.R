# paths
path.code <- as.character(ifelse(Sys.info()[1]=="Darwin", 
                                 "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/code/" ,
                                 "~/EBEN/code/"))
path.graph <- "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/graphs/"
path.res <- as.character(ifelse(Sys.info()[1]=="Darwin", 
                                "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/results/" ,
                                "~/EBEN/results/"))

### libraries
library(mvtnorm)
library(GRridge)
library(pROC)

### source grEBEN
source(paste(path.code, "grVBEM.R", sep=""))
source(paste(path.code, "mygrridge.R", sep=""))

## simulation 2
set.seed(456)
# set data characteristics
n <- 200
p <- 1000
G <- 5
pblock <- 20
rho <- 0.9
sigma <- matrix(rho, ncol=pblock, nrow=pblock); diag(sigma) <- 1
lambda <- 0.02
alpha <- 0.5
lambdag <- exp(seq(-1, 1, length.out=G))
m <- rep(1, n)
part.greben <- list(groups=rep(1:G, each=p/G))
part.grridge <- list(groups=CreatePartition(as.factor(part.greben$groups)))
ntest <- 1000

methods <- c("enet1", "enet2", "enet3")
nreps <- 50
auc2 <- list(vector(mode="list", length=nreps), 
             vector(mode="list", length=nreps), 
             vector(mode="list", length=nreps))
names(auc2) <- methods
briers2 <- list(vector(mode="list", length=nreps), 
                vector(mode="list", length=nreps), 
                vector(mode="list", length=nreps))
names(briers2) <- methods
varbeta2 <- matrix(nrow=nreps, ncol=G)
# the simulations
for(r in 1:nreps) {
  
  x <- do.call(cbind, replicate(p/pblock, rmvnorm(n, mean=rep(0, pblock), sigma=sigma), 
                                simplify=FALSE))
  beta <- as.numeric(sapply(1:G, function(g) {
    renbeta(p/G, 2*n*lambda*alpha*sqrt(lambdag[g]), n*lambda*(1 - alpha)*lambdag[g])}))
  prob <- as.numeric(exp(x %*% beta)/(1 + exp(x %*% beta)))
  y <- rbinom(n, 1, prob)
  
  xtest <- do.call(cbind, replicate(p/pblock, rmvnorm(ntest, mean=rep(0, pblock),
                                                      sigma=sigma), simplify=FALSE))
  probtest <- as.numeric(exp(xtest %*% beta)/(1 + exp(xtest %*% beta)))
  ytest <- rbinom(ntest, 1, probtest)
  
  fit1.greben <- cv.glmnet(x, y, family="binomial", standardize=FALSE, 
                           alpha=0.05)
  fit2.greben <- cv.glmnet(x, y, family="binomial", standardize=FALSE, 
                           alpha=0.5)
  fit3.greben <- cv.glmnet(x, y, family="binomial", standardize=FALSE,
                           alpha=0.95)
  
  auc1.enet <- sapply(fit1.greben$lambda, function(l) {
    pred <- as.numeric(predict(fit1.greben, xtest, s=l, type="response"))
    pROC::roc(ytest, pred)$auc})
  auc2.enet <- sapply(fit2.greben$lambda, function(l) {
    pred <- as.numeric(predict(fit2.greben, xtest, s=l, type="response"))
    pROC::roc(ytest, pred)$auc})
  auc3.enet <- sapply(fit3.greben$lambda, function(l) {
    pred <- as.numeric(predict(fit3.greben, xtest, s=l, type="response"))
    pROC::roc(ytest, pred)$auc})
  
  briers1.enet <- sapply(fit1.greben$lambda, function(l) {
    pred <- as.numeric(predict(fit1.greben, xtest, s=l, type="response"))
    1 - sum((ytest - pred)^2)/sum((ytest - mean(ytest))^2)})
  briers2.enet <- sapply(fit2.greben$lambda, function(l) {
    pred <- as.numeric(predict(fit2.greben, xtest, s=l, type="response"))
    1 - sum((ytest - pred)^2)/sum((ytest - mean(ytest))^2)})
  briers3.enet <- sapply(fit3.greben$lambda, function(l) {
    pred <- as.numeric(predict(fit3.greben, xtest, s=l, type="response"))
    1 - sum((ytest - pred)^2)/sum((ytest - mean(ytest))^2)})
  
  psel1.enet <- fit1.greben$nzero
  psel2.enet <- fit2.greben$nzero
  psel3.enet <- fit3.greben$nzero
  
  auc2[[1]][[r]] <- cbind(psel=psel1.enet, auc=auc1.enet)
  auc2[[2]][[r]] <- cbind(psel=psel2.enet, auc=auc2.enet)
  auc2[[3]][[r]] <- cbind(psel=psel3.enet, auc=auc3.enet)
  
  briers2[[1]][[r]] <- cbind(psel=psel1.enet, auc=briers1.enet)
  briers2[[2]][[r]] <- cbind(psel=psel2.enet, auc=briers2.enet)
  briers2[[3]][[r]] <- cbind(psel=psel3.enet, auc=briers3.enet)
  
  varbeta2[r, ] <- sapply(1:G, function(g) {var(beta[(p*(g - 1)/G + 1):(p*g/G)])})
  
  results2 <- list(auc=auc2, briers=briers2, varbeta=varbeta2)
  
  save(results2, file=paste(path.res, "tukey_test_res.Rdata", sep=""))

}

load(paste(path.res, "tukey_test_res.Rdata", sep=""))


cbind(results2$auc[[1]][[1]][c(1:50), ], results2$auc[[1]][[2]][c(1:50), ],
      results2$auc[[1]][[3]][c(1:50), ])

