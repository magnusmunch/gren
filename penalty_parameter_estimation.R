################################### preamble ###################################
# penalty parameter selection simulations                                      #
# version: 01                                                                  #
# author: Magnus M?nch                                                         #
# created: 01-12-2017                                                          #
# last edited: 01-12-2017                                                      #
################################################################################

##################################### notes ####################################
#                                                                              #
################################################################################

### paths
path.data <- ifelse(as.character(Sys.info()[1])!="Darwin", "~/EBEN/data/lung_cancer_methylation",
                    "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/data/lung_cancer_methylation/")
path.code <- ifelse(as.character(Sys.info()[1])!="Darwin", "~/EBEN/code/",
                    "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/code/")
path.res <- ifelse(as.character(Sys.info()[1])!="Darwin", "~/EBEN/results/",
                   "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/results/")
path.graph <- "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/graphs/"

### loading libraries
library(glmnet)
library(psych)
library(mvtnorm)

### simulations
set.seed(3002)
n <- c(50, 100)
p <- 1000
k <- c(20, 50, 100)
pblock <- 20
rho <- c(0, 0.7, 0.9)
settings <- expand.grid(n=n, p=p, k=k, pblock=pblock, rho=rho)

nreps <- 50
mat <- matrix(NA, ncol=5, nrow=nreps, dimnames=list(NULL, c("CV", "AIC", "BIC",
                                                            "HBIC", "GIC")))
jaccard <- kappa <- rep(list(list(lasso=mat, enet1=mat, enet2=mat, enet3=mat)), 
                        nrow(settings))
for(s in 1:nrow(settings)) {
  cat("setting", s, "\n")
  for(r in 1:nreps) {
    cat("rep", r, "\n")
    sigma <- matrix(settings$rho[s], ncol=settings$pblock[s], 
                    nrow=settings$pblock[s]); diag(sigma) <- 1
    x <- do.call(cbind, replicate(
      settings$p[s]/settings$pblock[s], rmvnorm(settings$n[s], 
                                                mean=rep(0, settings$pblock[s]), 
                                                sigma=sigma), simplify=FALSE))
    beta <- c(rnorm(settings$k[s]), rep(0, settings$p[s] - settings$k[s]))
    prob <- as.numeric(exp(x %*% beta)/(1 + exp(x %*% beta)))
    y <- rbinom(settings$n[s], 1, prob)
  
    fit.lasso <- cv.glmnet(x, y, family="binomial", alpha=1)
    fit.enet1 <- cv.glmnet(x, y, family="binomial", alpha=0.05)
    fit.enet2 <- cv.glmnet(x, y, family="binomial", alpha=0.5)
    fit.enet3 <- cv.glmnet(x, y, family="binomial", alpha=0.95)
  
    # lasso
    psel <- fit.lasso$glmnet.fit$df
    dev <- (fit.lasso$glmnet.fit$dev.ratio - 1)*fit.lasso$glmnet.fit$nulldev
    lam.aic.lasso <- fit.lasso$glmnet.fit$lambda[which.min(2*psel - dev)]
    lam.bic.lasso <- fit.lasso$glmnet.fit$lambda[
      which.min(log(fit.lasso$glmnet.fit$nobs)*psel - dev)]
    lam.hbic.lasso <- fit.lasso$glmnet.fit$lambda[
      which.min(log(settings$p[s])*psel - dev)]
    lam.gic.lasso <- fit.lasso$glmnet.fit$lambda[
      which.min(log(log(fit.lasso$glmnet.fit$nobs))*
                  log(settings$p[s])*psel - dev)]
  
    sel.cv.lasso <- which(as.numeric(coef(fit.lasso, "lambda.min"))[-1]!=0)
    sel.aic.lasso <- which(as.numeric(coef(fit.lasso, lam.aic.lasso))[-1]!=0)
    sel.bic.lasso <- which(as.numeric(coef(fit.lasso, lam.bic.lasso))[-1]!=0)
    sel.hbic.lasso <- which(as.numeric(coef(fit.lasso, lam.hbic.lasso))[-1]!=0)
    sel.gic.lasso <- which(as.numeric(coef(fit.lasso, lam.gic.lasso))[-1]!=0)
  
    jaccard[[s]]$lasso[r, 1] <- 
      length(intersect(sel.cv.lasso, 1:settings$k[s]))/
      length(union(sel.cv.lasso, 1:settings$k[s]))
    jaccard[[s]]$lasso[r, 2] <- 
      length(intersect(sel.aic.lasso, 1:settings$k[s]))/
      length(union(sel.aic.lasso, 1:settings$k[s]))
    jaccard[[s]]$lasso[r, 3] <- 
      length(intersect(sel.bic.lasso, 1:settings$k))/
      length(union(sel.bic.lasso, 1:settings$k))
    jaccard[[s]]$lasso[r, 4] <- 
      length(intersect(sel.hbic.lasso, 1:settings$k))/
      length(union(sel.cv.lasso, 1:settings$k))
    jaccard[[s]]$lasso[r, 5] <- 
      length(intersect(sel.gic.lasso, 1:settings$k))/
      length(union(sel.cv.lasso, 1:settings$k))
  
    kappa[[s]]$lasso[r, 1] <- cohen.kappa(cbind(
      beta!=0, replace(rep(0, settings$p[s]), sel.cv.lasso, 1)))$kappa
    kappa[[s]]$lasso[r, 2] <- cohen.kappa(cbind(
      beta!=0, replace(rep(0, settings$p[s]), sel.aic.lasso, 1)))$kappa
    kappa[[s]]$lasso[r, 3] <- cohen.kappa(cbind(
      beta!=0, replace(rep(0, settings$p[s]), sel.bic.lasso, 1)))$kappa
    kappa[[s]]$lasso[r, 4] <- cohen.kappa(cbind(
      beta!=0, replace(rep(0, settings$p[s]), sel.hbic.lasso, 1)))$kappa
    kappa[[s]]$lasso[r, 5] <- cohen.kappa(cbind(
      beta!=0, replace(rep(0, settings$p[s]), sel.gic.lasso, 1)))$kappa
    
    # enet1
    psel <- fit.enet1$glmnet.fit$df
    dev <- (fit.enet1$glmnet.fit$dev.ratio - 1)*fit.enet1$glmnet.fit$nulldev
    lam.aic.enet1 <- fit.enet1$glmnet.fit$lambda[which.min(2*psel - dev)]
    lam.bic.enet1 <- fit.enet1$glmnet.fit$lambda[
      which.min(log(fit.enet1$glmnet.fit$nobs)*psel - dev)]
    lam.hbic.enet1 <- fit.enet1$glmnet.fit$lambda[
      which.min(log(settings$p[s])*psel - dev)]
    lam.gic.enet1 <- fit.enet1$glmnet.fit$lambda[
      which.min(log(log(fit.enet1$glmnet.fit$nobs))*
                  log(settings$p[s])*psel - dev)]
    
    sel.cv.enet1 <- which(as.numeric(coef(fit.enet1, "lambda.min"))[-1]!=0)
    sel.aic.enet1 <- which(as.numeric(coef(fit.enet1, lam.aic.enet1))[-1]!=0)
    sel.bic.enet1 <- which(as.numeric(coef(fit.enet1, lam.bic.enet1))[-1]!=0)
    sel.hbic.enet1 <- which(as.numeric(coef(fit.enet1, lam.hbic.enet1))[-1]!=0)
    sel.gic.enet1 <- which(as.numeric(coef(fit.enet1, lam.gic.enet1))[-1]!=0)
    
    jaccard[[s]]$enet1[r, 1] <- 
      length(intersect(sel.cv.enet1, 1:settings$k[s]))/
      length(union(sel.cv.enet1, 1:settings$k[s]))
    jaccard[[s]]$enet1[r, 2] <- 
      length(intersect(sel.aic.enet1, 1:settings$k[s]))/
      length(union(sel.aic.enet1, 1:settings$k[s]))
    jaccard[[s]]$enet1[r, 3] <- 
      length(intersect(sel.bic.enet1, 1:settings$k))/
      length(union(sel.bic.enet1, 1:settings$k))
    jaccard[[s]]$enet1[r, 4] <- 
      length(intersect(sel.hbic.enet1, 1:settings$k))/
      length(union(sel.cv.enet1, 1:settings$k))
    jaccard[[s]]$enet1[r, 5] <- 
      length(intersect(sel.gic.enet1, 1:settings$k))/
      length(union(sel.cv.enet1, 1:settings$k))
    
    kappa[[s]]$enet1[r, 1] <- cohen.kappa(cbind(
      beta!=0, replace(rep(0, settings$p[s]), sel.cv.enet1, 1)))$kappa
    kappa[[s]]$enet1[r, 2] <- cohen.kappa(cbind(
      beta!=0, replace(rep(0, settings$p[s]), sel.aic.enet1, 1)))$kappa
    kappa[[s]]$enet1[r, 3] <- cohen.kappa(cbind(
      beta!=0, replace(rep(0, settings$p[s]), sel.bic.enet1, 1)))$kappa
    kappa[[s]]$enet1[r, 4] <- cohen.kappa(cbind(
      beta!=0, replace(rep(0, settings$p[s]), sel.hbic.enet1, 1)))$kappa
    kappa[[s]]$enet1[r, 5] <- cohen.kappa(cbind(
      beta!=0, replace(rep(0, settings$p[s]), sel.gic.enet1, 1)))$kappa
    
    # enet2
    psel <- fit.enet2$glmnet.fit$df
    dev <- (fit.enet2$glmnet.fit$dev.ratio - 1)*fit.enet2$glmnet.fit$nulldev
    lam.aic.enet2 <- fit.enet2$glmnet.fit$lambda[which.min(2*psel - dev)]
    lam.bic.enet2 <- fit.enet2$glmnet.fit$lambda[
      which.min(log(fit.enet2$glmnet.fit$nobs)*psel - dev)]
    lam.hbic.enet2 <- fit.enet2$glmnet.fit$lambda[
      which.min(log(settings$p[s])*psel - dev)]
    lam.gic.enet2 <- fit.enet2$glmnet.fit$lambda[
      which.min(log(log(fit.enet2$glmnet.fit$nobs))*
                  log(settings$p[s])*psel - dev)]
    
    sel.cv.enet2 <- which(as.numeric(coef(fit.enet2, "lambda.min"))[-1]!=0)
    sel.aic.enet2 <- which(as.numeric(coef(fit.enet2, lam.aic.enet2))[-1]!=0)
    sel.bic.enet2 <- which(as.numeric(coef(fit.enet2, lam.bic.enet2))[-1]!=0)
    sel.hbic.enet2 <- which(as.numeric(coef(fit.enet2, lam.hbic.enet2))[-1]!=0)
    sel.gic.enet2 <- which(as.numeric(coef(fit.enet2, lam.gic.enet2))[-1]!=0)
    
    jaccard[[s]]$enet2[r, 1] <- 
      length(intersect(sel.cv.enet2, 1:settings$k[s]))/
      length(union(sel.cv.enet2, 1:settings$k[s]))
    jaccard[[s]]$enet2[r, 2] <- 
      length(intersect(sel.aic.enet2, 1:settings$k[s]))/
      length(union(sel.aic.enet2, 1:settings$k[s]))
    jaccard[[s]]$enet2[r, 3] <- 
      length(intersect(sel.bic.enet2, 1:settings$k))/
      length(union(sel.bic.enet2, 1:settings$k))
    jaccard[[s]]$enet2[r, 4] <- 
      length(intersect(sel.hbic.enet2, 1:settings$k))/
      length(union(sel.cv.enet2, 1:settings$k))
    jaccard[[s]]$enet2[r, 5] <- 
      length(intersect(sel.gic.enet2, 1:settings$k))/
      length(union(sel.cv.enet2, 1:settings$k))
    
    kappa[[s]]$enet2[r, 1] <- cohen.kappa(cbind(
      beta!=0, replace(rep(0, settings$p[s]), sel.cv.enet2, 1)))$kappa
    kappa[[s]]$enet2[r, 2] <- cohen.kappa(cbind(
      beta!=0, replace(rep(0, settings$p[s]), sel.aic.enet2, 1)))$kappa
    kappa[[s]]$enet2[r, 3] <- cohen.kappa(cbind(
      beta!=0, replace(rep(0, settings$p[s]), sel.bic.enet2, 1)))$kappa
    kappa[[s]]$enet2[r, 4] <- cohen.kappa(cbind(
      beta!=0, replace(rep(0, settings$p[s]), sel.hbic.enet2, 1)))$kappa
    kappa[[s]]$enet2[r, 5] <- cohen.kappa(cbind(
      beta!=0, replace(rep(0, settings$p[s]), sel.gic.enet2, 1)))$kappa
    
    # enet3
    psel <- fit.enet3$glmnet.fit$df
    dev <- (fit.enet3$glmnet.fit$dev.ratio - 1)*fit.enet3$glmnet.fit$nulldev
    lam.aic.enet3 <- fit.enet3$glmnet.fit$lambda[which.min(2*psel - dev)]
    lam.bic.enet3 <- fit.enet3$glmnet.fit$lambda[
      which.min(log(fit.enet3$glmnet.fit$nobs)*psel - dev)]
    lam.hbic.enet3 <- fit.enet3$glmnet.fit$lambda[
      which.min(log(settings$p[s])*psel - dev)]
    lam.gic.enet3 <- fit.enet3$glmnet.fit$lambda[
      which.min(log(log(fit.enet3$glmnet.fit$nobs))*
                  log(settings$p[s])*psel - dev)]
    
    sel.cv.enet3 <- which(as.numeric(coef(fit.enet3, "lambda.min"))[-1]!=0)
    sel.aic.enet3 <- which(as.numeric(coef(fit.enet3, lam.aic.enet3))[-1]!=0)
    sel.bic.enet3 <- which(as.numeric(coef(fit.enet3, lam.bic.enet3))[-1]!=0)
    sel.hbic.enet3 <- which(as.numeric(coef(fit.enet3, lam.hbic.enet3))[-1]!=0)
    sel.gic.enet3 <- which(as.numeric(coef(fit.enet3, lam.gic.enet3))[-1]!=0)
    
    jaccard[[s]]$enet3[r, 1] <- 
      length(intersect(sel.cv.enet3, 1:settings$k[s]))/
      length(union(sel.cv.enet3, 1:settings$k[s]))
    jaccard[[s]]$enet3[r, 2] <- 
      length(intersect(sel.aic.enet3, 1:settings$k[s]))/
      length(union(sel.aic.enet3, 1:settings$k[s]))
    jaccard[[s]]$enet3[r, 3] <- 
      length(intersect(sel.bic.enet3, 1:settings$k))/
      length(union(sel.bic.enet3, 1:settings$k))
    jaccard[[s]]$enet3[r, 4] <- 
      length(intersect(sel.hbic.enet3, 1:settings$k))/
      length(union(sel.cv.enet3, 1:settings$k))
    jaccard[[s]]$enet3[r, 5] <- 
      length(intersect(sel.gic.enet3, 1:settings$k))/
      length(union(sel.cv.enet3, 1:settings$k))
    
    kappa[[s]]$enet3[r, 1] <- cohen.kappa(cbind(
      beta!=0, replace(rep(0, settings$p[s]), sel.cv.enet3, 1)))$kappa
    kappa[[s]]$enet3[r, 2] <- cohen.kappa(cbind(
      beta!=0, replace(rep(0, settings$p[s]), sel.aic.enet3, 1)))$kappa
    kappa[[s]]$enet3[r, 3] <- cohen.kappa(cbind(
      beta!=0, replace(rep(0, settings$p[s]), sel.bic.enet3, 1)))$kappa
    kappa[[s]]$enet3[r, 4] <- cohen.kappa(cbind(
      beta!=0, replace(rep(0, settings$p[s]), sel.hbic.enet3, 1)))$kappa
    kappa[[s]]$enet3[r, 5] <- cohen.kappa(cbind(
      beta!=0, replace(rep(0, settings$p[s]), sel.gic.enet3, 1)))$kappa
  }
}

results <- list(kappa, jaccard)
save(results, file=paste(path.res, "penalty_selection_sim1.Rdata", sep=""))

### graphs
load(paste(path.res, "penalty_selection_sim1.Rdata", sep=""))
for(s in 1:nrow(settings)) {
  png(paste(path.graph, "penalty_selection_sim1_rho", settings$rho[s], "_n",
            settings$n[s], "_k", settings$k[s], ".png", sep=""), units="in", 
      width=12, height=6, res=120)
  par(mfrow=c(2, 4))
  boxplot(kappa[[s]]$enet1, main=expression(alpha==0.05), ylab="Cohen's kappa")
  boxplot(kappa[[s]]$enet2, main=expression(alpha==0.5), ylab="Cohen's kappa")
  boxplot(kappa[[s]]$enet3, main=expression(alpha==0.95), ylab="Cohen's kappa")
  boxplot(kappa[[s]]$lasso, main=expression(alpha==1), ylab="Cohen's kappa")
  
  boxplot(jaccard[[s]]$enet1, main=expression(alpha==0.05), 
          ylab="Jaccard index")
  boxplot(jaccard[[s]]$enet2, main=expression(alpha==0.5), ylab="Jaccard index")
  boxplot(jaccard[[s]]$enet3, main=expression(alpha==0.95), 
          ylab="Jaccard index")
  boxplot(jaccard[[s]]$lasso, main=expression(alpha==1), ylab="Jaccard index")
  dev.off()
  
}



