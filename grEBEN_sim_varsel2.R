##############################  preamble  #############################
# Simulations to check variable selection                             #
# version: 01                                                         #
# author: Magnus Munch                                                #
# created: 12-09-2017                                                 #
# last edited: 12-09-2017                                             #
#######################################################################

###############################  notes  ###############################
#                                                                     #
#######################################################################

# paths
path.code <- as.character(ifelse(Sys.info()[1]=="Darwin","/Users/magnusmunch/Documents/PhD/EBEN/code/" ,"~/EBEN/code/"))
path.graph <- "/Users/magnusmunch/Documents/PhD/EBEN/graphs/"
path.res <- as.character(ifelse(Sys.info()[1]=="Darwin","/Users/magnusmunch/Documents/PhD/EBEN/results/" ,"~/EBEN/results/"))

### libraries
library(foreach)
library(psych)
library(doMC)
library(mvtnorm)
library(pROC)

### functions
# source function for variational Bayes
source(paste(path.code, "grVBEM.R", sep=""))

### setting seed for reproducibility
set.seed(1001)

# which methods are we comparing
methods <- c("enet", "enet+sel", "grEBEN+sel", "ridge", "GRridge", "GRridge+sel")

# dependent data characteristics
n <- 100
m <- 1
ntest <- 1000

# predictor variable characteristics
sigma2 <- 1
rhow <- 0.3
rhob <- 0.7

# model parameter characteristics
p <- 1000
G <- 10
f <- 1.3
beta.mean <- 0.1

# simulation and estimation characteristics
nreps <- 50
pselmin <- 10
pselmax <- 150
pselstep <- 10

# combining the simulation settings
set <- expand.grid(n=n, m=m, ntest=ntest, p=p, G=G, f=f, beta.mean=beta.mean, rhow=rhow, rhob=rhob,
                   sigma2=sigma2, pselmin=pselmin, pselmax=pselmax, pselstep=pselstep)

# loop over the settings
for(s in 1:nrow(set)) {
  
  beta <- as.numeric(sapply(1:set$G[s], function(g) {
    rep(c(0, c(0, set$beta.mean[s]*set$p[s]*set$f[s]^c(0:(set$G[s] - 2))/
                 sum(set$p[s]*set$f[s]^c(0:(set$G[s] - 2))*c(1:(set$G[s] - 1))/set$G[s]^2))[g]), 
        times=round(c(((1 - c(0:(set$G[s] - 1))/set$G[s])*set$p[s]/set$G[s])[g], 
                      (set$p[s]*c(0:(set$G[s] - 1))/set$G[s]^2)[g])))}))
  
  Sigma.group <- matrix(c(rep(rep(c(set$sigma2[s], set$rhow[s]), times=c(1, set$p[s]/set$G[s])), 
                              times=set$p[s]/set$G[s] - 1), set$sigma2[s]), 
                        ncol=set$p[s]/set$G[s], nrow=set$p[s]/set$G[s])
  Sigma.block <- as.matrix(bdiag(Sigma.group, Sigma.group))
  diag(Sigma.block[(set$p[s]/set$G[s] + 1):(1.5*set$p[s]/set$G[s]), 1:(set$p[s]/(2*set$G[s]))]) <- 
    diag(Sigma.block[1:(set$p[s]/(2*set$G[s])), (set$p[s]/set$G[s] + 1):(1.5*set$p[s]/set$G[s])]) <- 
    set$rhob[s]
  partitions1 <- list(groups=rep(1:set$G[s], each=set$p[s]/set$G[s]))
  partitions2 <- list(groups=CreatePartition(as.factor(partitions1$groups)))
  
  # loop over the repetitions of the simulation
  pselseq <- seq(pselmin, pselmax, by=pselstep)
  aucmat <- pselmat <- kappamat <- msemat <- precmat <- recmat <- f1mat <-
    matrix(NA, nrow=nreps*length(pselseq), ncol=length(methods) + 2)
  colnames(aucmat) <- colnames(pselmat) <- colnames(kappamat) <- colnames(msemat) <- colnames(precmat) <- 
    colnames(recmat) <- colnames(f1mat) <- c("rep", "psel", methods)
  for(crep in 1:nreps) {
    # print iteration info
    print(paste("setting ", s, ", repetition ", crep, sep=""))
    
    # create data
    x <- matrix(NA, ncol=set$p[s], nrow=set$n[s])
    for(b in 1:(set$G[s]/2)) {
      x.block <- rmvnorm(set$n[s], sigma=Sigma.block)
      x[, ((b - 1)*set$p[s]/set$G[s] + 1):(b*set$p[s]/set$G[s])] <- 
        x.block[, c(1:(set$p[s]/set$G[s]))]
      x[, c((set$p[s] - b*set$p[s]/set$G[s] + 1):(set$p[s] - (b - 1)*set$p[s]/set$G[s]))] <- 
        x.block[, c((set$p[s]/set$G[s] + 1):(2*set$p[s]/set$G[s]))]
    }
    y <- rbinom(set$n[s], rep(set$m[s], set$n[s]), as.numeric(exp(x %*% beta)/(1 + exp(x %*% beta))))
    
    xtest <- matrix(NA, ncol=set$p[s], nrow=set$ntest[s])
    for(b in 1:(set$G[s]/2)) {
      x.block <- rmvnorm(set$ntest[s], sigma=Sigma.block)
      xtest[, ((b - 1)*set$p[s]/set$G[s] + 1):(b*set$p[s]/set$G[s])] <- 
        x.block[, c(1:(set$p[s]/set$G[s]))]
      xtest[, c((set$p[s] - b*set$p[s]/set$G[s] + 1):(set$p[s] - (b - 1)*set$p[s]/set$G[s]))] <- 
        x.block[, c((set$p[s]/set$G[s] + 1):(2*set$p[s]/set$G[s]))]
    }
    ytest <- rbinom(set$ntest[s], rep(set$m[s], set$ntest[s]), as.numeric(exp(xtest %*% beta)/(1 + exp(xtest %*% beta))))
    
    
    # fit the methods that do not select variables
    fit.GRridge1 <- grridge(t(x), y, partitions=partitions2, monotone=FALSE, trace=FALSE)
    fit.ridge <- penalized(y, x, unpenalized=~1, lambda1=0, lambda2=fit.GRridge1$optl, 
                           model="logistic")
    fit.enet1 <- cv.pen(x, y, unpenalized=NULL, intercept=TRUE, psel=NULL)
    fit.enet2 <- glmnet(x, y, family="binomial", alpha=fit.enet1$alpha[which.min(fit.enet1$cvll)], 
                        lambda=fit.enet1$lambda[which.min(fit.enet1$cvll)], standardize=FALSE, intercept=TRUE)
    
    ### number of cores to use in parallel operations
    if(Sys.info()[1]=="Darwin") {
      registerDoMC(2)
    } else {
      registerDoMC(length(pselseq))
    }
    # loop over the different estimations of the models
    result <- foreach(csel=c(1:length(pselseq)), .combine=rbind) %dopar% {
      psel <- pselseq[csel]
      
      # variable selection methods
      fit.grEBEN <- grEBEN(x, y, rep(set$m[s], set$n[s]), unpenalized=NULL, intercept=TRUE, partitions=partitions1, lambda1=NULL, lambda2=NULL,
                           monotone=list(FALSE, FALSE), psel=psel, posterior=FALSE, ELBO=FALSE, eps=0.001, maxiter=500, trace=FALSE)
      fit.GRridge2 <- grridge(t(x), y, partitions=partitions2, monotone=FALSE, optl=fit.GRridge1$optl,
                              selectionEN=TRUE, maxsel=psel, trace=FALSE)
      test.GRridge <- mygrridge(t(x), y, partitions=partitions2, monotone=FALSE, optl=fit.GRridge1$optl,
                              selectionEN=TRUE, maxsel=psel, trace=FALSE)
      
      # calculate different metrics
      auc.enet <- pROC::roc(ytest, as.numeric(predict.glmnet(fit.enet2, xtest, type="response")))$auc
      auc.enetsel <- pROC::roc(ytest, predict.grEBEN(fit.grEBEN, xtest, unpenalized=NULL, type="nogroups"))$auc
      auc.grEBENsel <- pROC::roc(ytest, predict.grEBEN(fit.grEBEN, xtest, unpenalized=NULL, type="selection"))$auc
      auc.ridge <- pROC::roc(ytest, predict.grridge(fit.GRridge1, t(xtest))[, 1])$auc
      auc.GRridge <- pROC::roc(ytest, predict.grridge(fit.GRridge1, t(xtest))[, 2])$auc
      auc.GRridgesel <- pROC::roc(ytest, predict.grridge(fit.GRridge2, t(xtest))[, 3])$auc
      
      psel.enet <- fit.enet2$df
      psel.enetsel <- sum(fit.grEBEN$beta.nogroups[-1]!=0)
      psel.grEBENsel <- sum(fit.grEBEN$beta.sel[-1]!=0)
      psel.ridge <- set$p[s]
      psel.GRridge <- set$p[s]
      psel.GRridgesel <- length(fit.GRridge2$resEN$whichEN)
      
      kappa.enet <- cohen.kappa(cbind(as.numeric(beta!=0), as.numeric(coef(fit.enet2) > 0)[-1]))$kappa
      kappa.enetsel <- cohen.kappa(cbind(as.numeric(beta!=0), as.numeric(fit.grEBEN$beta.nogroups > 0)[-1]))$kappa
      kappa.grEBENsel <- cohen.kappa(cbind(as.numeric(beta!=0), as.numeric(fit.grEBEN$beta.sel > 0)[-1]))$kappa
      kappa.ridge <- 0
      kappa.GRridge <- 0
      kappa.GRridgesel <- cohen.kappa(cbind(as.numeric(beta!=0), as.numeric(c(1:set$p[s]) %in% fit.GRridge2$resEN$whichEN)))$kappa
      
      mse.enet <- mean((as.numeric(coef(fit.enet2))[-1] - beta)^2)
      mse.enetsel <- mean((fit.grEBEN$beta.nogroups[-1] - beta)^2)
      mse.grEBENsel <- mean((fit.grEBEN$beta.sel[-1] - beta)^2)
      mse.ridge <- mean((coef(fit.ridge)[-1] - beta)^2)
      mse.GRridge <- mean((fit.GRridge1$betas - beta)^2)
      mse.GRridgesel <- mean(c(fit.GRridge2$resEN$betasEN - beta[fit.GRridge2$resEN$whichEN], 
                               beta[-fit.GRridge2$resEN$whichEN])^2)
      
      prec.enet <- sum((coef(fit.enet2) > 0)[-1] & (beta!=0))/sum((coef(fit.enet2) > 0)[-1])
      prec.enetsel <- sum((fit.grEBEN$beta.nogroups > 0)[-1] & (beta!=0))/sum((fit.grEBEN$beta.nogroups > 0)[-1])
      prec.grEBENsel <- sum((fit.grEBEN$beta.sel > 0)[-1] & (beta!=0))/sum((fit.grEBEN$beta.sel > 0)[-1])
      prec.ridge <- sum(beta!=0)/set$p[s]
      prec.GRridge <- sum(beta!=0)/set$p[s]
      prec.GRridgesel <- sum((c(1:set$p[s]) %in% fit.GRridge2$resEN$whichEN) & (beta!=0))/
        sum(c(1:set$p[s]) %in% fit.GRridge2$resEN$whichEN)
      
      rec.enet <- sum((coef(fit.enet2) > 0)[-1] & (beta!=0))/sum(beta!=0)
      rec.enetsel <- sum((fit.grEBEN$beta.nogroups > 0)[-1] & (beta!=0))/sum(beta!=0)
      rec.grEBENsel <- sum((fit.grEBEN$beta.sel > 0)[-1] & (beta!=0))/sum(beta!=0)
      rec.ridge <- 1
      rec.GRridge <- 1
      rec.GRridgesel <- sum((c(1:set$p[s]) %in% fit.GRridge2$resEN$whichEN) & (beta!=0))/sum(beta!=0)
      
      f1.enet <- 2*prec.enet*rec.enet/(prec.enet + rec.enet)
      f1.enetsel <- 2*prec.enetsel*rec.enetsel/(prec.enetsel + rec.enetsel)
      f1.grEBENsel <- 2*prec.grEBENsel*rec.grEBENsel/(prec.grEBENsel + rec.grEBENsel)
      f1.ridge <- 2*prec.ridge*rec.ridge/(prec.ridge + rec.ridge)
      f1.GRridge <- 2*prec.GRridge*rec.GRridge/(prec.GRridge + rec.GRridge)
      f1.GRridgesel <- 2*prec.GRridgesel*rec.GRridgesel/(prec.GRridgesel + rec.GRridgesel)
      
      # return the vector of calculated metrics
      return(c(auc.enet, auc.enetsel, auc.grEBENsel,
               auc.ridge, auc.GRridge, auc.GRridgesel, psel.enet, psel.enetsel, psel.grEBENsel,
               psel.ridge, psel.GRridge, psel.GRridgesel, kappa.enet, kappa.enetsel, kappa.grEBENsel,
               kappa.ridge, kappa.GRridge, kappa.GRridgesel, mse.enet, mse.enetsel, mse.grEBENsel,
               mse.ridge, mse.GRridge, mse.GRridgesel, prec.enet, prec.enetsel, prec.grEBENsel,
               prec.ridge, prec.GRridge, prec.GRridgesel, rec.enet, rec.enetsel, rec.grEBENsel,
               rec.ridge, rec.GRridge, rec.GRridgesel, f1.enet, f1.enetsel, f1.grEBENsel,
               f1.ridge, f1.GRridge, f1.GRridgesel))
    }
    
    # assigning to matrices in parallel case
    aucmat[((crep - 1)*length(pselseq) + 1):(crep*length(pselseq)), ] <- cbind(crep, pselseq, result[, c(1:length(methods))])
    pselmat[((crep - 1)*length(pselseq) + 1):(crep*length(pselseq)), ] <- cbind(crep, pselseq, result[, c((length(methods) + 1):(2*length(methods)))])
    kappamat[((crep - 1)*length(pselseq) + 1):(crep*length(pselseq)), ] <- cbind(crep, pselseq, result[, c((2*length(methods) + 1):(3*length(methods)))])
    msemat[((crep - 1)*length(pselseq) + 1):(crep*length(pselseq)), ] <- cbind(crep, pselseq, result[, c((3*length(methods) + 1):(4*length(methods)))])
    precmat[((crep - 1)*length(pselseq) + 1):(crep*length(pselseq)), ] <- cbind(crep, pselseq, result[, c((4*length(methods) + 1):(5*length(methods)))])
    recmat[((crep - 1)*length(pselseq) + 1):(crep*length(pselseq)), ] <- cbind(crep, pselseq, result[, c((5*length(methods) + 1):(6*length(methods)))])
    f1mat[((crep - 1)*length(pselseq) + 1):(crep*length(pselseq)), ] <- cbind(crep, pselseq, result[, c((6*length(methods) + 1):(7*length(methods)))])
    
    # saving results at every repetition in case of errors
    assign(paste("res.setting", s, sep=""), list(auc=aucmat, psel=pselmat, kappa=kappamat, mse=msemat,
                                                 precision=precmat, recall=recmat, f1score=f1mat))
    save(list=paste("res.setting", s, sep=""), file=paste(path.res, "grEBEN_sim_varsel2_setting", s, ".Rdata", sep=""))
    
  }
  # remove results from setting to avoid memory overflow
  rm(list=paste("res.setting", s, sep=""))
  
}

# ### creating graphs
# # loading libraries for graphs
# library(ggplot2)
# library(gridExtra)
# 
# ### setting1 
# # loading the data and inspecting the variable selections
# load(paste(path.res, "grEBEN_sim_varsel1_setting1.Rdata", sep=""))
# unique(res.setting1$auc[, 1][!is.na(res.setting1$auc[, 1])])
# res.setting1$psel[!(res.setting1$psel[, 5] %in% pselseq) & !is.na(res.setting1$psel[, 5]), c(1, 2, 5)]
# pselseq <- sort(as.numeric(na.omit(unique(res.setting1$auc[, 2]))))
# methods <- colnames(res.setting1$auc)[-c(1:2)]
# res.setting1 <- lapply(res.setting1, function(met) {met[(res.setting1$psel[, 5] %in% pselseq) & !is.na(res.setting1$psel[, 5]), ]})
# 
# # creating the plotting tables
# plot.auc <- data.frame(psel=pselseq, t(sapply(pselseq, function(psel) {
#   apply(res.setting1$auc[, -c(1:2)][ifelse(is.na(res.setting1$auc[, 2]==psel), FALSE, res.setting1$auc[, 2]==psel), ], 2, function(vars) {
#     c(median(vars, na.rm=TRUE), mean(vars, na.rm=TRUE), sd(vars, na.rm=TRUE))})})))
# colnames(plot.auc) <- c("psel", as.vector(t(outer(methods, c("median", "mean", "sd"), paste, sep="."))))
# plot.kappa <- data.frame(psel=pselseq, t(sapply(pselseq, function(psel) {
#   apply(res.setting1$kappa[, -c(1:2)][ifelse(is.na(res.setting1$kappa[, 2]==psel), FALSE, res.setting1$kappa[, 2]==psel), ], 2, function(vars) {
#     c(median(vars, na.rm=TRUE), mean(vars, na.rm=TRUE), sd(vars, na.rm=TRUE))})})))
# colnames(plot.kappa) <- c("psel", as.vector(t(outer(methods, c("median", "mean", "sd"), paste, sep="."))))
# plot.mse <- data.frame(psel=pselseq, t(sapply(pselseq, function(psel) {
#   apply(res.setting1$mse[, -c(1:2)][ifelse(is.na(res.setting1$mse[, 2]==psel), FALSE, res.setting1$mse[, 2]==psel), ], 2, function(vars) {
#     c(median(vars, na.rm=TRUE), mean(vars, na.rm=TRUE), sd(vars, na.rm=TRUE))})})))
# colnames(plot.mse) <- c("psel", as.vector(t(outer(methods, c("median", "mean", "sd"), paste, sep="."))))
# plot.prec <- data.frame(psel=pselseq, t(sapply(pselseq, function(psel) {
#   apply(res.setting1$prec[, -c(1:2)][ifelse(is.na(res.setting1$prec[, 2]==psel), FALSE, res.setting1$prec[, 2]==psel), ], 2, function(vars) {
#     c(median(vars, na.rm=TRUE), mean(vars, na.rm=TRUE), sd(vars, na.rm=TRUE))})})))
# colnames(plot.prec) <- c("psel", as.vector(t(outer(methods, c("median", "mean", "sd"), paste, sep="."))))
# plot.rec <- data.frame(psel=pselseq, t(sapply(pselseq, function(psel) {
#   apply(res.setting1$rec[, -c(1:2)][ifelse(is.na(res.setting1$rec[, 2]==psel), FALSE, res.setting1$rec[, 2]==psel), ], 2, function(vars) {
#     c(median(vars, na.rm=TRUE), mean(vars, na.rm=TRUE), sd(vars, na.rm=TRUE))})})))
# colnames(plot.rec) <- c("psel", as.vector(t(outer(methods, c("median", "mean", "sd"), paste, sep="."))))
# plot.f1 <- data.frame(psel=pselseq, t(sapply(pselseq, function(psel) {
#   apply(res.setting1$f1[, -c(1:2)][ifelse(is.na(res.setting1$f1[, 2]==psel), FALSE, res.setting1$f1[, 2]==psel), ], 2, function(vars) {
#     c(median(vars, na.rm=TRUE), mean(vars, na.rm=TRUE), sd(vars, na.rm=TRUE))})})))
# colnames(plot.f1) <- c("psel", as.vector(t(outer(methods, c("median", "mean", "sd"), paste, sep="."))))
# 
# # plots of the medians over simulations
# png(paste(path.graph, "grEBEN_sim_varsel1_set1_median.png", sep=""), width=1200, height=720, res=90)
# par(mfrow=c(2, 3))
# plot(plot.auc[, 1], plot.auc[, 2], type="l", col=2, lty=1, ylim=range(plot.auc[, c(2, 5, 8, 11, 14, 17)]),
#      xlab="Number of selected variables", ylab="AUC", main="a)")
# lines(plot.auc[, 1], plot.auc[, 5], col=3, lty=1)
# lines(plot.auc[, 1], plot.auc[, 8], col=4, lty=1)
# lines(plot.auc[, 1], plot.auc[, 11], col=5, lty=1)
# lines(plot.auc[, 1], plot.auc[, 14], col=6, lty=1)
# lines(plot.auc[, 1], plot.auc[, 17], col=7, lty=1)
# abline(v=100, lty=3)
# legend("bottomright", legend=c(methods, "Active variables"), lty=c(rep(1, 6), 3), col=c(c(2:7), 1))
# 
# plot(plot.kappa[, 1], plot.kappa[, 2], type="l", col=2, lty=1, ylim=range(plot.kappa[, c(2, 5, 8, 11, 14, 17)]),
#      xlab="Number of selected variables", ylab="kappa", main="b)")
# lines(plot.kappa[, 1], plot.kappa[, 5], col=3, lty=1)
# lines(plot.kappa[, 1], plot.kappa[, 8], col=4, lty=1)
# lines(plot.kappa[, 1], plot.kappa[, 11], col=5, lty=1)
# lines(plot.kappa[, 1], plot.kappa[, 14], col=6, lty=1)
# lines(plot.kappa[, 1], plot.kappa[, 17], col=7, lty=1)
# abline(v=100, lty=3)
# 
# plot(plot.mse[, 1], plot.mse[, 2], type="l", col=2, lty=1, ylim=range(plot.mse[, c(2, 5, 8, 11, 14, 17)]),
#      xlab="Number of selected variables", ylab="mse", main="c)")
# lines(plot.mse[, 1], plot.mse[, 5], col=3, lty=1)
# lines(plot.mse[, 1], plot.mse[, 8], col=4, lty=1)
# lines(plot.mse[, 1], plot.mse[, 11], col=5, lty=1)
# lines(plot.mse[, 1], plot.mse[, 14], col=6, lty=1)
# lines(plot.mse[, 1], plot.mse[, 17], col=7, lty=1)
# abline(v=100, lty=3)
# 
# plot(plot.prec[, 1], plot.prec[, 2], type="l", col=2, lty=1, ylim=range(plot.prec[, c(2, 5, 8, 11, 14, 17)]),
#      xlab="Number of selected variables", ylab="precision", main="d)")
# lines(plot.prec[, 1], plot.prec[, 5], col=3, lty=1)
# lines(plot.prec[, 1], plot.prec[, 8], col=4, lty=1)
# lines(plot.prec[, 1], plot.prec[, 11], col=5, lty=1)
# lines(plot.prec[, 1], plot.prec[, 14], col=6, lty=1)
# lines(plot.prec[, 1], plot.prec[, 17], col=7, lty=1)
# abline(v=100, lty=3)
# 
# plot(plot.rec[, 1], plot.rec[, 2], type="l", col=2, lty=1, ylim=range(plot.rec[, c(2, 5, 8, 11, 14, 17)]),
#      xlab="Number of selected variables", ylab="recall", main="e)")
# lines(plot.rec[, 1], plot.rec[, 5], col=3, lty=1)
# lines(plot.rec[, 1], plot.rec[, 8], col=4, lty=1)
# lines(plot.rec[, 1], plot.rec[, 11], col=5, lty=1)
# lines(plot.rec[, 1], plot.rec[, 14], col=6, lty=1)
# lines(plot.rec[, 1], plot.rec[, 17], col=7, lty=1)
# abline(v=100, lty=3)
# 
# plot(plot.f1[, 1], plot.f1[, 2], type="l", col=2, lty=1, ylim=range(plot.f1[, c(2, 5, 8, 11, 14, 17)]),
#      xlab="Number of selected variables", ylab="f1", main="f)")
# lines(plot.f1[, 1], plot.f1[, 5], col=3, lty=1)
# lines(plot.f1[, 1], plot.f1[, 8], col=4, lty=1)
# lines(plot.f1[, 1], plot.f1[, 11], col=5, lty=1)
# lines(plot.f1[, 1], plot.f1[, 14], col=6, lty=1)
# lines(plot.f1[, 1], plot.f1[, 17], col=7, lty=1)
# abline(v=100, lty=3)
# dev.off()
# 
# # plots of the means plus error bars
# plot.auc2 <- reshape(plot.auc, varying=list(paste(methods, "median", sep="."), paste(methods, "mean", sep="."),
#                                             paste(methods, "sd", sep=".")), direction="long",
#                      times=methods, v.names=c("median", "mean", "sd"), idvar="psel")
# rownames(plot.auc2) <- NULL
# colnames(plot.auc2)[2] <- "method"
# plot.kappa2 <- reshape(plot.kappa, varying=list(paste(methods, "median", sep="."), paste(methods, "mean", sep="."),
#                                                 paste(methods, "sd", sep=".")), direction="long",
#                        times=methods, v.names=c("median", "mean", "sd"), idvar="psel")
# rownames(plot.kappa2) <- NULL
# colnames(plot.kappa2)[2] <- "method"
# plot.mse2 <- reshape(plot.mse, varying=list(paste(methods, "median", sep="."), paste(methods, "mean", sep="."),
#                                             paste(methods, "sd", sep=".")), direction="long",
#                      times=methods, v.names=c("median", "mean", "sd"), idvar="psel")
# rownames(plot.mse2) <- NULL
# colnames(plot.mse2)[2] <- "method"
# plot.prec2 <- reshape(plot.prec, varying=list(paste(methods, "median", sep="."), paste(methods, "mean", sep="."),
#                                               paste(methods, "sd", sep=".")), direction="long",
#                       times=methods, v.names=c("median", "mean", "sd"), idvar="psel")
# rownames(plot.prec2) <- NULL
# colnames(plot.prec2)[2] <- "method"
# plot.rec2 <- reshape(plot.rec, varying=list(paste(methods, "median", sep="."), paste(methods, "mean", sep="."),
#                                             paste(methods, "sd", sep=".")), direction="long",
#                      times=methods, v.names=c("median", "mean", "sd"), idvar="psel")
# rownames(plot.rec2) <- NULL
# colnames(plot.rec2)[2] <- "method"
# plot.f12 <- reshape(plot.f1, varying=list(paste(methods, "median", sep="."), paste(methods, "mean", sep="."),
#                                           paste(methods, "sd", sep=".")), direction="long",
#                     times=methods, v.names=c("median", "mean", "sd"), idvar="psel")
# rownames(plot.f12) <- NULL
# colnames(plot.f12)[2] <- "method"
# 
# mean.auc <- ggplot(plot.auc2, aes(x=psel, y=mean, colour=method)) +
#   geom_vline(xintercept=100, linetype="dotted", show.legend=TRUE) +
#   geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=1) + geom_line() + geom_point() +
#   xlab("Number of selected variables") + ylab("AUC") + ggtitle("a)") + theme_bw() +
#   theme(legend.position="none", plot.title=element_text(hjust=0.5))
# 
# mean.kappa <- ggplot(plot.kappa2, aes(x=psel, y=mean, colour=method)) +
#   geom_vline(xintercept=100, linetype="dotted", show.legend=TRUE) +
#   geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=1) + geom_line() + geom_point() +
#   xlab("Number of selected variables") + ylab("Cohen's kappa") + ggtitle("b)") + theme_bw() +
#   theme(legend.position="none", plot.title=element_text(hjust=0.5))
# 
# mean.mse <- ggplot(plot.mse2, aes(x=psel, y=mean, colour=method)) +
#   geom_vline(xintercept=100, linetype="dotted", show.legend=TRUE) +
#   geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=1) + geom_line() + geom_point() +
#   xlab("Number of selected variables") + ylab("MSE") + ggtitle("c)") + theme_bw() +
#   theme(legend.position="none", plot.title=element_text(hjust=0.5))
# 
# mean.prec <- ggplot(plot.prec2, aes(x=psel, y=mean, colour=method)) +
#   geom_vline(xintercept=100, linetype="dotted", show.legend=TRUE) +
#   geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=1) + geom_line() + geom_point() +
#   xlab("Number of selected variables") + ylab("Precision") + ggtitle("d)") + theme_bw() +
#   theme(legend.position="none", plot.title=element_text(hjust=0.5))
# 
# mean.rec <- ggplot(plot.rec2, aes(x=psel, y=mean, colour=method)) +
#   geom_vline(xintercept=100, linetype="dotted", show.legend=TRUE) +
#   geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=1) + geom_line() + geom_point() +
#   xlab("Number of selected variables") + ylab("Recall") + ggtitle("e)") + theme_bw() +
#   theme(legend.justification=c(0, 1), legend.position=c(0, 1), plot.title=element_text(hjust=0.5))
# 
# mean.f1 <- ggplot(plot.f12, aes(x=psel, y=mean, colour=method)) +
#   geom_vline(xintercept=100, linetype="dotted", show.legend=TRUE) +
#   geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=1) + geom_line() + geom_point() +
#   xlab("Number of selected variables") + ylab("F1 score") + ggtitle("f)") + theme_bw() +
#   theme(legend.position="none", plot.title=element_text(hjust=0.5))
# 
# png(paste(path.graph, "grEBEN_sim_varsel1_set1_mean.png", sep=""), width=1200, height=720, res=90)
# grid.arrange(mean.auc, mean.kappa, mean.mse, mean.prec, mean.rec, mean.f1, ncol=3)
# dev.off()
# 
# 
# ### setting2
# # loading the data and inspecting the variable selections
# load(paste(path.res, "grEBEN_sim_varsel1_setting2.Rdata", sep=""))
# unique(res.setting2$auc[, 1][!is.na(res.setting2$auc[, 1])])
# res.setting2$psel[!(res.setting2$psel[, 5] %in% pselseq) & !is.na(res.setting2$psel[, 5]), c(1, 2, 5)]
# pselseq <- sort(as.numeric(na.omit(unique(res.setting2$auc[, 2]))))
# methods <- colnames(res.setting2$auc)[-c(1:2)]
# res.setting2 <- lapply(res.setting2, function(met) {met[(res.setting2$psel[, 5] %in% pselseq) & !is.na(res.setting2$psel[, 5]), ]})
# 
# # creating the plotting tables
# plot.auc <- data.frame(psel=pselseq, t(sapply(pselseq, function(psel) {
#   apply(res.setting2$auc[, -c(1:2)][ifelse(is.na(res.setting2$auc[, 2]==psel), FALSE, res.setting2$auc[, 2]==psel), ], 2, function(vars) {
#     c(median(vars, na.rm=TRUE), mean(vars, na.rm=TRUE), sd(vars, na.rm=TRUE))})})))
# colnames(plot.auc) <- c("psel", as.vector(t(outer(methods, c("median", "mean", "sd"), paste, sep="."))))
# plot.kappa <- data.frame(psel=pselseq, t(sapply(pselseq, function(psel) {
#   apply(res.setting2$kappa[, -c(1:2)][ifelse(is.na(res.setting2$kappa[, 2]==psel), FALSE, res.setting2$kappa[, 2]==psel), ], 2, function(vars) {
#     c(median(vars, na.rm=TRUE), mean(vars, na.rm=TRUE), sd(vars, na.rm=TRUE))})})))
# colnames(plot.kappa) <- c("psel", as.vector(t(outer(methods, c("median", "mean", "sd"), paste, sep="."))))
# plot.mse <- data.frame(psel=pselseq, t(sapply(pselseq, function(psel) {
#   apply(res.setting2$mse[, -c(1:2)][ifelse(is.na(res.setting2$mse[, 2]==psel), FALSE, res.setting2$mse[, 2]==psel), ], 2, function(vars) {
#     c(median(vars, na.rm=TRUE), mean(vars, na.rm=TRUE), sd(vars, na.rm=TRUE))})})))
# colnames(plot.mse) <- c("psel", as.vector(t(outer(methods, c("median", "mean", "sd"), paste, sep="."))))
# plot.prec <- data.frame(psel=pselseq, t(sapply(pselseq, function(psel) {
#   apply(res.setting2$prec[, -c(1:2)][ifelse(is.na(res.setting2$prec[, 2]==psel), FALSE, res.setting2$prec[, 2]==psel), ], 2, function(vars) {
#     c(median(vars, na.rm=TRUE), mean(vars, na.rm=TRUE), sd(vars, na.rm=TRUE))})})))
# colnames(plot.prec) <- c("psel", as.vector(t(outer(methods, c("median", "mean", "sd"), paste, sep="."))))
# plot.rec <- data.frame(psel=pselseq, t(sapply(pselseq, function(psel) {
#   apply(res.setting2$rec[, -c(1:2)][ifelse(is.na(res.setting2$rec[, 2]==psel), FALSE, res.setting2$rec[, 2]==psel), ], 2, function(vars) {
#     c(median(vars, na.rm=TRUE), mean(vars, na.rm=TRUE), sd(vars, na.rm=TRUE))})})))
# colnames(plot.rec) <- c("psel", as.vector(t(outer(methods, c("median", "mean", "sd"), paste, sep="."))))
# plot.f1 <- data.frame(psel=pselseq, t(sapply(pselseq, function(psel) {
#   apply(res.setting2$f1[, -c(1:2)][ifelse(is.na(res.setting2$f1[, 2]==psel), FALSE, res.setting2$f1[, 2]==psel), ], 2, function(vars) {
#     c(median(vars, na.rm=TRUE), mean(vars, na.rm=TRUE), sd(vars, na.rm=TRUE))})})))
# colnames(plot.f1) <- c("psel", as.vector(t(outer(methods, c("median", "mean", "sd"), paste, sep="."))))
# 
# # plots of the medians over simulations
# png(paste(path.graph, "grEBEN_sim_varsel1_set2_median.png", sep=""), width=1200, height=720, res=90)
# par(mfrow=c(2, 3))
# plot(plot.auc[, 1], plot.auc[, 2], type="l", col=2, lty=1, ylim=range(plot.auc[, c(2, 5, 8, 11, 14, 17)]),
#      xlab="Number of selected variables", ylab="AUC", main="a)")
# lines(plot.auc[, 1], plot.auc[, 5], col=3, lty=1)
# lines(plot.auc[, 1], plot.auc[, 8], col=4, lty=1)
# lines(plot.auc[, 1], plot.auc[, 11], col=5, lty=1)
# lines(plot.auc[, 1], plot.auc[, 14], col=6, lty=1)
# lines(plot.auc[, 1], plot.auc[, 17], col=7, lty=1)
# abline(v=100, lty=3)
# legend("bottomright", legend=c(methods, "Active variables"), lty=c(rep(1, 6), 3), col=c(c(2:7), 1))
# 
# plot(plot.kappa[, 1], plot.kappa[, 2], type="l", col=2, lty=1, ylim=range(plot.kappa[, c(2, 5, 8, 11, 14, 17)]),
#      xlab="Number of selected variables", ylab="kappa", main="b)")
# lines(plot.kappa[, 1], plot.kappa[, 5], col=3, lty=1)
# lines(plot.kappa[, 1], plot.kappa[, 8], col=4, lty=1)
# lines(plot.kappa[, 1], plot.kappa[, 11], col=5, lty=1)
# lines(plot.kappa[, 1], plot.kappa[, 14], col=6, lty=1)
# lines(plot.kappa[, 1], plot.kappa[, 17], col=7, lty=1)
# abline(v=100, lty=3)
# 
# plot(plot.mse[, 1], plot.mse[, 2], type="l", col=2, lty=1, ylim=range(plot.mse[, c(2, 5, 8, 11, 14, 17)]),
#      xlab="Number of selected variables", ylab="mse", main="c)")
# lines(plot.mse[, 1], plot.mse[, 5], col=3, lty=1)
# lines(plot.mse[, 1], plot.mse[, 8], col=4, lty=1)
# lines(plot.mse[, 1], plot.mse[, 11], col=5, lty=1)
# lines(plot.mse[, 1], plot.mse[, 14], col=6, lty=1)
# lines(plot.mse[, 1], plot.mse[, 17], col=7, lty=1)
# abline(v=100, lty=3)
# 
# plot(plot.prec[, 1], plot.prec[, 2], type="l", col=2, lty=1, ylim=range(plot.prec[, c(2, 5, 8, 11, 14, 17)]),
#      xlab="Number of selected variables", ylab="precision", main="d)")
# lines(plot.prec[, 1], plot.prec[, 5], col=3, lty=1)
# lines(plot.prec[, 1], plot.prec[, 8], col=4, lty=1)
# lines(plot.prec[, 1], plot.prec[, 11], col=5, lty=1)
# lines(plot.prec[, 1], plot.prec[, 14], col=6, lty=1)
# lines(plot.prec[, 1], plot.prec[, 17], col=7, lty=1)
# abline(v=100, lty=3)
# 
# plot(plot.rec[, 1], plot.rec[, 2], type="l", col=2, lty=1, ylim=range(plot.rec[, c(2, 5, 8, 11, 14, 17)]),
#      xlab="Number of selected variables", ylab="recall", main="e)")
# lines(plot.rec[, 1], plot.rec[, 5], col=3, lty=1)
# lines(plot.rec[, 1], plot.rec[, 8], col=4, lty=1)
# lines(plot.rec[, 1], plot.rec[, 11], col=5, lty=1)
# lines(plot.rec[, 1], plot.rec[, 14], col=6, lty=1)
# lines(plot.rec[, 1], plot.rec[, 17], col=7, lty=1)
# abline(v=100, lty=3)
# 
# plot(plot.f1[, 1], plot.f1[, 2], type="l", col=2, lty=1, ylim=range(plot.f1[, c(2, 5, 8, 11, 14, 17)]),
#      xlab="Number of selected variables", ylab="f1", main="f)")
# lines(plot.f1[, 1], plot.f1[, 5], col=3, lty=1)
# lines(plot.f1[, 1], plot.f1[, 8], col=4, lty=1)
# lines(plot.f1[, 1], plot.f1[, 11], col=5, lty=1)
# lines(plot.f1[, 1], plot.f1[, 14], col=6, lty=1)
# lines(plot.f1[, 1], plot.f1[, 17], col=7, lty=1)
# abline(v=100, lty=3)
# dev.off()
# 
# # plots of the means plus error bars
# plot.auc2 <- reshape(plot.auc, varying=list(paste(methods, "median", sep="."), paste(methods, "mean", sep="."),
#                                             paste(methods, "sd", sep=".")), direction="long",
#                      times=methods, v.names=c("median", "mean", "sd"), idvar="psel")
# rownames(plot.auc2) <- NULL
# colnames(plot.auc2)[2] <- "method"
# plot.kappa2 <- reshape(plot.kappa, varying=list(paste(methods, "median", sep="."), paste(methods, "mean", sep="."),
#                                                 paste(methods, "sd", sep=".")), direction="long",
#                        times=methods, v.names=c("median", "mean", "sd"), idvar="psel")
# rownames(plot.kappa2) <- NULL
# colnames(plot.kappa2)[2] <- "method"
# plot.mse2 <- reshape(plot.mse, varying=list(paste(methods, "median", sep="."), paste(methods, "mean", sep="."),
#                                             paste(methods, "sd", sep=".")), direction="long",
#                      times=methods, v.names=c("median", "mean", "sd"), idvar="psel")
# rownames(plot.mse2) <- NULL
# colnames(plot.mse2)[2] <- "method"
# plot.prec2 <- reshape(plot.prec, varying=list(paste(methods, "median", sep="."), paste(methods, "mean", sep="."),
#                                               paste(methods, "sd", sep=".")), direction="long",
#                       times=methods, v.names=c("median", "mean", "sd"), idvar="psel")
# rownames(plot.prec2) <- NULL
# colnames(plot.prec2)[2] <- "method"
# plot.rec2 <- reshape(plot.rec, varying=list(paste(methods, "median", sep="."), paste(methods, "mean", sep="."),
#                                             paste(methods, "sd", sep=".")), direction="long",
#                      times=methods, v.names=c("median", "mean", "sd"), idvar="psel")
# rownames(plot.rec2) <- NULL
# colnames(plot.rec2)[2] <- "method"
# plot.f12 <- reshape(plot.f1, varying=list(paste(methods, "median", sep="."), paste(methods, "mean", sep="."),
#                                           paste(methods, "sd", sep=".")), direction="long",
#                     times=methods, v.names=c("median", "mean", "sd"), idvar="psel")
# rownames(plot.f12) <- NULL
# colnames(plot.f12)[2] <- "method"
# 
# mean.auc <- ggplot(plot.auc2, aes(x=psel, y=mean, colour=method)) +
#   geom_vline(xintercept=100, linetype="dotted", show.legend=TRUE) +
#   geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=1) + geom_line() + geom_point() +
#   xlab("Number of selected variables") + ylab("AUC") + ggtitle("a)") + theme_bw() +
#   theme(legend.position="none", plot.title=element_text(hjust=0.5))
# 
# mean.kappa <- ggplot(plot.kappa2, aes(x=psel, y=mean, colour=method)) +
#   geom_vline(xintercept=100, linetype="dotted", show.legend=TRUE) +
#   geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=1) + geom_line() + geom_point() +
#   xlab("Number of selected variables") + ylab("Cohen's kappa") + ggtitle("b)") + theme_bw() +
#   theme(legend.position="none", plot.title=element_text(hjust=0.5))
# 
# mean.mse <- ggplot(plot.mse2, aes(x=psel, y=mean, colour=method)) +
#   geom_vline(xintercept=100, linetype="dotted", show.legend=TRUE) +
#   geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=1) + geom_line() + geom_point() +
#   xlab("Number of selected variables") + ylab("MSE") + ggtitle("c)") + theme_bw() +
#   theme(legend.position="none", plot.title=element_text(hjust=0.5))
# 
# mean.prec <- ggplot(plot.prec2, aes(x=psel, y=mean, colour=method)) +
#   geom_vline(xintercept=100, linetype="dotted", show.legend=TRUE) +
#   geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=1) + geom_line() + geom_point() +
#   xlab("Number of selected variables") + ylab("Precision") + ggtitle("d)") + theme_bw() +
#   theme(legend.position="none", plot.title=element_text(hjust=0.5))
# 
# mean.rec <- ggplot(plot.rec2, aes(x=psel, y=mean, colour=method)) +
#   geom_vline(xintercept=100, linetype="dotted", show.legend=TRUE) +
#   geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=1) + geom_line() + geom_point() +
#   xlab("Number of selected variables") + ylab("Recall") + ggtitle("e)") + theme_bw() +
#   theme(legend.justification=c(0, 1), legend.position=c(0, 1), plot.title=element_text(hjust=0.5))
# 
# mean.f1 <- ggplot(plot.f12, aes(x=psel, y=mean, colour=method)) +
#   geom_vline(xintercept=100, linetype="dotted", show.legend=TRUE) +
#   geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=1) + geom_line() + geom_point() +
#   xlab("Number of selected variables") + ylab("F1 score") + ggtitle("f)") + theme_bw() +
#   theme(legend.position="none", plot.title=element_text(hjust=0.5))
# 
# png(paste(path.graph, "grEBEN_sim_varsel1_set2_mean.png", sep=""), width=1200, height=720, res=90)
# grid.arrange(mean.auc, mean.kappa, mean.mse, mean.prec, mean.rec, mean.f1, ncol=3)
# dev.off()
# 
# ### setting3
# # loading the data and inspecting the variable selections
# load(paste(path.res, "grEBEN_sim_varsel1_setting3.Rdata", sep=""))
# unique(res.setting3$auc[, 1][!is.na(res.setting3$auc[, 1])])
# res.setting3$psel[!(res.setting3$psel[, 5] %in% pselseq) & !is.na(res.setting3$psel[, 5]), c(1, 2, 5)]
# pselseq <- sort(as.numeric(na.omit(unique(res.setting3$auc[, 2]))))
# methods <- colnames(res.setting3$auc)[-c(1:2)]
# res.setting3 <- lapply(res.setting3, function(met) {met[(res.setting3$psel[, 5] %in% pselseq) & !is.na(res.setting3$psel[, 5]), ]})
# 
# # creating the plotting tables
# plot.auc <- data.frame(psel=pselseq, t(sapply(pselseq, function(psel) {
#   apply(res.setting3$auc[, -c(1:2)][ifelse(is.na(res.setting3$auc[, 2]==psel), FALSE, res.setting3$auc[, 2]==psel), ], 2, function(vars) {
#     c(median(vars, na.rm=TRUE), mean(vars, na.rm=TRUE), sd(vars, na.rm=TRUE))})})))
# colnames(plot.auc) <- c("psel", as.vector(t(outer(methods, c("median", "mean", "sd"), paste, sep="."))))
# plot.kappa <- data.frame(psel=pselseq, t(sapply(pselseq, function(psel) {
#   apply(res.setting3$kappa[, -c(1:2)][ifelse(is.na(res.setting3$kappa[, 2]==psel), FALSE, res.setting3$kappa[, 2]==psel), ], 2, function(vars) {
#     c(median(vars, na.rm=TRUE), mean(vars, na.rm=TRUE), sd(vars, na.rm=TRUE))})})))
# colnames(plot.kappa) <- c("psel", as.vector(t(outer(methods, c("median", "mean", "sd"), paste, sep="."))))
# plot.mse <- data.frame(psel=pselseq, t(sapply(pselseq, function(psel) {
#   apply(res.setting3$mse[, -c(1:2)][ifelse(is.na(res.setting3$mse[, 2]==psel), FALSE, res.setting3$mse[, 2]==psel), ], 2, function(vars) {
#     c(median(vars, na.rm=TRUE), mean(vars, na.rm=TRUE), sd(vars, na.rm=TRUE))})})))
# colnames(plot.mse) <- c("psel", as.vector(t(outer(methods, c("median", "mean", "sd"), paste, sep="."))))
# plot.prec <- data.frame(psel=pselseq, t(sapply(pselseq, function(psel) {
#   apply(res.setting3$prec[, -c(1:2)][ifelse(is.na(res.setting3$prec[, 2]==psel), FALSE, res.setting3$prec[, 2]==psel), ], 2, function(vars) {
#     c(median(vars, na.rm=TRUE), mean(vars, na.rm=TRUE), sd(vars, na.rm=TRUE))})})))
# colnames(plot.prec) <- c("psel", as.vector(t(outer(methods, c("median", "mean", "sd"), paste, sep="."))))
# plot.rec <- data.frame(psel=pselseq, t(sapply(pselseq, function(psel) {
#   apply(res.setting3$rec[, -c(1:2)][ifelse(is.na(res.setting3$rec[, 2]==psel), FALSE, res.setting3$rec[, 2]==psel), ], 2, function(vars) {
#     c(median(vars, na.rm=TRUE), mean(vars, na.rm=TRUE), sd(vars, na.rm=TRUE))})})))
# colnames(plot.rec) <- c("psel", as.vector(t(outer(methods, c("median", "mean", "sd"), paste, sep="."))))
# plot.f1 <- data.frame(psel=pselseq, t(sapply(pselseq, function(psel) {
#   apply(res.setting3$f1[, -c(1:2)][ifelse(is.na(res.setting3$f1[, 2]==psel), FALSE, res.setting3$f1[, 2]==psel), ], 2, function(vars) {
#     c(median(vars, na.rm=TRUE), mean(vars, na.rm=TRUE), sd(vars, na.rm=TRUE))})})))
# colnames(plot.f1) <- c("psel", as.vector(t(outer(methods, c("median", "mean", "sd"), paste, sep="."))))
# 
# # plots of the medians over simulations
# png(paste(path.graph, "grEBEN_sim_varsel1_set3_median.png", sep=""), width=1200, height=720, res=90)
# par(mfrow=c(2, 3))
# plot(plot.auc[, 1], plot.auc[, 2], type="l", col=2, lty=1, ylim=range(plot.auc[, c(2, 5, 8, 11, 14, 17)]),
#      xlab="Number of selected variables", ylab="AUC", main="a)")
# lines(plot.auc[, 1], plot.auc[, 5], col=3, lty=1)
# lines(plot.auc[, 1], plot.auc[, 8], col=4, lty=1)
# lines(plot.auc[, 1], plot.auc[, 11], col=5, lty=1)
# lines(plot.auc[, 1], plot.auc[, 14], col=6, lty=1)
# lines(plot.auc[, 1], plot.auc[, 17], col=7, lty=1)
# abline(v=100, lty=3)
# legend("bottomright", legend=c(methods, "Active variables"), lty=c(rep(1, 6), 3), col=c(c(2:7), 1))
# 
# plot(plot.kappa[, 1], plot.kappa[, 2], type="l", col=2, lty=1, ylim=range(plot.kappa[, c(2, 5, 8, 11, 14, 17)]),
#      xlab="Number of selected variables", ylab="kappa", main="b)")
# lines(plot.kappa[, 1], plot.kappa[, 5], col=3, lty=1)
# lines(plot.kappa[, 1], plot.kappa[, 8], col=4, lty=1)
# lines(plot.kappa[, 1], plot.kappa[, 11], col=5, lty=1)
# lines(plot.kappa[, 1], plot.kappa[, 14], col=6, lty=1)
# lines(plot.kappa[, 1], plot.kappa[, 17], col=7, lty=1)
# abline(v=100, lty=3)
# 
# plot(plot.mse[, 1], plot.mse[, 2], type="l", col=2, lty=1, ylim=range(plot.mse[, c(2, 5, 8, 11, 14, 17)]),
#      xlab="Number of selected variables", ylab="mse", main="c)")
# lines(plot.mse[, 1], plot.mse[, 5], col=3, lty=1)
# lines(plot.mse[, 1], plot.mse[, 8], col=4, lty=1)
# lines(plot.mse[, 1], plot.mse[, 11], col=5, lty=1)
# lines(plot.mse[, 1], plot.mse[, 14], col=6, lty=1)
# lines(plot.mse[, 1], plot.mse[, 17], col=7, lty=1)
# abline(v=100, lty=3)
# 
# plot(plot.prec[, 1], plot.prec[, 2], type="l", col=2, lty=1, ylim=range(plot.prec[, c(2, 5, 8, 11, 14, 17)]),
#      xlab="Number of selected variables", ylab="precision", main="d)")
# lines(plot.prec[, 1], plot.prec[, 5], col=3, lty=1)
# lines(plot.prec[, 1], plot.prec[, 8], col=4, lty=1)
# lines(plot.prec[, 1], plot.prec[, 11], col=5, lty=1)
# lines(plot.prec[, 1], plot.prec[, 14], col=6, lty=1)
# lines(plot.prec[, 1], plot.prec[, 17], col=7, lty=1)
# abline(v=100, lty=3)
# 
# plot(plot.rec[, 1], plot.rec[, 2], type="l", col=2, lty=1, ylim=range(plot.rec[, c(2, 5, 8, 11, 14, 17)]),
#      xlab="Number of selected variables", ylab="recall", main="e)")
# lines(plot.rec[, 1], plot.rec[, 5], col=3, lty=1)
# lines(plot.rec[, 1], plot.rec[, 8], col=4, lty=1)
# lines(plot.rec[, 1], plot.rec[, 11], col=5, lty=1)
# lines(plot.rec[, 1], plot.rec[, 14], col=6, lty=1)
# lines(plot.rec[, 1], plot.rec[, 17], col=7, lty=1)
# abline(v=100, lty=3)
# 
# plot(plot.f1[, 1], plot.f1[, 2], type="l", col=2, lty=1, ylim=range(plot.f1[, c(2, 5, 8, 11, 14, 17)]),
#      xlab="Number of selected variables", ylab="f1", main="f)")
# lines(plot.f1[, 1], plot.f1[, 5], col=3, lty=1)
# lines(plot.f1[, 1], plot.f1[, 8], col=4, lty=1)
# lines(plot.f1[, 1], plot.f1[, 11], col=5, lty=1)
# lines(plot.f1[, 1], plot.f1[, 14], col=6, lty=1)
# lines(plot.f1[, 1], plot.f1[, 17], col=7, lty=1)
# abline(v=100, lty=3)
# dev.off()
# 
# # plots of the means plus error bars
# plot.auc2 <- reshape(plot.auc, varying=list(paste(methods, "median", sep="."), paste(methods, "mean", sep="."),
#                                             paste(methods, "sd", sep=".")), direction="long",
#                      times=methods, v.names=c("median", "mean", "sd"), idvar="psel")
# rownames(plot.auc2) <- NULL
# colnames(plot.auc2)[2] <- "method"
# plot.kappa2 <- reshape(plot.kappa, varying=list(paste(methods, "median", sep="."), paste(methods, "mean", sep="."),
#                                                 paste(methods, "sd", sep=".")), direction="long",
#                        times=methods, v.names=c("median", "mean", "sd"), idvar="psel")
# rownames(plot.kappa2) <- NULL
# colnames(plot.kappa2)[2] <- "method"
# plot.mse2 <- reshape(plot.mse, varying=list(paste(methods, "median", sep="."), paste(methods, "mean", sep="."),
#                                             paste(methods, "sd", sep=".")), direction="long",
#                      times=methods, v.names=c("median", "mean", "sd"), idvar="psel")
# rownames(plot.mse2) <- NULL
# colnames(plot.mse2)[2] <- "method"
# plot.prec2 <- reshape(plot.prec, varying=list(paste(methods, "median", sep="."), paste(methods, "mean", sep="."),
#                                               paste(methods, "sd", sep=".")), direction="long",
#                       times=methods, v.names=c("median", "mean", "sd"), idvar="psel")
# rownames(plot.prec2) <- NULL
# colnames(plot.prec2)[2] <- "method"
# plot.rec2 <- reshape(plot.rec, varying=list(paste(methods, "median", sep="."), paste(methods, "mean", sep="."),
#                                             paste(methods, "sd", sep=".")), direction="long",
#                      times=methods, v.names=c("median", "mean", "sd"), idvar="psel")
# rownames(plot.rec2) <- NULL
# colnames(plot.rec2)[2] <- "method"
# plot.f12 <- reshape(plot.f1, varying=list(paste(methods, "median", sep="."), paste(methods, "mean", sep="."),
#                                           paste(methods, "sd", sep=".")), direction="long",
#                     times=methods, v.names=c("median", "mean", "sd"), idvar="psel")
# rownames(plot.f12) <- NULL
# colnames(plot.f12)[2] <- "method"
# 
# mean.auc <- ggplot(plot.auc2, aes(x=psel, y=mean, colour=method)) +
#   geom_vline(xintercept=100, linetype="dotted", show.legend=TRUE) +
#   geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=1) + geom_line() + geom_point() +
#   xlab("Number of selected variables") + ylab("AUC") + ggtitle("a)") + theme_bw() +
#   theme(legend.position="none", plot.title=element_text(hjust=0.5))
# 
# mean.kappa <- ggplot(plot.kappa2, aes(x=psel, y=mean, colour=method)) +
#   geom_vline(xintercept=100, linetype="dotted", show.legend=TRUE) +
#   geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=1) + geom_line() + geom_point() +
#   xlab("Number of selected variables") + ylab("Cohen's kappa") + ggtitle("b)") + theme_bw() +
#   theme(legend.position="none", plot.title=element_text(hjust=0.5))
# 
# mean.mse <- ggplot(plot.mse2, aes(x=psel, y=mean, colour=method)) +
#   geom_vline(xintercept=100, linetype="dotted", show.legend=TRUE) +
#   geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=1) + geom_line() + geom_point() +
#   xlab("Number of selected variables") + ylab("MSE") + ggtitle("c)") + theme_bw() +
#   theme(legend.position="none", plot.title=element_text(hjust=0.5))
# 
# mean.prec <- ggplot(plot.prec2, aes(x=psel, y=mean, colour=method)) +
#   geom_vline(xintercept=100, linetype="dotted", show.legend=TRUE) +
#   geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=1) + geom_line() + geom_point() +
#   xlab("Number of selected variables") + ylab("Precision") + ggtitle("d)") + theme_bw() +
#   theme(legend.position="none", plot.title=element_text(hjust=0.5))
# 
# mean.rec <- ggplot(plot.rec2, aes(x=psel, y=mean, colour=method)) +
#   geom_vline(xintercept=100, linetype="dotted", show.legend=TRUE) +
#   geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=1) + geom_line() + geom_point() +
#   xlab("Number of selected variables") + ylab("Recall") + ggtitle("e)") + theme_bw() +
#   theme(legend.justification=c(0, 1), legend.position=c(0, 1), plot.title=element_text(hjust=0.5))
# 
# mean.f1 <- ggplot(plot.f12, aes(x=psel, y=mean, colour=method)) +
#   geom_vline(xintercept=100, linetype="dotted", show.legend=TRUE) +
#   geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=1) + geom_line() + geom_point() +
#   xlab("Number of selected variables") + ylab("F1 score") + ggtitle("f)") + theme_bw() +
#   theme(legend.position="none", plot.title=element_text(hjust=0.5))
# 
# png(paste(path.graph, "grEBEN_sim_varsel1_set3_mean.png", sep=""), width=1200, height=720, res=90)
# grid.arrange(mean.auc, mean.kappa, mean.mse, mean.prec, mean.rec, mean.f1, ncol=3)
# dev.off()
# 
# ### setting4
# # loading the data and inspecting the variable selections
# load(paste(path.res, "grEBEN_sim_varsel1_setting4.Rdata", sep=""))
# unique(res.setting4$auc[, 1][!is.na(res.setting4$auc[, 1])])
# res.setting4$psel[!(res.setting4$psel[, 5] %in% pselseq) & !is.na(res.setting4$psel[, 5]), c(1, 2, 5)]
# pselseq <- sort(as.numeric(na.omit(unique(res.setting4$auc[, 2]))))
# methods <- colnames(res.setting4$auc)[-c(1:2)]
# res.setting4 <- lapply(res.setting4, function(met) {met[(res.setting4$psel[, 5] %in% pselseq) & !is.na(res.setting4$psel[, 5]), ]})
# 
# # creating the plotting tables
# plot.auc <- data.frame(psel=pselseq, t(sapply(pselseq, function(psel) {
#   apply(res.setting4$auc[, -c(1:2)][ifelse(is.na(res.setting4$auc[, 2]==psel), FALSE, res.setting4$auc[, 2]==psel), ], 2, function(vars) {
#     c(median(vars, na.rm=TRUE), mean(vars, na.rm=TRUE), sd(vars, na.rm=TRUE))})})))
# colnames(plot.auc) <- c("psel", as.vector(t(outer(methods, c("median", "mean", "sd"), paste, sep="."))))
# plot.kappa <- data.frame(psel=pselseq, t(sapply(pselseq, function(psel) {
#   apply(res.setting4$kappa[, -c(1:2)][ifelse(is.na(res.setting4$kappa[, 2]==psel), FALSE, res.setting4$kappa[, 2]==psel), ], 2, function(vars) {
#     c(median(vars, na.rm=TRUE), mean(vars, na.rm=TRUE), sd(vars, na.rm=TRUE))})})))
# colnames(plot.kappa) <- c("psel", as.vector(t(outer(methods, c("median", "mean", "sd"), paste, sep="."))))
# plot.mse <- data.frame(psel=pselseq, t(sapply(pselseq, function(psel) {
#   apply(res.setting4$mse[, -c(1:2)][ifelse(is.na(res.setting4$mse[, 2]==psel), FALSE, res.setting4$mse[, 2]==psel), ], 2, function(vars) {
#     c(median(vars, na.rm=TRUE), mean(vars, na.rm=TRUE), sd(vars, na.rm=TRUE))})})))
# colnames(plot.mse) <- c("psel", as.vector(t(outer(methods, c("median", "mean", "sd"), paste, sep="."))))
# plot.prec <- data.frame(psel=pselseq, t(sapply(pselseq, function(psel) {
#   apply(res.setting4$prec[, -c(1:2)][ifelse(is.na(res.setting4$prec[, 2]==psel), FALSE, res.setting4$prec[, 2]==psel), ], 2, function(vars) {
#     c(median(vars, na.rm=TRUE), mean(vars, na.rm=TRUE), sd(vars, na.rm=TRUE))})})))
# colnames(plot.prec) <- c("psel", as.vector(t(outer(methods, c("median", "mean", "sd"), paste, sep="."))))
# plot.rec <- data.frame(psel=pselseq, t(sapply(pselseq, function(psel) {
#   apply(res.setting4$rec[, -c(1:2)][ifelse(is.na(res.setting4$rec[, 2]==psel), FALSE, res.setting4$rec[, 2]==psel), ], 2, function(vars) {
#     c(median(vars, na.rm=TRUE), mean(vars, na.rm=TRUE), sd(vars, na.rm=TRUE))})})))
# colnames(plot.rec) <- c("psel", as.vector(t(outer(methods, c("median", "mean", "sd"), paste, sep="."))))
# plot.f1 <- data.frame(psel=pselseq, t(sapply(pselseq, function(psel) {
#   apply(res.setting4$f1[, -c(1:2)][ifelse(is.na(res.setting4$f1[, 2]==psel), FALSE, res.setting4$f1[, 2]==psel), ], 2, function(vars) {
#     c(median(vars, na.rm=TRUE), mean(vars, na.rm=TRUE), sd(vars, na.rm=TRUE))})})))
# colnames(plot.f1) <- c("psel", as.vector(t(outer(methods, c("median", "mean", "sd"), paste, sep="."))))
# 
# # plots of the medians over simulations
# png(paste(path.graph, "grEBEN_sim_varsel1_set4_median.png", sep=""), width=1200, height=720, res=90)
# par(mfrow=c(2, 3))
# plot(plot.auc[, 1], plot.auc[, 2], type="l", col=2, lty=1, ylim=range(plot.auc[, c(2, 5, 8, 11, 14, 17)]),
#      xlab="Number of selected variables", ylab="AUC", main="a)")
# lines(plot.auc[, 1], plot.auc[, 5], col=3, lty=1)
# lines(plot.auc[, 1], plot.auc[, 8], col=4, lty=1)
# lines(plot.auc[, 1], plot.auc[, 11], col=5, lty=1)
# lines(plot.auc[, 1], plot.auc[, 14], col=6, lty=1)
# lines(plot.auc[, 1], plot.auc[, 17], col=7, lty=1)
# abline(v=100, lty=3)
# legend("bottomright", legend=c(methods, "Active variables"), lty=c(rep(1, 6), 3), col=c(c(2:7), 1))
# 
# plot(plot.kappa[, 1], plot.kappa[, 2], type="l", col=2, lty=1, ylim=range(plot.kappa[, c(2, 5, 8, 11, 14, 17)]),
#      xlab="Number of selected variables", ylab="kappa", main="b)")
# lines(plot.kappa[, 1], plot.kappa[, 5], col=3, lty=1)
# lines(plot.kappa[, 1], plot.kappa[, 8], col=4, lty=1)
# lines(plot.kappa[, 1], plot.kappa[, 11], col=5, lty=1)
# lines(plot.kappa[, 1], plot.kappa[, 14], col=6, lty=1)
# lines(plot.kappa[, 1], plot.kappa[, 17], col=7, lty=1)
# abline(v=100, lty=3)
# 
# plot(plot.mse[, 1], plot.mse[, 2], type="l", col=2, lty=1, ylim=range(plot.mse[, c(2, 5, 8, 11, 14, 17)]),
#      xlab="Number of selected variables", ylab="mse", main="c)")
# lines(plot.mse[, 1], plot.mse[, 5], col=3, lty=1)
# lines(plot.mse[, 1], plot.mse[, 8], col=4, lty=1)
# lines(plot.mse[, 1], plot.mse[, 11], col=5, lty=1)
# lines(plot.mse[, 1], plot.mse[, 14], col=6, lty=1)
# lines(plot.mse[, 1], plot.mse[, 17], col=7, lty=1)
# abline(v=100, lty=3)
# 
# plot(plot.prec[, 1], plot.prec[, 2], type="l", col=2, lty=1, ylim=range(plot.prec[, c(2, 5, 8, 11, 14, 17)]),
#      xlab="Number of selected variables", ylab="precision", main="d)")
# lines(plot.prec[, 1], plot.prec[, 5], col=3, lty=1)
# lines(plot.prec[, 1], plot.prec[, 8], col=4, lty=1)
# lines(plot.prec[, 1], plot.prec[, 11], col=5, lty=1)
# lines(plot.prec[, 1], plot.prec[, 14], col=6, lty=1)
# lines(plot.prec[, 1], plot.prec[, 17], col=7, lty=1)
# abline(v=100, lty=3)
# 
# plot(plot.rec[, 1], plot.rec[, 2], type="l", col=2, lty=1, ylim=range(plot.rec[, c(2, 5, 8, 11, 14, 17)]),
#      xlab="Number of selected variables", ylab="recall", main="e)")
# lines(plot.rec[, 1], plot.rec[, 5], col=3, lty=1)
# lines(plot.rec[, 1], plot.rec[, 8], col=4, lty=1)
# lines(plot.rec[, 1], plot.rec[, 11], col=5, lty=1)
# lines(plot.rec[, 1], plot.rec[, 14], col=6, lty=1)
# lines(plot.rec[, 1], plot.rec[, 17], col=7, lty=1)
# abline(v=100, lty=3)
# 
# plot(plot.f1[, 1], plot.f1[, 2], type="l", col=2, lty=1, ylim=range(plot.f1[, c(2, 5, 8, 11, 14, 17)]),
#      xlab="Number of selected variables", ylab="f1", main="f)")
# lines(plot.f1[, 1], plot.f1[, 5], col=3, lty=1)
# lines(plot.f1[, 1], plot.f1[, 8], col=4, lty=1)
# lines(plot.f1[, 1], plot.f1[, 11], col=5, lty=1)
# lines(plot.f1[, 1], plot.f1[, 14], col=6, lty=1)
# lines(plot.f1[, 1], plot.f1[, 17], col=7, lty=1)
# abline(v=100, lty=3)
# dev.off()
# 
# # plots of the means plus error bars
# plot.auc2 <- reshape(plot.auc, varying=list(paste(methods, "median", sep="."), paste(methods, "mean", sep="."),
#                                             paste(methods, "sd", sep=".")), direction="long",
#                      times=methods, v.names=c("median", "mean", "sd"), idvar="psel")
# rownames(plot.auc2) <- NULL
# colnames(plot.auc2)[2] <- "method"
# plot.kappa2 <- reshape(plot.kappa, varying=list(paste(methods, "median", sep="."), paste(methods, "mean", sep="."),
#                                                 paste(methods, "sd", sep=".")), direction="long",
#                        times=methods, v.names=c("median", "mean", "sd"), idvar="psel")
# rownames(plot.kappa2) <- NULL
# colnames(plot.kappa2)[2] <- "method"
# plot.mse2 <- reshape(plot.mse, varying=list(paste(methods, "median", sep="."), paste(methods, "mean", sep="."),
#                                             paste(methods, "sd", sep=".")), direction="long",
#                      times=methods, v.names=c("median", "mean", "sd"), idvar="psel")
# rownames(plot.mse2) <- NULL
# colnames(plot.mse2)[2] <- "method"
# plot.prec2 <- reshape(plot.prec, varying=list(paste(methods, "median", sep="."), paste(methods, "mean", sep="."),
#                                               paste(methods, "sd", sep=".")), direction="long",
#                       times=methods, v.names=c("median", "mean", "sd"), idvar="psel")
# rownames(plot.prec2) <- NULL
# colnames(plot.prec2)[2] <- "method"
# plot.rec2 <- reshape(plot.rec, varying=list(paste(methods, "median", sep="."), paste(methods, "mean", sep="."),
#                                             paste(methods, "sd", sep=".")), direction="long",
#                      times=methods, v.names=c("median", "mean", "sd"), idvar="psel")
# rownames(plot.rec2) <- NULL
# colnames(plot.rec2)[2] <- "method"
# plot.f12 <- reshape(plot.f1, varying=list(paste(methods, "median", sep="."), paste(methods, "mean", sep="."),
#                                           paste(methods, "sd", sep=".")), direction="long",
#                     times=methods, v.names=c("median", "mean", "sd"), idvar="psel")
# rownames(plot.f12) <- NULL
# colnames(plot.f12)[2] <- "method"
# 
# mean.auc <- ggplot(plot.auc2, aes(x=psel, y=mean, colour=method)) +
#   geom_vline(xintercept=100, linetype="dotted", show.legend=TRUE) +
#   geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=1) + geom_line() + geom_point() +
#   xlab("Number of selected variables") + ylab("AUC") + ggtitle("a)") + theme_bw() +
#   theme(legend.position="none", plot.title=element_text(hjust=0.5))
# 
# mean.kappa <- ggplot(plot.kappa2, aes(x=psel, y=mean, colour=method)) +
#   geom_vline(xintercept=100, linetype="dotted", show.legend=TRUE) +
#   geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=1) + geom_line() + geom_point() +
#   xlab("Number of selected variables") + ylab("Cohen's kappa") + ggtitle("b)") + theme_bw() +
#   theme(legend.position="none", plot.title=element_text(hjust=0.5))
# 
# mean.mse <- ggplot(plot.mse2, aes(x=psel, y=mean, colour=method)) +
#   geom_vline(xintercept=100, linetype="dotted", show.legend=TRUE) +
#   geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=1) + geom_line() + geom_point() +
#   xlab("Number of selected variables") + ylab("MSE") + ggtitle("c)") + theme_bw() +
#   theme(legend.position="none", plot.title=element_text(hjust=0.5))
# 
# mean.prec <- ggplot(plot.prec2, aes(x=psel, y=mean, colour=method)) +
#   geom_vline(xintercept=100, linetype="dotted", show.legend=TRUE) +
#   geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=1) + geom_line() + geom_point() +
#   xlab("Number of selected variables") + ylab("Precision") + ggtitle("d)") + theme_bw() +
#   theme(legend.position="none", plot.title=element_text(hjust=0.5))
# 
# mean.rec <- ggplot(plot.rec2, aes(x=psel, y=mean, colour=method)) +
#   geom_vline(xintercept=100, linetype="dotted", show.legend=TRUE) +
#   geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=1) + geom_line() + geom_point() +
#   xlab("Number of selected variables") + ylab("Recall") + ggtitle("e)") + theme_bw() +
#   theme(legend.justification=c(0, 1), legend.position=c(0, 1), plot.title=element_text(hjust=0.5))
# 
# mean.f1 <- ggplot(plot.f12, aes(x=psel, y=mean, colour=method)) +
#   geom_vline(xintercept=100, linetype="dotted", show.legend=TRUE) +
#   geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=1) + geom_line() + geom_point() +
#   xlab("Number of selected variables") + ylab("F1 score") + ggtitle("f)") + theme_bw() +
#   theme(legend.position="none", plot.title=element_text(hjust=0.5))
# 
# png(paste(path.graph, "grEBEN_sim_varsel1_set4_mean.png", sep=""), width=1200, height=720, res=90)
# grid.arrange(mean.auc, mean.kappa, mean.mse, mean.prec, mean.rec, mean.f1, ncol=3)
# dev.off()


mygrridge <- function (highdimdata, response, partitions, unpenal = ~1, offset = NULL, 
                       method = "exactstable", niter = 10, monotone = NULL, optl = NULL, 
                       innfold = NULL, fixedfoldsinn = TRUE, selectionForward = FALSE, 
                       maxsel = 100, selectionEN = FALSE, stepsel = 1, cvlmarg = 1, 
                       savepredobj = "all", dataunpen = NULL, ord = 1:length(partitions), 
                       comparelasso = FALSE, optllasso = NULL, cvllasso = TRUE, 
                       compareEN = FALSE, compareunpenal = FALSE, trace = FALSE, 
                       modus = 1) 
{
  if (method == "adaptridge" | method == "exact") 
    niter <- 1
  if (class(partitions[[1]]) == "integer") {
    partitions = list(group = partitions)
  }
  nclass <- length(partitions)
  if (is.null(monotone)) 
    monotone <- rep(FALSE, nclass)
  if (length(monotone) != length(partitions)) {
    print(paste("ERROR: length 'monotone' unequal to length 'partitions' "))
    return(NULL)
  }
  partitions <- partitions[ord]
  monotone <- monotone[ord]
  nr <- nrow(highdimdata)
  for (ncl in 1:nclass) {
    indexset <- unlist(partitions[[ncl]])
    if (length(indexset) < nr) {
      print(paste("Warning: partition", ncl, "does not contain all row indices of the data"))
    }
    if (max(indexset) > nr | min(indexset) < 1) {
      print(paste("ERROR: partition", ncl, "contains an invalid index, e.g. larger than number of data rows"))
      return(NULL)
    }
  }
  overlap <- c()
  Wmat <- c()
  nfeattot <- c()
  for (ncl in 1:nclass) {
    indexset <- unlist(partitions[[ncl]])
    nfeatcl <- length(unique(indexset))
    nfeattot <- c(nfeattot, nfeatcl)
    if (length(indexset) > nfeatcl) {
      print(paste("Grouping", ncl, "contains overlapping groups"))
      overlap <- c(overlap, TRUE)
      whgroup <- partitions[[ncl]]
      nover <- rep(0, nr)
      for (k in 1:length(whgroup)) {
        wh <- whgroup[[k]]
        nover[wh] <- nover[wh] + 1
      }
      Wmat <- cbind(Wmat, sqrt(1/nover))
    }
    else {
      print(paste("Grouping", ncl, "contains mutually exclusive groups"))
      overlap <- c(overlap, FALSE)
      Wmat <- cbind(Wmat, rep(1, nr))
    }
  }
  arguments <- list(partitions = partitions, unpenal = unpenal, 
                    offset = offset, method = method, niter = niter, monotone = monotone, 
                    optl = optl, innfold = innfold, fixedfoldsinn = fixedfoldsinn, 
                    selectionForward = selectionForward, selectionEN = selectionEN, 
                    maxsel = maxsel, stepsel = stepsel, cvlmarg = cvlmarg, 
                    dataunpen = dataunpen, savepredobj = savepredobj, ord = ord, 
                    comparelasso = comparelasso, optllasso = optllasso, compareEN = compareEN, 
                    compareunpenal = compareunpenal, modus = modus)
  if (nr > 10000 & is.null(innfold)) 
    print("NOTE: consider setting innfold=10 to save computing time")
  nmp0 <- names(partitions)
  if (is.null(nmp0)) 
    nmp0 <- sapply(1:length(partitions), function(i) paste("Grouping", 
                                                           i))
  nmp0 <- sapply(1:length(partitions), function(i) {
    if (nmp0[i] == "") 
      return(paste("Grouping", i))
    else return(nmp0[i])
  })
  nmp <- c("NoGroups", "GroupRegul")
  nmpweight <- nmp
  if (selectionForward) 
    nmp <- c(nmp, "ForwSel")
  if (comparelasso) 
    nmp <- c(nmp, "lasso")
  if (selectionEN) 
    nmp <- c(nmp, "EN")
  if (compareunpenal) 
    nmp <- c(nmp, "modelunpen")
  if (class(response) == "factor") {
    nlevel <- length(levels(response))
    if (nlevel != 2) {
      print("Response is not binary, so not suitable for two-class classification.")
      return(NULL)
    }
    else {
      model = "logistic"
      print("Binary response, executing logistic ridge regression")
      lev <- levels(response)
      print(paste("Predicting probability on factor level", 
                  lev[2]))
    }
  }
  else {
    if (class(response) == "numeric" | class(response) == 
        "integer") {
      valresp <- sort(unique(response))
      if (length(valresp) == 2 & valresp[1] == 0 & valresp[2] == 
          1) {
        model = "logistic"
        print("Binary response, executing logistic ridge regression")
      }
      else {
        model = "linear"
        print("Numeric continuous response, executing linear ridge regression")
      }
    }
    else {
      if (class(response) == "Surv") {
        model = "survival"
        print("Survival response, executing cox ridge regression")
      }
      else {
        print("Non-valid response. Should be binary, numeric or survival.")
        return(NULL)
      }
    }
  }
  if ((unpenal != ~0) & (unpenal != ~1)) {
    if (is.null(dataunpen)) {
      print("If unpenal contains variables, data of \n                                  the unpenalized variables should be specified in the data slot!")
      return(NULL)
    }
  }
  nsam <- ncol(highdimdata)
  if (!is.null(offset)) {
    noffs <- length(offset)
    offsets <- "c("
    if (noffs == 1) {
      for (i in 1:(nsam - 1)) offsets <- paste(offsets, 
                                               offset, ",", sep = "")
    }
    else {
      for (i in 1:(nsam - 1)) offsets <- paste(offsets, 
                                               offset[i], ",", sep = "")
    }
    if (noffs == 1) 
      offsets <- paste(offsets, offset, ")", sep = "")
    else offsets <- paste(offsets, offset[nsam], ")", sep = "")
    if ((unpenal != ~0) & (unpenal != ~1)) {
      unpenal <- formula(paste(deparse(unpenal), "+ offset(", 
                               offsets, ")", sep = ""))
    }
    else {
      unpenal <- formula(paste("~", "offset(", offsets, 
                               ")", sep = ""))
    }
  }
  if (is.null(dataunpen)) 
    datapred <- data.frame(fake = rep(NA, ncol(highdimdata)))
  else datapred <- dataunpen
  nopen <- unpenal
  if (is.null(innfold)) 
    foldinit <- nsam
  else foldinit <- innfold
  pmt0 <- proc.time()
  optl0 <- optl
  if (is.null(optl)) {
    print("Finding lambda for initial ridge regression")
    if (fixedfoldsinn) 
      set.seed(346477)
    opt <- optL2(response, penalized = t(highdimdata), fold = foldinit, 
                 unpenalized = nopen, data = datapred, trace = trace)
    time1 <- proc.time() - pmt0
    print(opt$cv)
    print(paste("Computation time for cross-validating main penalty parameter:", 
                time1[3]))
    optl <- opt$lambda
    print(paste("lambda2", optl))
    arguments$optl <- optl
  }
  pmt <- proc.time()
  nsam <- ncol(highdimdata)
  nfeat <- nrow(highdimdata)
  XM0 <- t(highdimdata)
  response0 <- response
  if (!selectionForward) 
    cvlnstot <- rep(0, (nclass + 1))
  else cvlnstot <- rep(0, (nclass + 2))
  allpreds <- c()
  whsam <- 1:nsam
  responsemin <- response0
  pen0 <- penalized(responsemin, penalized = XM0, lambda2 = optl, 
                    unpenalized = nopen, data = cbind(XM0, datapred))
  nmunpen <- names(pen0@unpenalized)
  if (is.element("(Intercept)", nmunpen)) 
    addintercept <- TRUE
  else addintercept <- FALSE
  if (is.null(innfold)) {
    nf <- nrow(XM0)
  }
  else {
    if (!is.null(optl0)) {
      nf <- innfold
      if (fixedfoldsinn) 
        set.seed(346477)
    }
    else {
      nf <- opt$fold
    }
  }
  opt2 <- cvl(responsemin, penalized = XM0, fold = nf, lambda2 = optl, 
              unpenalized = nopen, data = datapred, trace = trace)
  nf <- opt2$fold
  cvln0 <- opt2$cvl
  cvlnprev <- cvln0
  penprev <- pen0
  pen <- pen0
  print(cvln0)
  XMw0 <- XM0
  XMw0prev <- XM0
  converged <- FALSE
  conv <- rep(FALSE, nclass)
  controlbound1 <- 1000
  controlbound2 <- controlbound3 <- 10
  almvecall <- rep(1, nfeat)
  lambdas <- lapply(partitions, function(cla) {
    ngroup <- length(cla)
    return(rep(1, ngroup))
  })
  lmvec <- lmvecprev <- array(1, nfeat)
  i <- 1
  while (!converged & i <= niter) {
    cl <- 1
    if (method == "adaptridge") 
      cl <- nclass
    while (cl <= nclass) {
      convcl <- conv[cl]
      if (!convcl) {
        whgr <- partitions[[cl]]
        lenggr <- unlist(lapply(whgr, length))
        ngroup1 <- length(whgr)
        names(lambdas[[cl]]) <- names(whgr)
        coeff <- penprev@penalized
        if (model == "survival") {
          preds <- predict(penprev, XMw0, data = datapred)
        }
        else {
          preds <- predict(penprev, XMw0, data = datapred)[1:nsam]
        }
        coeffsq <- coeff^2
        if (model == "logistic") {
          Wi <- sqrt(preds * (1 - preds))
          constlam <- 2
        }
        if (model == "linear") {
          Wi <- rep(1, length(preds))
          constlam <- 1
        }
        if (model == "survival") {
          resptime <- response[, 1]
          predsnew <- -log(sapply(1:nsam, function(k) survival(preds, 
                                                               time = resptime[k])[k]))
          Wi <- sqrt(predsnew)
          constlam <- 2
        }
        if (!is.null(dataunpen)) {
          mm <- model.matrix(nopen, dataunpen)
          XMW <- t(t(cbind(XMw0, 10^5 * mm)) %*% diag(Wi))
        }
        else {
          if (addintercept) 
            XMW <- t(t(cbind(XMw0, rep(10^5, nsam))) %*% 
                       diag(Wi))
          else XMW <- t(t(XMw0) %*% diag(Wi))
        }
        SVD <- svd(XMW)
        leftmat <- SVD$v %*% diag(1/((SVD$d)^2 + constlam * 
                                       optl)) %*% diag(SVD$d) %*% t(SVD$u)
        if (model == "linear") {
          Hatm <- XMW %*% leftmat
          df <- nsam - sum(diag(2 * Hatm - Hatm %*% t(Hatm)))
          VarRes <- sum((response - preds)^2)/df
          print(paste("Sigma^2 estimate:", VarRes))
          vars3 <- VarRes * rowSums(leftmat^2)
        }
        else {
          vars3 <- rowSums(leftmat^2)
        }
        which0 <- which(vars3 == 0)
        vars3[which0] <- 10^{
          -30
        }
        if (model == "linear") {
          mycoeff2svd <- (leftmat %*% response)^2
        }
        if (model == "logistic") {
          if (is.factor(response)) 
            respnum <- as.numeric(response) - 1
          else respnum <- response
          z <- matrix(log(preds/(1 - preds)) + (respnum - 
                                                  preds)/(preds * (1 - preds)), ncol = 1)
          if (modus == 1) 
            mycoeff2svd <- coeffsq
          if (modus == 2) 
            mycoeff2svd <- (leftmat %*% z)^2
        }
        if (model == "survival") {
          mycoeff2svd <- coeffsq
        }
        cii2 <- (rowSums(leftmat * t(XMW)))^2
        leftmat <- leftmat/sqrt(vars3)
        lowerleft <- 10^(-30)
        lefts2 <- function(group) {
          ind <- whgr[[group]]
          ngr <- lenggr[group]
          coefftau2 <- sum(sapply(mycoeff2svd[ind]/vars3[ind], 
                                  function(x) max(x, 1))) - length(ind)
          return(max(lowerleft, coefftau2/ngr))
        }
        leftside <- sapply(1:length(whgr), lefts2)
        ellarg0 <- length(leftside[leftside > lowerleft])/length(leftside)
        if (ellarg0 <= 0.5) {
          print(paste("Partition", nmp0[cl], "NOT ITERATED"))
          conv[cl] <- TRUE
          cvln1 <- cvlnprev
          XMw0 <- XMw0prev
          pen <- penprev
        }
        else {
          lefts2ran <- function(group, randomind) {
            ind <- whgr[[group]]
            ngr <- lenggr[group]
            coefftau2 <- sum(sapply(mycoeff2svd[1:nfeat][randomind[ind]]/vars3[1:nfeat][randomind[ind]], 
                                    function(x) max(x, 1))) - length(randomind[ind])
            return(max(lowerleft, coefftau2/ngr))
          }
          randomiz <- function(fakex) {
            randomind2 <- sample(1:nfeat)
            leftsideran <- sapply(1:length(whgr), lefts2ran, 
                                  randomind = randomind2)
            return(leftsideran)
          }
          nlefts <- 100
          leftsran <- sapply(1:nlefts, randomiz)
          means <- apply(leftsran, 1, mean)
          leftsrancen <- t(t(leftsran) - means)
          relerror <- sum(abs(leftsrancen))/(nlefts * 
                                               sum(abs(means)))
          if (cl == 1 & i == 1) 
            print(cvln0)
          print(paste("Relative error:", relerror))
          if (relerror >= 0.1) 
            print("WARNING: large relative error (>=0.1). Consider using larger groups of variable.")
          nadd <- ncol(XMW) - ncol(XMw0)
          rightmat = t(t(XMW) * c(Wmat[, cl], rep(1, 
                                                  nadd)))
          rightmats <- lapply(1:length(whgr), function(j) {
            rightj2 <- rightmat[, whgr[[j]]]
            rcp <- rightj2 %*% t(rightj2)
            return(rcp)
          })
          coefmatfast <- t(apply(matrix(1:length(whgr), 
                                        nrow = length(whgr)), 1, function(i) {
                                          lefti2 <- leftmat[whgr[[i]], ]
                                          lcp <- t(lefti2) %*% lefti2
                                          ckls <- sapply(1:length(whgr), function(j) {
                                            rcp <- rightmats[[j]]
                                            return(sum(lcp * rcp))
                                          })
                                          return(ckls)
                                        }))
          coefmatfast <- coefmatfast/lenggr
          CNfun <- function(lam, cfmmat = coefmatfast) {
            ng <- nrow(cfmmat)
            dmax <- max(diag(cfmmat))
            cfmlam <- (1 - lam) * cfmmat + lam * diag(dmax, 
                                                      nrow = ng)
            eigenvals <- eigen(cfmlam, only.values = TRUE)$values
            CN <- eigenvals[1]/eigenvals[ng]
            return(Re(CN))
          }
          lams <- seq(0, 1, by = 0.005)
          CNsRan <- sapply(lams, CNfun, cfmmat = coefmatfast)
          CNsRanre <- CNsRan * relerror
          if (relerror <= 0.1) {
            lam <- lams[which(CNsRanre <= 0.1)[1]]
          }
          else lam <- 1
          print(paste("Shrink Factor coefficient matrix", 
                      lam))
          cfmmat <- coefmatfast
          ng <- nrow(cfmmat)
          dmax <- max(diag(cfmmat))
          cfmlam <- (1 - lam) * cfmmat + lam * diag(dmax, 
                                                    nrow = ng)
          if (method == "exactstable") {
            soltau = solve(sum(cfmlam), sum(leftside))
            sol = solve(cfmlam, leftside)
            low <- soltau/controlbound1
            up = soltau * controlbound1
            parinint <- sapply(sol, function(x) min(max(low, 
                                                        x), up))
            minopt <- optim(par = parinint, fn = function(pars = c(parinint)) sum(leftside - 
                                                                                    cfmlam %*% pars)^2, method = "L-BFGS-B", 
                            lower = rep(low, ngroup1), upper = rep(up, 
                                                                   ngroup1))
            tausqest0 <- minopt$par
          }
          if (method == "exact") {
            soltau = solve(sum(coefmatfast), sum(leftside))
            sol = solve(coefmatfast, leftside)
            low <- soltau/controlbound2
            up = soltau * controlbound2
            parinint <- sapply(sol, function(x) min(max(low, 
                                                        x), up))
            minopt <- optim(par = parinint, fn = function(pars = c(parinint)) sum(leftside - 
                                                                                    cfmlam %*% pars)^2, method = "L-BFGS-B", 
                            lower = rep(low, ngroup1), upper = rep(up, 
                                                                   ngroup1))
            tausqest0 <- minopt$par
          }
          if (method == "stable") {
            soltau = solve(sum(coefmatfast), sum(leftside))
            solhyb <- sapply(1:ngroup1, function(i) {
              leftsidei <- leftside[i]
              rightsidei <- c(coefmatfast[i, i], sum(coefmatfast[i, 
                                                                 -i] * soltau))
              soli <- (leftsidei - rightsidei[2])/rightsidei[1]
              return(max(min(soltau * controlbound3, 
                             soli), soltau/controlbound3))
            })
          }
          if (method == "simple") {
            solsim = leftside
            tausqest <- solsim
            print("simple")
          }
          if (method == "exact") {
            tausqest <- tausqest0
            print("exact")
          }
          if (method == "exactstable") {
            tausqest <- tausqest0
            print("exactstable")
          }
          if (method == "stable") {
            print("stable")
            tausqest <- solhyb
          }
          if (method == "adaptridge") {
            print("adaptive ridge")
          }
          if (method == "stable" | method == "exact" | 
              method == "exactstable" | method == "simple") {
            lambdanoncal <- 1/tausqest
            if (monotone[cl]) {
              weigh = unlist(lapply(whgr, length))
              lambdamultnoncal <- pava(lambdanoncal, 
                                       w = weigh)
            }
            else lambdamultnoncal <- lambdanoncal
            tausqest <- 1/lambdamultnoncal
            nfeatcl <- nfeattot[cl]
            overl <- overlap[cl]
            if (!overl) {
              con3 <- sum(sapply(1:length(whgr), function(gr) {
                return(length(whgr[[gr]]) * tausqest[gr])
              }))
              tausqestcal <- nfeatcl/con3 * tausqest
              lambdamult <- 1/tausqestcal
              print(lambdamult)
              for (k in 1:length(whgr)) {
                wh <- whgr[[k]]
                XMw0[, wh] <- XMw0[, wh]/sqrt(lambdamult[k])
              }
            }
            else {
              tauk <- rep(0, nfeat)
              Wsq <- (Wmat[, cl])^2
              for (k in 1:length(whgr)) {
                wh <- whgr[[k]]
                tauk[wh] <- tauk[wh] + tausqest[k]
              }
              tauk <- tauk * Wsq
              whna <- which(is.na(tauk))
              if (length(whna) > 0) 
                con3 <- sum(tauk[-whna])
              else con3 <- sum(tauk)
              tausqestcal0 <- nfeatcl/con3 * tausqest
              lambdamult <- 1/tausqestcal0
              print(lambdamult)
              tausqestcal <- (nfeatcl/con3) * tauk
              lambdamultperk <- 1/tausqestcal
              lambdamultperk[whna] <- 1
              XMw0 <- t(t(XMw0)/sqrt(lambdamultperk))
            }
          }
          else {
            tausqest <- coeffsq
            con3 <- sum(tausqest)
            tausqestcal <- nfeat/con3 * tausqest
            lambdamult <- 1/tausqestcal
            XMw0 <- t(t(XMw0)/sqrt(lambdamult))
          }
          opt2w <- cvl(responsemin, penalized = XMw0, 
                       fold = nf, lambda2 = optl, unpenalized = nopen, 
                       data = datapred, trace = trace)
          cvln1 <- opt2w$cvl
          print(cvln1)
          if ((cvln1 - cvlnprev)/abs(cvlnprev) > 1/100 | 
              ((cvln1 - cvlnprev)/abs(cvlnprev) >= 0 & 
               i == 1)) {
            pen <- penalized(responsemin, penalized = XMw0, 
                             lambda2 = optl, unpenalized = nopen, data = datapred)
            if (niter > 1) {
              if (!overl) {
                for (group in 1:ngroup1) {
                  ind <- whgr[[group]]
                  lmvec[ind] <- lmvec[ind] * lambdamult[group]
                }
              }
              else {
                lmvec <- lmvec * lambdamultperk
              }
            }
            lambdas[[cl]] <- lambdas[[cl]] * lambdamult
            cvlnprev <- cvln1
            penprev <- pen
            XMw0prev <- XMw0
            print(paste("Partition", nmp0[cl], "improved results"))
          }
          else {
            if (niter > 1) 
              print(paste("Partition", nmp0[cl], "CONVERGED after", 
                          i, "iterations"))
            else print(paste("Partition", nmp0[cl], "did not improve results"))
            conv[cl] <- TRUE
            cvln1 <- cvlnprev
            XMw0 <- XMw0prev
            pen <- penprev
          }
        }
      }
      cl <- cl + 1
    }
    if (sum(conv) == nclass) {
      converged <- TRUE
      if (niter > 1) 
        print(paste("All partitions CONVERGED after", 
                    i, "iterations"))
    }
    i <- i + 1
  }
  if (niter == 0) {
    pen <- pen0
    XMw0 <- XM0
    cvln1 <- cvln0
    soltau <- NULL
  }
  if (model == "survival") {
    pred0 <- predict(pen0, XM0, unpenalized = nopen, data = datapred)
    predw <- predict(pen, XMw0, unpenalized = nopen, data = datapred)
  }
  else {
    pred0 <- predict(pen0, XM0, unpenalized = nopen, data = datapred)[1:nsam]
    predw <- predict(pen, XMw0, unpenalized = nopen, data = datapred)[1:nsam]
  }
  predshere <- cbind(pred0, predw)
  cvlnssam <- c(cvln0, cvln1)
  lmvecall <- lmvec
  almvecall <- cbind(almvecall, lmvecall)
  predobj <- c(pen0, pen)
  allpreds <- predshere
  whichsel <- NULL
  betassel <- NULL
  npr <- length(predobj)
  pred2 <- predobj[[npr]]
  lambs <- lmvecall
  oldbeta <- pred2@penalized
  newbeta <- oldbeta/sqrt(lambs)
  time2 <- proc.time() - pmt
  print(paste("Computation time for adaptive weigthing:", time2[3]))
  if (selectionForward) {
    print("Start posthoc variable selection")
    pmt <- proc.time()
    ord <- order(abs(newbeta), decreasing = TRUE)
    sequ <- seq(0, maxsel, by = stepsel)
    cvlsels <- c()
    for (nsel in sequ) {
      whsel <- ord[1:nsel]
      datwsel <- XMw0[, whsel]
      pensel <- penalized(responsemin, datwsel, lambda2 = optl, 
                          unpenalized = nopen, data = datapred, trace = FALSE)
      optsel <- cvl(responsemin, datwsel, fold = nf, lambda2 = optl, 
                    unpenalized = nopen, data = datapred, trace = trace)
      cvlsel <- optsel$cvl
      cvlsels <- c(cvlsels, cvlsel)
    }
    whbest <- which.max(cvlsels)
    cvmax <- cvlsels[whbest]
    nsel2 <- sequ[(which((cvlsels - cvmax) >= cvlmarg/100 * 
                           (cvmax)))[1]]
    whsel2 <- ord[1:nsel2]
    whichsel <- whsel2
    betassel <- newbeta[whsel2]
    datwsel2 <- XMw0[, whsel2]
    pensel2 <- penalized(responsemin, datwsel2, lambda2 = optl, 
                         unpenalized = nopen, data = datapred)
    optsel2 <- cvl(responsemin, datwsel2, fold = nf, lambda2 = optl, 
                   unpenalized = nopen, data = datapred, trace = trace)
    cvlsel2 <- optsel2$cvl
    if (model == "survival") {
      predsel <- predict(pensel2, penalized = XMw0[, whsel2, 
                                                   drop = FALSE], unpenalized = nopen, data = datapred)
    }
    else {
      predsel <- predict(pensel2, penalized = XMw0[, whsel2, 
                                                   drop = FALSE], unpenalized = nopen, data = datapred)[1:nsam]
    }
    allpreds <- cbind(allpreds, predsel)
    predobj <- c(predobj, pensel2)
    cvlnssam <- c(cvlnssam, cvlsel2)
    print(paste("Number of selected markers by forward selection:", 
                nsel2))
    time3 <- proc.time() - pmt
    print(paste("Computation time for feature selection by forward selection:", 
                time3[3]))
  }
  if (selectionForward) 
    cvlssel <- data.frame(nsel = sequ, cvl = cvlsels)
  else cvlssel <- NULL
  cvlnstot <- cvlnssam
  reslasso <- NULL
  if (comparelasso) {
    print("Starting lasso")
    if (is.null(optllasso)) {
      print("Finding lambda for lasso regression")
      opt <- optL1(response, penalized = t(highdimdata), 
                   fold = nf, unpenalized = nopen, data = datapred, 
                   trace = trace)
      print(opt$cv)
      optllasso <- opt$lambda
      print(paste("lambda1", optllasso))
      arguments$optllasso <- optllasso
      cvllasso <- opt$cv
    }
    else {
      cvliklasso <- if (cvllasso) 
        try(cvl(response, penalized = t(highdimdata), 
                lambda1 = optllasso, fold = nf, unpenalized = nopen, 
                data = datapred, trace = trace))
      if (class(cvliklasso) == "try-error" | !cvllasso) 
        cvllasso <- NA
      else cvllasso <- cvliklasso$cvl
    }
    cvlnstot <- c(cvlnstot, cvllasso)
    penlasso <- penalized(response, penalized = t(highdimdata), 
                          lambda1 = optllasso, unpenalized = nopen, data = cbind(XM0, 
                                                                                 datapred))
    whichlasso <- which(penlasso@penalized != 0)
    betaslasso <- penlasso@penalized[whichlasso]
    predobj <- c(predobj, penlasso)
    reslasso <- list(cvllasso = cvllasso, whichlasso = whichlasso, 
                     betaslasso = betaslasso)
  }
  resEN <- NULL
  if (selectionEN) {
    print("Variable selection by elastic net started...")
    if (selectionForward) 
      maxsel <- length(whsel2)
    fsel <- function(lam1, maxselec = maxsel, lam2) {
      if (lam1 == 0) 
        return(nfeat - maxselec)
      else {
        penselEN <- penalized(responsemin, XMw0, lambda1 = lam1, 
                              lambda2 = lam2, unpenalized = nopen, data = datapred, 
                              trace = FALSE, maxiter = 100)
        coef <- penselEN@penalized
        return(length(coef[coef != 0]) - maxselec)
      }
    }
    lam1 <- uniroot(fsel, interval = c(0, optl * 100), maxiter = 50, 
                    lam2 = optl)$root
    penselEN0 <- penalized(responsemin, XMw0, lambda1 = lam1, 
                           lambda2 = optl, unpenalized = nopen, data = datapred, 
                           trace = FALSE, maxiter = 100)
    coefEN0 <- penselEN0@penalized
    whichEN <- which(coefEN0 != 0)
    penselEN <- penalized(responsemin, XMw0[, whichEN, drop = FALSE], 
                          lambda2 = optl, unpenalized = nopen, data = datapred, 
                          trace = FALSE, maxiter = 100)
    coefEN <- penselEN@penalized
    predobj <- c(predobj, penselEN)
    resEN <- list(whichEN = whichEN, betasEN = coefEN)
  }
  if (compareunpenal) {
    if (model == "survival") {
      print("Starting unpenalized Cox-model")
      bogus <- matrix(rnorm(nsam), ncol = 1)
      print(dim(datapred))
      penlambdas0 <- penalized(response, penalized = bogus, 
                               unpenalized = nopen, lambda1 = 0, lambda2 = 10^8, 
                               data = datapred)
      predobj <- c(predobj, penlambdas0)
    }
    else {
      if (model == "logistic") 
        famglm <- "binomial"
      else famglm <- "gaussian"
      print("Starting unpenalized glm")
      form <- formula(paste("response", "~", as.character(unpenal)[2]))
      modelglm <- glm(form, family = famglm, data = dataunpen)
      predobj <- c(predobj, list(modelglm))
    }
  }
  printlam <- function(lambs) {
    if (length(lambs) <= 10) 
      return(lambs)
    else return(summary(lambs))
  }
  suml <- lapply(lambdas, printlam)
  print("Final lambda multipliers (summary):")
  print(suml)
  print(paste("CVLs", cvlnstot))
  timetot <- proc.time() - pmt0
  print(paste("Total computation time:", timetot[3]))
  names(predobj) <- nmp
  colnames(almvecall) <- nmpweight
  if (savepredobj == "last") {
    predobj <- predobj[length(predobj)]
    almvecall <- matrix(lmvecall, ncol = 1)
  }
  if (savepredobj == "none") 
    predobj <- NULL
  cat("\n")
  if (model == "linear") 
    print("PLEASE NOTE THAT WE HIGHLY RECOMMEND TO USE grridgelin INSTEAD OF THIS FUNCTION FOR LINEAR RESPONSE")
  return(list(true = response, cvls = cvlnstot, lambdamults = lambdas, 
              optl = optl, lambdamultvec = almvecall, predobj = predobj, 
              betas = newbeta, whichsel = whichsel, cvlssel = cvlssel, 
              reslasso = reslasso, resEN = resEN, model = model, arguments = arguments, 
              allpreds = allpreds))
}


