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
# SETTING 1 REP 15
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

# source version of grridge that works
source(paste(path.code, "mygrridge.R", sep=""))

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
q <- 0.9
beta.mean <- 0.1

# simulation and estimation characteristics
nreps <- 50
pselmin <- 10
pselmax <- 150
pselstep <- 10

# combining the simulation settings
set <- expand.grid(n=n, m=m, ntest=ntest, p=p, G=G, f=f, q=q, beta.mean=beta.mean, rhow=rhow, 
                   rhob=rhob, sigma2=sigma2, pselmin=pselmin, pselmax=pselmax, pselstep=pselstep)

# loop over the settings
for(s in 1:nrow(set)) {
  
  betag <- 0.5*set$beta.mean[s]*(set$f[s] - 1)^2*(set$G[s]^2 - set$G[s])*set$f[s]^c(0:(set$G[s] - 1))/
    ((1 - set$q[s])*(set$G[s]*set$f[s]^(set$G[s] + 1) - set$f[s]^(set$G[s] + 1) - 
                       set$G[s]*set$f[s]^set$G[s] + set$f[s]))
  pactiveg <- 2*(1 - set$q[s])*c(0:(set$G[s] - 1))*set$p[s]/(set$G[s]^2 - set$G[s])
  beta <- as.numeric(sapply(1:set$G[s], function(g) {
    rep(c(0, betag[g]), times=round(c(set$p[s]/set$G[s] - pactiveg[g], pactiveg[g])))}))

  Sigma.group <- matrix(c(rep(rep(c(set$sigma2[s], set$rhow[s]), times=c(1, set$p[s]/set$G[s])),
                              times=set$p[s]/set$G[s] - 1), set$sigma2[s]),
                        ncol=set$p[s]/set$G[s], nrow=set$p[s]/set$G[s])
  Sigma.block <- as.matrix(bdiag(Sigma.group, Sigma.group))
  diag(Sigma.block[(set$p[s]/set$G[s] + 1):(2*set$p[s]/set$G[s]), 1:(set$p[s]/set$G[s])]) <-
    diag(Sigma.block[1:(set$p[s]/set$G[s]), (set$p[s]/set$G[s] + 1):(2*set$p[s]/set$G[s])]) <-
    rep(c(set$rhob[s], 0), times=0.5*set$p[s]/set$G[s])
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
    # if(Sys.info()[1]=="Darwin") {
    #   registerDoMC(2)
    # } else {
    #   registerDoMC(length(pselseq))
    # }
    # # loop over the different estimations of the models
    # result <- foreach(csel=c(1:length(pselseq)), .combine=rbind) %dopar% {
    for(csel in 1:length(pselseq)) {
      psel <- pselseq[csel]

      # variable selection methods
      fit.grEBEN <- grEBEN(x, y, rep(set$m[s], set$n[s]), unpenalized=NULL, intercept=TRUE, partitions=partitions1, lambda1=NULL, lambda2=NULL,
                           monotone=list(FALSE, FALSE), psel=psel, posterior=FALSE, ELBO=FALSE, eps=0.001, maxiter=500, trace=FALSE)
      fit.GRridge2 <- mygrridge(t(x), y, partitions=partitions2, monotone=FALSE, optl=fit.GRridge1$optl,
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

      # # return the vector of calculated metrics
      # return(c(auc.enet, auc.enetsel, auc.grEBENsel,
      #          auc.ridge, auc.GRridge, auc.GRridgesel, psel.enet, psel.enetsel, psel.grEBENsel,
      #          psel.ridge, psel.GRridge, psel.GRridgesel, kappa.enet, kappa.enetsel, kappa.grEBENsel,
      #          kappa.ridge, kappa.GRridge, kappa.GRridgesel, mse.enet, mse.enetsel, mse.grEBENsel,
      #          mse.ridge, mse.GRridge, mse.GRridgesel, prec.enet, prec.enetsel, prec.grEBENsel,
      #          prec.ridge, prec.GRridge, prec.GRridgesel, rec.enet, rec.enetsel, rec.grEBENsel,
      #          rec.ridge, rec.GRridge, rec.GRridgesel, f1.enet, f1.enetsel, f1.grEBENsel,
      #          f1.ridge, f1.GRridge, f1.GRridgesel))
      
      # assigning to matrices in parallel case
      aucmat[(crep - 1)*nreps + csel, ] <- c(crep, psel, auc.enet, auc.enetsel, auc.grEBENsel, auc.ridge, auc.GRridge, auc.GRridgesel)
      pselmat[(crep - 1)*nreps + csel, ] <- c(crep, psel, psel.enet, psel.enetsel, psel.grEBENsel, psel.ridge, psel.GRridge, psel.GRridgesel)
      kappamat[(crep - 1)*nreps + csel, ] <- c(crep, psel, kappa.enet, kappa.enetsel, kappa.grEBENsel, kappa.ridge, kappa.GRridge, kappa.GRridgesel)
      msemat[(crep - 1)*nreps + csel, ] <- c(crep, psel, mse.enet, mse.enetsel, mse.grEBENsel, mse.ridge, mse.GRridge, mse.GRridgesel)
      precmat[(crep - 1)*nreps + csel, ] <- c(crep, psel, prec.enet, prec.enetsel, prec.grEBENsel, prec.ridge, prec.GRridge, prec.GRridgesel)
      recmat[(crep - 1)*nreps + csel, ] <- c(crep, psel, rec.enet, rec.enetsel, rec.grEBENsel, rec.ridge, rec.GRridge, rec.GRridgesel)
      f1mat[(crep - 1)*nreps + csel, ] <- c(crep, psel, f1.enet, f1.enetsel, f1.grEBENsel, f1.ridge, f1.GRridge, f1.GRridgesel)
      
    }

    # # assigning to matrices in parallel case
    # aucmat[((crep - 1)*length(pselseq) + 1):(crep*length(pselseq)), ] <- cbind(crep, pselseq, result[, c(1:length(methods))])
    # pselmat[((crep - 1)*length(pselseq) + 1):(crep*length(pselseq)), ] <- cbind(crep, pselseq, result[, c((length(methods) + 1):(2*length(methods)))])
    # kappamat[((crep - 1)*length(pselseq) + 1):(crep*length(pselseq)), ] <- cbind(crep, pselseq, result[, c((2*length(methods) + 1):(3*length(methods)))])
    # msemat[((crep - 1)*length(pselseq) + 1):(crep*length(pselseq)), ] <- cbind(crep, pselseq, result[, c((3*length(methods) + 1):(4*length(methods)))])
    # precmat[((crep - 1)*length(pselseq) + 1):(crep*length(pselseq)), ] <- cbind(crep, pselseq, result[, c((4*length(methods) + 1):(5*length(methods)))])
    # recmat[((crep - 1)*length(pselseq) + 1):(crep*length(pselseq)), ] <- cbind(crep, pselseq, result[, c((5*length(methods) + 1):(6*length(methods)))])
    # f1mat[((crep - 1)*length(pselseq) + 1):(crep*length(pselseq)), ] <- cbind(crep, pselseq, result[, c((6*length(methods) + 1):(7*length(methods)))])

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
# load(paste(path.res, "grEBEN_sim_varsel2_setting1.Rdata", sep=""))
# pselseq <- sort(as.numeric(na.omit(unique(res.setting1$auc[, 2]))))
# methods <- colnames(res.setting1$auc)[-c(1:2)]
# unique(res.setting1$auc[, 1][!is.na(res.setting1$auc[, 1])])
# res.setting1$psel[!(res.setting1$psel[, 5] %in% pselseq) & !is.na(res.setting1$psel[, 5]), c(1, 2, 5)]
# rle(sort(res.setting1$psel[!(res.setting1$psel[, 5] %in% pselseq) & !is.na(res.setting1$psel[, 5]), 2]))
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
# png(paste(path.graph, "grEBEN_sim_varsel2_set1_median.png", sep=""), width=1200, height=720, res=90)
# par(mfrow=c(2, 3))
# plot(plot.auc[, 1], plot.auc[, 2], type="l", col=2, lty=1, ylim=range(plot.auc[, c(2, 5, 8, 11, 14, 17)]),
#      xlab="Number of selected variables", ylab="AUC", main="a)")
# lines(plot.auc[, 1], plot.auc[, 5], col=3, lty=1)
# lines(plot.auc[, 1], plot.auc[, 8], col=4, lty=1)
# lines(plot.auc[, 1], plot.auc[, 11], col=5, lty=1)
# lines(plot.auc[, 1], plot.auc[, 14], col=6, lty=1)
# lines(plot.auc[, 1], plot.auc[, 17], col=7, lty=1)
# abline(v=sum(set$p[1]/set$G[1] - set$p[1]*c(0:(set$G[1] - 1))/set$G[1]^2), lty=3)
# legend("bottomright", legend=c(methods, "Active variables"), lty=c(rep(1, 6), 3), col=c(c(2:7), 1))
# 
# plot(plot.kappa[, 1], plot.kappa[, 2], type="l", col=2, lty=1, ylim=range(plot.kappa[, c(2, 5, 8, 11, 14, 17)]),
#      xlab="Number of selected variables", ylab="kappa", main="b)")
# lines(plot.kappa[, 1], plot.kappa[, 5], col=3, lty=1)
# lines(plot.kappa[, 1], plot.kappa[, 8], col=4, lty=1)
# lines(plot.kappa[, 1], plot.kappa[, 11], col=5, lty=1)
# lines(plot.kappa[, 1], plot.kappa[, 14], col=6, lty=1)
# lines(plot.kappa[, 1], plot.kappa[, 17], col=7, lty=1)
# abline(v=sum(set$p[1]/set$G[1] - set$p[1]*c(0:(set$G[1] - 1))/set$G[1]^2), lty=3)
# 
# plot(plot.mse[, 1], plot.mse[, 2], type="l", col=2, lty=1, ylim=range(plot.mse[, c(2, 5, 8, 11, 14, 17)]),
#      xlab="Number of selected variables", ylab="mse", main="c)")
# lines(plot.mse[, 1], plot.mse[, 5], col=3, lty=1)
# lines(plot.mse[, 1], plot.mse[, 8], col=4, lty=1)
# lines(plot.mse[, 1], plot.mse[, 11], col=5, lty=1)
# lines(plot.mse[, 1], plot.mse[, 14], col=6, lty=1)
# lines(plot.mse[, 1], plot.mse[, 17], col=7, lty=1)
# abline(v=sum(set$p[1]/set$G[1] - set$p[1]*c(0:(set$G[1] - 1))/set$G[1]^2), lty=3)
# 
# plot(plot.prec[, 1], plot.prec[, 2], type="l", col=2, lty=1, ylim=range(plot.prec[, c(2, 5, 8, 11, 14, 17)]),
#      xlab="Number of selected variables", ylab="precision", main="d)")
# lines(plot.prec[, 1], plot.prec[, 5], col=3, lty=1)
# lines(plot.prec[, 1], plot.prec[, 8], col=4, lty=1)
# lines(plot.prec[, 1], plot.prec[, 11], col=5, lty=1)
# lines(plot.prec[, 1], plot.prec[, 14], col=6, lty=1)
# lines(plot.prec[, 1], plot.prec[, 17], col=7, lty=1)
# abline(v=sum(set$p[1]/set$G[1] - set$p[1]*c(0:(set$G[1] - 1))/set$G[1]^2), lty=3)
# 
# plot(plot.rec[, 1], plot.rec[, 2], type="l", col=2, lty=1, ylim=range(plot.rec[, c(2, 5, 8, 11, 14, 17)]),
#      xlab="Number of selected variables", ylab="recall", main="e)")
# lines(plot.rec[, 1], plot.rec[, 5], col=3, lty=1)
# lines(plot.rec[, 1], plot.rec[, 8], col=4, lty=1)
# lines(plot.rec[, 1], plot.rec[, 11], col=5, lty=1)
# lines(plot.rec[, 1], plot.rec[, 14], col=6, lty=1)
# lines(plot.rec[, 1], plot.rec[, 17], col=7, lty=1)
# abline(v=sum(set$p[1]/set$G[1] - set$p[1]*c(0:(set$G[1] - 1))/set$G[1]^2), lty=3)
# 
# plot(plot.f1[, 1], plot.f1[, 2], type="l", col=2, lty=1, ylim=range(plot.f1[, c(2, 5, 8, 11, 14, 17)]),
#      xlab="Number of selected variables", ylab="f1", main="f)")
# lines(plot.f1[, 1], plot.f1[, 5], col=3, lty=1)
# lines(plot.f1[, 1], plot.f1[, 8], col=4, lty=1)
# lines(plot.f1[, 1], plot.f1[, 11], col=5, lty=1)
# lines(plot.f1[, 1], plot.f1[, 14], col=6, lty=1)
# lines(plot.f1[, 1], plot.f1[, 17], col=7, lty=1)
# abline(v=sum(set$p[1]/set$G[1] - set$p[1]*c(0:(set$G[1] - 1))/set$G[1]^2), lty=3)
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
#   geom_vline(xintercept=sum(set$p[1]/set$G[1] - set$p[1]*c(0:(set$G[1] - 1))/set$G[1]^2), linetype="dotted", show.legend=TRUE) +
#   geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=1) + geom_line() + geom_point() +
#   xlab("Number of selected variables") + ylab("AUC") + ggtitle("a)") + theme_bw() +
#   theme(legend.position="none", plot.title=element_text(hjust=0.5))
# 
# mean.kappa <- ggplot(plot.kappa2, aes(x=psel, y=mean, colour=method)) +
#   geom_vline(xintercept=sum(set$p[1]/set$G[1] - set$p[1]*c(0:(set$G[1] - 1))/set$G[1]^2), linetype="dotted", show.legend=TRUE) +
#   geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=1) + geom_line() + geom_point() +
#   xlab("Number of selected variables") + ylab("Cohen's kappa") + ggtitle("b)") + theme_bw() +
#   theme(legend.position="none", plot.title=element_text(hjust=0.5))
# 
# mean.mse <- ggplot(plot.mse2, aes(x=psel, y=mean, colour=method)) +
#   geom_vline(xintercept=sum(set$p[1]/set$G[1] - set$p[1]*c(0:(set$G[1] - 1))/set$G[1]^2), linetype="dotted", show.legend=TRUE) +
#   geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=1) + geom_line() + geom_point() +
#   xlab("Number of selected variables") + ylab("MSE") + ggtitle("c)") + theme_bw() +
#   theme(legend.position="none", plot.title=element_text(hjust=0.5))
# 
# mean.prec <- ggplot(plot.prec2, aes(x=psel, y=mean, colour=method)) +
#   geom_vline(xintercept=sum(set$p[1]/set$G[1] - set$p[1]*c(0:(set$G[1] - 1))/set$G[1]^2), linetype="dotted", show.legend=TRUE) +
#   geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=1) + geom_line() + geom_point() +
#   xlab("Number of selected variables") + ylab("Precision") + ggtitle("d)") + theme_bw() +
#   theme(legend.position="none", plot.title=element_text(hjust=0.5))
# 
# mean.rec <- ggplot(plot.rec2, aes(x=psel, y=mean, colour=method)) +
#   geom_vline(xintercept=sum(set$p[1]/set$G[1] - set$p[1]*c(0:(set$G[1] - 1))/set$G[1]^2), linetype="dotted", show.legend=TRUE) +
#   geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=1) + geom_line() + geom_point() +
#   xlab("Number of selected variables") + ylab("Recall") + ggtitle("e)") + theme_bw() +
#   theme(legend.justification=c(0, 1), legend.position=c(0, 1), plot.title=element_text(hjust=0.5))
# 
# mean.f1 <- ggplot(plot.f12, aes(x=psel, y=mean, colour=method)) +
#   geom_vline(xintercept=set$p[1] - choose(set$G[1], 2)*set$p[1]/set$G[1]^2, linetype="dotted", show.legend=TRUE) +
#   geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=1) + geom_line() + geom_point() +
#   xlab("Number of selected variables") + ylab("F1 score") + ggtitle("f)") + theme_bw() +
#   theme(legend.position="none", plot.title=element_text(hjust=0.5))
# 
# png(paste(path.graph, "grEBEN_sim_varsel2_set1_mean.png", sep=""), width=1200, height=720, res=90)
# grid.arrange(mean.auc, mean.kappa, mean.mse, mean.prec, mean.rec, mean.f1, ncol=3)
# dev.off()




