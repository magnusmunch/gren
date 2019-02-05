##############################  preamble  #############################
# Putri's miRNA data                                                  #
# version: 01                                                         #
# author: Magnus M?nch                                                #
# created: 06-09-2017                                                 #
# last edited: 06-09-2017                                             #
#######################################################################

###############################  notes  ###############################
#                                                                     #
#######################################################################

# paths
path.code <- as.character(ifelse(Sys.info()[1]=="Darwin","/Users/magnusmunch/Documents/PhD/EBEN/code/" ,"~/EBEN/code/"))
path.graph <- "/Users/magnusmunch/Documents/PhD/EBEN/graphs/"
path.res <- as.character(ifelse(Sys.info()[1]=="Darwin","/Users/magnusmunch/Documents/PhD/EBEN/results/" ,"~/EBEN/results/"))
path.data <- as.character(ifelse(Sys.info()[1]=="Darwin","/Users/magnusmunch/Documents/PhD/EBEN/data/" ,"~/EBEN/data/"))

### libraries
library(GRridge)
library(edgeR)
library(pROC)
library(reshape)

### functions
# source function for variational Bayes
source(paste(path.code, "grVBEM.R", sep=""))

# version of grridge that doesn't give errors
source(paste(path.code, "mygrridge.R", sep=""))

# p-values calculation
logreguni <- function(dat, resp) {
  #dat <- mirblst[1,];resp <- respblood
  glm1 <- glm(resp ~ 1, family="binomial")
  glm2 <- glm(resp ~ 1 + dat, family="binomial")
  pval <- anova(glm1, glm2, test="Chisq")$Pr[2]
  return(pval)  
}

### setting seed for reproducibility
set.seed(1001)

### loading data
load(paste(path.data, "normIsoProstate.Rdata", sep=""))
load(paste(path.data, "ann2.Rdata", sep=""))
mirsfamily <- read.table(paste(path.data, "miR_Family_Info_edit.txt", sep=""), header=TRUE, 
                         fill=TRUE, as.is=c(1, 2, 4, 5, 7))

### data manipulation
# response
resp <- c(rep(0, 25), rep(1, 36)) # 0 is normal, 1 is prostate cancer
n <- length(resp)

# mirsData (combine normalized-isomirs)
mirs <- as.character(ann2$Name)
datMirNorm <- aggregate(normIso$datCount, by=list(mirs), sum) # p=879, n=62
datMirTrans <- sqrt(datMirNorm[, -1] + 3/8) - sqrt(3/8)
datMirStd <- t(scale(t(datMirTrans), center=TRUE, scale=TRUE))
rownames(datMirTrans) <- rownames(datMirStd) <- datMirNorm[, 1]
data <- t(datMirStd)

# conservation status of miRNAs, 0: not conserved, 1: conserved accross vertebrates, 2: conserved
mirsfamily$Family.Conservation.2 <- (mirsfamily$Family.Conservation. > 0)*mirsfamily$Family.Conservation.

# partitions
abundance <- rowSums(datMirNorm[, -1])
sds <- apply(datMirTrans, 1, sd)
conservation <- mirsfamily$Family.Conservation.2[match(colnames(data), mirsfamily$MiRBase.ID)] + 1

parAbund1 <- CreatePartition(abundance, ngroup=5, uniform=TRUE, decreasing=TRUE)
parSds1 <- CreatePartition(sds, decreasing=TRUE, uniform=TRUE, ngroup=5)
parCons1 <- CreatePartition(as.factor(conservation))

parAbund2 <- rep(1:length(parAbund1), times=unlist(lapply(parAbund1, length)))[order(unlist(parAbund1))]
parSds2 <- rep(1:length(parSds1), times=unlist(lapply(parSds1, length)))[order(unlist(parSds1))]
parCons2 <- conservation

partitions1.1 <- list(abundance=lapply(parAbund1, function(group) {
  ifelse(group > (ncol(data) - sum(is.na(abundance) | is.na(sds) | is.na(conservation))), 
         group - sum(is.na(abundance) | is.na(sds) | is.na(conservation)), group)}),
  sds=lapply(parSds1, function(group) {
    ifelse(group > (ncol(data) - sum(is.na(abundance) | is.na(sds) | is.na(conservation))), group - sum(is.na(abundance) | is.na(sds) | is.na(conservation)), group)}), 
  conservation=lapply(parCons1, function(group) {
    ifelse(group > (ncol(data) - sum(is.na(abundance) | is.na(sds) | is.na(conservation))), group - sum(is.na(abundance) | is.na(sds) | is.na(conservation)), group)}))
partitions1.2 <- list(abundance=lapply(parAbund1, function(group) {
  ifelse(group > (ncol(data) - sum(is.na(abundance))), group - sum(is.na(abundance)), group)}))
partitions1.3 <- list(conservation=lapply(parCons1, function(group) {
  ifelse(group > (ncol(data) - sum(is.na(conservation))), group - sum(is.na(conservation)), group)}))

partitions2.1 <- list(abundance=parAbund2[!(is.na(abundance) | is.na(sds) | is.na(conservation))],
                      sds=parSds2[!(is.na(abundance) | is.na(sds) | is.na(conservation))],
                      conservation=parCons2[!(is.na(abundance) | is.na(sds) | is.na(conservation))])
partitions2.2 <- list(abundance=parAbund2[!is.na(abundance)])
partitions2.3 <- list(conservation=parCons2[!is.na(conservation)])

# ### fit models
# fit1.GRridge <- grridge(t(data[, !(is.na(abundance) | is.na(sds) | is.na(conservation))]), resp, partitions1.1,
#                         monotone=c(TRUE, TRUE, FALSE))
# fit2.GRridge <- grridge(t(data[, !is.na(abundance)]), resp, partitions1.2, monotone=TRUE)
# fit3.GRridge <- grridge(t(data[, !is.na(conservation)]), resp, partitions1.3, monotone=FALSE)
# fit4.GRridge <- grridge(t(data[, !(is.na(abundance) | is.na(sds) | is.na(conservation))]), resp, partitions1.1,
#                         monotone=c(FALSE, FALSE, FALSE))
# fit5.GRridge <- grridge(t(data[, !is.na(abundance)]), resp, partitions1.2, monotone=TRUE,
#                         selectionEN=TRUE, maxsel=10)
# 
# fit1.grEBEN <- grEBEN(data[, !(is.na(abundance) | is.na(sds) | is.na(conservation))], resp, rep(1, n),
#                       unpenalized=NULL, intercept=TRUE, partitions=partitions2.1, lambda1=NULL,
#                       lambda2=NULL, monotone=list(monotone=c(TRUE, TRUE, FALSE), decreasing=c(FALSE, FALSE, FALSE)),
#                       psel=NULL, posterior=FALSE, ELBO=TRUE, eps=0.001, maxiter=1000, trace=TRUE)
# fit2.grEBEN <- grEBEN(data[, !is.na(abundance)], resp, rep(1, n), unpenalized=NULL, intercept=TRUE,
#                       partitions=partitions2.2, lambda1=NULL, lambda2=NULL,
#                       monotone=list(monotone=FALSE, decreasing=FALSE),
#                       psel=NULL, posterior=FALSE, ELBO=TRUE, eps=0.001, maxiter=500, trace=TRUE)
# fit3.grEBEN <- grEBEN(data[, !is.na(conservation)], resp, rep(1, n), unpenalized=NULL, intercept=TRUE,
#                       partitions=partitions2.3, lambda1=NULL, lambda2=NULL,
#                       monotone=list(monotone=FALSE, decreasing=FALSE),
#                       psel=NULL, posterior=FALSE, ELBO=TRUE, eps=0.001, maxiter=500, trace=TRUE)
# fit4.grEBEN <- grEBEN(data[, !(is.na(abundance) | is.na(sds) | is.na(conservation))], resp, rep(1, n),
#                       unpenalized=NULL, intercept=TRUE, partitions=partitions2.1, lambda1=NULL,
#                       lambda2=NULL, monotone=list(monotone=c(FALSE, FALSE, FALSE), decreasing=c(FALSE, FALSE, FALSE)),
#                       psel=NULL, posterior=FALSE, ELBO=TRUE, eps=0.001, maxiter=500, trace=TRUE)
# fit5.grEBEN <- grEBEN(data[, !is.na(abundance)], resp, rep(1, n), unpenalized=NULL, intercept=TRUE,
#                       partitions=partitions2.2, lambda1=NULL, lambda2=NULL,
#                       monotone=list(monotone=FALSE, decreasing=FALSE),
#                       psel=10, posterior=FALSE, ELBO=TRUE, eps=0.001, maxiter=500, trace=TRUE)
# 
# fit.enet <- glmnet(data, resp, family="binomial", intercept=TRUE, standardize=FALSE, lambda=0.00000001)
# fit.lasso <- glmnet(data, resp, family="binomial", intercept=TRUE, standardize=FALSE, alpha=1, dfmax=10)

# # save the models
# save(fit1.GRridge, fit2.GRridge, fit3.GRridge, fit4.GRridge, fit1.grEBEN, fit2.grEBEN, fit3.grEBEN,
#      fit4.grEBEN, fit.enet, file=paste(path.res, "grEBEN_mirna_putri_fitted1.Rdata", sep=""))


# # cross-validate predictive performance for only conservation status
# set.seed(1001)
# data2 <- data[, !is.na(conservation)]
# p <- ncol(data2)
# n <- nrow(data2)
# glob.pen1 <- cv.pen(data2, resp, unpenalized=NULL, intercept=TRUE, psel=NULL)
# glob.pen2 <- cv.glmnet(data2, resp, family="binomial", alpha=0, standardize=FALSE, intercept=TRUE, 
#                        nfolds=n, grouped=FALSE)
# 
# lambda1 <- glob.pen1$lambda1bayes[which.min(glob.pen1$cvll)]
# lambda2 <- glob.pen1$lambda2bayes[which.min(glob.pen1$cvll)]
# 
# nfolds <- n
# rest <- n %% nfolds
# foldid <- sample(rep(1:nfolds, times=round(c(rep(n %/% nfolds + as.numeric(rest!=0), times=rest), 
#                                              rep(n %/% nfolds, times=nfolds - rest)))))
# 
# pred <- matrix(NA, ncol=5, nrow=n)
# colnames(pred) <- c("GRridge", "grEBEN", "grEBEN+VB","ridge", "enet")
# for(k in 1:nfolds) {
#   print(paste("Fold", k, sep=" "))
#   xtrain <- data2[foldid!=k, ]
#   xtest <- data2[foldid==k, ]
#   ytrain <- resp[foldid!=k]
#   ytest <- resp[foldid==k]
#   ntrain <- length(ytrain)
#   ntest <- length(ytest)
# 
#   fit.grEBEN <- grEBEN(xtrain, ytrain, rep(1, ntrain), partitions=partitions2.3, lambda1=lambda1, 
#                        lambda2=lambda2, ELBO=FALSE, trace=FALSE)
#   fit.GRridge <- grridge(t(xtrain), ytrain, partitions1.3, optl=glob.pen2$lambda.min*0.5*n, 
#                          trace=FALSE)
# 
#   pred[foldid==k, 1] <- predict(fit.GRridge$predobj$GroupRegul, matrix(xtest, ncol=p))
#   pred[foldid==k, 2] <- predict.grEBEN(fit.grEBEN, matrix(xtest, ncol=p), type="penalized")
#   pred[foldid==k, 3] <- predict.grEBEN(fit.grEBEN, matrix(xtest, ncol=p), type="VB")
#   pred[foldid==k, 4] <- predict(fit.GRridge$predobj$NoGroups, matrix(xtest, ncol=p))
#   pred[foldid==k, 5] <- predict.grEBEN(fit.grEBEN, matrix(xtest, ncol=p), type="nogroups")
# 
# }
# 
# auc <- sapply(1:ncol(pred), function(method) {return(pROC::roc(resp, pred[, method])$auc)})
# names(auc) <- colnames(pred)
# setting1 <- list(pred=pred, auc=auc)
# save(setting1, file=paste(path.res, "grEBEN_mirna_putri_setting1.Rdata", sep=""))



# # cross-validate predictive performance for only abundance
# set.seed(1001)
# data2 <- data[, !is.na(abundance)]
# p <- ncol(data2)
# n <- nrow(data2)
# glob.pen1 <- cv.pen(data2, resp, unpenalized=NULL, intercept=TRUE, psel=NULL)
# glob.pen2 <- cv.glmnet(data2, resp, family="binomial", alpha=0, standardize=FALSE, intercept=TRUE,
#                        nfolds=n, grouped=FALSE)
# 
# lambda1 <- glob.pen1$lambda1bayes[which.min(glob.pen1$cvll)]
# lambda2 <- glob.pen1$lambda2bayes[which.min(glob.pen1$cvll)]
# 
# nfolds <- n
# rest <- n %% nfolds
# foldid <- sample(rep(1:nfolds, times=round(c(rep(n %/% nfolds + as.numeric(rest!=0), times=rest),
#                                              rep(n %/% nfolds, times=nfolds - rest)))))
# 
# pred <- matrix(NA, ncol=5, nrow=n)
# methods <- c("GRridge", "grEBEN", "grEBEN+VB","ridge", "enet")
# pred <- matrix(NA, ncol=length(methods), nrow=n)
# for(k in 1:nfolds) {
# 
#   print(paste("Fold", k, sep=" "))
#   xtrain <- data2[foldid!=k, ]
#   xtest <- data2[foldid==k, ]
#   ytrain <- resp[foldid!=k]
#   ytest <- resp[foldid==k]
#   ntrain <- length(ytrain)
#   ntest <- length(ytest)
# 
#   fit.grEBEN <- grEBEN(xtrain, ytrain, rep(1, ntrain), partitions=partitions2.2, lambda1=lambda1,
#                        lambda2=lambda2, ELBO=FALSE, trace=FALSE)
#   fit.GRridge <- grridge(t(xtrain), ytrain, partitions1.2, optl=glob.pen2$lambda.min*0.5*n,
#                          trace=FALSE)
# 
#   pred[foldid==k, 1] <- predict(fit.GRridge$predobj$GroupRegul, matrix(xtest, ncol=p))
#   pred[foldid==k, 2] <- predict.grEBEN(fit.grEBEN, matrix(xtest, ncol=p), type="penalized")
#   pred[foldid==k, 3] <- predict.grEBEN(fit.grEBEN, matrix(xtest, ncol=p), type="VB")
#   pred[foldid==k, 4] <- predict(fit.GRridge$predobj$NoGroups, matrix(xtest, ncol=p))
#   pred[foldid==k, 5] <- predict.grEBEN(fit.grEBEN, matrix(xtest, ncol=p), type="nogroups")
# 
# }
# auc <- sapply(1:length(methods), function(method) {return(pROC::roc(resp, pred[, method])$auc)})
# names(auc) <- colnames(pred) <- methods
# setting2 <- list(pred=pred, auc=auc)
# save(setting2, file=paste(path.res, "grEBEN_mirna_putri_setting2.Rdata", sep=""))





# # cross-validate predictive performance for only abundance with variable selection
# set.seed(1001)
# data2 <- data[, !is.na(abundance)]
# p <- ncol(data2)
# n <- nrow(data2)
# 
# pselseq <- c(5, 10, 15, 20, 30, 40, 50, 75, 100)
# nfolds <- n
# rest <- n %% nfolds
# foldid <- sample(rep(1:nfolds, times=round(c(rep(n %/% nfolds + as.numeric(rest!=0), times=rest),
#                                              rep(n %/% nfolds, times=nfolds - rest)))))
# 
# methods <- c("GRridge", "GRridge+sel", "grEBEN", "grEBEN+VB", "grEBEN+sel", "ridge", "enet")
# pred <- matrix(NA, ncol=length(methods), nrow=n*length(pselseq))
# selmat <- matrix(NA, ncol=length(methods), nrow=nfolds*length(pselseq))
# auc <- cbind(pselseq, matrix(NA, ncol=length(methods), nrow=length(pselseq)))
# colnames(auc)[-1] <- colnames(selmat) <- colnames(pred) <- methods
# for(csel in 1:length(pselseq)) {
#   psel <- pselseq[csel]
# 
#   glob.pen1 <- cv.pen(data2, resp, unpenalized=NULL, intercept=TRUE, psel=psel)
#   glob.pen2 <- mygrridge(t(data2), resp, partitions1.2, trace=FALSE, selectionEN=TRUE, maxsel=psel)
# 
#   lambda1 <- glob.pen1$lambda1bayes[which.min(glob.pen1$cvll)]
#   lambda2 <- glob.pen1$lambda2bayes[which.min(glob.pen1$cvll)]
# 
#   for(k in 1:nfolds) {
# 
#     print(paste("Psel ", psel, ", Fold ", k, sep=""))
#     xtrain <- data2[foldid!=k, ]
#     xtest <- data2[foldid==k, ]
#     ytrain <- resp[foldid!=k]
#     ytest <- resp[foldid==k]
#     ntrain <- length(ytrain)
#     ntest <- length(ytest)
# 
#     fit.grEBEN <- grEBEN(xtrain, ytrain, rep(1, ntrain), partitions=partitions2.2, lambda1=lambda1,
#                          lambda2=lambda2, ELBO=FALSE, trace=FALSE, psel=psel)
#     fit.GRridge <- mygrridge(t(xtrain), ytrain, partitions1.2, optl=glob.pen2$optl, trace=FALSE,
#                              selectionEN=TRUE, maxsel=psel)
# 
#     pred[(csel - 1)*n + which(foldid==k), 1] <- predict.grridge(fit.GRridge, t(matrix(xtest, ncol=p)))[, 2]
#     pred[(csel - 1)*n + which(foldid==k), 2] <- predict.grridge(fit.GRridge, t(matrix(xtest, ncol=p)))[, 3]
#     pred[(csel - 1)*n + which(foldid==k), 3] <- predict.grEBEN(fit.grEBEN, matrix(xtest, ncol=p), type="penalized")
#     pred[(csel - 1)*n + which(foldid==k), 4] <- predict.grEBEN(fit.grEBEN, matrix(xtest, ncol=p), type="VB")
#     pred[(csel - 1)*n + which(foldid==k), 5] <- predict.grEBEN(fit.grEBEN, matrix(xtest, ncol=p), type="selection")
#     pred[(csel - 1)*n + which(foldid==k), 6] <- predict.grridge(fit.GRridge, t(matrix(xtest, ncol=p)))[, 1]
#     pred[(csel - 1)*n + which(foldid==k), 7] <- predict.grEBEN(fit.grEBEN, matrix(xtest, ncol=p), type="nogroups")
# 
#     selmat[(csel - 1)*nfolds + k, 1] <- p
#     selmat[(csel - 1)*nfolds + k, 2] <- length(fit.GRridge$resEN$whichEN)
#     selmat[(csel - 1)*nfolds + k, 3] <- sum(fit.grEBEN$beta[-1]!=0)
#     selmat[(csel - 1)*nfolds + k, 4] <- p
#     selmat[(csel - 1)*nfolds + k, 5] <- sum(fit.grEBEN$beta.sel[-1]!=0)
#     selmat[(csel - 1)*nfolds + k, 6] <- p
#     selmat[(csel - 1)*nfolds + k, 7] <- sum(fit.grEBEN$beta.nogroups[-1]!=0)
# 
# 
#   }
# 
#   auc[csel, -1] <- apply(pred[c(((csel - 1)*n + 1):(csel*n)), ], 2, function(cpred) {
#     pROC::roc(resp, cpred, na.rm=TRUE)$auc})
#   setting3 <- list(pred=pred, auc=auc, psel=selmat)
#   save(setting3, file=paste(path.res, "grEBEN_mirna_putri_setting3.Rdata", sep=""))
# 
# }





load(paste(path.res, "grEBEN_mirna_putri_fitted1.Rdata", sep=""))

# create partitions based on the regular elastic net and pvalues
pvalues <- apply(data, 2, logreguni, resp=resp)
benet <- as.numeric(coef(fit.enet, s="lambda.min"))[-1]
parPval1 <- CreatePartition(pvalues, ngroup=5, uniform=TRUE)
parPval2 <- rep(1:length(parPval1), times=unlist(lapply(parPval1, length)))[order(unlist(parPval1))]
parEst1 <- CreatePartition(benet, ngroup=5, uniform=TRUE, decreasing=TRUE)
parEst2 <- rep(1:length(parEst1), times=unlist(lapply(parEst1, length)))[order(unlist(parEst1))]

# plot p-values and elastic net estimates
plot(log(pvalues[abs(benet) < 50]) ~ benet[abs(benet) < 50])
boxplot(log(pvalues) ~ parEst2) # against elastic net estimates
boxplot(log(pvalues) ~ parAbund2) # against abundances
boxplot(benet[abs(benet) < 50] ~ parCons2[abs(benet) < 50]) # against conservation status
boxplot(log(abs(benet[abs(benet) < 50]) + 1) ~ parCons2[abs(benet) < 50])

# plot variance of elastic net estimates against estimated lambdag
var.enest <- sapply(1:length(unique(parCons2)[!is.na(unique(parCons2))]), function(cons) {
  var(as.numeric(coef(fit.enet, s="lambda.min"))[-1][parCons2==cons], na.rm=TRUE)})
plot(var.enest, fit3.grEBEN$lambdag$conservation[, fit3.grEBEN$nouteriter + 1],
     ylab=expression(lambda[g]), xlab=expression(paste(hat(V), "(", hat(beta)[EN], ")", sep="")))
lambdag <- exp(-log(var.enest) + sum(rle(sort(conservation))$lengths*log(var.enest))/sum(!is.na(conservation)))
plot(lambdag, fit3.grEBEN$lambdag$conservation[, fit3.grEBEN$nouteriter + 1])
plot(lambdag, fit3.GRridge$lambdamults$conservation)
plot(fit3.grEBEN$lambdag$conservation[, fit3.grEBEN$nouteriter + 1],
     fit3.GRridge$lambdamults$conservation)

# first model
colvec <- heat.colors(5)
plot(fit1.grEBEN$lambdag$abundance[1, ], type="l", ylim=range(fit1.grEBEN$lambdag$abundance), col=colvec[1], lwd=2,
     ylab=expression(paste(lambda[g], "'")), xlab="Iteration")
for(i in 2:nrow(fit1.grEBEN$lambdag$abundance)) {
  lines(fit1.grEBEN$lambdag$abundance[i, ], col=colvec[i], lwd=2)
}
legend("topright", legend=c("large", "small") ,lty=1, lwd=3, col=range(colvec), title="Abundance")

plot(fit1.grEBEN$lambdag$sds[1, ], type="l", ylim=range(fit1.grEBEN$lambdag$sds), col=colvec[1], lwd=2,
     ylab=expression(paste(lambda[g], "'")), xlab="Iteration")
for(i in 2:nrow(fit1.grEBEN$lambdag$sds)) {
  lines(fit1.grEBEN$lambdag$sds[i, ], col=colvec[i], lwd=2)
}
legend("topleft", legend=c("large", "small") ,lty=1, lwd=3, col=range(colvec), title="Standard deviation")

plot(fit1.grEBEN$lambdag$conservation[1, ], type="l", ylim=range(fit1.grEBEN$lambdag$conservation),
     ylab=expression(paste(lambda[g], "'")), xlab="Iteration")
lines(fit1.grEBEN$lambdag$conservation[2, ], col=2)
lines(fit1.grEBEN$lambdag$conservation[3, ], col=3)
legend("topleft", legend=c("not conserved", "in vertebrates", "conserved") ,lty=1, col=c(1:3))

barplot(rbind(fit1.GRridge$lambdamults$abundance,
              fit1.grEBEN$lambdag$abundance[, fit1.grEBEN$nouteriter + 1]), beside=TRUE,
        xlab="", axisnames=FALSE,
        args.legend=list(x="topleft", fill=c(gray.colors(2), 0, 0), lty=c(NA, NA, 2, 2),
                         border=c(rep(1, 2), 0, 0), merge=TRUE, seg.len=1),
        ylab=expression(paste(lambda[g], "'")),
        legend.text=c("GRridge", "grEBEN", "ridge", "enet"))
title(xlab="Abundances in decreasing order", line=1)
abline(h=1, lty=2)

barplot(rbind(fit1.GRridge$lambdamults$sds,
              fit1.grEBEN$lambdag$sds[, fit1.grEBEN$nouteriter + 1]), beside=TRUE,
        xlab="", axisnames=FALSE,
        args.legend=list(x="topleft", fill=c(gray.colors(2), 0, 0), lty=c(NA, NA, 2, 2),
                         border=c(rep(1, 2), 0, 0), merge=TRUE, seg.len=1),
        ylab=expression(paste(lambda[g], "'")),
        legend.text=c("GRridge", "grEBEN", "ridge", "enet"))
title(xlab="Standard deviations in decreasing order", line=1)
abline(h=1, lty=2)

barplot(rbind(fit1.GRridge$lambdamults$conservation,
              fit1.grEBEN$lambdag$conservation[, fit1.grEBEN$nouteriter + 1]), beside=TRUE,
        xlab="", names.arg=c("not \n conserved", "conserved \n in mammals", "broadly \n conserved"),
        args.legend=list(x="topleft", fill=c(gray.colors(2), 0, 0), lty=c(NA, NA, 2, 2),
                         border=c(rep(1, 2), 0, 0), merge=TRUE, seg.len=1),
        ylab=expression(paste(lambda[g], "'")),
        legend.text=c("GRridge", "grEBEN", "ridge", "enet"))
abline(h=1, lty=2)

### Only abundance
load(paste(path.res, "grEBEN_mirna_putri_setting2.Rdata", sep=""))
colvec <- heat.colors(5)
plot(fit2.grEBEN$lambdag$abundance[1, ], type="l", ylim=range(fit2.grEBEN$lambdag$abundance), col=colvec[1], lwd=2)
for(i in 2:nrow(fit2.grEBEN$lambdag$abundance)) {
  lines(fit2.grEBEN$lambdag$abundance[i, ], col=colvec[i], lwd=2)
}
legend("topleft", legend=c("large", "small") ,lty=1, lwd=3, col=range(colvec), title="Abundance")

png(paste(path.graph, "grEBEN_mirna_putri_bar2.png", sep=""), units="in", width=6, height=5.33, res=200)
barplot(rbind(fit2.GRridge$lambdamults$abundance,
              fit2.grEBEN$lambdag$abundance[, fit2.grEBEN$nouteriter + 1]), beside=TRUE,
        xlab="", axisnames=FALSE,
        args.legend=list(x="topleft", fill=c(gray.colors(2), 0, 0), lty=c(NA, NA, 2, 2),
                         border=c(rep(1, 2), 0, 0), merge=TRUE, seg.len=1),
        ylab=expression(paste(lambda[g], "'")),
        legend.text=paste(c("GRridge", "grEBEN", "ridge", "enet"), ", AUC=",
                          round(setting2$auc[-3], 2), sep=""))
title(xlab="Abundances in decreasing order", line=1)
abline(h=1, lty=2)
dev.off()

brier <- apply(setting2$pred, 2, function(x) {mean((resp - x)^2)})
par(mfrow=c(2, 3))
plot(sort(setting2$pred[, 1]), col=resp[order(setting2$pred[, 1])] + 1, main="GRridge",
     ylab="Predicted probability")
legend("topleft", legend=paste("Brier score:", round(brier[1], 2)))

plot(sort(setting2$pred[, 2]), col=resp[order(setting2$pred[, 2])] + 1, main="grEBEN",
     ylab="Predicted probability")
legend("topleft", legend=paste("Brier score:", round(brier[2], 2)))

plot(sort(setting2$pred[, 3]), col=resp[order(setting2$pred[, 3])] + 1, main="grEBEN VB",
     ylab="Predicted probability")
legend("topleft", legend=paste("Brier score:", round(brier[3], 2)))

plot(sort(setting2$pred[, 4]), col=resp[order(setting2$pred[, 4])] + 1, main="Ridge",
     ylab="Predicted probability")
legend("topleft", legend=paste("Brier score:", round(brier[4], 2)))

plot(sort(setting2$pred[, 5]), col=resp[order(setting2$pred[, 5])] + 1, main="Elastic net",
     ylab="Predicted probability")
legend("topleft", legend=paste("Brier score:", round(brier[5], 2)))
par(mfrow=c(1, 1))

# with selection
load(paste(path.res, "grEBEN_mirna_putri_setting3.Rdata", sep=""))
setting3$brier <- cbind(psel=setting3$auc[, 1], t(sapply(1:nrow(setting3$auc), function(csel) {
  apply(setting3$pred[c(((csel - 1)*n + 1):(csel*n)), ], 2, function(x) {mean((resp - x)^2)})})))
png(paste(path.graph, "grEBEN_mirna_putri_select2.png", sep=""), units="in", width=12, height=5.33, res=200)
par(mfrow=c(1, 2))
plot(setting3$auc[!is.na(setting3$auc[, 2]), 1], setting3$auc[!is.na(setting3$auc[, 2]), 2], col=2, 
     ylim=range(setting3$auc[, -1], na.rm=TRUE), type="l", xlab="Number of selected variables",
     ylab="AUC", main="a)")
lines(setting3$auc[!is.na(setting3$auc[, 3]), 1], setting3$auc[!is.na(setting3$auc[, 3]), 3], col=3)
lines(setting3$auc[!is.na(setting3$auc[, 4]), 1], setting3$auc[!is.na(setting3$auc[, 4]), 4], col=4)
lines(setting3$auc[!is.na(setting3$auc[, 6]), 1], setting3$auc[!is.na(setting3$auc[, 6]), 6], col=5)
lines(setting3$auc[!is.na(setting3$auc[, 7]), 1], setting3$auc[!is.na(setting3$auc[, 7]), 7], col=6)
lines(setting3$auc[!is.na(setting3$auc[, 8]), 1], setting3$auc[!is.na(setting3$auc[, 8]), 8], col=7)
legend("bottomright", legend=colnames(setting3$auc)[-c(1, 5)], lty=1, col=2:7)

plot(setting3$mse[!is.na(setting3$mse[, 2]), 1], setting3$mse[!is.na(setting3$mse[, 2]), 2], col=2, 
     ylim=range(setting3$mse[, -1], na.rm=TRUE), type="l", xlab="Number of selected variables",
     ylab="Mean brier score", main="b)")
lines(setting3$mse[!is.na(setting3$mse[, 3]), 1], setting3$mse[!is.na(setting3$mse[, 3]), 3], col=3)
lines(setting3$mse[!is.na(setting3$mse[, 4]), 1], setting3$mse[!is.na(setting3$mse[, 4]), 4], col=4)
lines(setting3$mse[!is.na(setting3$mse[, 6]), 1], setting3$mse[!is.na(setting3$mse[, 6]), 6], col=5)
lines(setting3$mse[!is.na(setting3$mse[, 7]), 1], setting3$mse[!is.na(setting3$mse[, 7]), 7], col=6)
lines(setting3$mse[!is.na(setting3$mse[, 8]), 1], setting3$mse[!is.na(setting3$mse[, 8]), 8], col=7)
dev.off()

### Only conservation status
load(paste(path.res, "grEBEN_mirna_putri_setting1.Rdata", sep=""))
plot(fit3.grEBEN$lambdag$conservation[1, ], type="l", ylim=range(fit3.grEBEN$lambdag$conservation))
lines(fit3.grEBEN$lambdag$conservation[2, ], col=2)
lines(fit3.grEBEN$lambdag$conservation[3, ], col=3)
legend("topright", legend=c("not conserved", "vertebrates", "conserved"), lty=1, col=1:3)

png(paste(path.graph, "grEBEN_mirna_putri_bar3.png", sep=""), units="in", width=6, height=5.33, res=200)
barplot(rbind(fit3.GRridge$lambdamults$conservation,
              fit3.grEBEN$lambdag$conservation[, fit3.grEBEN$nouteriter + 1]), beside=TRUE,
        xlab="", names.arg=c("not \n conserved", "conserved \n in mammals", "broadly \n conserved"),
        ylab=expression(paste(lambda[g], "'")),
        args.legend=list(x="topright", fill=c(gray.colors(2), 0, 0), lty=c(NA, NA, 2, 2),
                         border=c(rep(1, 2), 0, 0), merge=TRUE, seg.len=1),
        legend.text=paste(c("GRridge", "grEBEN", "ridge", "enet"), ", AUC=",
                          round(setting1$auc[-3], 2), sep=""))
abline(h=1, lty=2)
dev.off()

# fourth model
colvec <- heat.colors(5)
plot(fit4.grEBEN$lambdag$abundance[1, ], type="l", ylim=range(fit4.grEBEN$lambdag$abundance), col=colvec[1], lwd=2,
     ylab=expression(paste(lambda[g], "'")), xlab="Iteration")
for(i in 2:nrow(fit4.grEBEN$lambdag$abundance)) {
  lines(fit4.grEBEN$lambdag$abundance[i, ], col=colvec[i], lwd=2)
}
legend("topright", legend=c("large", "small") ,lty=1, lwd=3, col=range(colvec), title="Abundance")

plot(fit4.grEBEN$lambdag$sds[1, ], type="l", ylim=range(fit4.grEBEN$lambdag$sds), col=colvec[1], lwd=2,
     ylab=expression(paste(lambda[g], "'")), xlab="Iteration")
for(i in 2:nrow(fit4.grEBEN$lambdag$sds)) {
  lines(fit4.grEBEN$lambdag$sds[i, ], col=colvec[i], lwd=2)
}
legend("topleft", legend=c("large", "small") ,lty=1, lwd=3, col=range(colvec), title="Standard deviation")

plot(fit4.grEBEN$lambdag$conservation[1, ], type="l", ylim=range(fit4.grEBEN$lambdag$conservation),
     ylab=expression(paste(lambda[g], "'")), xlab="Iteration")
lines(fit4.grEBEN$lambdag$conservation[2, ], col=2)
lines(fit4.grEBEN$lambdag$conservation[3, ], col=3)
legend("topleft", legend=c("not conserved", "in vertebrates", "conserved") ,lty=1, col=c(1:3))

barplot(rbind(fit4.GRridge$lambdamults$abundance,
              fit4.grEBEN$lambdag$abundance[, fit4.grEBEN$nouteriter + 1]), beside=TRUE,
        xlab="", axisnames=FALSE,
        args.legend=list(x="topleft", fill=c(gray.colors(2), 0, 0), lty=c(NA, NA, 2, 2),
                         border=c(rep(1, 2), 0, 0), merge=TRUE, seg.len=1),
        ylab=expression(paste(lambda[g], "'")),
        legend.text=c("GRridge", "grEBEN", "ridge", "enet"))
title(xlab="Abundances in decreasing order", line=1)
abline(h=1, lty=2)

barplot(rbind(fit4.GRridge$lambdamults$sds,
              fit4.grEBEN$lambdag$sds[, fit4.grEBEN$nouteriter + 1]), beside=TRUE,
        xlab="", axisnames=FALSE,
        args.legend=list(x="topleft", fill=c(gray.colors(2), 0, 0), lty=c(NA, NA, 2, 2),
                         border=c(rep(1, 2), 0, 0), merge=TRUE, seg.len=1),
        ylab=expression(paste(lambda[g], "'")),
        legend.text=c("GRridge", "grEBEN", "ridge", "enet"))
title(xlab="Standard deviations in decreasing order", line=1)
abline(h=1, lty=2)

barplot(rbind(fit4.GRridge$lambdamults$conservation,
              fit4.grEBEN$lambdag$conservation[, fit4.grEBEN$nouteriter + 1]), beside=TRUE,
        xlab="", names.arg=c("not \n conserved", "conserved \n in mammals", "broadly \n conserved"),
        args.legend=list(x="topleft", fill=c(gray.colors(2), 0, 0), lty=c(NA, NA, 2, 2),
                         border=c(rep(1, 2), 0, 0), merge=TRUE, seg.len=1),
        ylab=expression(paste(lambda[g], "'")),
        legend.text=c("GRridge", "grEBEN", "ridge", "enet"))
abline(h=1, lty=2)

# first versus third for conservation status
barplot(rbind(fit1.grEBEN$lambdag$conservation[, fit1.grEBEN$nouteriter + 1],
              fit3.grEBEN$lambdag$conservation[, fit3.grEBEN$nouteriter + 1]), beside=TRUE,
        xlab="", names.arg=c("not \n conserved", "conserved \n in mammals", "broadly \n conserved"),
        args.legend=list(x="topleft", fill=c(gray.colors(2), 0, 0), lty=c(NA, NA, 2, 2),
                         border=c(rep(1, 2), 0, 0), merge=TRUE, seg.len=1),
        ylab=expression(paste(lambda[g], "'")),
        legend.text=c("All partitions included", "Only conservation", "ridge", "enet"))
abline(h=1, lty=2)

# checking monotonicity
barplot(rbind(fit1.grEBEN$lambdag$abundance[, fit1.grEBEN$nouteriter + 1],
              fit4.grEBEN$lambdag$abundance[, fit4.grEBEN$nouteriter + 1]), beside=TRUE,
        xlab="", axisnames=FALSE,
        args.legend=list(x="topleft", fill=c(gray.colors(2), 0, 0), lty=c(NA, NA, 2, 2),
                         border=c(rep(1, 2), 0, 0), merge=TRUE, seg.len=1),
        ylab=expression(paste(lambda[g], "'")),
        legend.text=c("Monotonicity enforced", "Not enforced", "ridge", "enet"))
title(xlab="Abundances in decreasing order", line=1)
abline(h=1, lty=2)

barplot(rbind(fit1.grEBEN$lambdag$sds[, fit1.grEBEN$nouteriter + 1],
              fit4.grEBEN$lambdag$sds[, fit4.grEBEN$nouteriter + 1]), beside=TRUE,
        xlab="", axisnames=FALSE,
        args.legend=list(x="topleft", fill=c(gray.colors(2), 0, 0), lty=c(NA, NA, 2, 2),
                         border=c(rep(1, 2), 0, 0), merge=TRUE, seg.len=1),
        ylab=expression(paste(lambda[g], "'")),
        legend.text=c("Monotonicity enforced", "Not enforced", "ridge", "enet"))
title(xlab="Standard deviations in decreasing order", line=1)
abline(h=1, lty=2)
