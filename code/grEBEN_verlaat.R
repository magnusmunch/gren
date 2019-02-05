##############################  preamble  #############################
# grMCEM on real data                                                 #
# version: 01                                                         #
# author: Magnus M?nch                                                #
# created: 21-03-2017                                                 #
# last edited: 21-03-2017                                             #
#######################################################################

###############################  notes  ###############################
# 21-03-2017: Results are used in presentation for channel meeting    #
#             2017                                                    #
#######################################################################

##########################  data description  #########################
# A deep sequencing analysis on small non-coding ribonucleic acid     #
# (miRNAseq) was performed on 56 samples (24 women with high-grade    #
# cervical intraepithelial neoplasia (CIN3) and 32 healthy women) for #
# the purpose of finding relevant screening markers for cervical      #
# cancer screening. The next generation sequencing analysis resulted  #
# in 2,576 transcripts. The data was normalized and pre-processed,    #
# rendering 772 transcripts. More detail description of the data sets #
# and the preprocessing step are available on the supplementary       #
# material of this following publication "Better diagnostic           #
# signatures from RNAseq data through use of auxiliary co-data",      #
# Bioinformatics (2017).                                              #
#######################################################################

### paths
path.code <- as.character(ifelse(Sys.info()[1]=="Darwin","/Users/magnusmunch/Documents/PhD/EBEN/code/" ,"~/EBEN/code/"))
path.res <- as.character(ifelse(Sys.info()[1]=="Darwin","/Users/magnusmunch/Documents/PhD/EBEN/results/" ,"~/EBEN/results/"))
path.data <- as.character(ifelse(Sys.info()[1]=="Darwin","/Users/magnusmunch/Documents/PhD/EBEN/data/" ,"~/EBEN/data/"))
path.graph <- "/Users/magnusmunch/Documents/PhD/EBEN/graphs/"

### libraries
library(glmnet)
library(penalized)
library(GRridge)
library(pROC)

### functions
# source grENVB functions
source(paste(path.code, "grVBEM.R", sep=""))

# function to cross-validate elastic net
cv.en <- function(x, y, intercept=intercept) {
  fit.pen <- cv.pen(x, y, intercept=intercept)
  fit.en <- glmnet(x, y, family="binomial", alpha=fit.pen$alpha[which.min(fit.pen$cvll)],
                   lambda=fit.pen$lambda[which.min(fit.pen$cvll)])
  return(fit.en)
}

### loading data
data(dataVerlaat)

### partitions
annotation <- as.numeric(CpGann)
pvalues <- CreatePartition(pvalFarkas, decreasing=FALSE, mingr=10, ngroup=10)
pvalues <- rep(1:length(pvalues), times=unlist(lapply(pvalues, length)))[order(unlist(pvalues))]

### using annotation as grouping information
partitions <- list(annotation=annotation[order(annotation)])
x <- apply(t(datcenVerlaat), 2, function(x) {(x - mean(x))/sd(x)})[, order(annotation)]
y <- respVerlaat
n <- nrow(x)
p <- ncol(x)
m <- rep(1, n)

# # cross validation of AUC
# set.seed(5001)
# nfolds <- 10
# rest <- n %% nfolds
# foldsize <- c(rep(n %/% nfolds + as.numeric(rest!=0), times=rest),
#               rep(n %/% nfolds, times=nfolds - rest))
# foldid <- sample(rep(1:nfolds, times=foldsize))
#
# cv.ridge <- cv.glmnet(x, y, family="binomial", alpha=0, standardize=FALSE, intercept=TRUE)
# cv.en <- cv.pen(x, y, intercept=TRUE)
# # cv.SGL <- cvSGL(list(x=x, y=y), index=groups, type="logit", standardize=FALSE)
#
# lambda2ridge <- 0.5*n*cv.ridge$lambda.min
# lambdaridge <- cv.ridge$lambda.min
# alphaglmnet <- cv.en$alpha[which.min(cv.en$cvll)]
# lambdaglmnet <- cv.en$lambda[which.min(cv.en$cvll)]
# lambda1gren <- 2*n*lambdaglmnet*alphaglmnet
# lambda2gren <- n*lambdaglmnet*(1 - alphaglmnet)
# # lambdaSGL <- cv.SGL$lambdas[which.min(cv.SGL$lldiff)]
#
# pred1.ridge <- pred1.en <- pred1.grEBEN <- pred1.grridge <- numeric(n) # pred1.SGL <- numeric(n)
# for(k in 1:nfolds) {
#   print(paste("Fold", k, sep=" "))
#   xtrain <- x[foldid!=k, ]
#   xtest <- x[foldid==k, ]
#   ytrain <- y[foldid!=k]
#   ytest <- y[foldid==k]
#   mtrain <- m[foldid!=k]
#   mtest <- m[foldid==k]
#   ntrain <- length(ytrain)
#   ntest <- length(ytest)
#
#   # grVBEM
#   fit.ridge <- glmnet(xtrain, ytrain, family="binomial", alpha=0, lambda=lambdaridge, standardize=FALSE,
#                       intercept=TRUE)
#   fit.en <- glmnet(xtrain, ytrain, family="binomial", alpha=alphaglmnet, lambda=lambdaglmnet, standardize=FALSE,
#                    intercept=TRUE)
#   fit.grEBEN <- grVBEM2(xtrain, ytrain, mtrain, partitions, lambda1=lambda1gren, lambda2=lambda2gren,
#                         intercept=TRUE, posterior=TRUE, eps=0.001, maxiter=500, trace=TRUE)
#   fit.grridge <- grridge(t(xtrain), ytrain, list(group=CreatePartition(as.factor(partitions[[1]]))), unpenal=~1,
#                          optl=lambda2ridge, trace=FALSE)
#   # fit.SGL <- SGL(list(x=xtrain, y=ytrain), index=groups, type="logit", standardize=FALSE)
#
#   pred1.ridge[foldid==k] <- predict(fit.ridge, xtest, type="response")
#   pred1.en[foldid==k] <- predict(fit.en, xtest, type="response")
#   pred1.grEBEN[foldid==k] <- as.numeric(exp(cbind(1, xtest) %*% fit.grEBEN$mu)/
#                                            (1 + exp(cbind(1, xtest) %*% fit.grEBEN$mu)))
#   pred1.grridge[foldid==k] <- predict(fit.grridge$predobj$GroupRegul, xtest)
#   # pred1.SGL[foldid==k] <- as.numeric(exp(cbind(1, matrix(xtest, nrow=1)) %*%
#   #                                          rbind(fit.SGL$fit$intercepts, fit.SGL$fit$beta)
#   #                                        [, which.min(fit.SGL$lldiff)])/
#   #                                      (1 + exp(cbind(1, matrix(xtest, nrow=1)) %*%
#   #                                                 rbind(fit.SGL$fit$intercepts, fit.SGL$fit$beta)
#   #                                               [, which.min(fit.SGL$lldiff)])))
#
# }
#
# auc1.ridge <- pROC::roc(y, pred1.ridge)$auc
# auc1.en <- pROC::roc(y, pred1.en)$auc
# auc1.grEBEN <- pROC::roc(y, pred1.grEBEN)$auc
# auc1.grridge <- pROC::roc(y, pred1.grridge)$auc
# # auc1.SGL <- pROC::roc(y, pred1.SGL)$auc
#
# pred1 <- cbind(ridge=pred1.ridge, en=pred1.en, grEBEN=pred1.grEBEN, grridge=pred1.grridge)# , SGL=pred1.SGL)
# auc1 <- c(ridge=auc1.ridge, en=auc1.en, grEBEN=auc1.grEBEN, grridge=auc1.grridge)# , SGL=auc1.SGL)
# anno1 <- list(pred=pred1, auc=auc1)
#
# save(anno1, file=paste(path.res, "grEBEN_verlaat_anno_res1.Rdata", sep=""))

set.seed(2017)
fit1.grEBEN <- grVBEM2(x, y, m, partitions, lambda1=NULL, lambda2=NULL, intercept=TRUE,
                       posterior=FALSE, eps=0.001, maxiter=300, trace=TRUE)
fit1.grridge <- grridge(t(x), y, list(annotation=CreatePartition(as.factor(partitions$pvalues))),
                        unpenal=~1)


### using pvalues as grouping information
partitions <- list(pvalues=pvalues[order(pvalues)])
x <- apply(t(datcenVerlaat), 2, function(x) {(x - mean(x))/sd(x)})[, order(pvalues)]
y <- respVerlaat
n <- nrow(x)
p <- ncol(x)
m <- rep(1, n)

# # cross validation of AUC
# set.seed(5001)
# nfolds <- 10
# rest <- n %% nfolds
# foldsize <- c(rep(n %/% nfolds + as.numeric(rest!=0), times=rest),
#               rep(n %/% nfolds, times=nfolds - rest))
# foldid <- sample(rep(1:nfolds, times=foldsize))
# 
# cv.ridge <- cv.glmnet(x, y, family="binomial", alpha=0, standardize=FALSE, intercept=TRUE)
# cv.en <- cv.pen(x, y, intercept=TRUE)
# # cv.SGL <- cvSGL(list(x=x, y=y), index=groups, type="logit", standardize=FALSE)
# 
# lambda2ridge <- 0.5*n*cv.ridge$lambda.min
# lambdaridge <- cv.ridge$lambda.min
# alphaglmnet <- cv.en$alpha[which.min(cv.en$cvll)]
# lambdaglmnet <- cv.en$lambda[which.min(cv.en$cvll)]
# lambda1gren <- 2*n*lambdaglmnet*alphaglmnet
# lambda2gren <- n*lambdaglmnet*(1 - alphaglmnet)
# # lambdaSGL <- cv.SGL$lambdas[which.min(cv.SGL$lldiff)]
# 
# pred1.ridge <- pred1.en <- pred1.grEBEN <- pred1.grridge <- numeric(n) # pred1.SGL <- numeric(n)
# for(k in 1:nfolds) {
#   print(paste("Fold", k, sep=" "))
#   xtrain <- x[foldid!=k, ]
#   xtest <- x[foldid==k, ]
#   ytrain <- y[foldid!=k]
#   ytest <- y[foldid==k]
#   mtrain <- m[foldid!=k]
#   mtest <- m[foldid==k]
#   ntrain <- length(ytrain)
#   ntest <- length(ytest)
# 
#   # grVBEM
#   fit.ridge <- glmnet(xtrain, ytrain, family="binomial", alpha=0, lambda=lambdaridge, standardize=FALSE,
#                       intercept=TRUE)
#   fit.en <- glmnet(xtrain, ytrain, family="binomial", alpha=alphaglmnet, lambda=lambdaglmnet, standardize=FALSE,
#                    intercept=TRUE)
#   fit.grEBEN <- grVBEM2(xtrain, ytrain, mtrain, partitions, lambda1=lambda1gren, lambda2=lambda2gren,
#                         intercept=TRUE, posterior=TRUE, eps=0.001, maxiter=500, trace=TRUE)
#   fit.grridge <- grridge(t(xtrain), ytrain, list(group=CreatePartition(as.factor(partitions[[1]]))), unpenal=~1,
#                          optl=lambda2ridge, trace=FALSE)
#   # fit.SGL <- SGL(list(x=xtrain, y=ytrain), index=groups, type="logit", standardize=FALSE)
# 
#   pred1.ridge[foldid==k] <- predict(fit.ridge, xtest, type="response")
#   pred1.en[foldid==k] <- predict(fit.en, xtest, type="response")
#   pred1.grEBEN[foldid==k] <- as.numeric(exp(cbind(1, xtest) %*% fit.grEBEN$mu)/
#                                            (1 + exp(cbind(1, xtest) %*% fit.grEBEN$mu)))
#   pred1.grridge[foldid==k] <- predict(fit.grridge$predobj$GroupRegul, xtest)
#   # pred1.SGL[foldid==k] <- as.numeric(exp(cbind(1, matrix(xtest, nrow=1)) %*%
#   #                                          rbind(fit.SGL$fit$intercepts, fit.SGL$fit$beta)
#   #                                        [, which.min(fit.SGL$lldiff)])/
#   #                                      (1 + exp(cbind(1, matrix(xtest, nrow=1)) %*%
#   #                                                 rbind(fit.SGL$fit$intercepts, fit.SGL$fit$beta)
#   #                                               [, which.min(fit.SGL$lldiff)])))
# 
# }
# 
# auc1.ridge <- pROC::roc(y, pred1.ridge)$auc
# auc1.en <- pROC::roc(y, pred1.en)$auc
# auc1.grEBEN <- pROC::roc(y, pred1.grEBEN)$auc
# auc1.grridge <- pROC::roc(y, pred1.grridge)$auc
# # auc1.SGL <- pROC::roc(y, pred1.SGL)$auc
# 
# pred1 <- cbind(ridge=pred1.ridge, en=pred1.en, grEBEN=pred1.grEBEN, grridge=pred1.grridge)# , SGL=pred1.SGL)
# auc1 <- c(ridge=auc1.ridge, en=auc1.en, grEBEN=auc1.grEBEN, grridge=auc1.grridge)# , SGL=auc1.SGL)
# pvalues1 <- list(pred=pred1, auc=auc1)
# 
# save(pvalues1, file=paste(path.res, "grEBEN_verlaat_pvalues_res1.Rdata", sep=""))

set.seed(2018)
fit2.grEBEN <- grVBEM2(x, y, m, partitions, lambda1=NULL, lambda2=NULL, intercept=TRUE,
                       posterior=FALSE, eps=0.001, maxiter=300, trace=TRUE)
fit2.grridge <- grridge(t(x), y, list(pvalues=CreatePartition(as.factor(partitions$pvalues))),
                        unpenal=~1)


# ### saving fitted models
# save(fit1.grEBEN, fit1.grridge, fit2.grEBEN, fit2.grridge, file=paste(path.res, "grEBEN_verlaat_anno_fitted.Rdata", sep=""))

### plots
load(paste(path.res, "grEBEN_verlaat_anno_fitted.Rdata", sep=""))

# annotation
load(paste(path.res, "grEBEN_verlaat_anno_res1.Rdata", sep=""))
xlabels1 <- c("CpG \nisland", "North \nShore", "South \nShore", "North \nShelf", "South \nShelf", "Distant")
png(paste(path.graph, "grEBEN_verlaat_bar1.png", sep=""), units="in", width=7, height=5, res=200)
barplot(rbind(fit1.grridge$lambdamults[[1]],
              fit1.grEBEN$lambdag[[1]][, fit1.grEBEN$nouteriter + 1]), beside=TRUE,
        names.arg=xlabels1, ylab=expression(paste(lambda[g], "'")),
        args.legend=list(x="topleft", fill=c(gray.colors(2), 0, 0), lty=c(NA, NA, 2, 2),
                         border=c(rep(1, 2), 0, 0), merge=TRUE, seg.len=1),
        legend.text=paste(c("GRridge", "grEBEN", "ridge", "enet"), ", AUC=",
                          round(anno1$auc[c(4, 3, 1, 2)], 2), sep=""))
abline(h=1, lty=2)
dev.off()

cols <- gray.colors(4)
png(paste(path.graph, "grEBEN_verlaat_roc1.png", sep=""), units="in", width=5, height=5, res=200)
plot(roc(y, anno1$pred[, 4]), asp=1, col=cols[1], lty=1)
plot(roc(y, anno1$pred[, 3]), add=TRUE, col=cols[2], lty=2)
plot(roc(y, anno1$pred[, 1]), add=TRUE, col=cols[3], lty=1)
plot(roc(y, anno1$pred[, 2]), add=TRUE, col=cols[4], lty=2)
legend("bottomright", legend=paste(c("GRridge", "grEBEN", "ridge", "enet"), ", AUC=",
                                   round(anno1$auc[c(4, 3, 1, 2)], 2), sep=""),
       lty=c(1, 2, 1, 2), col=cols)
dev.off()

# pvalues
load(paste(path.res, "grEBEN_verlaat_pvalues_res1.Rdata", sep=""))
png(paste(path.graph, "grEBEN_verlaat_bar2.png", sep=""), units="in", width=6, height=5, res=200)
barplot(rbind(fit2.grridge$lambdamults$pvalues,
              fit2.grEBEN$lambdag$pvalues[, fit2.grEBEN$nouteriter + 1]), beside=TRUE,
        xlab="", axisnames=FALSE,
        args.legend=list(x="top", fill=c(gray.colors(2), 0, 0), lty=c(NA, NA, 2, 2),
                         border=c(rep(1, 2), 0, 0), merge=TRUE, seg.len=1),
        ylab=expression(paste(lambda[g], "'")),
        legend.text=paste(c("GRridge", "grEBEN", "ridge", "enet"), ", AUC=",
                          round(pvalues1$auc[c(4, 3, 1, 2)], 2), sep=""))
title(xlab="P-values in increasing order", line=1)
abline(h=1, lty=2)
dev.off()

cols <- gray.colors(4)
png(paste(path.graph, "grEBEN_verlaat_roc2.png", sep=""), units="in", width=5, height=5, res=200)
plot(roc(y, pvalues1$pred[, 4]), asp=1, col=cols[1], lty=1)
plot(roc(y, pvalues1$pred[, 3]), add=TRUE, col=cols[2], lty=2)
plot(roc(y, pvalues1$pred[, 1]), add=TRUE, col=cols[3], lty=1)
plot(roc(y, pvalues1$pred[, 2]), add=TRUE, col=cols[4], lty=2)
legend("bottomright", legend=paste(c("GRridge", "grEBEN", "ridge", "enet"), ", AUC=",
                                   round(pvalues1$auc[c(4, 3, 1, 2)], 2), sep=""),
       lty=c(1, 2, 1, 2), col=cols)
dev.off()


