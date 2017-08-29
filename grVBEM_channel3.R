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

### loading data
load(paste(path.data, "mirsData.RData", sep=""))

# ### creating partitions
parCons <- CreatePartition(mirsData$conservation) # using conservation status as grouping
parAbund <- CreatePartition(rowSums(mirsData$countData), mingr=25, ngroup=10, decreasing=TRUE)




### using conservation status as grouping information
groups <- rep(1:length(parCons), unlist(lapply(parCons, length)))
x <- apply(t(as.matrix(mirsData$transformedData)), 2,
           function(x) {(x - mean(x))/sd(x)})[, unlist(parCons)]
y <- as.numeric(mirsData$response) - 1
partitions1 <- list(conservation=groups)
n <- nrow(x)
p <- ncol(x)
m <- rep(1, n)

# cross validation of AUC
set.seed(1001)
nfolds <- 10
rest <- n %% nfolds
foldsize <- c(rep(n %/% nfolds + as.numeric(rest!=0), times=rest),
              rep(n %/% nfolds, times=nfolds - rest))
foldid <- sample(rep(1:nfolds, times=foldsize))

cv.ridge <- cv.glmnet(x, y, family="binomial", alpha=0, standardize=FALSE, intercept=TRUE)
cv.en <- cv.pen(x, y, intercept=TRUE)

lambdaridge <- cv.ridge$lambda.min
lambda2ridge <- cv.ridge$lambda.min*0.5*n
lambda1bayes <- cv.en$lambda1bayes[which.min(cv.en$cvll)]
lambda2bayes <- cv.en$lambda2bayes[which.min(cv.en$cvll)]
lambdaenet <- cv.en$lambda[which.min(cv.en$cvll)]
alphaenet <- cv.en$alpha[which.min(cv.en$cvll)]

# pred1.ridge <- pred1.en <- pred1.grEBEN <- pred1.grridge <- numeric(n)
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
#   fit.grridge <- grridge(t(xtrain), ytrain, list(group=CreatePartition(as.factor(groups))), unpenal=~1,
#                          optl=lambda2ridge, trace=FALSE)
#   fit.en <- glmnet(xtrain, ytrain, family="binomial", alpha=alphaenet, lambda=lambdaenet, standardize=FALSE,
#                    intercept=TRUE)
#   fit.grEBEN <- grVBEM2(xtrain, ytrain, m=mtrain, partitions=partitions1, lambda1=lambda1bayes, 
#                         lambda2=lambda2bayes, intercept=TRUE, monotone=FALSE, iso=NULL,
#                         posterior=TRUE, eps=0.001, maxiter=500, trace=TRUE, alphastart=NULL)
#   
#   pred1.ridge[foldid==k] <- as.numeric(predict(fit.ridge, xtest, type="response"))
#   pred1.en[foldid==k] <- as.numeric(predict(fit.en, xtest, type="response"))
#   pred1.grEBEN[foldid==k] <- as.numeric(exp(cbind(1, xtest) %*% fit.grEBEN$mu)/
#                                            (1 + exp(cbind(1, xtest) %*% fit.grEBEN$mu)))
#   pred1.grridge[foldid==k] <- as.numeric(predict(fit.grridge$predobj$GroupRegul, xtest))
# 
# }
# 
# auc1.ridge <- pROC::roc(y, pred1.ridge)$auc
# auc1.en <- pROC::roc(y, pred1.en)$auc
# auc1.grEBEN <- pROC::roc(y, pred1.grEBEN)$auc
# auc1.grridge <- pROC::roc(y, pred1.grridge)$auc
# 
# pred1 <- cbind(grridge=pred1.grridge, grEBEN=pred1.grEBEN, ridge=pred1.ridge, en=pred1.en)
# auc1 <- c(GRridge=auc1.grridge, grEBEN=auc1.grEBEN, ridge=auc1.ridge, enet=auc1.en)
# cons1 <- list(pred=pred1, auc=auc1)
# 
# save(cons1, file=paste(path.res, "grVBEM_channel2_cons_res1.2.Rdata", sep=""))

# co-data
set.seed(1001)
# conservation status as co-data
fit1.grEBEN <- grVBEM2(x, y, m, partitions1, lambda1=lambda1bayes, lambda2=lambda2bayes, 
                       intercept=TRUE, monotone=FALSE, iso=NULL, posterior=TRUE, eps=0.001,
                       maxiter=500, trace=TRUE, alphastart=NULL)
fit1.grridge <- grridge(t(x), y, list(conservation=CreatePartition(as.factor(groups))), 
                        unpenal=~1)


### abundance on top of conservation status
groups1 <- rep(1:length(parCons), times=unlist(lapply(parCons, length)))[order(unlist(parCons))]
groups2 <- rep(1:length(parAbund), times=unlist(lapply(parAbund, length)))[order(unlist(parAbund))]
x <- apply(t(as.matrix(mirsData$transformedData)), 2, function(x) {
  (x - mean(x))/sd(x)})[, order(groups1, groups2)]
partitions1 <- list(conservation=groups1[order(groups1, groups2)],
                    abundance=groups2[order(groups1, groups2)])
y <- as.numeric(mirsData$response) - 1
n <- nrow(x)
p <- ncol(x)
m <- rep(1, n)

# cross validation of AUC
set.seed(1002)
nfolds <- 10
rest <- n %% nfolds
foldsize <- c(rep(n %/% nfolds + as.numeric(rest!=0), times=rest),
              rep(n %/% nfolds, times=nfolds - rest))
foldid <- sample(rep(1:nfolds, times=foldsize))

cv.ridge <- cv.glmnet(x, y, family="binomial", alpha=0, standardize=FALSE, intercept=TRUE)
cv.en <- cv.pen(x, y, intercept=TRUE)

lambdaridge <- cv.ridge$lambda.min
lambda2ridge <- cv.ridge$lambda.min*0.5*n
lambda1bayes <- cv.en$lambda1bayes[which.min(cv.en$cvll)]
lambda2bayes <- cv.en$lambda2bayes[which.min(cv.en$cvll)]
lambdaenet <- cv.en$lambda[which.min(cv.en$cvll)]
alphaenet <- cv.en$alpha[which.min(cv.en$cvll)]

# pred1.ridge <- pred1.en <- pred1.grEBEN <- pred1.grridge <- numeric(n)
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
#   fit.grridge <- grridge(t(xtrain), ytrain, list(conservation=CreatePartition(as.factor(partitions1$conservation)),
#                                                  abundance=CreatePartition(as.factor(partitions1$abundance))),
#                            unpenal=~1, optl=lambda2ridge, trace=FALSE, monotone=c(FALSE, TRUE))  
#   fit.en <- glmnet(xtrain, ytrain, family="binomial", alpha=alphaenet, lambda=lambdaenet, standardize=FALSE,
#                      intercept=TRUE)
#   fit.grEBEN <- grVBEM2(xtrain, ytrain, m=mtrain, partitions=partitions1, lambda1=lambda1bayes,
#                         lambda2=lambda2bayes, intercept=TRUE, monotone=FALSE, 
#                         iso=list(conservation=FALSE, abundance=TRUE), posterior=TRUE, eps=0.001, 
#                         maxiter=500, trace=TRUE, alphastart=NULL)
#   
#   pred1.ridge[foldid==k] <- as.numeric(predict(fit.ridge, xtest, type="response"))
#   pred1.en[foldid==k] <- as.numeric(predict(fit.en, xtest, type="response"))
#   pred1.grEBEN[foldid==k] <- as.numeric(exp(cbind(1, xtest) %*% fit.grEBEN$mu)/
#                                            (1 + exp(cbind(1, xtest) %*% fit.grEBEN$mu)))
#   pred1.grridge[foldid==k] <- as.numeric(predict(fit.grridge$predobj$GroupRegul, xtest))
# 
# }
# 
# auc1.ridge <- pROC::roc(y, pred1.ridge)$auc
# auc1.en <- pROC::roc(y, pred1.en)$auc
# auc1.grEBEN <- pROC::roc(y, pred1.grEBEN)$auc
# auc1.grridge <- pROC::roc(y, pred1.grridge)$auc
# 
# pred1 <- cbind(grridge=pred1.grridge, grEBEN=pred1.grEBEN, ridge=pred1.ridge, enet=pred1.en)
# auc1 <- c(GRridge=auc1.grridge, grEBEN=auc1.grEBEN, ridge=auc1.ridge, enet=auc1.en)
# cons2 <- list(pred=pred1, auc=auc1)
# 
# save(cons2, file=paste(path.res, "grVBEM_channel2_cons_res2.2.Rdata", sep=""))

# fitting the models to the full data
set.seed(1002)
fit2.grEBEN <- grVBEM2(x, y, m, partitions1, lambda1=lambda1bayes,lambda2=lambda2bayes, 
                       intercept=TRUE, monotone=FALSE, iso=list(conservation=FALSE, abundance=TRUE), 
                       posterior=TRUE, eps=0.001, maxiter=300, trace=TRUE, alphastart=NULL)
fit2.grridge <- grridge(t(x), y, list(conservation=CreatePartition(as.factor(partitions1$conservation)),
                                        abundance=CreatePartition(as.factor(partitions1$abundance))), 
                          monotone=c(FALSE, TRUE), unpenal=~1, optl=lambda2ridge)




### using only abundance
groups <- rep(1:length(parAbund), unlist(lapply(parAbund, length)))
x <- apply(t(as.matrix(mirsData$transformedData)), 2,
           function(x) {(x - mean(x))/sd(x)})[, unlist(parAbund)]
y <- as.numeric(mirsData$response) - 1
n <- nrow(x)
p <- ncol(x)
m <- rep(1, n)
partitions1 <- list(abundance=groups)

# cross validation of AUC
set.seed(1003)
nfolds <- 10
rest <- n %% nfolds
foldsize <- c(rep(n %/% nfolds + as.numeric(rest!=0), times=rest),
              rep(n %/% nfolds, times=nfolds - rest))
foldid <- sample(rep(1:nfolds, times=foldsize))

cv.ridge <- cv.glmnet(x, y, family="binomial", alpha=0, standardize=FALSE, intercept=TRUE)
cv.en <- cv.pen(x, y, intercept=TRUE)

lambdaridge <- cv.ridge$lambda.min
lambda2ridge <- cv.ridge$lambda.min*0.5*n
lambda1bayes <- cv.en$lambda1bayes[which.min(cv.en$cvll)]
lambda2bayes <- cv.en$lambda2bayes[which.min(cv.en$cvll)]
lambdaenet <- cv.en$lambda[which.min(cv.en$cvll)]
alphaenet <- cv.en$alpha[which.min(cv.en$cvll)]

# pred1.ridge <- pred1.en <- pred1.grEBEN <- pred1.grridge <- numeric(n)
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
#   # fitting the models to k-1 folds
#   fit.ridge <- glmnet(xtrain, ytrain, family="binomial", alpha=0, lambda=lambdaridge, 
#                       standardize=FALSE, intercept=TRUE)
#   fit.en <- glmnet(xtrain, ytrain, family="binomial", alpha=alphaenet, lambda=lambdaenet, 
#                    standardize=FALSE, intercept=TRUE)
#   fit.grEBEN <- grVBEM2(xtrain, ytrain, m=mtrain, partitions=partitions1, lambda1=lambda1bayes,
#                         lambda2=lambda2bayes, intercept=TRUE, monotone=FALSE, 
#                         iso=list(abundance=TRUE), posterior=TRUE, eps=0.001, maxiter=500, 
#                         trace=TRUE, alphastart=NULL)
#   fit.grridge <- grridge(t(xtrain), ytrain, list(abundance=CreatePartition(as.factor(partitions1$abundance))),
#                          unpenal=~1, optl=lambda2ridge, monotone=TRUE, trace=FALSE)
#   
#   # make predictions on test data
#   pred1.ridge[foldid==k] <- as.numeric(predict(fit.ridge, xtest, type="response"))
#   pred1.en[foldid==k] <- as.numeric(predict(fit.en, xtest, type="response"))
#   pred1.grEBEN[foldid==k] <- as.numeric(exp(cbind(1, xtest) %*% fit.grEBEN$mu)/
#                                           (1 + exp(cbind(1, xtest) %*% fit.grEBEN$mu)))
#   pred1.grridge[foldid==k] <- as.numeric(predict(fit.grridge$predobj$GroupRegul, xtest))
# 
# }
# 
# auc1.ridge <- pROC::roc(y, pred1.ridge)$auc
# auc1.en <- pROC::roc(y, pred1.en)$auc
# auc1.grEBEN <- pROC::roc(y, pred1.grEBEN)$auc
# auc1.grridge <- pROC::roc(y, pred1.grridge)$auc
# 
# pred1 <- cbind(grridge=pred1.grridge, grEBEN=pred1.grEBEN, ridge=pred1.ridge, enet=pred1.en)
# auc1 <- c(GRridge=auc1.grridge, grEBEN=auc1.grEBEN, ridge=auc1.ridge, enet=auc1.en)
# cons3 <- list(pred=pred1, auc=auc1)
# 
# save(cons3, file=paste(path.res, "grVBEM_channel2_cons_res3.2.Rdata", sep=""))

# fitting the model
set.seed(1003)
fit3.grEBEN <- grVBEM2(x, y, m, partitions1, lambda1=lambda1bayes, 
                       lambda2=lambda2bayes, intercept=TRUE, monotone=FALSE, 
                       iso=list(abundance=TRUE), posterior=FALSE, eps=0.001, maxiter=300, 
                       trace=TRUE, alphastart=NULL)
fit3.grridge <- grridge(t(x), y, list(abundance=CreatePartition(as.factor(groups))), unpenal=~1,
                        monotone=TRUE)

# ### saving fitted models
# save(fit1.grEBEN, fit1.grridge, fit2.grEBEN, fit2.grridge, fit3.grEBEN, fit3.grridge,
#      file=paste(path.res, "grVBEM_channel2_cons_fitted2.Rdata", sep=""))



### plots
load(paste(path.res, "grVBEM_channel2_cons_fitted2.Rdata", sep=""))

# conservation status of genes
load(paste(path.res, "grVBEM_channel2_cons_res1.2.Rdata", sep=""))
xlabels1 <- c("not \n conserved", "conserved \n in mammals", "broadly \n conserved")
png(paste(path.graph, "grEBEN_channel_bar1.2.png", sep=""), units="in", width=4.5, height=4, res=200)
barplot(rbind(fit1.grridge$lambdamults[[1]],
              fit1.grEBEN$lambdag[[1]][, fit1.grEBEN$nouteriter + 1]), beside=TRUE,
        names.arg=xlabels1, ylab=expression(paste(lambda[g], "'")),
        args.legend=list(x="topright", fill=c(gray.colors(2), 0, 0), lty=c(NA, NA, 2, 2),
                         border=c(rep(1, 2), 0, 0), merge=TRUE, seg.len=1),
        legend.text=paste(c("GRridge", "grEBEN", "ridge", "enet"), ", AUC=",
                          round(cons1$auc, 2), sep=""))
abline(h=1, lty=2)
dev.off()

cols <- gray.colors(4)
png(paste(path.graph, "grEBEN_channel_roc1.2.png", sep=""), units="in", width=5, height=5, res=200)
plot(roc(y, cons1$pred[, 1]), asp=1, col=cols[1], lty=1)
plot(roc(y, cons1$pred[, 2]), add=TRUE, col=cols[2], lty=2)
plot(roc(y, cons1$pred[, 3]), add=TRUE, col=cols[3], lty=1)
plot(roc(y, cons1$pred[, 4]), add=TRUE, col=cols[4], lty=2)
legend("bottomright", legend=paste(c("GRridge", "grEBEN", "ridge", "enet"), ", AUC=",
                                   round(cons1$auc, 2), sep=""),
       lty=c(1, 2, 1, 2), col=cols)
dev.off()


# conservation and abundance
xlabels1 <- c("not \n conserved", "conserved \n in mammals", "broadly \n conserved")
load(paste(path.res, "grVBEM_channel2_cons_res2.2.Rdata", sep=""))
png(paste(path.graph, "grEBEN_channel_bar2.2.png", sep=""), units="in", width=9, height=5.5, res=200)
par(mfrow=c(1, 2))
barplot(rbind(fit2.grridge$lambdamults$abundance,
              fit2.grEBEN$lambdag$abundance[, fit2.grEBEN$nouteriter + 1]), beside=TRUE,
        xlab="", axisnames=FALSE,
        ylab=expression(paste(lambda[g], "'")))
# args.legend=list(x="topleft", fill=c(gray.colors(2), 0, 0), lty=c(NA, NA, 2, 2),
#                  border=c(rep(1, 2), 0, 0), merge=TRUE, seg.len=1),
# legend.text=paste(c("GRridge", "grEBEN", "ridge", "enet"), ", AUC=",
#                   round(cons2$auc[c(4, 3, 1, 2)], 2), sep=""))
title(xlab="Abundances in decreasing order", line=1)
abline(h=1, lty=2)

barplot(rbind(fit2.grridge$lambdamults$conservation,
              fit2.grEBEN$lambdag$conservation[, fit2.grEBEN$nouteriter + 1]), beside=TRUE,
        names.arg=xlabels1, ylab=expression(paste(lambda[g], "'")),
        args.legend=list(x="topright", fill=c(gray.colors(2), 0, 0), lty=c(NA, NA, 2, 2),
                         border=c(rep(1, 2), 0, 0), merge=TRUE, seg.len=1),
        legend.text=paste(c("GRridge", "grEBEN", "ridge", "enet"), ", AUC=",
                          round(cons2$auc, 2), sep=""))
abline(h=1, lty=2)
dev.off()

cols <- gray.colors(4)
png(paste(path.graph, "grEBEN_channel_roc2.2.png", sep=""), units="in", width=5, height=5, res=200)
plot(roc(y, cons2$pred[, 1]), asp=1, col=cols[1], lty=1)
plot(roc(y, cons2$pred[, 2]), add=TRUE, col=cols[2], lty=2)
plot(roc(y, cons2$pred[, 3]), add=TRUE, col=cols[3], lty=1)
plot(roc(y, cons2$pred[, 4]), add=TRUE, col=cols[4], lty=2)
legend("bottomright", legend=paste(c("GRridge", "grEBEN", "ridge", "enet"), ", AUC=",
                                   round(cons2$auc, 2), sep=""),
       lty=c(1, 2, 1, 2), col=cols)
dev.off()


# only abundance
load(paste(path.res, "grVBEM_channel2_cons_res3.2.Rdata", sep=""))
png(paste(path.graph, "grEBEN_channel_bar3.2.png", sep=""), units="in", width=4.5, height=4, res=200)
barplot(rbind(fit3.grridge$lambdamults$abundance,
              fit3.grEBEN$lambdag$abundance[, fit3.grEBEN$nouteriter + 1]), beside=TRUE,
        xlab="", axisnames=FALSE,
        args.legend=list(x="topleft", fill=c(gray.colors(2), 0, 0), lty=c(NA, NA, 2, 2),
                         border=c(rep(1, 2), 0, 0), merge=TRUE, seg.len=1),
        ylab=expression(paste(lambda[g], "'")),
        legend.text=paste(c("GRridge", "grEBEN", "ridge", "enet"), ", AUC=",
                          round(cons3$auc, 2), sep=""))
title(xlab="Abundances in decreasing order", line=1)
abline(h=1, lty=2)
dev.off()

cols <- gray.colors(4)
png(paste(path.graph, "grEBEN_channel_roc3.2.png", sep=""), units="in", width=5, height=5, res=200)
plot(roc(y, cons3$pred[, 1]), asp=1, col=cols[1], lty=1)
plot(roc(y, cons3$pred[, 2]), add=TRUE, col=cols[2], lty=2)
plot(roc(y, cons3$pred[, 3]), add=TRUE, col=cols[3], lty=1)
plot(roc(y, cons3$pred[, 4]), add=TRUE, col=cols[4], lty=2)
legend("bottomright", legend=paste(c("GRridge", "grEBEN", "ridge", "enet"), ", AUC=",
                                   round(cons3$auc, 2), sep=""),
       lty=c(1, 2, 1, 2), col=cols)
dev.off()
