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

# # cross validation of AUC
# set.seed(1001)
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
#   fit.grEBEN <- grVBEM(xtrain, ytrain, m=mtrain, partitions=groups, lambda1=lambda1gren, lambda2=lambda2gren,
#                        intercept=TRUE, eps=0.001, maxiter=500, trace=TRUE, QNacc=FALSE)
#   fit.grridge <- grridge(t(xtrain), ytrain, list(group=CreatePartition(as.factor(groups))), unpenal=~1,
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
# cons1 <- list(pred=pred1, auc=auc1)
#
# save(cons1, file=paste(path.res, "grVBEM_channel2_cons_res1.Rdata", sep=""))

# co-data
set.seed(1001)
# conservation status as co-data
# fit1.grEBEN <- grVBEM(x, y, m, partitions1, lambda1=NULL, lambda2=NULL, intercept=TRUE, eps=0.001,
#                       maxiter=300, trace=TRUE)
fit1.grEBEN <- grVBEM(x, y, m, partitions1, lambda1=NULL, lambda2=NULL, intercept=TRUE, eps=0.001,
                      maxiter=300, trace=TRUE)
fit1.grridge <- grridge(t(x), y, list(conservation=CreatePartition(as.factor(groups))), unpenal=~1)



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

# # cross validation of AUC
# set.seed(1002)
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
#   fit.grEBEN <- grVBEM(xtrain, ytrain, m=mtrain, partitions=partitions1, lambda1=lambda1gren, lambda2=lambda2gren,
#                        intercept=TRUE, eps=0.001, maxiter=500, trace=TRUE)
#   fit.grridge <- grridge(t(xtrain), ytrain, list(conservation=CreatePartition(as.factor(partitions1$conservation)),
#                                                  abundance=CreatePartition(as.factor(partitions1$abundance))),
#                          unpenal=~1, optl=lambda2ridge, trace=FALSE)
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
# cons2 <- list(pred=pred1, auc=auc1)
#
# save(cons2, file=paste(path.res, "grVBEM_channel2_cons_res2.Rdata", sep=""))

set.seed(1002)
fit2.grEBEN <- grVBEM2(x, y, m, partitions1, lambda1=NULL, lambda2=NULL, intercept=TRUE,
                       monotone=FALSE, posterior=FALSE, eps=0.001, maxiter=300, trace=TRUE)
fit2.grridge <- grridge(t(x), y, list(conservation=CreatePartition(as.factor(partitions1$conservation)),
                                      abundance=CreatePartition(as.factor(partitions1$abundance))), unpenal=~1)

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
# set.seed(1003)
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
#   fit.grEBEN <- grVBEM(xtrain, ytrain, m=mtrain, partitions=partitions1, lambda1=lambda1gren, lambda2=lambda2gren,
#                        intercept=TRUE, eps=0.001, maxiter=500, trace=TRUE, QNacc=FALSE)
#   fit.grridge <- grridge(t(xtrain), ytrain, list(abundance=CreatePartition(as.factor(partitions1$abundance))),
#                          unpenal=~1, optl=lambda2ridge, trace=FALSE)
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
# cons3 <- list(pred=pred1, auc=auc1)
#
# save(cons3, file=paste(path.res, "grVBEM_channel2_cons_res3.Rdata", sep=""))

# fitting the model
set.seed(1003)
fit3.grEBEN <- grVBEM2(x, y, m, partitions1, lambda1=NULL, lambda2=NULL, intercept=TRUE,
                       monotone=FALSE, posterior=FALSE, eps=0.001, maxiter=300, trace=TRUE)
fit3.grridge <- grridge(t(x), y, list(abundance=CreatePartition(as.factor(groups))), unpenal=~1,
                        monotone=FALSE)
plot(fit3.grEBEN$lambdag$abundance[, ncol(fit3.grEBEN$lambdag$abundance)], fit3.grridge$lambdamults$abundance)

# ### random groups
# x <- apply(t(as.matrix(mirsData$transformedData)), 2,
#            function(x) {(x - mean(x))/sd(x)})[, unlist(parAbund)]
# y <- as.numeric(mirsData$response) - 1
# n <- nrow(x)
# p <- ncol(x)
# m <- rep(1, n)
# G <- 10
# set.seed(1004)
# groups <- sample(rep(1:G, times=c(rep(p %/% G + 1, p %% G), rep(p %/% G, G - p %% G))))
# partitions1 <- list(random=groups)
#
# # cross validation of AUC
# set.seed(1004)
# nfolds <- n
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
# mat.mult <- matrix(nrow=nfolds, ncol=G*2)
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
#   fit.grEBEN <- grVBEM(xtrain, ytrain, m=mtrain, partitions=partitions1, lambda1=lambda1gren, lambda2=lambda2gren,
#                        intercept=TRUE, eps=0.001, maxiter=500, trace=TRUE)
#   fit.grridge <- grridge(t(xtrain), ytrain, list(random=CreatePartition(as.factor(partitions1$random))),
#                          unpenal=~1, optl=lambda2ridge, trace=FALSE)
#   # fit.SGL <- SGL(list(x=xtrain, y=ytrain), index=groups, type="logit", standardize=FALSE)
#
#   pred1.ridge[foldid==k] <- predict(fit.ridge, matrix(xtest, nrow=1), type="response")
#   pred1.en[foldid==k] <- predict(fit.en, matrix(xtest, nrow=1), type="response")
#   pred1.grEBEN[foldid==k] <- as.numeric(exp(matrix(c(1, xtest), nrow=1) %*% fit.grEBEN$mu)/
#                                            (1 + exp(matrix(c(1, xtest), nrow=1) %*% fit.grEBEN$mu)))
#   pred1.grridge[foldid==k] <- predict(fit.grridge$predobj$GroupRegul, matrix(xtest, nrow=1))
#   # pred1.SGL[foldid==k] <- as.numeric(exp(cbind(1, matrix(xtest, nrow=1)) %*%
#   #                                          rbind(fit.SGL$fit$intercepts, fit.SGL$fit$beta)
#   #                                        [, which.min(fit.SGL$lldiff)])/
#   #                                      (1 + exp(cbind(1, matrix(xtest, nrow=1)) %*%
#   #                                                 rbind(fit.SGL$fit$intercepts, fit.SGL$fit$beta)
#   #                                               [, which.min(fit.SGL$lldiff)])))
#
#   mat.mult[k, ] <- c(fit.grridge$lambdamults$random,
#                      fit.grEBEN$lambdag$random[, fit.grEBEN$nouteriter + 1])
#
# }
#
# auc1.ridge <- pROC::roc(y, pred1.ridge)$auc
# auc1.en <- pROC::roc(y, pred1.en)$auc
# auc1.grEBEN <- pROC::roc(y, pred1.grEBEN)$auc
# auc1.grridge <- pROC::roc(y, pred1.grridge)$auc
# # auc1.SGL <- pROC::roc(y, pred1.SGL)$auc
#
# colnames(mat.mult) <- paste(rep(c("grridge", "grEBEN"), each=G), 1:G, sep="")
# pred1 <- cbind(ridge=pred1.ridge, en=pred1.en, grEBEN=pred1.grEBEN, grridge=pred1.grridge)# , SGL=pred1.SGL)
# auc1 <- c(ridge=auc1.ridge, en=auc1.en, grEBEN=auc1.grEBEN, grridge=auc1.grridge)# , SGL=auc1.SGL)
# cons4.2 <- list(pred=pred1, auc=auc1, lambdag=mat.mult)
#
# save(cons4.2, file=paste(path.res, "grVBEM_channel2_cons_res4.2.Rdata", sep=""))

set.seed(1004)
fit4.grEBEN <- grVBEM(x, y, m, partitions1, lambda1=NULL, lambda2=NULL, intercept=TRUE, eps=0.001,
                      maxiter=300, trace=TRUE)
fit4.grridge <- grridge(t(x), y, list(random=CreatePartition(as.factor(groups))), unpenal=~1)


### random groups 100 times
x <- apply(t(as.matrix(mirsData$transformedData)), 2,
           function(x) {(x - mean(x))/sd(x)})[, unlist(parAbund)]
y <- as.numeric(mirsData$response) - 1
n <- nrow(x)
p <- ncol(x)
m <- rep(1, n)
G <- 10

# # cross validation of AUC
# set.seed(1005)
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
# nreps <- 100
# mat.mult <- matrix(nrow=nreps, ncol=G*2)
# for(k in 1:nreps) {
#   print(paste("Repeat", k, sep=" "))
#   
#   groups <- sample(rep(1:G, times=c(rep(p %/% G + 1, p %% G), rep(p %/% G, G - p %% G))))
#   partitions1 <- list(random=groups)
# 
#   # grVBEM
#   fit.grEBEN <- grVBEM(x, y, m=m, partitions=partitions1, lambda1=lambda1gren, lambda2=lambda2gren,
#                        intercept=TRUE, eps=0.001, maxiter=500, trace=TRUE)
#   fit.grridge <- grridge(t(x), y, list(random=CreatePartition(as.factor(partitions1$random))),
#                          unpenal=~1, optl=lambda2ridge, trace=FALSE)
# 
#   mat.mult[k, ] <- c(fit.grridge$lambdamults$random,
#                      fit.grEBEN$lambdag$random[, fit.grEBEN$nouteriter + 1])
# 
# }
# 
# colnames(mat.mult) <- paste(rep(c("grridge", "grEBEN"), each=G), 1:G, sep="")
# cons5 <- list(lambdag=mat.mult)
# 
# save(cons5, file=paste(path.res, "grVBEM_channel2_cons_res5.Rdata", sep=""))




# ### saving fitted models
# save(fit1.grEBEN, fit1.grridge, fit2.grEBEN, fit2.grridge, fit3.grEBEN, fit3.grridge,
#      fit4.grEBEN, fit4.grridge, file=paste(path.res, "grVBEM_channel2_cons_fitted.Rdata", sep=""))



### plots
load(paste(path.res, "grVBEM_channel2_cons_fitted.Rdata", sep=""))

# conservation status of genes
load(paste(path.res, "grVBEM_channel2_cons_res1.Rdata", sep=""))
xlabels1 <- c("not \n conserved", "conserved \n in mammals", "broadly \n conserved")
png(paste(path.graph, "grEBEN_channel_bar1.png", sep=""), units="in", width=4.5, height=4, res=200)
barplot(rbind(fit1.grridge$lambdamults[[1]],
              fit1.grEBEN$lambdag[[1]][, fit1.grEBEN$nouteriter + 1]), beside=TRUE,
        names.arg=xlabels1, ylab=expression(paste(lambda[g], "'")),
        args.legend=list(x="topright", fill=c(gray.colors(2), 0, 0), lty=c(NA, NA, 2, 2),
                         border=c(rep(1, 2), 0, 0), merge=TRUE, seg.len=1),
        legend.text=paste(c("GRridge", "grEBEN", "ridge", "enet"), ", AUC=",
                          round(cons1$auc[c(4, 3, 1, 2)], 2), sep=""))
abline(h=1, lty=2)
dev.off()

cols <- gray.colors(4)
png(paste(path.graph, "grEBEN_channel_roc1.png", sep=""), units="in", width=5, height=5, res=200)
plot(roc(y, cons1$pred[, 4]), asp=1, col=cols[1], lty=1)
plot(roc(y, cons1$pred[, 3]), add=TRUE, col=cols[2], lty=2)
plot(roc(y, cons1$pred[, 1]), add=TRUE, col=cols[3], lty=1)
plot(roc(y, cons1$pred[, 2]), add=TRUE, col=cols[4], lty=2)
legend("bottomright", legend=paste(c("GRridge", "grEBEN", "ridge", "enet"), ", AUC=",
                                   round(cons1$auc[c(4, 3, 1, 2)], 2), sep=""),
       lty=c(1, 2, 1, 2), col=cols)
dev.off()


# conservation and abundance
xlabels1 <- c("not \n conserved", "conserved \n in mammals", "broadly \n conserved")
load(paste(path.res, "grVBEM_channel2_cons_res2.Rdata", sep=""))
png(paste(path.graph, "grEBEN_channel_bar2.png", sep=""), units="in", width=9, height=5.5, res=200)
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
                          round(cons2$auc[c(4, 3, 1, 2)], 2), sep=""))
abline(h=1, lty=2)
dev.off()

cols <- gray.colors(4)
png(paste(path.graph, "grEBEN_channel_roc2.png", sep=""), units="in", width=5, height=5, res=200)
plot(roc(y, cons2$pred[, 4]), asp=1, col=cols[1], lty=1)
plot(roc(y, cons2$pred[, 3]), add=TRUE, col=cols[2], lty=2)
plot(roc(y, cons2$pred[, 1]), add=TRUE, col=cols[3], lty=1)
plot(roc(y, cons2$pred[, 2]), add=TRUE, col=cols[4], lty=2)
legend("bottomright", legend=paste(c("GRridge", "grEBEN", "ridge", "enet"), ", AUC=",
                                   round(cons2$auc[c(4, 3, 1, 2)], 2), sep=""),
       lty=c(1, 2, 1, 2), col=cols)
dev.off()


# only abundance
load(paste(path.res, "grVBEM_channel2_cons_res3.Rdata", sep=""))
png(paste(path.graph, "grEBEN_channel_bar3.png", sep=""), units="in", width=4.5, height=4, res=200)
barplot(rbind(fit3.grridge$lambdamults$abundance,
              fit3.grEBEN$lambdag$abundance[, fit3.grEBEN$nouteriter + 1]), beside=TRUE,
        xlab="", axisnames=FALSE,
        args.legend=list(x="topleft", fill=c(gray.colors(2), 0, 0), lty=c(NA, NA, 2, 2),
                         border=c(rep(1, 2), 0, 0), merge=TRUE, seg.len=1),
        ylab=expression(paste(lambda[g], "'")),
        legend.text=paste(c("GRridge", "grEBEN", "ridge", "enet"), ", AUC=",
                          round(cons3$auc[c(4, 3, 1, 2)], 2), sep=""))
title(xlab="Abundances in decreasing order", line=1)
abline(h=1, lty=2)
dev.off()

cols <- gray.colors(4)
png(paste(path.graph, "grEBEN_channel_roc3.png", sep=""), units="in", width=5, height=5, res=200)
plot(roc(y, cons3$pred[, 4]), asp=1, col=cols[1], lty=1)
plot(roc(y, cons3$pred[, 3]), add=TRUE, col=cols[2], lty=2)
plot(roc(y, cons3$pred[, 1]), add=TRUE, col=cols[3], lty=1)
plot(roc(y, cons3$pred[, 2]), add=TRUE, col=cols[4], lty=2)
legend("bottomright", legend=paste(c("GRridge", "grEBEN", "ridge", "enet"), ", AUC=",
                                   round(cons3$auc[c(4, 3, 1, 2)], 2), sep=""),
       lty=c(1, 2, 1, 2), col=cols)
dev.off()

# random groups
load(paste(path.res, "grVBEM_channel2_cons_res4.Rdata", sep=""))
png(paste(path.graph, "grEBEN_channel_bar4.png", sep=""), units="in", width=4.5, height=4, res=200)
barplot(rbind(fit4.grridge$lambdamults$random,
              fit4.grEBEN$lambdag$random[, fit4.grEBEN$nouteriter + 1]), beside=TRUE,
        xlab="", axisnames=FALSE,
        args.legend=list(x="topright", fill=c(gray.colors(2), 0, 0), lty=c(NA, NA, 2, 2),
                         border=c(rep(1, 2), 0, 0), merge=TRUE, seg.len=1),
        ylab=expression(paste(lambda[g], "'")),
        legend.text=paste(c("GRridge", "grEBEN", "ridge", "enet"), ", AUC=",
                          round(cons4$auc[c(4, 3, 1, 2)], 2), sep=""))
title(xlab="Random groups", line=1)
abline(h=1, lty=2)
dev.off()

# random groups 2
load(paste(path.res, "grVBEM_channel2_cons_res4.2.Rdata", sep=""))
png(paste(path.graph, "grEBEN_channel_bar4.2.png", sep=""), units="in", width=4.5, height=4, res=200)
barplot(rbind(fit4.grridge$lambdamults$random,
              fit4.grEBEN$lambdag$random[, fit4.grEBEN$nouteriter + 1]), beside=TRUE,
        xlab="", axisnames=FALSE,
        args.legend=list(x="topright", fill=c(gray.colors(2), 0, 0), lty=c(NA, NA, 2, 2),
                         border=c(rep(1, 2), 0, 0), merge=TRUE, seg.len=1),
        ylab=expression(paste(lambda[g], "'")),
        legend.text=paste(c("GRridge", "grEBEN", "ridge", "enet"), ", AUC=",
                          round(cons4.2$auc[c(4, 3, 1, 2)], 2), sep=""))
title(xlab="Random groups", line=1)
abline(h=1, lty=2)
dev.off()

# random groups permuted 100 times
load(paste(path.res, "grVBEM_channel2_cons_res5.Rdata", sep=""))
png(paste(path.graph, "grEBEN_channel_boxplot2.png", sep=""), units="in", width=4.5, height=4, res=400)
par(fig=c(0, 1, 0, 1))
boxplot(cbind(c(cons5$lambdag[, c(1:10)]), c(cons5$lambdag[, c(11:20)])),
        col=gray.colors(2), names=c("GRridge", "grEBEN"))
par(fig=c(0.40, 0.95, 0.19, 0.94), new=T)
boxplot(cbind(c(cons5$lambdag[, c(1:10)]), c(cons5$lambdag[, c(11:20)])), ylim=c(0, 6), 
        col=gray.colors(2), xaxt="n")
abline(h=1, lty=2)
par(fig=c(0, 1, 0, 1))
dev.off()

# random 2 and random groups permuted 100 times
load(paste(path.res, "grVBEM_channel2_cons_res4.Rdata", sep=""))
load(paste(path.res, "grVBEM_channel2_cons_res5.Rdata", sep=""))
png(paste(path.graph, "grEBEN_channel_boxplot1.png", sep=""), units="in", width=9, height=5.5, res=400)
par(fig=c(0, 0.52, 0, 1), new=FALSE)
barplot(rbind(fit4.grridge$lambdamults$random,
              fit4.grEBEN$lambdag$random[, fit4.grEBEN$nouteriter + 1]), beside=TRUE,
        xlab="", axisnames=FALSE,
        ylab=expression(paste(lambda[g], "'")), main="a)")
title(xlab="Random groups", line=1)
abline(h=1, lty=2)
par(fig=c(0.48, 1, 0, 1), new=TRUE)
# boxplot(cbind(c(cons5$lambdag[, c(1:10)]), c(cons5$lambdag[, c(11:20)])),
#         col=gray.colors(2), names=c("GRridge", "grEBEN"), main="b)")
# par(fig=c(0.65, 0.98, 0.21, 0.97), new=TRUE)
boxplot(cbind(c(cons5$lambdag[, c(1:10)]), c(cons5$lambdag[, c(11:20)])), ylim=c(0, 6), 
        col=gray.colors(2), main="b)", names=c("GRridge", "grEBEN"))
abline(h=1, lty=2)
legend("topright", legend=c("GRridge", "grEBEN"), fill=c(gray.colors(2), 0, 0), lty=c(NA, NA, 2, 2),
       border=c(rep(1, 2), 0, 0), merge=TRUE, seg.len=1)
# par(fig=c(0, 1, 0, 1))
dev.off()

