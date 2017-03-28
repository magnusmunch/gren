##############################  preamble  #############################
# grMCEM on real data                                                 #
# version: 01                                                         #
# author: Magnus Münch                                                #
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
path.rcode <- "~/EBEN/code/"
path.res <- "~/EBEN/results/"
path.data <- "~/EBEN/data/"

### libraries
library(glmnet)
library(penalized)
library(GRridge)
library(SGL)
library(e1071)

### functions
# source grENVB functions
source(paste(path.rcode, "grVBEM.R", sep=""))

# function to cross-validate elastic net
cv.en <- function(x, y, intercept=intercept) {
  fit.pen <- cv.pen(x, y, intercept=intercept)
  fit.en <- glmnet(x, y, family="binomial", alpha=fit.pen$alpha[which.min(fit.pen$cvll)],
                   lambda=fit.pen$lambda[which.min(fit.pen$cvll)])
  return(fit.en)
}

########## using conservation status as grouping information
# using the group information
load(paste(path.data, "mirsData.RData", sep=""))
parCons <- CreatePartition(mirsData$conservation) # using conservation status as grouping
set.seed(789)
groups <- rep(1:length(parCons), unlist(lapply(parCons, length)))
x <- apply(t(as.matrix(mirsData$transformedData))[, unlist(parCons)], 2, 
           function(x) {(x - mean(x))/sd(x)})       
y <- as.numeric(mirsData$response) - 1
n <- nrow(x)
p <- ncol(x)
m <- rep(1, n)
nfolds <- n
rest <- n %% nfolds
foldsize <- c(rep(n %/% nfolds + as.numeric(rest!=0), times=rest),
              rep(n %/% nfolds, times=nfolds - rest))
foldid <- sample(rep(1:nfolds, times=foldsize))

pred1.ridge <- pred1.lasso <- pred1.en <- pred1.svm <- pred1.grVBEM <- pred1.gren <- 
  pred1.grridge <- pred1.SGL <- numeric(n)
for(k in 12:nfolds) {
  print(paste("Fold", k, sep=" "))
  xtrain <- x[foldid!=k, ]
  xtest <- x[foldid==k, ]
  ytrain <- y[foldid!=k]
  ytest <- y[foldid==k]
  mtrain <- m[foldid!=k]
  mtest <- m[foldid==k]
  ntrain <- length(ytrain)
  ntest <- length(ytest)
  
  # grVBEM
  fit.ridge <- cv.glmnet(xtrain, ytrain, family="binomial", alpha=0, standardize=FALSE, intercept=TRUE)
  fit.lasso <- cv.glmnet(xtrain, ytrain, family="binomial", alpha=1, standardize=FALSE, intercept=TRUE)
  fit.en <- cv.en(xtrain, ytrain, intercept=TRUE)
  fit.svm <- svm(xtrain, ytrain, type="C-classification", probability=TRUE)
  fit.grVBEM <- grVBEM(xtrain, ytrain, m=mtrain, groups=groups, lambda1=NULL, 
                       lambda2=NULL, intercept=TRUE, eps=0.001, maxiter=500)
  fit.gren <- penalized(ytrain, xtrain, unpenalized=~1, model="logistic",
                        lambda1=0.5*fit.grVBEM$lambda1*sqrt(rep(fit.grVBEM$lambdag[, fit.grVBEM$nouteriter + 1], times=rle(groups)$lengths)),
                        lambda2=0.5*fit.grVBEM$lambda2*rep(fit.grVBEM$lambdag[, fit.grVBEM$nouteriter + 1], times=rle(groups)$lengths))
  fit.grridge <- grridge(t(xtrain), ytrain, list(group=CreatePartition(as.factor(groups))), unpenal=~1, trace=FALSE)
  fit.SGL <- cvSGL(list(x=xtrain, y=ytrain), index=groups, type="logit", standardize=FALSE)
  
  pred1.ridge[foldid==k] <- predict(fit.ridge, matrix(xtest, nrow=1), type="response")
  pred1.lasso[foldid==k] <- predict(fit.lasso, matrix(xtest, nrow=1), type="response")
  pred1.en[foldid==k] <- predict(fit.en, matrix(xtest, nrow=1), type="response")
  pred1.svm[foldid==k] <- attributes(predict(fit.svm, matrix(xtest, nrow=1), probability=TRUE))$probabilities[2]
  pred1.grVBEM[foldid==k] <- as.numeric(exp(cbind(1, matrix(xtest, nrow=1)) %*% fit.grVBEM$mu)/
                                           (1 + exp(cbind(1, matrix(xtest, nrow=1)) %*% fit.grVBEM$mu)))
  pred1.gren[foldid==k] <- predict(fit.gren, matrix(xtest, nrow=1))
  pred1.grridge[foldid==k] <- predict(fit.grridge$predobj$GroupRegul, matrix(xtest, nrow=1))
  pred1.SGL[foldid==k] <- as.numeric(exp(cbind(1, matrix(xtest, nrow=1)) %*% 
                                           rbind(fit.SGL$fit$intercepts, fit.SGL$fit$beta)
                                         [, which.min(fit.SGL$lldiff)])/
                                       (1 + exp(cbind(1, matrix(xtest, nrow=1)) %*% 
                                                  rbind(fit.SGL$fit$intercepts, fit.SGL$fit$beta)
                                                [, which.min(fit.SGL$lldiff)])))
  
}

auc1.ridge <- pROC::roc(y, pred1.ridge)$auc
auc1.lasso <- pROC::roc(y, pred1.lasso)$auc
auc1.en <- pROC::roc(y, pred1.en)$auc
auc1.svm <- pROC::roc(y, pred1.svm)$auc
auc1.grVBEM <- pROC::roc(y, pred1.grVBEM)$auc
auc1.gren <- pROC::roc(y, pred1.gren)$auc
auc1.grridge <- pROC::roc(y, pred1.grridge)$auc
auc1.SGL <- pROC::roc(y, pred1.SGL)$auc

pred1 <- cbind(ridge=pred1.ridge, lasso=pred1.lasso, en=pred1.en, svm=pred1.svm, grVBEM=pred1.grVBEM,
               gren=pred1.gren, grridge=pred1.grridge, SGL=pred1.SGL)
auc1 <- c(ridge=auc1.ridge, lasso=auc1.lasso, en=auc1.en, svm=auc1.svm, grVBEM=auc1.grVBEM,
          gren=auc1.gren, grridge=auc1.grridge, SGL=auc1.SGL)
cons1 <- list(pred=pred1, auc=auc1)

save(cons1, file=paste(path.res, "grVBEM_channel_cons_res1.Rdata", sep=""))

# # data 4 (permuted data (no info))
# load(paste(path.data, "mirsData.RData", sep=""))
# parCons <- CreatePartition(mirsData$conservation) # using conservation status as grouping
# set.seed(123)
# n <- ncol(mirsData$transformedData)
# p <- nrow(mirsData$transformedData)
# groups <- rep(1:length(parCons), unlist(lapply(parCons, length)))
# x <- apply(t(as.matrix(mirsData$transformedData))[, unlist(parCons)], 2, 
#            function(x) {(x - mean(x))/sd(x)})[, sample(1:p)]       
# y <- as.numeric(mirsData$response) - 1
# m <- rep(1, n)
# fit.optL2 <- optL2(y, x, unpenalized=~1, lambda1=0, model="logistic")
# b0 <- coef(fit.optL2$fullfit)
# pred.b0 <- as.numeric(exp(cbind(1, x) %*% b0)/(1 + exp(cbind(1, x) %*% b0)))
# W <- diag(sqrt(pred.b0*(1 - pred.b0)))
# Xw <- W %*% cbind(1, x)
# invmat <- solve(t(Xw) %*% Xw + diag(c(0, rep(2*fit.optL2$lambda, p))))
# sigma0 <- invmat %*% t(Xw) %*% Xw %*% invmat
# nfolds <- 10
# rest <- n %% nfolds
# foldsize <- c(rep(n %/% nfolds + as.numeric(rest!=0), times=rest),
#               rep(n %/% nfolds, times=nfolds - rest))
# foldid <- sample(rep(1:nfolds, times=foldsize))
# 
# pred20.grVBEM <- numeric(n)
# for(k in 1:nfolds) {
#   print(paste("Fold", k, sep=" "))
#   xtrain <- x[foldid!=k, ]
#   xtest <- x[foldid==k, ]
#   ytrain <- y[foldid!=k]
#   mtrain <- rep(1, times=sum(foldid!=k))
#   
#   fit.grVBEM <- grVBEM(xtrain, ytrain, m=mtrain, groups=groups, lambda1=NULL, 
#                        lambda2=NULL, sigma0=sigma0, intercept=TRUE, 
#                        eps=0.001, maxiter=500)
#   best <- fit.grVBEM$mu
#   pred20.grVBEM[foldid==k] <- as.numeric(exp(cbind(1, xtest) %*% best)/
#                                            (1 + exp(cbind(1, xtest) %*% best)))
#   
# }
# 
# # test20.grVBEM <- cv.grVBEM(x, y, m, groups, lambda1=NULL, lambda2=NULL, sigma0, 
# #                            intercept=TRUE, eps=0.001, maxiter=500, nfolds=nfolds,
# #                            foldid=foldid)
# test20a.grridge <- grridge(t(x), y, list(group=CreatePartition(as.factor(groups))), unpenal=~1)
# test20b.grridge <- grridgeCV(test20a.grridge, t(x), y, outerfold=foldid, fixedfolds=TRUE)
# 
# auc20.grVBEM <- pROC::roc(y, pred20.grVBEM)$auc
# auc20.grridge <- pROC::roc(y, test20b.grridge[, 3])$auc
# 

load("C:/Users/Magnus/Documents/phd/ENVB/results/grVBEM_channel_cons_res1.Rdata")
cons1$auc








