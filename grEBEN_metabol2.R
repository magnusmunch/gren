##############################  preamble  #############################
# grEBEN on real data                                                 #
# version: 01                                                         #
# author: Magnus Munch                                                #
# created: 04-04-2017                                                 #
# last edited: 04-04-2017                                             #
#######################################################################

###############################  notes  ###############################
#######################################################################

##########################  data description  #########################
# 248 samples of cases of Alzheimers and controls with measurements   #
# on 230 metabolites.                                                 #
#######################################################################

### paths
path.data <- as.character(ifelse(Sys.info()[1]=="Darwin","/Users/magnusmunch/Documents/PhD/EBEN/data/" ,"~/EBEN/data/"))
path.code <- as.character(ifelse(Sys.info()[1]=="Darwin","/Users/magnusmunch/Documents/PhD/EBEN/code/" ,"~/EBEN/code/"))
path.res <- as.character(ifelse(Sys.info()[1]=="Darwin","/Users/magnusmunch/Documents/PhD/EBEN/results/" ,"~/EBEN/results/"))

### libraries
library(glmnet)
library(penalized)
library(GRridge)
library(pROC)
library(Biobase)

### loading data
load(paste(path.data, "ESetMbolCSFPR2.Rdata", sep=""))

### loading functions 
source(paste(path.code, "grVBEM.R", sep=""))

# cross-validate to find penalty parameters
cv.pen <- function(x, y, intercept, psel) {
  
  n <- nrow(x)
  seq.alpha <- seq(0.01, 0.99, length.out=50)
  seq.lambda <- seq.df <- seq.cvm <- numeric(length(seq.alpha))
  
  if(is.null(psel)) {
    for(a in 1:length(seq.alpha)) {
      fit.cv <- cv.glmnet(x, y, family="binomial", alpha=seq.alpha[a], standardize=FALSE,
                          intercept=intercept)
      ind <- which(fit.cv$lambda==fit.cv$lambda.min)
      seq.lambda[a] <- fit.cv$lambda.min
      seq.df[a] <- fit.cv$nzero[ind]
      seq.cvm[a] <- fit.cv$cvm[ind]
    }
  } else {
    fit.glmnet <- glmnet(x, y, family="binomial", alpha=seq.alpha[a], standardize=FALSE, 
                         intercept=intercept, dfmax=psel, nlambda=floor(5000*20/ncol(x)))
    fit.lambda <- tail(fit.glmnet$lambda, n=2L)
    if(length(fit.lambda)==1) {
      fit.lambda <- c(fit.lambda - 0.1, fit.lambda)
    }
    fit.cv <- cv.glmnet(x, y, family="binomial", alpha=seq.alpha[a], intercept=intercept, 
                        lambda=fit.lambda)
    seq.lambda[a] <- fit.cv$lambda[2]
    seq.df[a] <- fit.cv$nzero[2]
    seq.cvm[a] <- fit.cv$cvm[2]
  }
  lambda1 <- 2*n*seq.alpha*seq.lambda
  #lambda1 <- n*seq.alpha*seq.lam
  lambda2 <- n*(1 - seq.alpha)*seq.lambda
  
  out <- list(lambda1=lambda1, lambda2=lambda2, alpha=seq.alpha, lambda=seq.lambda, cvll=seq.cvm)
  return(out)
}

### data analysis
# phenotype data
y <- as.numeric(pData(ESetMbolCSFPR2)$D_diag_name) - 1

# metabolite data
x <- apply(t(exprs(ESetMbolCSFPR2)), 2, function(x) {(x - mean(x))/sd(x)})

# co-data (two partitions)
partitions1 <- list(platform=replace(as.numeric(fData(ESetMbolCSFPR2)$PlatformCode),
                                     as.numeric(fData(ESetMbolCSFPR2)$PlatformCode)==5, 4))
partitions2 <- list(platform=as.numeric(fData(ESetMbolCSFPR2)$PlatformCode))

# more info
n <- nrow(x)
p <- ncol(x)
m <- y*y + (y - 1)*(y - 1)
nfolds <- 10
rest <- n %% nfolds
foldsize <- c(rep(n %/% nfolds + as.numeric(rest!=0), times=rest),
              rep(n %/% nfolds, times=nfolds - rest))
set.seed(1002)
foldid <- sample(rep(1:nfolds, times=foldsize))

# cross validation of AUC
cv.ridge <- cv.glmnet(x, y, family="binomial", alpha=0, standardize=FALSE, intercept=TRUE)
cv.en <- cv.pen(x, y, intercept=TRUE, psel=NULL)
glmnet(x, y, family="binomial", alpha=cv.en$alpha[which.min(cv.en$cvll)],
       lambda=cv.en$lambda[which.min(cv.en$cvll)], standardize=FALSE, intercept=TRUE)

lambda2ridge <- 0.5*n*cv.ridge$lambda.min
lambdaridge <- cv.ridge$lambda.min
alphaglmnet <- cv.en$alpha[which.min(cv.en$cvll)]
lambdaglmnet <- cv.en$lambda[which.min(cv.en$cvll)]
lambda1gren <- 2*n*lambdaglmnet*alphaglmnet
lambda2gren <- n*lambdaglmnet*(1 - alphaglmnet)

pred.ridge <- pred.en <- pred1.grEBEN <- pred1.grridge <- 
  pred2.grEBEN <- pred2.grridge <- numeric(n) # pred1.SGL <- numeric(n)
for(k in 1:nfolds) {
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
  fit.ridge <- glmnet(xtrain, ytrain, family="binomial", alpha=0, lambda=lambdaridge, standardize=FALSE, 
                      intercept=TRUE)
  fit.en <- glmnet(xtrain, ytrain, family="binomial", alpha=alphaglmnet, lambda=lambdaglmnet, standardize=FALSE,
                   intercept=TRUE)
  fit1.grEBEN <- grVBEM(xtrain, ytrain, mtrain, partitions1, lambda1=lambda1gren, lambda2=lambda2gren, 
                        intercept=TRUE, eps=0.001, maxiter=500, trace=TRUE, QNacc=FALSE)
  fit2.grEBEN <- grVBEM(xtrain, ytrain, mtrain, partitions2, lambda1=lambda1gren, lambda2=lambda2gren, 
                        intercept=TRUE, eps=0.001, maxiter=500, trace=TRUE, QNacc=FALSE)
  fit1.grridge <- grridge(t(xtrain), ytrain, list(platform=CreatePartition(as.factor(partitions1[[1]]))), 
                          unpenal=~1, optl=lambda2ridge, trace=FALSE)
  fit2.grridge <- grridge(t(xtrain), ytrain, list(platform=CreatePartition(as.factor(partitions2[[1]]))), 
                          unpenal=~1, optl=lambda2ridge, trace=FALSE)
  # fit.SGL <- SGL(list(x=xtrain, y=ytrain), index=groups, type="logit", standardize=FALSE)
  
  pred.ridge[foldid==k] <- predict(fit.ridge, xtest, type="response")
  pred.en[foldid==k] <- predict(fit.en, xtest, type="response")
  pred1.grEBEN[foldid==k] <- as.numeric(exp(cbind(1, xtest) %*% fit1.grEBEN$mu)/
                                          (1 + exp(cbind(1, xtest) %*% fit1.grEBEN$mu)))
  pred2.grEBEN[foldid==k] <- as.numeric(exp(cbind(1, xtest) %*% fit2.grEBEN$mu)/
                                          (1 + exp(cbind(1, xtest) %*% fit2.grEBEN$mu)))
  pred1.grridge[foldid==k] <- predict(fit1.grridge$predobj$GroupRegul, xtest)
  pred2.grridge[foldid==k] <- predict(fit2.grridge$predobj$GroupRegul, xtest)
  # pred1.SGL[foldid==k] <- as.numeric(exp(cbind(1, matrix(xtest, nrow=1)) %*% 
  #                                          rbind(fit.SGL$fit$intercepts, fit.SGL$fit$beta)
  #                                        [, which.min(fit.SGL$lldiff)])/
  #                                      (1 + exp(cbind(1, matrix(xtest, nrow=1)) %*% 
  #                                                 rbind(fit.SGL$fit$intercepts, fit.SGL$fit$beta)
  #                                               [, which.min(fit.SGL$lldiff)])))
  
}

auc.ridge <- pROC::roc(y, pred.ridge)$auc
auc.en <- pROC::roc(y, pred.en)$auc
auc1.grEBEN <- pROC::roc(y, pred1.grEBEN)$auc
auc2.grEBEN <- pROC::roc(y, pred2.grEBEN)$auc
auc1.grridge <- pROC::roc(y, pred1.grridge)$auc
auc2.grridge <- pROC::roc(y, pred2.grridge)$auc
# auc1.SGL <- pROC::roc(y, pred1.SGL)$auc

pred <- cbind(ridge=pred.ridge, en=pred.en, grEBEN1=pred1.grEBEN, grEBEN2=pred2.grEBEN,
              grridge1=pred1.grridge, grridge2=pred2.grridge)# , SGL=pred1.SGL)
auc <- c(ridge=auc.ridge, en=auc.en, grEBEN1=auc1.grEBEN, grEBEN2=auc2.grEBEN, 
         grridge1=auc1.grridge, grridge2=auc2.grridge)# , SGL=auc1.SGL)
metabol <- list(pred=pred1, auc=auc1)

save(metabol, file=paste(path.res, "grEBEN_metabol_res1.Rdata", sep=""))



# # co-data
# # combining the two oxidative stress groups
# partitions1 <- list(platform=replace(as.numeric(fData(ESetMbolCSFPR2)$PlatformCode),
#                                      as.numeric(fData(ESetMbolCSFPR2)$PlatformCode)==5, 4))
# fit1.grEBEN <- grVBEM(x, y, m, partitions1, lambda1=NULL, lambda2=NULL, intercept=TRUE, eps=0.001, 
#                       maxiter=300, trace=TRUE, QNacc=FALSE)
# fit1.grridge <- grridge(t(x), y, list(platform=CreatePartition(as.factor(partitions1[[1]]))), 
#                         unpenal=~1)
# 
# # using the two oxidative stress levels seperately
# partitions2 <- list(platform=as.numeric(fData(ESetMbolCSFPR2)$PlatformCode))
# fit2.grEBEN <- grVBEM(x, y, m, partitions2, lambda1=NULL, lambda2=NULL, intercept=TRUE, eps=0.001, 
#                       maxiter=300, trace=TRUE, QNacc=FALSE)
# fit2.grridge <- grridge(t(x), y, list(platform=CreatePartition(as.factor(partitions2[[1]]))), 
#                         unpenal=~1)
# 
# ### plots
# # two oxidative stress groups combined
# xlabels1 <- levels(fData(ESetMbolCSFPR2)$Platform)[-5]
# xlabels1[4] <- "Oxidative Stress"
# barplot(rbind(fit1.grridge$lambdamults[[1]], 
#               fit1.grEBEN$lambdag[[1]][, fit1.grEBEN$nouteriter + 1]), beside=TRUE,
#         names.arg=xlabels1, legend.text=c("GRridge", "grEBEN"),
#         args.legend=list(x="topright"))
# abline(h=1, lty=2)
# 
# # separate oxidative stress groups
# xlabels2 <- gsub("\\s+Stress", " \n Stress", levels(fData(ESetMbolCSFPR2)$Platform))
# barplot(rbind(fit2.grridge$lambdamults[[1]], 
#               fit2.grEBEN$lambdag[[1]][, fit2.grEBEN$nouteriter + 1]), beside=TRUE,
#         names.arg=xlabels2, legend.text=c("GRridge", "grEBEN"),
#         args.legend=list(x="topleft"))
# abline(h=1, lty=2)
