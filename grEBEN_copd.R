##########################  data description  #########################
# Data:                                                               #
# http://www.ebi.ac.uk/metabolights/MTBLS358                          #
#                                                                     #
# Reference:                                                          #
# Alterations in Serum Polyunsaturated Fatty Acids and Eicosanoids in #
# Patients with Mild to Moderate Chronic Obstructive Pulmonary        #
# Disease (COPD),                                                     #
#######################################################################

### paths
path.code <- as.character(ifelse(Sys.info()[1]=="Darwin","/Users/magnusmunch/Documents/PhD/EBEN/code/" ,"~/EBEN/code/"))
path.res <- as.character(ifelse(Sys.info()[1]=="Darwin","/Users/magnusmunch/Documents/PhD/EBEN/results/" ,"~/EBEN/results/"))
path.data <- as.character(ifelse(Sys.info()[1]=="Darwin","/Users/magnusmunch/Documents/PhD/EBEN/data/copd_MTBLS358_20170411_140414/" ,"~/EBEN/data/MTBLS358_20170411_140414/"))
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

### reading in data
data <- read.table(paste(path.data, "s_Study.txt", sep=""), header=TRUE)
CER <- read.table(paste(path.data, "m_CER_mass_spectrometry_v4.maf", sep=""), header=TRUE)
EICO <- read.table(paste(path.data, "m_EICO_mass_spectrometry_v4.maf", sep=""), header=TRUE)
SHOT <- read.table(paste(path.data, "m_SHOT_mass_spectrometry_v4.maf", sep=""), header=TRUE)
TAG <- read.table(paste(path.data, "m_TAG_mass_spectrometry_v4.maf", sep=""), header=TRUE)

### data formatting
# clinical variables
status <- data$Factor.Value.Study.Group # FS: former smokers, CS: current smokers, NS: never-smokers
id <- as.character(data$Sample.Name)

# metabolomics data
CER2 <- data.frame(id=as.character(names(CER)[c(22:176)]), t(CER[, c(22:176)]))
CER2[sapply(CER2, is.factor)] <- lapply(CER2[sapply(CER2, is.factor)], as.character)
colnames(CER2)[-1] <- as.character(CER[, 5])

EICO2 <- data.frame(id=as.character(names(EICO)[c(22:176)]), t(EICO[, c(22:176)]))
EICO2[sapply(EICO2, is.factor)] <- lapply(EICO2[sapply(EICO2, is.factor)], as.character)
colnames(EICO2)[-1] <- as.character(EICO[, 5])

SHOT2 <- data.frame(id=as.character(names(SHOT)[c(22:176)]), t(SHOT[, c(22:176)]))
SHOT2[sapply(SHOT2, is.factor)] <- lapply(SHOT2[sapply(SHOT2, is.factor)], as.character)
colnames(SHOT2)[-1] <- as.character(SHOT[, 5])

TAG2 <- data.frame(id=as.character(names(TAG)[c(22:176)]), t(TAG[, c(22:176)]))
TAG2[sapply(TAG2, is.factor)] <- lapply(TAG2[sapply(TAG2, is.factor)], as.character)
colnames(TAG2)[-1] <- as.character(TAG[, 5])

metabol <- data.frame(CER2[, -1], EICO2[, -1], SHOT2[, -1], TAG2[, -1])

# data preparation
# remove sample/feature if more than 10% missing
# impute half of the minimum of a feature
metabol2 <- metabol[, apply(metabol, 2, function(x) {sum(is.na(x))})/nrow(metabol) < 0.10]
indsel1 <- which(apply(metabol, 2, function(x) {sum(is.na(x))})/nrow(metabol) < 0.10)
metabol3 <- metabol2[apply(metabol2, 1, function(x) {sum(is.na(x))})/ncol(metabol) < 0.1, ]
indsel2 <- which(apply(metabol2, 1, function(x) {sum(is.na(x))})/ncol(metabol) < 0.1)
metabolprep <- apply(metabol3, 2, function(x) {
  imput <- min(x, na.rm=TRUE)/2; ifelse(is.na(x), imput, x)})
statusprep <- status[indsel2]

### data analysis
y <- ifelse(as.numeric(statusprep) > 1, 0, as.numeric(statusprep))
x <- apply(metabolprep, 2, function(x) {(x - mean(x))/sd(x)})
n <- nrow(x)
p <- ncol(x)
m <- rep(1, n)

# co-data
groups1 <- rep(1:4, times=c(ncol(CER2[, -1]), ncol(EICO2[, -1]), ncol(SHOT2[, -1]), 
                            ncol(TAG2[, -1])))[indsel1]

# # cross validation of AUC
# set.seed(3001)
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
#   fit.grEBEN <- grVBEM(xtrain, ytrain, m=mtrain, partitions=groups1, lambda1=lambda1gren, lambda2=lambda2gren,
#                        intercept=TRUE, eps=0.001, maxiter=500, trace=TRUE)
#   fit.grridge <- grridge(t(xtrain), ytrain, list(group=CreatePartition(as.factor(groups1))), unpenal=~1,
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
# copd1 <- list(pred=pred1, auc=auc1)
# 
# save(copd1, file=paste(path.res, "grVBEM_copd_res1.Rdata", sep=""))

fit1.grEBEN <- grVBEM(x, y, m, list(platform=groups1), lambda1=NULL, lambda2=NULL, intercept=TRUE, eps=0.001,
                      maxiter=300, trace=TRUE)
fit1.grridge <- grridge(t(x), y, list(assay=CreatePartition(as.factor(groups1))), unpenal=~1)


### plots
# with platform as co-data
load(paste(path.res, "grVBEM_copd_res1.Rdata", sep=""))
xlabels1 <- c("CER", "EICO", "SHOT", "TAG")
barplot(rbind(fit1.grridge$lambdamults$platform,
              fit1.grEBEN$lambdag$platform[, fit1.grEBEN$nouteriter + 1]), beside=TRUE,
        args.legend=list(x="topleft", fill=c(col=gray.colors(2), rep("#FFFFFF", 2)),
                         border=c(rep(c("black", "white"), each=2))), 
        ylab=expression(paste(lambda[g], "'")),
        legend.text=paste(c("GRridge", "grEBEN", "ridge", "enet"), ", AUC=",
                          round(copd1$auc, 2)[c(4, 3, 1, 2)], sep=""))
abline(h=1, lty=2)
