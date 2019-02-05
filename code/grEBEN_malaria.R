##########################  data description  #########################
# Data:                                                               #
# http://www.ebi.ac.uk/metabolights/MTBLS315                          #
#                                                                     #
# Reference:                                                          #
# Towards improving point-of-care diagnosis of non-malaria febrile    #
# illness: a metabolomics approach                                    #
#######################################################################

### paths
path.code <- as.character(ifelse(Sys.info()[1]=="Darwin","/Users/magnusmunch/Documents/PhD/EBEN/code/" ,"~/EBEN/code/"))
path.res <- as.character(ifelse(Sys.info()[1]=="Darwin","/Users/magnusmunch/Documents/PhD/EBEN/results/" ,"~/EBEN/results/"))
path.data <- as.character(ifelse(Sys.info()[1]=="Darwin","/Users/magnusmunch/Documents/PhD/EBEN/data/malaria_MTBLS315_20170412_084021/" ,"~/EBEN/data/MTBLS358_20170411_140414/"))
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
data <- read.table(paste(path.data, "s_NMFI and BSI diagnosis.txt", sep=""), header=TRUE)
GC <- read.table(paste(path.data, "m_GC_nmfi_and_bsi_diagnosis_v2_maf.tsv", sep=""), header=TRUE)
LC <- read.table(paste(path.data, "m_LC_nmfi_and_bsi_diagnosis_v2_maf.tsv", sep=""), header=TRUE)
UPLC_NEG <- read.table(paste(path.data, "m_UPLC_NEG_nmfi_and_bsi_diagnosis_v2_maf.tsv", sep=""), header=TRUE)
UPLC_POS <- read.table(paste(path.data, "m_UPLC_POS_nmfi_and_bsi_diagnosis_v2_maf.tsv", sep=""), header=TRUE)

### data formatting
# clinical variables
status <- data$Factor.Value.patient.group. 
id <- as.numeric(substring(data$Source.Name, 9))
id2 <- gsub("^.*A|^.*B0", "", data$Sample.Name)

# metabolomics data
GC2 <- data.frame(id=id, t(GC[, c(22:82)])[order(substring(names(GC)[22:82], 2)), ])
GC2[sapply(GC2, is.factor)] <- lapply(GC2[sapply(GC2, is.factor)], as.character)
colnames(GC2)[-1] <- as.character(GC[, 5])

LC2 <- data.frame(id=id, t(LC[22:82])[
  order(gsub("_", "", substr(gsub("^.*A|^.*B0", "", names(LC)[22:82]), 1, 3))), ])
LC2[sapply(LC2, is.factor)] <- lapply(LC2[sapply(LC2, is.factor)], as.character)
colnames(LC2)[-1] <- as.character(LC[, 5])

UPLC_NEG2 <- data.frame(id=id, t(UPLC_NEG[, c(22:82)])[
  order(as.numeric(substr(names(UPLC_NEG)[22:82], 4, 5))), ])
UPLC_NEG2[sapply(UPLC_NEG2, is.factor)] <- lapply(UPLC_NEG2[sapply(UPLC_NEG2, is.factor)], 
                                                  as.character)
colnames(UPLC_NEG2)[-1] <- as.character(UPLC_NEG[, 5])

UPLC_POS2 <- data.frame(id=id, t(UPLC_POS[, c(22:82)])[
  order(as.numeric(substr(names(UPLC_POS)[22:82], 4, 5))), ])
UPLC_POS2[sapply(UPLC_POS2, is.factor)] <- lapply(UPLC_POS2[sapply(UPLC_NEG2, is.factor)], 
                                                  as.character)
colnames(UPLC_POS2)[-1] <- as.character(UPLC_POS[, 5])

metabol <- data.frame(GC2[, -1], LC2[, -1], UPLC_NEG2[, -1], UPLC_POS2[, -1])
metabolun <- unname(metabol)


# data preparation
# remove sample/feature if more than 10% missing
# impute half of the minimum of a feature
metabol2 <- metabolun[, apply(metabolun, 2, function(x) {sum(is.na(x))})/nrow(metabolun) < 0.10]
indsel1 <- which(apply(metabolun, 2, function(x) {sum(is.na(x))})/nrow(metabolun) < 0.10)
metabol3 <- metabol2[apply(metabol2, 1, function(x) {sum(is.na(x))})/ncol(metabolun) < 0.1, ]
indsel2 <- which(apply(metabol2, 1, function(x) {sum(is.na(x))})/ncol(metabolun) < 0.1)
metabolprep <- apply(metabol3, 2, function(x) {
  imput <- min(x, na.rm=TRUE)/2; ifelse(is.na(x), imput, x)})
statusprep <- status[indsel2]

### data analysis
y <- ifelse(as.numeric(statusprep)!=2, 0, 1)
x <- apply(metabolprep, 2, function(x) {(x - mean(x))/sd(x)})
n <- nrow(x)
p <- ncol(x)
m <- rep(1, n)

# co-data
assay <- rep(1:4, times=c(ncol(GC2[, -1]), ncol(LC2[, -1]), ncol(UPLC_NEG2[, -1]), 
                          ncol(UPLC_POS2[, -1])))[indsel1]

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

fit1.grEBEN <- grVBEM(x, y, m, list(assay=assay), lambda1=NULL, lambda2=NULL, intercept=TRUE, eps=0.001,
                      maxiter=300, trace=TRUE)
fit1.grridge <- grridge(t(x), y, list(conservation=CreatePartition(as.factor(assay))), unpenal=~1)
fit1.2 <- grridgeCV(fit1.grridge, t(x), y)
pROC::roc(y, fit1.2$NoGroups)$auc
pROC::roc(y, fit1.2$GroupRegul)$auc




