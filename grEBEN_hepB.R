##########################  data description  #########################
# Data:                                                               #
# http://www.ebi.ac.uk/metabolights/MTBLS253                          #
# http://www.ebi.ac.uk/metabolights/MTBLS279                          #
# http://www.ebi.ac.uk/metabolights/MTBLS280                          #
#                                                                     #
# Reference:                                                          #
# Metabolic characterization of the natural progression of chronic    #
# Hepatitis B                                                         #
#######################################################################

### paths
path.code <- as.character(ifelse(Sys.info()[1]=="Darwin","/Users/magnusmunch/Documents/PhD/EBEN/code/" ,"~/EBEN/code/"))
path.res <- as.character(ifelse(Sys.info()[1]=="Darwin","/Users/magnusmunch/Documents/PhD/EBEN/results/" ,"~/EBEN/results/"))
path.data <- as.character(ifelse(Sys.info()[1]=="Darwin","/Users/magnusmunch/Documents/PhD/EBEN/data/hepB/" ,"~/EBEN/data/hepB/"))
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
data_oxy <- read.table(paste(path.data, "s_mtbls253.txt", sep=""), header=TRUE)
data_ami_acy <- read.table(paste(path.data, "s_mtbls280.txt", sep=""), header=TRUE)
data_lip <- read.table(paste(path.data, "s_mtbls279.txt", sep=""), header=TRUE)

met_ami <- read.table(paste(path.data, "m_amine_mass_spectrometry_v2_maf.tsv", sep=""), dec=",", sep="\t", header=TRUE)
met_POS_lip <- read.table(paste(path.data, "m_POS_lipid_mass_spectrometry_v2_maf.tsv", sep=""), dec=",", sep="\t", header=TRUE)
met_NEG_lip <- read.table(paste(path.data, "m_NEG_lipid_mass_spectrometry_v2_maf.tsv", sep=""), dec=",", sep="\t", header=TRUE)
met_oxy <- read.table(paste(path.data, "m_oxylipin_mass_spectrometry_v2_maf.tsv", sep=""), dec=",", sep="\t", header=TRUE)
met_acy <- read.table(paste(path.data, "m_acyl_carnitine_mass_spectrometry_v2_maf.tsv", sep=""), dec=",", sep="\t", header=TRUE)

### data formatting
# clinical stuff
status1_oxy <- data_oxy$Factor.Value.control
status1_oxy <- factor(c(status1_oxy[1:37], NA, status1_oxy[38:66], NA, 
                      status1_oxy[67:length(status1_oxy)]), levels=c(1, 2), 
                      labels=c("False", "True"))
status2_oxy <- data_oxy$Factor.Value.phase.
status2_oxy <- factor(c(status2_oxy[1:37], NA, status2_oxy[38:66], NA, 
                        status2_oxy[67:length(status2_oxy)]), levels=c(1, 2, 3, 4, 5), 
                      labels=c("", "ENEG", "IA", "IC", "IT"))
id_oxy <- data_oxy$Sample.Name
id_oxy <- c(id_oxy[1:37], 78, id_oxy[38:66], 76, id_oxy[67:length(id_oxy)])
# oxy is missing patients 76 and 78

status1_ami_acy <- data_ami_acy$Factor.Value.control
status2_ami_acy <- data_ami_acy$Factor.Value.phase.
id_ami_acy <- data_ami_acy$Sample.Name

status1_lip <- data_lip$Factor.Value.control
status2_lip <- data_lip$Factor.Value.phase.
id_lip <- data_lip$Sample.Name
# lip and ami contain the same samples

# metabolomics data
met_ami2 <- t(met_ami[, c(22:109)])
met_ami3 <- apply(met_ami2, 2, function(x) {
  as.numeric(ifelse(substr(x, 1, 2)=="0,", gsub(",", ".", x), gsub(",", "", x)))})
met_ami4 <- data.frame(id=as.numeric(substring(names(met_ami)[22:109], 2)), met_ami3)
colnames(met_ami4)[-1] <- as.character(met_ami[, 5])

met_POS_lip2 <- t(met_POS_lip[, c(22:109)])
met_POS_lip3 <- apply(met_POS_lip2, 2, function(x) {
  as.numeric(ifelse(substr(x, 1, 2)=="0,", gsub(",", ".", x), gsub(",", "", x)))})
met_POS_lip4 <- data.frame(id=as.numeric(substring(names(met_POS_lip)[22:109], 2)), met_POS_lip3)
colnames(met_POS_lip4)[-1] <- as.character(met_POS_lip[, 5])

met_NEG_lip2 <- t(met_NEG_lip[, c(22:109)])
met_NEG_lip3 <- apply(met_NEG_lip2, 2, function(x) {
  as.numeric(ifelse(substr(x, 1, 2)=="0,", gsub(",", ".", x), gsub(",", "", x)))})
met_NEG_lip4 <- data.frame(id=as.numeric(substring(names(met_NEG_lip)[22:109], 2)), met_NEG_lip3)
colnames(met_NEG_lip4)[-1] <- as.character(met_NEG_lip[, 5])

met_oxy2 <- t(met_oxy[, c(22:107)])
met_oxy3 <- apply(met_oxy2, 2, function(x) {
  as.numeric(ifelse(substr(x, 1, 2)=="0,", gsub(",", ".", x), gsub(",", "", x)))})
met_oxy4 <- data.frame(id=as.numeric(substring(names(met_oxy)[22:107], 2)), met_oxy3)
colnames(met_oxy4)[-1] <- as.character(met_oxy[, 5])
met_oxy4 <- rbind(met_oxy4[c(1:37), ], c(78, rep(NA, ncol(met_oxy4) - 1)), 
                  met_oxy4[c(38:66), ], c(76, rep(NA, ncol(met_oxy4) - 1)),
                  met_oxy4[c(67:nrow(met_oxy4)), ])

met_acy2 <- t(met_acy[, c(22:109)])
met_acy3 <- apply(met_acy2, 2, function(x) {
  as.numeric(ifelse(substr(x, 1, 2)=="0,", gsub(",", ".", x), gsub(",", "", x)))})
met_acy4 <- data.frame(id=as.numeric(substring(names(met_acy)[22:109], 2)), met_acy3)
colnames(met_acy4)[-1] <- as.character(met_acy[, 5])

metabol <- data.frame(met_ami4[, -1], met_POS_lip4[, -1], met_NEG_lip4[, -1], 
                      met_oxy4[, -1], met_acy4[, -1])


# data preparation
# remove sample/feature if more than 10% missing
# impute half of the minimum of a feature
metabol2 <- metabol[, apply(metabol, 2, function(x) {sum(is.na(x))})/nrow(metabol) < 0.10]
indsel1 <- which(apply(metabol, 2, function(x) {sum(is.na(x))})/nrow(metabol) < 0.10)
metabol3 <- metabol2[apply(metabol2, 1, function(x) {sum(is.na(x))})/ncol(metabol) < 0.1, ]
indsel2 <- which(apply(metabol2, 1, function(x) {sum(is.na(x))})/ncol(metabol) < 0.1)
metabolprep <- apply(metabol3, 2, function(x) {
  imput <- min(x, na.rm=TRUE)/2; ifelse(is.na(x), imput, x)})
statusprep <- status1_ami_acy[indsel2]

### data analysis
y <- ifelse(as.numeric(statusprep)==2, 0, 1)
x <- apply(metabolprep, 2, function(x) {(x - mean(x))/sd(x)})
n <- nrow(x)
p <- ncol(x)
m <- rep(1, n)

# co-data
assay <- rep(1:5, times=c(ncol(met_ami4[, -1]), ncol(met_POS_lip4[, -1]), ncol(met_NEG_lip4[, -1]), 
                          ncol(met_oxy4[, -1]), ncol(met_acy4[, -1])))[indsel1]

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
#   fit.grEBEN <- grVBEM(xtrain, ytrain, m=mtrain, partitions=assay, lambda1=lambda1gren, lambda2=lambda2gren,
#                        intercept=TRUE, eps=0.001, maxiter=500, trace=TRUE)
#   fit.grridge <- grridge(t(xtrain), ytrain, list(group=CreatePartition(as.factor(assay))), unpenal=~1,
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
# hepB1 <- list(pred=pred1, auc=auc1)
# 
# save(hepB1, file=paste(path.res, "grVBEM_hepB_res1.Rdata", sep=""))

set.seed(3001)
fit1.grEBEN <- grVBEM(x, y, m, list(assay=assay), lambda1=NULL, lambda2=NULL, intercept=TRUE, eps=0.001,
                      maxiter=300, trace=TRUE)
fit1.grridge <- grridge(t(x), y, list(assay=CreatePartition(as.factor(assay))), unpenal=~1)



### plots
# with platform as co-data
load(paste(path.res, "grVBEM_hepB_res1.Rdata", sep=""))
xlabels1 <- c("Biogenic amine", "Positive lipid", "Negative lipid", "Oxylipins", "Acyl-carnitines")
par(fig=c(0, 1, 0, 1))
barplot(rbind(fit1.grridge$lambdamults$assay,
              fit1.grEBEN$lambdag$assay[, fit1.grEBEN$nouteriter + 1]), beside=TRUE,
        names.arg=xlabels1, 
        args.legend=list(x="topleft", fill=c(gray.colors(2), 0, 0), lty=c(NA, NA, 2, 2),
                         border=c(rep(1, 2), 0, 0), merge=TRUE, seg.len=1),
        ylab=expression(paste(lambda[g], "'")),
        legend.text=paste(c("GRridge", "grEBEN", "ridge", "enet"), ", AUC=",
                          round(hepB1$auc, 2)[c(4, 3, 1, 2)], sep=""))
par(fig=c(0.1, 0.62, 0.08, 0.73), new=T)
barplot(rbind(fit1.grridge$lambdamults$assay,
              fit1.grEBEN$lambdag$assay[, fit1.grEBEN$nouteriter + 1]), beside=TRUE,
        ylim=c(0, max(c(fit1.grridge$lambdamults$assay[1:3], 
                        fit1.grEBEN$lambdag$assay[, fit1.grEBEN$nouteriter + 1])) + 1),
        xpd=FALSE, names.arg=rep(NA, 5))
box()
abline(h=1, lty=2)
par(fig=c(0, 1, 0, 1))

gray.colors(2)
