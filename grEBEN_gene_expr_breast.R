################################### preamble ###################################
# breast cancer data from R-package cancerdata                                 #
# version: 01                                                                  #
# author: Magnus M?nch                                                         #
# created: 01-11-2017                                                          #
# last edited: 01-11-2017                                                      #
# references:                                                                  #
# - non-filtered: van 't Veer LJ et al. (2002), Gene expression profiling      #
#   predicts clinical outcome of breast cancer, Nature 415:530-536.            #
# - filtered: Michiels S, Koscielny S, Hill C (2005), Prediction of cancer     #
#   outcome with microarrays: a multiple random validation strategy,           #
#   Lancet 365:488-492.                                                        #
################################################################################

##################################### notes ####################################
#                                                                              #
################################################################################

### paths
path.data <- ifelse(as.character(Sys.info()[1])!="Darwin", "~/EBEN/data/",
                    "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/data/")
path.code <- ifelse(as.character(Sys.info()[1])=="Darwin", "~/EBEN/code/",
                    "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/code/")
path.res <- ifelse(as.character(Sys.info()[1])=="Darwin", "~/EBEN/results/",
                   "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/results/")
path.graph <- "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/graphs/"

### loading libraries
library(cancerdata)
library(impute)
library(glmnet)
library(pROC)

### loading data
data(VEER) 
data(VEER1) # filtered on genes version

### data manipulation
resp.breast <- VEER$class 
data1.expr <- exprs(VEER)
data2.expr <- data1.expr[
  apply(data1.expr, 1, function(x) {mean(is.na(x)) < 0.1}), ]
data3.expr <- impute.knn(data2.expr)$data # impute data using knn clustering
# DM is distant metastases within 5 years after diagnosis of breast cancer
# NODM is no metastases within 5 years
set.seed(4005)
id.train <- c(sample(which(VEER$class=="DM"), floor(sum(VEER$class=="DM")/2)), 
              sample(which(VEER$class!="DM"), floor(sum(VEER$class!="DM")/2)))
train.breast <- as.numeric(VEER$class=="DM")[id.train]
train.expr <- t(data3.expr)[id.train, ]
test.breast <- as.numeric(VEER$class=="DM")[-id.train]
test.expr <- t(data3.expr)[-id.train, ]

### standardise data
train.norm <- apply(train.expr, 2, function(x) {(x - mean(x))/sd(x)})
test.norm <- apply(test.expr, 2, function(x) {(x - mean(x))/sd(x)})

### data exploration
ntrain <- nrow(train.norm)

### fitting models
fit1.ridge <- cv.glmnet(train.norm, train.breast, family="binomial", alpha=0, 
                        nfolds=ntrain, grouped=FALSE, standardize=FALSE)
fit1.lasso <- cv.glmnet(train.norm, train.breast, family="binomial", alpha=1, 
                        nfolds=ntrain, grouped=FALSE, standardize=FALSE)
fit1.enet <- cv.glmnet(train.norm, train.breast, family="binomial", alpha=0.5, 
                       nfolds=ntrain, grouped=FALSE, standardize=FALSE)

### testing performance
methods <- c("ridge", "lasso", "enet")
pred1 <- cbind(predict(fit1.ridge, test.norm, "lambda.min", type="response"), 
               predict(fit1.lasso, test.norm, "lambda.min", type="response"),
               predict(fit1.enet, test.norm, "lambda.min", type="response"))
colnames(pred1) <- methods
auc1 <- setNames(apply(pred1, 2, function(pred) {
  pROC::roc(test.breast, pred)$auc}), methods)
psel1 <- setNames(c(fit1.ridge$nzero[fit1.ridge$lambda==fit1.ridge$lambda.min],
                    fit1.lasso$nzero[fit1.lasso$lambda==fit1.lasso$lambda.min],
                    fit1.enet$nzero[fit1.enet$lambda==fit1.enet$lambda.min]),
                  methods)
brier1 <- setNames(apply(pred1, 2, function(pred) {
  mean((test.breast - pred)^2)}), methods)
brierskill1 <- 1 - brier1/mean((test.breast - mean(test.breast))^2)

colnames(train.norm)




library("biomaRt")
fData(VEER)$symbol
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- as.character(fData(VEER)$symbol)[as.character(fData(VEER)$symbol)!=""]
test <- getBM(filters)
G_list <- getBM(filters="hgnc_symbol", attributes=c("ensembl_gene_id", "hgnc_symbol"),
                values=genes, mart=mart)


merge(df,G_list,by.x="gene",by.y="ensembl_peptide_id")
genes <- df$genes


df<-df[,-4]

as.character(fData(VEER)$symbol)[as.character(fData(VEER)$symbol)!=""]