################################### preamble ###################################
# Golub dataset (from original enet paper)                                     #
# version: 01                                                                  #
# author: Magnus M?nch                                                         #
# created: 31-10-2017                                                          #
# last edited: 31-10-2017                                                      #
################################################################################

##################################### notes ####################################

### paths
path.data <- ifelse(as.character(Sys.info()[1])!="Darwin", "~/EBEN/data/leukaemia_gene-expression_golub/",
                    "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/data/leukaemia_gene-expression_golub/")
path.code <- ifelse(as.character(Sys.info()[1])=="Darwin", "~/EBEN/code/",
                    "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/code/")
path.res <- ifelse(as.character(Sys.info()[1])=="Darwin", "~/EBEN/results/",
                   "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/results/")
path.graph <- "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/graphs/"

### libraries
library(golubEsets)
library(glmnet)
library(pROC)

### reading in data
data(golubTrain)
data(golubTest)

### data manipulation
train.data <- t(golubTrain@exprs)
test.data <- t(golubTest@exprs)
train.leuk <- as.numeric(golubTrain@phenoData@pData$ALL.AML=="ALL")
test.leuk <- as.numeric(golubTest@phenoData@pData$ALL.AML=="ALL")

### standardise data
train.norm <- apply(train.data, 2, function(x) {(x - mean(x))/sd(x)})
test.norm <- apply(test.data, 2, function(x) {(x - mean(x))/sd(x)})

### data exploration
ntrain <- nrow(train.norm)
ptrain <- ncol(train.norm)
prop.ALL.train <- mean(train.leuk)

### fitting models
fit1.ridge <- cv.glmnet(train.norm, train.leuk, family="binomial", alpha=0, 
                        nfolds=ntrain, grouped=FALSE, standardize=FALSE)
fit1.lasso <- cv.glmnet(train.norm, train.leuk, family="binomial", alpha=1, 
                        nfolds=ntrain, grouped=FALSE, standardize=FALSE, 
                        pmax=100)
fit1.enet <- cv.glmnet(train.norm, train.leuk, family="binomial", alpha=0.5, 
                       nfolds=ntrain, grouped=FALSE, standardize=FALSE, 
                       pmax=100)

### testing performance
methods <- c("ridge", "lasso", "enet")
pred1 <- cbind(predict(fit1.ridge, test.norm, "lambda.min", type="response"), 
               predict(fit1.lasso, test.norm, "lambda.min", type="response"),
               predict(fit1.enet, test.norm, "lambda.min", type="response"))
colnames(pred1) <- methods
auc1 <- setNames(apply(pred1, 2, function(pred) {
  pROC::roc(test.leuk, pred)$auc}), methods)
psel1 <- setNames(c(fit1.ridge$nzero[fit1.ridge$lambda==fit1.ridge$lambda.min],
                    fit1.lasso$nzero[fit1.lasso$lambda==fit1.lasso$lambda.min],
                    fit1.enet$nzero[fit1.enet$lambda==fit1.enet$lambda.min]),
                  methods)
brier1 <- setNames(apply(pred1, 2, function(pred) {
  mean((test.leuk - pred)^2)}), methods)
brierskill1 <- 1 - brier1/mean((test.leuk - mean(test.leuk))^2)

### savind results
results1 <- list(auc=auc1, brier=brier1, brierskill=brierskill1, psel=psel1,
                 fitted=list(fit1.ridge, fit1.lasso, fit1.enet))
save(results1, file=paste(path.res, "grEBEN_gene_expr_leukaemia_res1.Rdata", sep=""))








library(colonCA)
data.expr <- exprs(colonCA)
resp.colon <- colonCA$class # n is normal, t is colon tumour