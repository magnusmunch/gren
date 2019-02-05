### libraries
library(TCGAbiolinks)
library(glmnet)
library(pROC)

### paths
path.data <- as.character(ifelse(Sys.info()[1]=="Darwin","/Users/magnusmunch/Documents/PhD/EBEN/data/" ,
                                 "~/EBEN/data/"))

### TCGA projects
projects <- getGDCprojects()

### miRNA and liver cancer data
# miRNAs
query.mirna <- GDCquery(project=projects$project_id[1], legacy=TRUE,
                        data.category="Gene expression", data.type="miRNA gene quantification", 
                        sample.type=c("Primary solid Tumor", "Solid Tissue Normal"), file.type="mirna")
rle(sort(as.character(query.mirna$results[[1]]$tissue.definition))) # check the number of samples
GDCdownload(query.mirna, directory=path.data)
rawdata.mirna <- GDCprepare(query.mirna, directory=path.data)
counts.mirna <- matrix(t(rawdata.mirna[, substr(colnames(rawdata.mirna), 1, 10)=="read_count"]), 
                       ncol=nrow(rawdata.mirna), 
                       dimnames=list(substr(colnames(rawdata.mirna)[substr(colnames(rawdata.mirna), 1, 10)=="read_count"], 12, 1000000L),
                                     rawdata.mirna$miRNA_ID))
select.mirna <- counts.mirna[, colSums(counts.mirna) >= 100]
trans.mirna <- log(select.mirna + 1)
norm.mirna <- apply(trans.mirna, 2, function(x) {(x - mean(x))/sd(x)})

# response: 1 is Primary solid tumor, 0 is solid tissue normal
resp.liver <- as.numeric(query.mirna$results[[1]]$tissue.definition
                         [match(query.mirna$results[[1]]$cases, rownames(counts.mirna))]==
                           "Primary solid Tumor")

# cross-validate prediction
set.seed(2017)
n <- nrow(norm.mirna)
p <- ncol(norm.mirna)
nfolds <- n
rest <- n %% nfolds
foldid <- sample(rep(1:nfolds, times=round(c(rep(n %/% nfolds + as.numeric(rest!=0), times=rest),
                                             rep(n %/% nfolds, times=nfolds - rest)))))

lambda1 <- cv.glmnet(norm.mirna, resp.liver, family="binomial", standardize=FALSE, intercept=TRUE, 
                     alpha=1, nfolds=nfolds)$lambda.min # 0.002317547
lambda2 <- cv.glmnet(norm.mirna, resp.liver, family="binomial", standardize=FALSE, intercept=TRUE, 
                     alpha=0, nfolds=nfolds)$lambda.min # 2.317547

methods <- c("lasso", "ridge")
pred <- matrix(NA, nrow=n, ncol=length(methods), dimnames=list(NULL, methods))
for(k in 1:nfolds) {
  xtrain <- norm.mirna[foldid!=k, ]
  ytrain <- resp.liver[foldid!=k]
  xtest <- norm.mirna[foldid==k, ]
  
  fit.lasso <- glmnet(xtrain, ytrain, family="binomial", standardize=FALSE, intercept=TRUE, 
                      alpha=1, lambda=lambda1) 
  fit.ridge <- glmnet(xtrain, ytrain, family="binomial", standardize=FALSE, intercept=TRUE, 
                      alpha=0, lambda=lambda2) 
  
  pred[foldid==k, 1] <- as.numeric(predict(fit.lasso, matrix(xtest, ncol=p), type="response"))
  pred[foldid==k, 2] <- as.numeric(predict(fit.ridge, matrix(xtest, ncol=p), type="response"))
}

auc <- sapply(1:ncol(pred), function(method) {pROC::roc(resp.liver, pred[, method])$auc})
brier <- c(mean((resp.liver - pred[, 1])^2), mean((resp.liver - pred[, 2])^2))

### miRNA and liver cancer data
# miRNAs
query.met <- GDCquery(project=projects$project_id[34], data.category="DNA Methylation",
                      legacy=FALSE, sample.type=c("Primary solid Tumor"))
query.clin <- GDCquery_clinic(project=projects$project_id[34], type="clinical")
View(query.clin)
rle(sort(as.character(query.clin$progression_or_recurrence))) # check the number of samples


GDCdownload(query.mirna, directory=path.data)
rawdata.mirna <- GDCprepare(query.mirna, directory=path.data)

