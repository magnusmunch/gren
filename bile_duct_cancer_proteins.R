# library(rpx)
# library(mzR)
library(MSnbase)
library(glmnet)
library(pROC)
library(GRridge)

path.data <- as.character(ifelse(Sys.info()[1]=="Darwin","/Users/magnusmunch/Documents/PhD/EBEN/data/" ,
                                 "~/EBEN/data/"))

# # downloading files
# px <- PXDataset("PXD000534")
# pxurl(px)
# files <- pxfiles(px)
# 
# owd <- getwd()
# dir.create(paste(path.data, "bile_duct_cancer_proteomics/", sep=""), showWarnings = FALSE)
# setwd(paste(path.data, "bile_duct_cancer_proteomics/", sep=""))
# test1 <- pxget(px, files[33])
# setwd(owd)
# 
# test2 <- openMSfile(test1)
# test3 <- peaks(test2)
# 
# owd <- getwd()
# dir.create(paste(path.data, "bile_duct_cancer_proteomics/", sep=""), showWarnings = FALSE)
# setwd(paste(path.data, "bile_duct_cancer_proteomics/", sep=""))
# mzf <- pxget(px, files[substr(files, 1, 5)=="PRIDE"])
# setwd(owd)
# ms.duct <- openMSfile(paste(path.data, "bile_duct_cancer_proteomics/", files[33], sep=""))
# counts.duct <- peaksCount(ms.duct)
# peaks.duct <- peaks(ms.duct)

### reading in data
files.mztab <- list.files(paste(path.data, "bile_duct_cancer_proteomics/generated", sep=""), 
                          pattern="*.mztab", full.names=TRUE)
n <- length(files.mztab)

resp.duct <- numeric(n)
proteins <- numeric(0)
for(i in 1:n) {
  assign(paste("object.duct.subject", i, sep=""), readMzTabData(files.mztab[i], what="PRT"))
  assign(paste("rawdata.duct.subject", i, sep=""), fData(get(paste("object.duct.subject", i, sep=""))))
  assign(paste("prot.duct.subject", i, sep=""), get(paste("rawdata.duct.subject", i, sep=""))$accession)
  proteins <- c(proteins, get(paste("prot.duct.subject", i, sep="")))
  resp.duct[i] <- as.numeric(substr(experimentData(get(paste("object.duct.subject", i, sep="")))@
                                      other$mzTab$title, 1, 5)=="Tumor") # 0 is normal, 1 is tumor
}
proteins <- unique(proteins)
subject.ids <- sapply(paste("object.duct.subject", c(1:n), sep=""), function(x) {
  experimentData(get(x))@other$mzTab[[3]]}, USE.NAMES=FALSE)
p <- length(proteins)

### data manipulations
counts.duct <- matrix(t(sapply(paste("rawdata.duct.subject", c(1:n), sep=""), function(x) {
  counts <- get(x)$num_peptides_distinct_ms_run.1.[match(proteins, get(x)$accession)]
  replace(counts, is.na(counts), 0)})), ncol=p, dimnames=list(subject.ids, proteins))
select.duct <- counts.duct[, (colSums(counts.duct) >= 20) & (apply(counts.duct, 2, sd)!=0)]
trans.duct <- sqrt(select.duct)
norm.duct <- apply(trans.duct, 2, function(x) {(x - mean(x))/sd(x)})


parAbund <- CreatePartition(colSums(select.duct), uniform=TRUE, ngroup=5)

fit.grridge <- grridge(t(norm.duct), resp.duct, list(abundance=parAbund))

# cross-validate prediction
set.seed(2017)
n <- nrow(norm.duct)
p <- ncol(norm.duct)
nfolds <- n
rest <- n %% nfolds
foldid <- sample(rep(1:nfolds, times=round(c(rep(n %/% nfolds + as.numeric(rest!=0), times=rest),
                                             rep(n %/% nfolds, times=nfolds - rest)))))

fit.cv1 <- cv.glmnet(norm.duct, resp.duct, family="binomial", standardize=FALSE, intercept=TRUE, 
                     alpha=1, nfolds=nfolds, grouped=FALSE)
lambda1 <- fit.cv1$lambda.min
fit.cv2 <- cv.glmnet(norm.duct, resp.duct, family="binomial", standardize=FALSE, intercept=TRUE, 
                     alpha=0, nfolds=nfolds, grouped=FALSE)
lambda2 <- fit.cv2$lambda.min

methods <- c("lasso", "ridge")
pred <- matrix(NA, nrow=n, ncol=length(methods), dimnames=list(subject.ids, methods))
for(k in 1:nfolds) {
  xtrain <- norm.duct[foldid!=k, ]
  ytrain <- resp.duct[foldid!=k]
  xtest <- norm.duct[foldid==k, ]
  
  fit.lasso <- glmnet(xtrain, ytrain, family="binomial", standardize=FALSE, intercept=TRUE, 
                      alpha=1, lambda=lambda1) 
  fit.ridge <- glmnet(xtrain, ytrain, family="binomial", standardize=FALSE, intercept=TRUE, 
                      alpha=0, lambda=lambda2) 
  
  pred[foldid==k, 1] <- as.numeric(predict(fit.lasso, matrix(xtest, ncol=p), type="response"))
  pred[foldid==k, 2] <- as.numeric(predict(fit.ridge, matrix(xtest, ncol=p), type="response"))
}

auc <- sapply(1:ncol(pred), function(method) {pROC::roc(resp.duct, pred[, method])$auc})
brier <- apply(pred, 2, function(x) {mean((resp.duct - x)^2)})


proteins

pred

fit.cv1$nzero[which.min(fit.cv1$cvm)]



