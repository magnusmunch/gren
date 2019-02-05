### paths
path.code <- as.character(ifelse(
  Sys.info()[1]=="Darwin", "/Users/magnusmunch/Documents/PhD/EBEN/code/" ,"~/EBEN/code/"))
path.data <- path.data <- as.character(ifelse(
  Sys.info()[1]=="Darwin", "/Users/magnusmunch/Documents/PhD/EBEN/data/" ,"~/EBEN/data/"))

### libraries
library(pROC)
library(glmnet)
library(GRridge)

# source function for variational Bayes
source(paste(path.code, "grVBEM.R", sep=""))

### reading in the data
id.mirna <- read.table(paste(path.data, "ms_mirna/A-MEXP-1820.adf.txt", sep=""), 
                       header=TRUE, stringsAsFactors=FALSE, skip=19, sep="\t", blank.lines.skip=TRUE)
object.mirna <- read.table(paste(path.data, "ms_mirna/miRNA_profile_nonorm_nobkgd_noscaled_ArrayExpress_DataFile.txt", sep=""), 
                           header=TRUE, stringsAsFactors=FALSE, skip=7, sep="\t", blank.lines.skip=TRUE)
object.ms <- read.table(paste(path.data, "ms_mirna/E-MTAB-359.sdrf.txt", sep=""), 
                        header=TRUE, stringsAsFactors=FALSE, sep="\t", blank.lines.skip=TRUE)
mirsfamily <- read.table(paste(path.data, "miR_Family_Info_edit.txt", sep=""), header=TRUE, 
                         fill=TRUE, as.is=c(1, 2, 4, 5, 7))

### data manipulation
# checking the correlations between different measurements of the same subjects
cor.meas <- sapply(c(1:length(unique(object.ms$Source.Name))), function(s) {
  cor(object.mirna[, which((substr(colnames(object.mirna), 2, 13) %in% object.ms$Hybridization.Name[
    object.ms$Source.Name==unique(object.ms$Source.Name)[s]]) &
      (substring(colnames(object.mirna), 19)=="Signal"))])})
rankcor.meas <- sapply(c(1:length(unique(object.ms$Source.Name))), function(s) {
  cor(object.mirna[, which((substr(colnames(object.mirna), 2, 13) %in% object.ms$Hybridization.Name[
    object.ms$Source.Name==unique(object.ms$Source.Name)[s]]) &
      (substring(colnames(object.mirna), 19)=="Signal"))], method="spearman")})

# taking only miRNAs with at least one detection limit p-value < 0.01 in samples or controls
detlim <- apply(object.mirna[, substring(colnames(object.mirna), 15)=="Detection.Pval" &
                               substr(colnames(object.mirna), 2, 13) %in% 
                               object.ms$Hybridization.Name[
                                 object.ms$Factor.Value..DISEASESTATE.!="control"]], 1, function(var) {
                                   any(var < 0.01)}) | 
  apply(object.mirna[, substring(colnames(object.mirna), 15)=="Detection.Pval" &
                       substr(colnames(object.mirna), 2, 13) %in% object.ms$Hybridization.Name[
                         object.ms$Factor.Value..DISEASESTATE.=="control"]], 1, function(var) {
                           any(var < 0.01)})

# averaging if there multiple measurements for one subject
rawdata1.mirna <- object.mirna[, c(1, which(substring(colnames(object.mirna), 19)=="Signal"))]
rawdata2.mirna <- cbind(sapply(c(1:length(unique(object.ms$Source.Name))), function(s) {
  rowMeans(rawdata1.mirna[, which(substr(colnames(rawdata1.mirna), 2, 13) %in% object.ms$Hybridization.Name[
    object.ms$Source.Name==unique(object.ms$Source.Name)[s]])])}))
colnames(rawdata2.mirna) <- unique(object.ms$Source.Name)
rawdata3.mirna <- data.frame(PROBE_ID=rawdata1.mirna$PROBE_ID, rawdata2.mirna)

# creating the data table with the miRNA names
data.mirna <- t(rawdata3.mirna[, -1])
colnames(data.mirna) <- id.mirna$Reporter.Database.Entry.mirbase.[
  match(rawdata3.mirna$PROBE_ID, id.mirna$Reporter.Name)]
trans.mirna <- sqrt(data.mirna)
norm.mirna <- apply(trans.mirna, 2, function(var) {(var - mean(var))/sd(var)})

# response data
data.ms <- object.ms[!duplicated(object.ms$Source.Name), ]
resp.ms <- as.numeric(data.ms$Factor.Value..DISEASESTATE.[
  match(rownames(norm.mirna), data.ms$Source.Name)]!= "control")

### co-data
# conservation status of miRNAs, 0: not conserved, 1: conserved accross vertebrates, 2: conserved
mirsfamily$Family.Conservation.2 <- (mirsfamily$Family.Conservation. > 0)*
  mirsfamily$Family.Conservation.
conservation1 <- mirsfamily$Family.Conservation.2[match(colnames(norm.mirna), mirsfamily$MiRBase.ID)]
conservation2 <- mirsfamily$Family.Conservation.2[
  match(gsub("T", "U", id.mirna$Reporter.Sequence[
    match(colnames(norm.mirna), id.mirna$Reporter.Database.Entry.mirbase.)]),
    mirsfamily$Mature.sequence)]
conservation <- ifelse(is.na(conservation1), conservation2, conservation1) + 1
abundance <- colSums(data.mirna)
cf <- colMeans(data.mirna)/apply(data.mirna, 2, sd)
sds <- apply(data.mirna, 2, sd)

partitions1 <- list(conservation=CreatePartition(as.factor(ifelse(is.na(conservation), 4, 
                                                                  conservation))))
partitions2 <- list(conservation=ifelse(is.na(conservation), 4, conservation))
partitions3 <- list(conservation=CreatePartition(as.factor(ifelse(is.na(conservation), 4, 
                                                                 conservation))),
                    abundance=CreatePartition(abundance, ngroup=5, uniform=TRUE))
partitions4 <- list(conservation=ifelse(is.na(conservation), 4, conservation),
                    abundance=rep(1:length(partitions3$abundance), 
                                  times=lapply(partitions3$abundance, length))[
                                    order(unlist(partitions3$abundance))])
partitions5 <- list(conservation=CreatePartition(as.factor(ifelse(is.na(conservation), 4, 
                                                                  conservation))),
                    cf=CreatePartition(cf, ngroup=5, uniform=TRUE))
partitions6 <- list(conservation=ifelse(is.na(conservation), 4, conservation),
                    cf=rep(1:length(partitions5$cf), times=lapply(partitions5$cf, length))[
                             order(unlist(partitions5$cf))])
partitions7 <- list(conservation=CreatePartition(as.factor(ifelse(is.na(conservation), 4, 
                                                                  conservation))),
                    sds=CreatePartition(sds, ngroup=5, uniform=TRUE))
partitions8 <- list(conservation=ifelse(is.na(conservation), 4, conservation),
                    sds=rep(1:length(partitions7$sds), times=lapply(partitions7$sds, length))[
                      order(unlist(partitions7$sds))])

### with sqrt scaled data and conservation status
# fit models
enet.pen <- cv.pen(norm.mirna, resp.ms, unpenalized=NULL, intercept=TRUE, psel=NULL)
fit.ridge <- cv.glmnet(norm.mirna, resp.ms, family="binomial", alpha=0, grouped=FALSE,
                       standardize=FALSE, nfolds=length(resp.ms))
fit.lasso <- cv.glmnet(norm.mirna, resp.ms, family="binomial", alpha=1, grouped=FALSE,
                       standardize=FALSE, nfolds=length(resp.ms))
fit.enet <- glmnet(norm.mirna, resp.ms, family="binomial", standardize=FALSE,
                   alpha=enet.pen$alpha[which.min(enet.pen$cvll)],
                   lambda=enet.pen$lambda[which.min(enet.pen$cvll)])
fit.grridge <- grridge(t(norm.mirna), resp.ms, partitions1)
fit.greben <- grEBEN(norm.mirna, resp.ms, rep(1, length(resp.ms)), partitions=partitions2,
                     lambda1=enet.pen$lambda1bayes[which.min(enet.pen$cvll)],
                     lambda2=enet.pen$lambda2bayes[which.min(enet.pen$cvll)])

# cross-validate models
set.seed(1001)
p <- ncol(norm.mirna)
n <- nrow(norm.mirna)

nfolds <- n
rest <- n %% nfolds
foldid <- sample(rep(1:nfolds, times=round(c(rep(n %/% nfolds + as.numeric(rest!=0), times=rest),
                                             rep(n %/% nfolds, times=nfolds - rest)))))

methods <- c("ridge", "lasso", "enet", "GRridge", "grEBEN+enet", "grEBEN+VB")
pred1 <- matrix(NA, ncol=length(methods), nrow=n, dimnames=list(NULL, methods))
for(k in 1:nfolds) {
  cat(paste("Fold ", k, "\n"))
  xtrain <- norm.mirna[foldid!=k, ]
  xtest <- norm.mirna[foldid==k, ]
  ytrain <- resp.ms[foldid!=k]
  ytest <- resp.ms[foldid==k]
  ntrain <- length(ytrain)
  ntest <- length(ytest)
  
  enet.pen <- cv.pen(xtrain, ytrain, unpenalized=NULL, intercept=TRUE, psel=NULL)
  fitcv.ridge <- cv.glmnet(xtrain, ytrain, family="binomial", alpha=0, grouped=FALSE, 
                           standardize=FALSE, nfolds=ntrain)
  fitcv.lasso <- cv.glmnet(xtrain, ytrain, family="binomial", alpha=1, grouped=FALSE,
                           standardize=FALSE, nfolds=ntrain)
  fitcv.enet <- glmnet(xtrain, ytrain, family="binomial", standardize=FALSE,
                       alpha=enet.pen$alpha[which.min(enet.pen$cvll)],
                       lambda=enet.pen$lambda[which.min(enet.pen$cvll)])
  fitcv.grridge <- grridge(t(xtrain), ytrain, partitions1)
  fitcv.greben <- grEBEN(xtrain, ytrain, rep(1, length(ytrain)), partitions=partitions2,
                         lambda1=enet.pen$lambda1bayes[which.min(enet.pen$cvll)],
                         lambda2=enet.pen$lambda2bayes[which.min(enet.pen$cvll)])
  
  pred1[foldid==k, 1] <- predict(fitcv.ridge, matrix(xtest, ncol=p), s="lambda.min")
  pred1[foldid==k, 2] <- predict(fitcv.lasso, matrix(xtest, ncol=p), s="lambda.min")
  pred1[foldid==k, 3] <- predict(fitcv.enet, matrix(xtest, ncol=p))
  pred1[foldid==k, 4] <- predict.grridge(fitcv.grridge, t(matrix(xtest, ncol=p)))[, 2]
  pred1[foldid==k, 5] <- predict.grEBEN(fitcv.greben, matrix(xtest, ncol=p), type="penalized")
  pred1[foldid==k, 6] <- predict.grEBEN(fitcv.greben, matrix(xtest, ncol=p), type="VB")
  
}

auc1 <- setNames(sapply(1:ncol(pred1), function(method) {
  return(pROC::roc(resp.ms, pred1[, method])$auc)}), methods)
brier1 <- setNames(apply(pred1, 2, function(cpred) {mean((resp.ms - cpred)^2)}), methods)
brierskill1 <- 1 - brier1/mean((resp.ms - mean(resp.ms))^2)
nzero <- setNames(c(p, fit.lasso$nzero[fit.lasso$lambda==fit.lasso$lambda.min], fit.enet$df,
                     p, sum(fit.greben$beta[-1]!=0), p), methods)

tab1 <- cbind(auc=round(auc1, 2), brier=round(brierskill1, 2), psel=nzero)
results1 <- list(auc=auc1, brier=brier1, brierskill=brierskill1, nzero=nzero, table=tab1)
save(results1, file=paste(path.res, "grEBEN_mirna_ms_res1.Rdata", sep=""))





### without transforming the data
trans.mirna <- data.mirna
norm.mirna <- apply(trans.mirna, 2, function(var) {(var - mean(var))/sd(var)})

# fitting models
enet2.pen <- cv.pen(norm.mirna, resp.ms, unpenalized=NULL, intercept=TRUE, psel=NULL)
fit2.ridge <- cv.glmnet(norm.mirna, resp.ms, family="binomial", alpha=0, grouped=FALSE,
                       standardize=FALSE, nfolds=length(resp.ms))
fit2.lasso <- cv.glmnet(norm.mirna, resp.ms, family="binomial", alpha=1, grouped=FALSE,
                       standardize=FALSE, nfolds=length(resp.ms))
fit2.enet <- glmnet(norm.mirna, resp.ms, family="binomial", standardize=FALSE,
                   alpha=enet.pen$alpha[which.min(enet.pen$cvll)],
                   lambda=enet.pen$lambda[which.min(enet.pen$cvll)])
fit2.grridge <- grridge(t(norm.mirna), resp.ms, partitions1)
fit2.greben <- grEBEN(norm.mirna, resp.ms, rep(1, length(resp.ms)), partitions=partitions2,
                     lambda1=enet.pen$lambda1bayes[which.min(enet.pen$cvll)],
                     lambda2=enet.pen$lambda2bayes[which.min(enet.pen$cvll)])

# cross-validate
set.seed(1001)
p <- ncol(norm.mirna)
n <- nrow(norm.mirna)

nfolds <- n
rest <- n %% nfolds
foldid <- sample(rep(1:nfolds, times=round(c(rep(n %/% nfolds + as.numeric(rest!=0), times=rest),
                                             rep(n %/% nfolds, times=nfolds - rest)))))

methods <- c("ridge", "lasso", "enet", "GRridge", "grEBEN+enet", "grEBEN+VB")
pred2 <- matrix(NA, ncol=length(methods), nrow=n, dimnames=list(NULL, methods))
for(k in 1:nfolds) {
  cat(paste("Fold ", k, "\n"))
  xtrain <- norm.mirna[foldid!=k, ]
  xtest <- norm.mirna[foldid==k, ]
  ytrain <- resp.ms[foldid!=k]
  ytest <- resp.ms[foldid==k]
  ntrain <- length(ytrain)
  ntest <- length(ytest)
  
  enet.pen <- cv.pen(xtrain, ytrain, unpenalized=NULL, intercept=TRUE, psel=NULL)
  fitcv.ridge <- cv.glmnet(xtrain, ytrain, family="binomial", alpha=0, grouped=FALSE, 
                           standardize=FALSE, nfolds=ntrain)
  fitcv.lasso <- cv.glmnet(xtrain, ytrain, family="binomial", alpha=1, grouped=FALSE,
                           standardize=FALSE, nfolds=ntrain)
  fitcv.enet <- glmnet(xtrain, ytrain, family="binomial", standardize=FALSE,
                       alpha=enet.pen$alpha[which.min(enet.pen$cvll)],
                       lambda=enet.pen$lambda[which.min(enet.pen$cvll)])
  fitcv.grridge <- grridge(t(xtrain), ytrain, partitions1)
  fitcv.greben <- grEBEN(xtrain, ytrain, rep(1, length(ytrain)), partitions=partitions2,
                         lambda1=enet.pen$lambda1bayes[which.min(enet.pen$cvll)],
                         lambda2=enet.pen$lambda2bayes[which.min(enet.pen$cvll)])
  
  pred2[foldid==k, 1] <- predict(fitcv.ridge, matrix(xtest, ncol=p), s="lambda.min")
  pred2[foldid==k, 2] <- predict(fitcv.lasso, matrix(xtest, ncol=p), s="lambda.min")
  pred2[foldid==k, 3] <- predict(fitcv.enet, matrix(xtest, ncol=p))
  pred2[foldid==k, 4] <- predict.grridge(fitcv.grridge, t(matrix(xtest, ncol=p)))[, 2]
  pred2[foldid==k, 5] <- predict.grEBEN(fitcv.greben, matrix(xtest, ncol=p), type="penalized")
  pred2[foldid==k, 6] <- predict.grEBEN(fitcv.greben, matrix(xtest, ncol=p), type="VB")
  
}

auc2 <- setNames(sapply(1:ncol(pred2), function(method) {
  return(pROC::roc(resp.ms, pred2[, method])$auc)}), methods)
brier2 <- setNames(apply(pred2, 2, function(cpred) {mean((resp.ms - cpred)^2)}), methods)
brierskill2 <- 1 - brier2/mean((resp.ms - mean(resp.ms))^2)
nzero <- setNames(c(p, fit.lasso$nzero[fit.lasso$lambda==fit.lasso$lambda.min], fit.enet$df,
                    p, sum(fit.greben$beta[-1]!=0), p), methods)

tab2 <- cbind(auc=round(auc2, 2), brier=round(brierskill2, 2), psel=nzero)
results2 <- list(auc=auc2, brier=brier2, brierskill=brierskill2, nzero=nzero, table=tab2)
save(results2, file=paste(path.res, "grEBEN_mirna_ms_res2.Rdata", sep=""))




### log-transforming the data
trans.mirna <- log(data.mirna)
norm.mirna <- apply(trans.mirna, 2, function(var) {(var - mean(var))/sd(var)})

# fitting models
enet3.pen <- cv.pen(norm.mirna, resp.ms, unpenalized=NULL, intercept=TRUE, psel=NULL)
fit3.ridge <- cv.glmnet(norm.mirna, resp.ms, family="binomial", alpha=0, grouped=FALSE,
                        standardize=FALSE, nfolds=length(resp.ms))
fit3.lasso <- cv.glmnet(norm.mirna, resp.ms, family="binomial", alpha=1, grouped=FALSE,
                        standardize=FALSE, nfolds=length(resp.ms))
fit3.enet <- glmnet(norm.mirna, resp.ms, family="binomial", standardize=FALSE,
                    alpha=enet.pen$alpha[which.min(enet.pen$cvll)],
                    lambda=enet.pen$lambda[which.min(enet.pen$cvll)])
fit3.grridge <- grridge(t(norm.mirna), resp.ms, partitions1)
fit3.greben <- grEBEN(norm.mirna, resp.ms, rep(1, length(resp.ms)), partitions=partitions2,
                      lambda1=enet.pen$lambda1bayes[which.min(enet.pen$cvll)],
                      lambda2=enet.pen$lambda2bayes[which.min(enet.pen$cvll)])

# cross-validate
set.seed(1001)
p <- ncol(norm.mirna)
n <- nrow(norm.mirna)

nfolds <- n
rest <- n %% nfolds
foldid <- sample(rep(1:nfolds, times=round(c(rep(n %/% nfolds + as.numeric(rest!=0), times=rest),
                                             rep(n %/% nfolds, times=nfolds - rest)))))

methods <- c("ridge", "lasso", "enet", "GRridge", "grEBEN+enet", "grEBEN+VB")
pred3 <- matrix(NA, ncol=length(methods), nrow=n, dimnames=list(NULL, methods))
for(k in 1:nfolds) {
  cat(paste("Fold ", k, "\n"))
  xtrain <- norm.mirna[foldid!=k, ]
  xtest <- norm.mirna[foldid==k, ]
  ytrain <- resp.ms[foldid!=k]
  ytest <- resp.ms[foldid==k]
  ntrain <- length(ytrain)
  ntest <- length(ytest)
  
  enet.pen <- cv.pen(xtrain, ytrain, unpenalized=NULL, intercept=TRUE, psel=NULL)
  fitcv.ridge <- cv.glmnet(xtrain, ytrain, family="binomial", alpha=0, grouped=FALSE, 
                           standardize=FALSE, nfolds=ntrain)
  fitcv.lasso <- cv.glmnet(xtrain, ytrain, family="binomial", alpha=1, grouped=FALSE,
                           standardize=FALSE, nfolds=ntrain)
  fitcv.enet <- glmnet(xtrain, ytrain, family="binomial", standardize=FALSE,
                       alpha=enet.pen$alpha[which.min(enet.pen$cvll)],
                       lambda=enet.pen$lambda[which.min(enet.pen$cvll)])
  fitcv.grridge <- grridge(t(xtrain), ytrain, partitions1)
  fitcv.greben <- grEBEN(xtrain, ytrain, rep(1, length(ytrain)), partitions=partitions2,
                         lambda1=enet.pen$lambda1bayes[which.min(enet.pen$cvll)],
                         lambda2=enet.pen$lambda2bayes[which.min(enet.pen$cvll)])
  
  pred3[foldid==k, 1] <- predict(fitcv.ridge, matrix(xtest, ncol=p), s="lambda.min")
  pred3[foldid==k, 2] <- predict(fitcv.lasso, matrix(xtest, ncol=p), s="lambda.min")
  pred3[foldid==k, 3] <- predict(fitcv.enet, matrix(xtest, ncol=p))
  pred3[foldid==k, 4] <- predict.grridge(fitcv.grridge, t(matrix(xtest, ncol=p)))[, 2]
  pred3[foldid==k, 5] <- predict.grEBEN(fitcv.greben, matrix(xtest, ncol=p), type="penalized")
  pred3[foldid==k, 6] <- predict.grEBEN(fitcv.greben, matrix(xtest, ncol=p), type="VB")
  
}

auc3 <- setNames(sapply(1:ncol(pred3), function(method) {
  return(pROC::roc(resp.ms, pred3[, method])$auc)}), methods)
brier3 <- setNames(apply(pred3, 2, function(cpred) {mean((resp.ms - cpred)^2)}), methods)
brierskill3 <- 1 - brier3/mean((resp.ms - mean(resp.ms))^2)
nzero <- setNames(c(p, fit.lasso$nzero[fit.lasso$lambda==fit.lasso$lambda.min], fit.enet$df,
                    p, sum(fit.greben$beta[-1]!=0), p), methods)

tab3 <- cbind(auc=round(auc3, 2), brier=round(brierskill3, 2), psel=nzero)
results3 <- list(auc=auc3, brier=brier3, brierskill=brierskill3, nzero=nzero, table=tab3)
save(results3, file=paste(path.res, "grEBEN_mirna_ms_res3.Rdata", sep=""))




### including abundances and conservation status as co-data
# fitting the models
fit4.grridge <- grridge(t(norm.mirna), resp.ms, partitions3)
fit4.greben <- grEBEN(norm.mirna, resp.ms, rep(1, length(resp.ms)), partitions=partitions4,
                      lambda1=enet.pen$lambda1bayes[which.min(enet.pen$cvll)],
                      lambda2=enet.pen$lambda2bayes[which.min(enet.pen$cvll)])

set.seed(1001)
p <- ncol(norm.mirna)
n <- nrow(norm.mirna)

nfolds <- n
rest <- n %% nfolds
foldid <- sample(rep(1:nfolds, times=round(c(rep(n %/% nfolds + as.numeric(rest!=0), times=rest),
                                             rep(n %/% nfolds, times=nfolds - rest)))))

methods <- c("ridge", "lasso", "enet", "GRridge", "grEBEN+enet", "grEBEN+VB")
pred4 <- matrix(NA, ncol=length(methods), nrow=n, dimnames=list(NULL, methods))
for(k in 1:nfolds) {
  cat(paste("Fold ", k, "\n"))
  xtrain <- norm.mirna[foldid!=k, ]
  xtest <- norm.mirna[foldid==k, ]
  ytrain <- resp.ms[foldid!=k]
  ytest <- resp.ms[foldid==k]
  ntrain <- length(ytrain)
  ntest <- length(ytest)
  
  enet.pen <- cv.pen(xtrain, ytrain, unpenalized=NULL, intercept=TRUE, psel=NULL)
  fitcv.ridge <- cv.glmnet(xtrain, ytrain, family="binomial", alpha=0, grouped=FALSE, 
                           standardize=FALSE, nfolds=ntrain)
  fitcv.lasso <- cv.glmnet(xtrain, ytrain, family="binomial", alpha=1, grouped=FALSE,
                           standardize=FALSE, nfolds=ntrain)
  fitcv.enet <- glmnet(xtrain, ytrain, family="binomial", standardize=FALSE,
                       alpha=enet.pen$alpha[which.min(enet.pen$cvll)],
                       lambda=enet.pen$lambda[which.min(enet.pen$cvll)])
  fitcv.grridge <- grridge(t(xtrain), ytrain, partitions3)
  fitcv.greben <- grEBEN(xtrain, ytrain, rep(1, length(ytrain)), partitions=partitions4,
                         lambda1=enet.pen$lambda1bayes[which.min(enet.pen$cvll)],
                         lambda2=enet.pen$lambda2bayes[which.min(enet.pen$cvll)])
  
  pred4[foldid==k, 1] <- predict(fitcv.ridge, matrix(xtest, ncol=p), s="lambda.min")
  pred4[foldid==k, 2] <- predict(fitcv.lasso, matrix(xtest, ncol=p), s="lambda.min")
  pred4[foldid==k, 3] <- predict(fitcv.enet, matrix(xtest, ncol=p))
  pred4[foldid==k, 4] <- predict.grridge(fitcv.grridge, t(matrix(xtest, ncol=p)))[, 2]
  pred4[foldid==k, 5] <- predict.grEBEN(fitcv.greben, matrix(xtest, ncol=p), type="penalized")
  pred4[foldid==k, 6] <- predict.grEBEN(fitcv.greben, matrix(xtest, ncol=p), type="VB")
  
}

auc4 <- setNames(sapply(1:ncol(pred4), function(method) {
  return(pROC::roc(resp.ms, pred4[, method])$auc)}), methods)
brier4 <- setNames(apply(pred4, 2, function(cpred) {mean((resp.ms - cpred)^2)}), methods)
brierskill4 <- 1 - brier4/mean((resp.ms - mean(resp.ms))^2)
nzero <- setNames(c(p, fit.lasso$nzero[fit.lasso$lambda==fit.lasso$lambda.min], fit.enet$df,
                    p, sum(fit2.greben$beta[-1]!=0), p), methods)

tab4 <- cbind(auc=round(auc4, 2), brier=round(brierskill4, 2), psel=nzero)
results4 <- list(auc=auc4, brier=brier4, brierskill=brierskill4, nzero=nzero, table=tab4)
save(results4, file=paste(path.res, "grEBEN_mirna_ms_res4.Rdata", sep=""))



### including abundances and conservation status as monotone co-data
# fitting the models
fit5.grridge <- grridge(t(norm.mirna), resp.ms, partitions3, monotone=c(FALSE, TRUE))
fit5.greben <- grEBEN(norm.mirna, resp.ms, rep(1, length(resp.ms)), partitions=partitions4,
                      lambda1=enet.pen$lambda1bayes[which.min(enet.pen$cvll)],
                      lambda2=enet.pen$lambda2bayes[which.min(enet.pen$cvll)],
                      monotone=list(monotone=c(FALSE, TRUE), decreasing=c(FALSE, FALSE)))




### including conservatino status and coefficient of variation
fit6.grridge <- grridge(t(norm.mirna), resp.ms, partitions5, monotone=c(FALSE, TRUE))
fit6.greben <- grEBEN(norm.mirna, resp.ms, rep(1, length(resp.ms)), partitions=partitions6,
                      lambda1=enet.pen$lambda1bayes[which.min(enet.pen$cvll)],
                      lambda2=enet.pen$lambda2bayes[which.min(enet.pen$cvll)],
                      monotone=list(monotone=c(FALSE, TRUE), decreasing=c(FALSE, FALSE)))



### including conservatino status and standard deviation
fit6.grridge <- grridge(t(norm.mirna), resp.ms, partitions7, monotone=c(FALSE, FALSE))
fit6.greben <- grEBEN(norm.mirna, resp.ms, rep(1, length(resp.ms)), partitions=partitions8,
                      lambda1=enet.pen$lambda1bayes[which.min(enet.pen$cvll)],
                      lambda2=enet.pen$lambda2bayes[which.min(enet.pen$cvll)],
                      monotone=list(monotone=c(FALSE, FALSE), decreasing=c(FALSE, FALSE)))



### graphs
tcor <- cor(norm.mirna)

hist(tcor[lower.tri(tcor)], xlab="Correlation", freq=FALSE, 
     main="Data correlation")

quartz()
par(mfrow=c(length(unique(partitions2$conservation)), length(unique(partitions2$conservation))))
count <- 0
for(i in c(1:length(unique(partitions2$conservation)))) {
  for(j in c(1:length(unique(partitions2$conservation)))) {
    if(i < j) {
      plot.new()
    } else {
      count <- count + 1
      hist(tcor[partitions2$conservation==i, partitions2$conservation==j], freq=FALSE,
           xlab="Correlation", main=paste(letters[count], ")", sep=""))
    }
  }
}
  


plot(fit2.grridge$lambdamults$abundance, fit2.greben$lambdag$abundance[, fit2.greben$nouteriter + 1])
plot(fit2.grridge$lambdamults$conservation, fit2.greben$lambdag$conservation[, fit2.greben$nouteriter + 1])



barplot(rbind(fit2.grridge$lambdamults$abundance, 
              fit2.greben$lambdag$abundance[, fit2.greben$nouteriter + 1]),
        beside=TRUE)




