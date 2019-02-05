################################### preamble ###################################
# lung cancer data from Putri                                                  #
# version: 01                                                                  #
# author: Magnus M?nch                                                         #
# created: 27-11-2017                                                          #
# last edited: 27-11-2017                                                      #
################################################################################

##################################### notes ####################################
#                                                                              #
################################################################################

### paths
path.data <- ifelse(as.character(Sys.info()[1])!="Darwin", "~/EBEN/data/lung_cancer_methylation",
                    "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/data/lung_cancer_methylation/")
path.code <- ifelse(as.character(Sys.info()[1])!="Darwin", "~/EBEN/code/",
                    "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/code/")
path.res <- ifelse(as.character(Sys.info()[1])!="Darwin", "~/EBEN/results/",
                   "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/results/")
path.graph <- "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/graphs/"

### loading libraries
library(GRridge)
library(glmnet)
library(pROC)
library(limma)

### loading function
source(paste(path.code, "grVBEM.R", sep=""))

### loading data
load(paste(path.data, "datSputum.Rdata", sep=""))

### data preparation
resp <- as.numeric(pd$Sample_Group=="LungCancer")
rownames(pd) <- pd$Sample_Name
rownames(AnnoFile) <- AnnoFile$Name
eset <- ExpressionSet(betaPreProc, as(pd, "AnnotatedDataFrame"), 
                      as(AnnoFile, "AnnotatedDataFrame"))
design <- matrix(as.numeric(c(pData(eset)$Sample_Group=="control", 
                              pData(eset)$Sample_Group=="LungCancer")), 
                 ncol=2, dimnames=list(rownames(pData(eset)), 
                                       c("control", "LungCancer")))
contrast <- makeContrasts(LungCancer - control, levels=design)

### data exploration
dim(betaPreProc)

# calculate t- and Wilcoxon-test p-values
pval.t <- apply(betaPreProc, 1, function(dat) {
  t.test(dat[resp==1], dat[resp==0])$p.value})

pval.w <- apply(betaPreProc, 1, function(dat) {
  wilcox.test(dat[resp==1], dat[resp==0])$p.value})

fit.limma <- lmFit(eset, design)
fit.limma <- contrasts.fit(fit.limma, contrast)
fit.limma <- eBayes(fit.limma)
tab.limma <- toptable(fit.limma, coef=1, number=nrow(fit.limma))
pval.limma <- tab.limma$P.Value
logfc <- tab.limma$logFC

# data selections
smeth <- betaPreProc[tab.limma$logFC > 0, ]
smeth.norm <- apply(betaPreProc, 1, function(x) {(x - mean(x))/sd(x)})

# fitting models
fit.lasso <- cv.glmnet(t(betaPreProc), resp, family="binomial", 
                       standardize=FALSE, alpha=1)
sfit.lasso <- cv.glmnet(smeth.norm, resp, family="binomial",
                        standardize=FALSE, alpha=1)
fit.ridge <- cv.glmnet(t(betaPreProc), resp, family="binomial",
                       standardize=FALSE, alpha=0)
sfit.ridge <- cv.glmnet(smeth.norm, resp, family="binomial",
                        standardize=FALSE, alpha=0)
fit.enet1 <- cv.glmnet(t(betaPreProc), resp, family="binomial",
                       standardize=FALSE, alpha=0.05)
fit.enet2 <- cv.glmnet(t(betaPreProc), resp, family="binomial",
                       standardize=FALSE, alpha=0.5)
fit.enet3 <- cv.glmnet(t(betaPreProc), resp, family="binomial",
                       standardize=FALSE, alpha=0.95)
sfit.enet1 <- cv.glmnet(smeth.norm, resp, family="binomial",
                        standardize=FALSE, alpha=0.05)
sfit.enet2 <- cv.glmnet(smeth.norm, resp, family="binomial",
                        standardize=FALSE, alpha=0.5)
sfit.enet3 <- cv.glmnet(smeth.norm, resp, family="binomial",
                        standardize=FALSE, alpha=0.95)
est.lasso <- as.numeric(coef(fit.lasso, s="lambda.min"))[-1]
est.ridge <- as.numeric(coef(fit.ridge, s=min(fit.lasso$lambda)))[-1]
sest.lasso <- as.numeric(coef(sfit.lasso, s="lambda.min"))[-1]
est.enet1 <- as.numeric(coef(fit.enet1, s="lambda.min"))[-1]
est.enet2 <- as.numeric(coef(fit.enet2, s="lambda.min"))[-1]
est.enet3 <- as.numeric(coef(fit.enet3, s="lambda.min"))[-1]

# cross-validating performance
pred.lasso <- pred.ridge <- pred.enet1 <- pred.enet2 <- pred.enet3 <-
  spred.lasso <- spred.ridge <- spred.enet1 <- spred.enet2 <- spred.enet3 <-
  numeric(length(resp))
for(i in 1:length(resp)) {
  xtrain <- t(betaPreProc)[-i, ]
  xtest <- t(betaPreProc)[i, ]
  sxtrain <- smeth.norm[-i, ]
  sxtest <- smeth.norm[i, ]
  ytrain <- resp[-i]
  cv.ridge <- glmnet(xtrain, ytrain, family="binomial", standardize=FALSE,
                     alpha=0, lambda=fit.ridge$lambda.min)
  cv.lasso <- glmnet(xtrain, ytrain, family="binomial", standardize=FALSE,
                     alpha=1, lambda=fit.lasso$lambda.min)
  cv.enet1 <- glmnet(xtrain, ytrain, family="binomial", standardize=FALSE,
                     alpha=0.05, lambda=fit.enet1$lambda.min)
  cv.enet2 <- glmnet(xtrain, ytrain, family="binomial", standardize=FALSE,
                     alpha=0.5, lambda=fit.enet2$lambda.min)
  cv.enet3 <- glmnet(xtrain, ytrain, family="binomial", standardize=FALSE,
                     alpha=0.95, lambda=fit.enet3$lambda.min)
  scv.ridge <- glmnet(sxtrain, ytrain, family="binomial", standardize=FALSE,
                      alpha=0, lambda=sfit.ridge$lambda.min)
  scv.lasso <- glmnet(sxtrain, ytrain, family="binomial", standardize=FALSE,
                      alpha=1, lambda=sfit.lasso$lambda.min)
  scv.enet1 <- glmnet(sxtrain, ytrain, family="binomial", standardize=FALSE,
                      alpha=0.05, lambda=sfit.enet1$lambda.min)
  scv.enet2 <- glmnet(sxtrain, ytrain, family="binomial", standardize=FALSE,
                      alpha=0.5, lambda=sfit.enet2$lambda.min)
  scv.enet3 <- glmnet(sxtrain, ytrain, family="binomial", standardize=FALSE,
                      alpha=0.95, lambda=sfit.enet3$lambda.min)
  pred.lasso[i] <- as.numeric(predict(cv.lasso, matrix(xtest, nrow=1),
                                      type="response"))
  pred.ridge[i] <- as.numeric(predict(cv.ridge, matrix(xtest, nrow=1),
                                      type="response"))
  pred.enet1[i] <- as.numeric(predict(cv.enet1, matrix(xtest, nrow=1),
                                      type="response"))
  pred.enet2[i] <- as.numeric(predict(cv.enet2, matrix(xtest, nrow=1),
                                      type="response"))
  pred.enet3[i] <- as.numeric(predict(cv.enet3, matrix(xtest, nrow=1),
                                      type="response"))
  spred.lasso[i] <- as.numeric(predict(scv.lasso, matrix(sxtest, nrow=1),
                                      type="response"))
  spred.ridge[i] <- as.numeric(predict(scv.ridge, matrix(sxtest, nrow=1),
                                      type="response"))
  spred.enet1[i] <- as.numeric(predict(scv.enet1, matrix(sxtest, nrow=1),
                                      type="response"))
  spred.enet2[i] <- as.numeric(predict(scv.enet2, matrix(sxtest, nrow=1),
                                      type="response"))
  spred.enet3[i] <- as.numeric(predict(scv.enet3, matrix(sxtest, nrow=1),
                                      type="response"))
}


auc.ridge <- pROC::roc(resp, pred.ridge)$auc
auc.lasso <- pROC::roc(resp, pred.lasso)$auc
auc.enet1 <- pROC::roc(resp, pred.enet1)$auc
auc.enet2 <- pROC::roc(resp, pred.enet2)$auc
auc.enet3 <- pROC::roc(resp, pred.enet3)$auc
sauc.ridge <- pROC::roc(resp, spred.ridge)$auc
sauc.lasso <- pROC::roc(resp, spred.lasso)$auc
sauc.enet1 <- pROC::roc(resp, spred.enet1)$auc
sauc.enet2 <- pROC::roc(resp, spred.enet2)$auc
sauc.enet3 <- pROC::roc(resp, spred.enet3)$auc

fit.ridge$nzero[fit.ridge$lambda==fit.ridge$lambda.min]
fit.lasso$nzero[fit.lasso$lambda==fit.lasso$lambda.min]
fit.enet1$nzero[fit.enet1$lambda==fit.enet1$lambda.min]
fit.enet2$nzero[fit.enet2$lambda==fit.enet2$lambda.min]
fit.enet3$nzero[fit.enet3$lambda==fit.enet3$lambda.min]
sfit.ridge$nzero[sfit.ridge$lambda==sfit.ridge$lambda.min]
sfit.lasso$nzero[sfit.lasso$lambda==sfit.lasso$lambda.min]
sfit.enet1$nzero[sfit.enet1$lambda==sfit.enet1$lambda.min]
sfit.enet2$nzero[sfit.enet2$lambda==sfit.enet2$lambda.min]
sfit.enet3$nzero[sfit.enet3$lambda==sfit.enet3$lambda.min]

varr.cpg <- tapply(est.ridge, AnnoFile$Relation_to_UCSC_CpG_Island, function(g) {
  var(g)})[c(1, 3, 5, 2, 4)]
sizes.cpg <- rle(sort(AnnoFile$Relation_to_UCSC_CpG_Island))$lengths[c(1, 3, 5, 2, 4)]

varr.chr <- tapply(est.ridge, as.character(AnnoFile$CHR), function(g) {
  var(g)})[-c(24, 25)]
varr.chr <- varr.chr[order(as.numeric(names(varr.chr)))]
sizes.chr <- rle(sort(as.character(AnnoFile$CHR)))
sizes.chr <- sizes.chr$lengths[order(as.numeric(sizes.chr$values))]

varl.cpg <- tapply(est.lasso, AnnoFile$Relation_to_UCSC_CpG_Island, function(g) {
  var(g)})[c(1, 3, 5, 2, 4)]

varl.chr <- tapply(est.lasso, as.character(AnnoFile$CHR), function(g) {
  var(g)})[-c(24, 25)]
varl.chr <- varl.chr[order(as.numeric(names(varl.chr)))]

plot(varr.cpg)
plot(varr.chr)
plot(varl.cpg)
plot(varl.chr)

boxplot(pval.t ~ as.character(AnnoFile$Relation_to_UCSC_CpG_Island))
boxplot(pval.w ~ as.character(AnnoFile$Relation_to_UCSC_CpG_Island))
boxplot(abs(est.ridge) ~ as.character(AnnoFile$Relation_to_UCSC_CpG_Island))
boxplot(abs(est.lasso) ~ as.character(AnnoFile$Relation_to_UCSC_CpG_Island))
boxplot(abs(est.enet1) ~ as.character(AnnoFile$Relation_to_UCSC_CpG_Island))
boxplot(abs(est.enet2) ~ as.character(AnnoFile$Relation_to_UCSC_CpG_Island))
boxplot(abs(est.enet3) ~ as.character(AnnoFile$Relation_to_UCSC_CpG_Island))

boxplot(pval.t ~ as.character(AnnoFile$CHR))
boxplot(pval.w ~ as.character(AnnoFile$CHR))
boxplot(abs(est.ridge) ~ as.character(AnnoFile$CHR))
boxplot(abs(est.lasso) ~ as.character(AnnoFile$CHR))

rle(sort(AnnoFile$Relation_to_UCSC_CpG_Island[which(est.lasso2!=0)]))
rle(sort(as.character(AnnoFile$CHR[which(est.lasso2!=0)])))$lengths

cpg.grridge <- CreatePartition(as.factor(AnnoFile$Relation_to_UCSC_CpG_Island))
chr.grridge <- CreatePartition(as.factor(as.character(AnnoFile$CHR)))
part.grridge <- list(cpg=cpg.grridge, chr=chr.grridge)
fit1.grridge <- grridge(betaPreProc, resp, list(cpg=cpg.grridge, 
                                                chr=chr.grridge),
                        optl=28.8439273148345)#lambda=28.8439273148345
fit2.grridge <- grridge(betaPreProc, resp, list(chr=chr.grridge),
                        optl=28.8439273148345)
fit3.grridge <- grridge(betaPreProc, resp, list(cpg=cpg.grridge),
                       optl=28.8439273148345)
fit4.grridge <- grridge(t(smeth.norm), resp, part.grridge, 
                        innfold=10)#lambda=132.842714898256


# randomly selecting features for estimation of penalties
cpg.greben <- list(cpg=as.numeric(as.factor(
  AnnoFile$Relation_to_UCSC_CpG_Island)))
chr.greben <- list(chr=as.numeric(as.factor(as.character(AnnoFile$CHR))))
cpg.G <- length(unique(cpg.greben$cpg))
chr.G <- length(unique(chr.greben$chr))
cpg.samp <- as.numeric(sapply(1:cpg.G, function(g) {
  sample(which(cpg.greben$cpg==g), 10000/cpg.G)}))
chr.samp <- as.numeric(sapply(1:chr.G, function(g) {
  sample(which(chr.greben$chr==g), 10000/chr.G)}))
rcpg.meth <- t(betaPreProc)[, cpg.samp]
rchr.meth <- t(betaPreProc)[, chr.samp]
rcpg.greben <- list(cpg=cpg.greben$cpg[cpg.samp])
rchr.greben <- list(chr=chr.greben$chr[chr.samp])
fit.greben1 <- grEBEN3(rcpg.meth, resp, rep(1, length(resp)), 
                       partitions=rcpg.greben, alpha=0.05, psel=TRUE)
fit.greben2 <- grEBEN3(rchr.meth, resp, rep(1, length(resp)), 
                       partitions=rchr.greben, alpha=0.05, psel=TRUE)

smeth <- t(betaPreProc[tab.limma$logFC > 0.5, ])
scpg.greben <- list(cpg=cpg.greben$cpg[tab.limma$logFC > 0.5])
fit.greben3 <- grEBEN3(smeth, resp, rep(1, length(resp)), 
                       partitions=scpg.greben, alpha=0.5, psel=TRUE)

