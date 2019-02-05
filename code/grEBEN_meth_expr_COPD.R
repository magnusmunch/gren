# field	example	SQL type	info	description
# bin	585	smallint(6)	range	Indexing field to speed chromosome range queries.
# chrom	chr1	varchar(255)	values	Reference sequence chromosome or scaffold
# chromStart	28735	int(10) unsigned	range	Start position in chromosome
# chromEnd	29737	int(10) unsigned	range	End position in chromosome
# name	CpG: 111	varchar(255)	values	CpG Island
# length	1002	int(10) unsigned	range	Island Length
# cpgNum	111	int(10) unsigned	range	Number of CpGs in island
# gcNum	731	int(10) unsigned	range	Number of C and G in island
# perCpg	22.2	float	range	Percentage of island that is CpG
# perGc	73	float	range	Percentage of island that is C or G
# obsExp	0.85	float	range	Ratio of observed(cpgNum) to expected(numC*numG/length) CpG in island

### paths
path.data <- ifelse(as.character(Sys.info()[1])!="Darwin", "~/EBEN/data/",
                    "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/data/")
path.code <- ifelse(as.character(Sys.info()[1])!="Darwin", "~/EBEN/code/",
                    "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/code/")
path.res <- ifelse(as.character(Sys.info()[1])!="Darwin", "~/EBEN/results/",
                   "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/results/")
path.graph <- "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/graphs/"

### loading libraries
library(GRridge)
library(glmnet)
library(pROC)
library(COPDSexualDimorphism.data)
library(GenomicRanges)

### loading function
source(paste(path.code, "grVBEM.R", sep=""))

### reading data
data(lgrc.meta)
data(lgrc.expr.meta)
data(lgrc.expr)
data(lgrc.genes)
data(lgrc.methp)

cpg_locations <- read.table(paste(path.data, "cpg_locations.txt", sep=""))
cpg_locations <- cbind(cpg_locations[, 1:4], 
                       paste(cpg_locations[, 5], cpg_locations[, 6]),
                       cpg_locations[, 7:12])
colnames(cpg_locations) <- c("bin", "chrom", "chromStart", "chromEnd", "name",
                             "length", "cpgNum", "gcNum", "perCpg", "perGc", 
                             "obsExp")
# convert 0, to 1-based start
cpg_islands <- makeGRangesFromDataFrame(cpg_locations,
                                        seqnames.field="chrom", 
                                        start.field="chromStart", 
                                        end.field="chromEnd",
                                        starts.in.df.are.0based=TRUE)
cpg_islands$name <- "island"

###############################################################
#             Extract CpG island shores
###############################################################
# extract the shore defined by 2000 bp upstream of cpg islands
shore1 <- flank(cpg_islands, 2000)
# extract the shore defined by 2000 bp downstream of cpg islands
shore2 <- flank(cpg_islands, 2000, FALSE)
# perform intersection and combine the shores where they overlap
shore1_2 <- reduce(c(shore1, shore2))
# extract the features (ranges) that are present in shores only and not in cpg_islands (ie., shores not overlapping islands)
cpgi_shores <- setdiff(shore1_2, cpg_islands)
cpgi_shores$name <- "shore"

###############################################################
#             Extract CpG island shelves
###############################################################
# extract the shore defined by 4000 bp upstream of cpg islands
shelves1 <- flank(cpg_islands, 4000)
# extract the shore defined by 2000 bp downstream of cpg islands
shelves2 <- flank(cpg_islands, 4000, FALSE)
# perform intersection and combine the shelves where they overlap
shelves1_2 <- reduce(c(shelves1, shelves2))
# create a set of ranges consisting CpG Islands, Shores
island_shores <- c(cpg_islands, cpgi_shores)
# extract the features (ranges) that are present in shelves only and not in cpg_islands  or shores(ie., shelves not overlapping islands or shores)
cpgi_shelves=setdiff(shelves1_2, island_shores)
cpgi_shelves$name <- "shelf"

cpg_elements_gr <- GRangesList("islands"=cpg_islands, 
                               "shores"=cpgi_shores, 
                               "shelves"=cpgi_shelves)

names(methp)
site_meth <- methp[, c("chr", "start", "end")]
site_ranges <- makeGRangesFromDataFrame(site_meth,
                                        seqnames.field="chr", 
                                        start.field="start", 
                                        end.field="end",
                                        starts.in.df.are.0based=TRUE)
site_ranges$name <- "sites"

islands <- unique(to(findOverlaps(cpg_islands, site_ranges, type="any", 
                                  select="all", ignore.strand=FALSE)))
shores <- unique(to(findOverlaps(cpgi_shores, site_ranges, type="any", 
                                 select="all", ignore.strand=FALSE)))
shelves <- unique(to(findOverlaps(cpgi_shelves, site_ranges, type="any", 
                                  select="all", ignore.strand=FALSE)))

shores <- shores[!(shores %in% islands)]
shelves <- shelves[!((shelves %in% islands) | (shelves %in% shores))]
distants <- c(1:nrow(methp))[!(c(1:nrow(methp)) %in% 
                                 c(islands, shores, shelves))]


cpg_dist <- rep(c("island", "shore", "shelf", "distant"), 
                times=c(length(islands), length(shores), length(shelves), 
                        length(distants)))[order(c(islands, shores, 
                                                   shelves, distants))]
part.greben <- list(cpg=match(cpg_dist, 
                              c("island", "shore", "shelf", "distant")))
part.grridge <- list(cpg=CreatePartition(as.factor(cpg_dist)))


resp.meth <- as.numeric(meta$diagmaj[
  na.omit(match(colnames(methp), meta$tissueid))]=="2-COPD/Emphysema")
sel.meth <- t(methp[, colnames(methp) %in% meta$tissueid])
norm.meth <- apply(sel.meth, 2, function(x) {(x - mean(x))/sd(x)})

set.seed(8911)
ind.test <- c(sample(which(resp.meth==0), floor(sum(resp.meth==0)*0.25)),
              sample(which(resp.meth==1), floor(sum(resp.meth==1)*0.25)))
xtrain <- norm.meth[-ind.test, ]
ytrain <- resp.meth[-ind.test]
xtest <- norm.meth[ind.test, ]
ytest <- resp.meth[ind.test]

fit.ridge <- cv.glmnet(xtrain, ytrain, family="binomial", alpha=0, 
                       standardize=FALSE)
fit.lasso <- cv.glmnet(xtrain, ytrain, family="binomial", alpha=1, 
                       standardize=FALSE)
fit.enet <- cv.glmnet(xtrain, ytrain, family="binomial", alpha=0.5,
                      standardize=FALSE)
fit.grridge <- grridge(t(xtrain), ytrain, partitions=part.grridge, 
                       optl=0.5*length(ytrain)*fit.ridge$lambda.min)
fit.greben <- grEBEN3(xtrain, ytrain, rep(1, length(ytrain)), 
                      lambda=fit.enet$lambda.min, partitions=part.greben, 
                      alpha=0.5, psel=TRUE)

pred.ridge <- predict(fit.ridge, xtest, s="lambda.min", type="response")
pred.lasso <- predict(fit.lasso, xtest, s="lambda.min", type="response")
pred.enet <- predict(fit.enet, xtest, s="lambda.min", type="response")
pred.grridge <- predict.grridge(fit.grridge, t(xtest))
auc.ridge <- pROC::roc(ytest, as.numeric(pred.ridge))$auc
auc.lasso <- pROC::roc(ytest, as.numeric(pred.lasso))$auc
auc.enet <- pROC::roc(ytest, as.numeric(pred.enet))$auc
auc.grridge <- pROC::roc(ytest, pred.grridge[, 2])$auc

boxplot(abs(as.numeric(coef(fit.ridge, "lambda.min"))[-1]) ~ cpg_dist)





library(RTCGA.clinical)
library(RTCGA.methylation)
library(impute)

data(COAD.methylation)
data(COAD.clinical)

meth <- COAD.methylation
clin <- COAD.clinical

meth$bcr_patient_barcode <- substr(as.character(meth$bcr_patient_barcode), 1, 12)
meth <- t(sapply(split(meth[, -1], meth$bcr_patient_barcode), function(x) {
  apply(x, 2, median, na.rm=TRUE)}))

barcode.clinical <- toupper(clin$patient.bcr_patient_barcode)[
  !is.na(clin$patient.drugs.drug.measure_of_response)]
barcode.meth <- rownames(meth)
barcodes <- barcode.clinical[barcode.clinical %in% barcode.meth]

data <- meth[barcodes, ]
data <- data[, apply(data, 2, function(x) {mean(is.na(x))}) <= 0.4]
data <- t(impute.knn(t(data))$data)
norm <- apply(data, 2, function(x) {
  (x - mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE)})
resp <- as.numeric(clin$patient.drugs.drug.measure_of_response[
  which(toupper(clin$patient.bcr_patient_barcode) %in% barcodes)] %in% 
  c("complete response"))


fit.lasso <- cv.glmnet(norm, resp, family="binomial", standardize=FALSE,
                       alpha=1)
fit.ridge <- cv.glmnet(norm, resp, family="binomial", standardize=FALSE,
                       alpha=0)
methods <- c("lasso", "ridge")
set.seed(2017)
n <- nrow(norm)
p <- ncol(norm)
nfolds <- n
rest <- n %% nfolds
foldid <- sample(rep(1:nfolds, times=round(
  c(rep(n %/% nfolds + as.numeric(rest!=0), times=rest),
    rep(n %/% nfolds, times=nfolds - rest)))))
pred <- matrix(NA, nrow=n, ncol=length(methods), dimnames=list(NULL, methods))
for(k in 1:nfolds) {
  xtrain <- norm[foldid!=k, ]
  ytrain <- resp[foldid!=k]
  xtest <- norm[foldid==k, ]
  cv.lasso <- glmnet(xtrain, ytrain, family="binomial", standardize=FALSE, 
                     alpha=1, lambda=fit.lasso$lambda.min) 
  cv.ridge <- glmnet(xtrain, ytrain, family="binomial", standardize=FALSE, 
                     alpha=0, lambda=fit.ridge$lambda.min) 
  
  pred[foldid==k, 1] <- as.numeric(predict(fit.lasso, matrix(xtest, ncol=p), 
                                           type="response"))
  pred[foldid==k, 2] <- as.numeric(predict(fit.ridge, matrix(xtest, ncol=p), 
                                           type="response"))
}

pROC::roc(resp, pred[, 1])$auc
pROC::roc(resp, pred[, 2])$auc

#### mRNA

library("RTCGA.rnaseq")
data(COAD.rnaseq)

rnaseq <- COAD.rnaseq

rnaseq$bcr_patient_barcode <- substr(as.character(rnaseq$bcr_patient_barcode), 1, 12)
rnaseq <- t(sapply(split(rnaseq[, -1], rnaseq$bcr_patient_barcode), function(x) {
  apply(x, 2, median, na.rm=TRUE)}))

barcode.rnaseq <- rownames(rnaseq)
barcodes <- barcode.clinical[barcode.clinical %in% barcode.rnaseq]

data <- rnaseq[barcodes, ]
data <- data[, apply(data, 2, function(x) {mean(is.na(x))}) <= 0.4 &
               apply(data, 2, sd)!=0]
trans <- sqrt(data + 3/8) - sqrt(3/8)
norm <- apply(trans, 2, function(x) {
  (x - mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE)})
resp <- as.numeric(clin$patient.drugs.drug.measure_of_response[
  which(toupper(clin$patient.bcr_patient_barcode) %in% barcodes)] %in% 
    c("complete response"))

fit.lasso <- cv.glmnet(norm, resp, family="binomial", standardize=FALSE,
                       alpha=1)
fit.ridge <- cv.glmnet(norm, resp, family="binomial", standardize=FALSE,
                       alpha=0)
methods <- c("lasso", "ridge")
set.seed(2017)
n <- nrow(norm)
p <- ncol(norm)
nfolds <- n
rest <- n %% nfolds
foldid <- sample(rep(1:nfolds, times=round(
  c(rep(n %/% nfolds + as.numeric(rest!=0), times=rest),
    rep(n %/% nfolds, times=nfolds - rest)))))
pred <- matrix(NA, nrow=n, ncol=length(methods), dimnames=list(NULL, methods))
for(k in 1:nfolds) {
  xtrain <- norm[foldid!=k, ]
  ytrain <- resp[foldid!=k]
  xtest <- norm[foldid==k, ]
  cv.lasso <- cv.glmnet(xtrain, ytrain, family="binomial", standardize=FALSE, 
                        alpha=1)#, lambda=fit.lasso$lambda.min) 
  cv.ridge <- cv.glmnet(xtrain, ytrain, family="binomial", standardize=FALSE,
                        alpha=0)#, lambda=fit.ridge$lambda.min) 
  
  pred[foldid==k, 1] <- as.numeric(predict(cv.lasso, matrix(xtest, ncol=p), 
                                           s="lambda.min", type="response"))
  pred[foldid==k, 2] <- as.numeric(predict(cv.ridge, matrix(xtest, ncol=p), 
                                           s="lambda.min", type="response"))
}

pROC::roc(resp, pred[, 1])$auc
pROC::roc(resp, pred[, 2])$auc

1 - sum((resp - pred[, 1])^2)/sum((resp - mean(resp))^2)
1 - sum((resp - pred[, 2])^2)/sum((resp - mean(resp))^2)





