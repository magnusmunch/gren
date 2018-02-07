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
path.code <- ifelse(as.character(Sys.info()[1])!="Darwin", "~/EBEN/code/",
                    "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/code/")
path.res <- ifelse(as.character(Sys.info()[1])!="Darwin", "~/EBEN/results/",
                   "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/results/")
path.graph <- "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/graphs/"

### loading libraries
library(cancerdata)
library(impute)
library(GRridge)
library(glmnet)
library(pROC)
library(org.Hs.eg.db)
library(biomaRt)
library(VarfromPDB)

### loading data
data(VEER) 
data(VEER1) # filtered on genes version
# localPDB(paste(path.data, "localPDB", sep=""))
info.clinvar <- extract_clinvar("breast", paste(path.data, "localPDB", sep=""))
info.orphanet <- extract_genes_orphanet("breast", paste(path.data, "localPDB", 
                                                        sep=""))
info.omim <- extract_omim("breast", "kIWps7RTToWh8LhHlnIdFg", paste(path.data, "localPDB", sep=""))
# info.pubmed <- extract_pubmed("pubmed AND gene", "breast", 
#                         paste(path.data, "localPDB", sep=""))
info.uniprot <- extract_uniprot("breast", paste(path.data, "localPDB", sep=""))
ass.genes <- setdiff(unique(c(as.character(info.clinvar$gene2dis$AssociatedGenes), 
                              as.character(info.clinvar$gene2dis$RelatedGenes),
                              info.orphanet[, 3], 
                              info.uniprot[[1]]$GeneSymbol)), "")
ass.genelist <- getBM(filters="hgnc_symbol", values=ass.genes, mart=mart,
                      attributes=c("ensembl_gene_id", "hgnc_symbol", "embl"))


test <- read.table



### mapping gene names
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
filters <- listFilters(mart)
alias <- as.data.frame(org.Hs.egALIAS2EG)
genbank <- as.data.frame(org.Hs.egACCNUM)



c("ensembl_gene_id", "hgnc_symbol") 0.2857726

# different sources of gene names
genes1 <- trimws(as.character(fData(VEER)$symbol), "both")
genes1[genes1==""] <- NA
genelist1 <- getBM(values=genes1, mart=mart,
                   attributes=c("ensembl_gene_id", "hgnc_symbol"))
hgnc1 <- genelist1$hgnc_symbol[match(genes1, genelist1$hgnc_symbol)]
mean(!is.na(hgnc1))

genes2 <- nomen$Genbank.Accession.Number[
  match(trimws(rownames(exprs(VEER)), "both"), nomen$Phil.Green.Contig)]
genes2 <- as.numeric(genbank$gene_id[match(genes2, genbank$accession)])
genelist2 <- getBM(values=genes2, mart=mart,
                   attributes=c("entrezgene", "ensembl_gene_id", "hgnc_symbol"))
hgnc2 <- ifelse(!is.na(hgnc1), hgnc1, genelist2$hgnc_symbol[
  match(genes2, genelist2$entrezgene)])
mean(!is.na(hgnc2))

# they are accession numbers
# not entrezgene 
# maybe protein id
# not ensemble_gene_id
# maybe embl

genes3 <- alias$gene_id[match(genes1, alias$alias_symbol)]
genelist3 <- getBM(filters="entrezgene", values=genes3, mart=mart,
                   attributes=c("ensembl_gene_id", "hgnc_symbol", "entrezgene"))
hgnc3 <- genelist3$hgnc_symbol[match(genes3, genelist3$entrezgene)]

listAttributes(mart)$description[1001:length(listAttributes(mart)$description)]
genes2
genelist2$embl
all(is.na(match(genes2, genelist2$embl)))

str(genelist2)


hgnc2 <- ifelse(!is.na(hgnc1), hgnc1, )
genelist2$hgnc_symbol[match(genes2, genelist2$refseq_mrna)]


genes1 <- trimws(as.character(fData(VEER)$symbol), "both")
genelist1 <- getBM(values=genes1, mart=mart,
                   attributes=c("ensembl_gene_id", "hgnc_symbol"))
hgnc1 <- genelist1$hgnc_symbol[match(genes1, genelist1$hgnc_symbol)]
mean(!is.na(hgnc1))

genes2 <- trimws(rownames(exprs(VEER)), "both")
genelist2 <- getBM(filters="embl", values=genes2, mart=mart,
                   attributes=c("ensembl_gene_id", "hgnc_symbol", "embl"))
hgnc2 <- genelist2$hgnc_symbol[match(genes2, genelist2$embl)]

genes3 <- alias$gene_id[match(genes2, alias$alias_symbol)]
genelist3 <- getBM(filters="entrezgene", values=genes3, mart=mart,
                   attributes=c("ensembl_gene_id", "hgnc_symbol", "entrezgene"))
hgnc3 <- genelist3$hgnc_symbol[match(genes3, genelist3$entrezgene)]

genes4 <- trimws(rownames(exprs(VEER)), "both")
genelist4 <- getBM(filters=c("refseq_mrna"), values=genes4, mart=mart,
                   attributes=c("ensembl_gene_id", "hgnc_symbol", "embl", 
                                "refseq_mrna", "chromosome_name"))
hgnc4 <- genelist4$hgnc_symbol[match(genes4, genelist4$refseq_mrna)]

nomen <- read.table("/Users/magnusmunch/Downloads/ArrayNomenclature_contig_accession.csv",
                    header=TRUE, sep=",", stringsAsFactors=FALSE)
genes5 <- nomen$Genbank.Accession.Number[
  match(trimws(rownames(exprs(VEER)), "both"), nomen$Phil.Green.Contig)]
genelist5 <- getBM(values=genes5, mart=mart,
                   attributes = c("protein_id", "embl", "hgnc_symbol"))
str(genelist5)
listFilters(mart)$name
any(genes5 %in% genelist5$embl)
any(genes5 %in% genelist5$protein_id)
genes4
  
  
cbind(rownames(exprs(VEER)), trimws(as.character(fData(VEER)$symbol), "both"),
      )


hgnc2 <- ifelse(is.na(hgnc2), hgnc4, hgnc2)
cbind(trimws(rownames(exprs(VEER)), "both"), hgnc1, hgnc2, hgnc4)
dim(exprs(VEER))

# combining them
hgnc <- apply(cbind(hgnc1, hgnc2, hgnc3), 1, function(x) {
  temp <- rle(sort(x, decreasing=TRUE))
  if(length(x[!is.na(x)])==1) {
    out <- x[!is.na(x)]
  } else if(any(temp$lengths>=2)) {
    out <- as.character(temp$values[temp$lengths>=2])
  } else if(length(x[!is.na(x)])>=2) {
    out <- ifelse(!is.na(x[1]), x[1], ifelse(!is.na(x[2]), x[2], x[3]))
  } else {
    out <- NA
  }
  return(out)})

### data manipulation
resp.breast <- VEER$class 
data1.expr <- exprs(VEER)
which.keep <- which(apply(data1.expr, 1, function(x) {mean(is.na(x)) < 0.1}))
data2.expr <- data1.expr[which.keep, ]
data3.expr <- t(impute.knn(data2.expr)$data) # impute data using knn clustering
# DM is distant metastases within 5 years after diagnosis of breast cancer
# NODM is no metastases within 5 years
set.seed(4005)
id.train <- c(sample(which(VEER$class=="DM"), floor(sum(VEER$class=="DM")/2)), 
              sample(which(VEER$class!="DM"), floor(sum(VEER$class!="DM")/2)))
train.breast <- as.numeric(VEER$class=="DM")[id.train]
train.expr <- data3.expr[id.train, ]
test.breast <- as.numeric(VEER$class=="DM")[-id.train]
test.expr <- data3.expr[-id.train, ]

### standardise data
train.norm <- apply(train.expr, 2, function(x) {(x - mean(x))/sd(x)})
test.norm <- apply(test.expr, 2, function(x) {(x - mean(x))/sd(x)})

### co-data
ass.hgnc1 <- as.numeric(hgnc[which.keep] %in% 
                          unique(ass.genelist$hgnc_symbol)) + 1 
part.ass.grridge1 <- CreatePartition(as.factor(ass.hgnc1))

sds1 <- apply(train.expr, 2, sd)
part.sds.grridge1 <- CreatePartition(sds1, decreasing=TRUE, uniform=TRUE, 
                                     ngroup=10)
part.sds.greben1 <- rep(1:length(part.sds.grridge1), 
                        times=lapply(part.sds.grridge1, length))[
                          order(unlist(part.sds.grridge1))]
part.grridge1 <- list(sds=part.sds.grridge1, association=part.ass.grridge1)
part.greben1 <- list(sds=part.sds.greben1)

### data exploration
p <- ncol(train.norm)
ntrain <- nrow(train.norm)

### fitting models
fit1.ridge <- cv.glmnet(train.norm, train.breast, family="binomial", alpha=0, 
                        nfolds=ntrain, grouped=FALSE, standardize=FALSE)
fit1.lasso <- cv.glmnet(train.norm, train.breast, family="binomial", alpha=1, 
                        nfolds=ntrain, grouped=FALSE, standardize=FALSE)
fit1.enet <- cv.glmnet(train.norm, train.breast, family="binomial", alpha=0.5, 
                       nfolds=ntrain, grouped=FALSE, standardize=FALSE)
fit1.grridge <- grridge(t(train.norm), train.breast, part.grridge1)

### testing performance
methods <- c("ridge", "lasso", "enet", "GRridge")
pred1 <- cbind(predict(fit1.ridge, test.norm, "lambda.min", type="response"), 
               predict(fit1.lasso, test.norm, "lambda.min", type="response"),
               predict(fit1.enet, test.norm, "lambda.min", type="response"),
               predict.grridge(fit1.grridge, t(test.norm))[, 2])
colnames(pred1) <- methods
auc1 <- setNames(apply(pred1, 2, function(pred) {
  pROC::roc(test.breast, pred)$auc}), methods)
psel1 <- setNames(c(fit1.ridge$nzero[fit1.ridge$lambda==fit1.ridge$lambda.min],
                    fit1.lasso$nzero[fit1.lasso$lambda==fit1.lasso$lambda.min],
                    fit1.enet$nzero[fit1.enet$lambda==fit1.enet$lambda.min], p),
                  methods)
brier1 <- setNames(apply(pred1, 2, function(pred) {
  mean((test.breast - pred)^2)}), methods)
brierskill1 <- 1 - brier1/mean((test.breast - mean(test.breast))^2)

### saving results
results1 <- list(auc=auc1, brier=brier1, brierskill=brierskill1, psel=psel1,
                 fitted=list(fit1.ridge, fit1.lasso, fit1.enet, fit1.grridge))
save(results1, file=paste(path.res, "grEBEN_gene_expr_breast_res1.RData", sep=""))



### Only using known genes
# selecting the genes
sel1.expr <- data1.expr[!is.na(hgnc), ]
rownames(sel1.expr) <- hgnc[!is.na(hgnc)]
sel2.expr <- sel1.expr[
  apply(sel1.expr, 1, function(x) {mean(is.na(x)) < 0.1}), ]
sel3.expr <- t(impute.knn(sel2.expr)$data) 

train.sel.expr <- sel3.expr[id.train, ]
test.sel.expr <- sel3.expr[-id.train, ]

### co-data
ass.hgnc2 <- as.numeric(colnames(sel3.expr) %in% 
                          unique(ass.genelist$hgnc_symbol)) + 1 

part.ass.grridge2 <- CreatePartition(as.factor(ass.hgnc2))
part.ass.greben2 <- rep(1:length(part.ass.grridge2), 
                                 times=lapply(part.ass.grridge2, length))[
                                   order(unlist(part.ass.grridge2))]

part.grridge2 <- list(association=part.ass.grridge2)
part.greben2 <- list(association=part.ass.greben2)

### standardise data
train.sel.norm <- apply(train.sel.expr, 2, function(x) {(x - mean(x))/sd(x)})
test.sel.norm <- apply(test.sel.expr, 2, function(x) {(x - mean(x))/sd(x)})
 
### data exploration
p <- ncol(train.sel.norm)

### fit models
fit2.ridge <- cv.glmnet(train.sel.norm, train.breast, family="binomial", 
                        alpha=0, nfolds=ntrain, grouped=FALSE, 
                        standardize=FALSE)
fit2.lasso <- cv.glmnet(train.sel.norm, train.breast, family="binomial", 
                        alpha=1, nfolds=ntrain, grouped=FALSE, 
                        standardize=FALSE)
fit2.enet <- cv.glmnet(train.sel.norm, train.breast, family="binomial", 
                       alpha=0.5, nfolds=ntrain, grouped=FALSE, 
                       standardize=FALSE)
fit2.grridge <- grridge(t(train.sel.norm), train.breast, part.grridge2)

### testing performance
methods <- c("ridge", "lasso", "enet", "grridge")
pred2 <- cbind(predict(fit2.ridge, test.sel.norm, "lambda.min", 
                       type="response"), 
               predict(fit2.lasso, test.sel.norm, "lambda.min", 
                       type="response"),
               predict(fit2.enet, test.sel.norm, "lambda.min", 
                       type="response"),
               predict.grridge(fit2.grridge, t(test.sel.norm))[, 2])
colnames(pred2) <- methods
auc2 <- setNames(apply(pred2, 2, function(pred) {
  pROC::roc(test.breast, pred)$auc}), methods)
psel2 <- setNames(c(fit2.ridge$nzero[fit2.ridge$lambda==fit2.ridge$lambda.min],
                    fit2.lasso$nzero[fit2.lasso$lambda==fit2.lasso$lambda.min],
                    fit2.enet$nzero[fit2.enet$lambda==fit2.enet$lambda.min], p),
                  methods)
brier2 <- setNames(apply(pred2, 2, function(pred) {
  mean((test.breast - pred)^2)}), methods)
brierskill2 <- 1 - brier2/mean((test.breast - mean(test.breast))^2)

### saving results
results2 <- list(auc=auc2, brier=brier2, brierskill=brierskill2, psel=psel2,
                 fitted=list(fit2.ridge, fit2.lasso, fit2.enet, fit2.grridge))
save(results2, file=paste(path.res, "grEBEN_gene_expr_breast_res2.RData", sep=""))









### VEER1 data
### data manipulation
resp.filt.breast <- VEER1$class 
data1.filt.expr <- exprs(VEER1)
which.keep <- which(apply(data1.expr, 1, function(x) {mean(is.na(x)) < 0.1}))
data2.expr <- data1.expr[which.keep, ]
data3.expr <- t(impute.knn(data2.expr)$data) # impute data using knn clustering
# DM is distant metastases within 5 years after diagnosis of breast cancer
# NODM is no metastases within 5 years
set.seed(4005)
id.train <- c(sample(which(VEER$class=="DM"), floor(sum(VEER$class=="DM")/2)), 
              sample(which(VEER$class!="DM"), floor(sum(VEER$class!="DM")/2)))
train.breast <- as.numeric(VEER$class=="DM")[id.train]
train.expr <- data3.expr[id.train, ]
test.breast <- as.numeric(VEER$class=="DM")[-id.train]
test.expr <- data3.expr[-id.train, ]

rownames(data1.filt.expr)


cbind(rownames(data1.expr), as.character(fData(VEER)$symbol))
sum(as.character(fData(VEER)$symbol)!="")


"affy Hu25K"
reffilters <- c("refseq_mrna", "refseq_mrna_predicted")
test1 <- getBM(filters="refseq_mrna", values=rownames(data1.expr), mart=mart,
               attributes=c("ensembl_gene_id", "hgnc_symbol", "embl", 
                            "refseq_mrna", "chromosome_name"))

mean(apply(cbind(test1$hgnc_symbol[match(rownames(data1.expr), test1$refseq_mrna)],
      ifelse(as.character(fData(VEER)$symbol)=="", NA, 
             as.character(fData(VEER)$symbol))), 1, function(x) {any(!is.na(x))}))

mean(!is.na(test1$hgnc_symbol[match(rownames(data1.expr), test1$refseq_mrna)]))




dim(test1)
unique(test1$hgnc_symbol)



attr <- listAttributes(mart)
attr$name



