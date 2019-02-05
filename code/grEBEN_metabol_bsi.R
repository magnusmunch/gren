### paths
path.data <- as.character("/Users/magnusmunch/Documents/PhD/EBEN/data/bsi_MTBLS315_20170922_111318/")

### libraries
library(glmnet)
library(pROC)
library(GRridge)

### reading in data
object.bsi <- read.table(paste(path.data, "s_NMFI and BSI diagnosis.txt", sep=""), 
                       header=TRUE, stringsAsFactors=FALSE, fill=TRUE)
object.gc <- read.table(paste(path.data, "m_GC_nmfi_and_bsi_diagnosis_v2_maf.tsv", sep=""), 
                      header=TRUE, stringsAsFactors=FALSE, fill=TRUE)
object.gc <- object.gc[!is.na(object.gc$metabolite_identification) & 
                         object.gc$metabolite_identification!=" " &
                         !duplicated(object.gc), ]
object.lc <- read.table(paste(path.data, "m_LC_nmfi_and_bsi_diagnosis_v2_maf.tsv", sep=""), 
                      header=TRUE, stringsAsFactors=FALSE, fill=TRUE, sep="\t", quote="", 
                      na.strings="")
object.lc <- object.lc[!is.na(object.lc$metabolite_identification) & 
                         object.lc$metabolite_identification!=" " &
                         !duplicated(object.lc), ]
object.uplc_neg <- read.table(paste(path.data, "m_UPLC_NEG_nmfi_and_bsi_diagnosis_v2_maf.tsv", sep=""), 
                              header=TRUE, stringsAsFactors=FALSE, fill=TRUE, sep="\t", quote="", 
                              na.strings="")
object.uplc_neg <- object.uplc_neg[!is.na(object.uplc_neg$metabolite_identification) & 
                                     object.uplc_neg$metabolite_identification!=" " & 
                                     !duplicated(object.uplc_neg), ]
object.uplc_pos <- read.table(paste(path.data, "m_UPLC_POS_nmfi_and_bsi_diagnosis_v2_maf.tsv", sep=""), 
                            header=TRUE, stringsAsFactors=FALSE, fill=TRUE, sep="\t", quote="", 
                            na.strings="")
object.uplc_pos <- object.uplc_pos[!is.na(object.uplc_pos$metabolite_identification) & 
                                     object.uplc_pos$metabolite_identification!=" " &
                                     !duplicated(object.uplc_pos), ]

meta.uplc_neg <- read.table(paste(path.data, "a_UPLC_NEG_nmfi_and_bsi_diagnosis.txt", sep=""), 
                      header=TRUE, stringsAsFactors=FALSE, fill=TRUE)
meta.uplc_pos <- read.table(paste(path.data, "a_UPLC_POS_nmfi_and_bsi_diagnosis.txt", sep=""), 
                            header=TRUE, stringsAsFactors=FALSE, fill=TRUE)

### data manipulation
# matching sample names
samp.uplc_neg <- meta.uplc_neg$Sample.Name[match(sub("X", "", colnames(object.uplc_neg)[
  22:ncol(object.uplc_neg)]), meta.uplc_neg$MS.Assay.Name)]
samp.uplc_pos <- meta.uplc_pos$Sample.Name[match(sub("X", "", colnames(object.uplc_pos)[
  22:ncol(object.uplc_pos)]), meta.uplc_pos$MS.Assay.Name)]
samp.gc <- ifelse(nchar(sub("X", "", colnames(object.gc)[22:ncol(object.gc)]))==2, 
                  paste("MURB0", sub("X", "", colnames(object.gc)[22:ncol(object.gc)]), sep=""),
                  paste("MCMA", sub("X", "", colnames(object.gc)[22:ncol(object.gc)]), sep=""))
samp.lc <- sub("[^.]*.[^.]*.", "", sub("_([^_]*)$", "", colnames(object.lc)[22:ncol(object.lc)]))

# creating data object
object.metabol <- rbind(as.matrix(object.uplc_neg[, 22:ncol(object.uplc_neg)][, order(samp.uplc_neg)]),
                        as.matrix(object.uplc_pos[, 22:ncol(object.uplc_pos)][, order(samp.uplc_pos)]),
                        as.matrix(object.gc[, 22:ncol(object.gc)][, order(samp.gc)]),
                        as.matrix(object.lc[, 22:ncol(object.lc)][, order(samp.lc)]))
rawdata.metabol <- matrix(t(object.metabol), nrow=nrow(object.bsi),
                          dimnames=list(sort(samp.uplc_neg), 
                                        c(object.uplc_neg$metabolite_identification, 
                                          object.uplc_pos$metabolite_identification,
                                          object.gc$metabolite_identification, 
                                          object.lc$metabolite_identification)))

# unique identifiers of the metabolites
id.metabol <- c(paste(object.uplc_neg$metabolite_identification, object.uplc_neg$mass_to_charge, 
                      object.uplc_neg$fragmentation, object.uplc_neg$modifications,
                      object.uplc_neg$charge, object.uplc_neg$retention_time),
                paste(object.uplc_pos$metabolite_identification, object.uplc_pos$mass_to_charge, 
                      object.uplc_pos$fragmentation, object.uplc_pos$modifications,
                      object.uplc_pos$charge, object.uplc_pos$retention_time),
                paste(object.gc$metabolite_identification, object.gc$mass_to_charge, 
                      object.gc$fragmentation, object.gc$modifications,
                      object.gc$charge, object.gc$retention_time),
                paste(object.lc$metabolite_identification, object.lc$mass_to_charge, 
                      object.lc$fragmentation, object.lc$modifications,
                      object.lc$charge, object.lc$retention_time))

# removing variables and samples with more than 10% missing
data.metabol <- rawdata.metabol[, !apply(rawdata.metabol, 2, function(metabol) {
  (sum(is.na(metabol))/length(metabol)) > 0.1})]
data.metabol <- data.metabol[!apply(data.metabol, 1, function(sample) {
  (sum(is.na(sample))/length(sample)) > 0.1}), ]
id.metabol <- id.metabol[!apply(rawdata.metabol, 2, function(metabol) {
  (sum(is.na(metabol))/length(metabol)) > 0.1})]

# impute half of the minimum for the missing values
impute.metabol <- apply(data.metabol, 2, function(metabol) {
  imput <- min(metabol, na.rm=TRUE)/2; ifelse(is.na(metabol), imput, metabol)})

# selecting metabolites that show variation
select.metabol <- impute.metabol[, apply(impute.metabol, 2, sd)!=0]
id.metabol <- id.metabol[apply(impute.metabol, 2, sd)!=0]

# removing duplicate variables
whichdupli.metabol <- sapply(which(!duplicated(t(select.metabol))), function(um) {
  as.numeric(which(apply(select.metabol, 2, function(metabol) {all(metabol==select.metabol[, um])})))})
names.metabol <- lapply(whichdupli.metabol, function(umetabol) {colnames(select.metabol)[umetabol]})
id.metabol <- lapply(whichdupli.metabol, function(umetabol) {id.metabol[umetabol]})
dupli.metabol <- select.metabol[, which(!duplicated(t(select.metabol)))]

# transforming and standarizing data
trans.metabol <- sqrt(dupli.metabol)
norm.metabol <- apply(trans.metabol, 2, function(x) {(x - mean(x))/sd(x)})

# co-data
sds <- apply(trans.metabol, 2, sd)
parSds <- CreatePartition(sds, ngroup=5, uniform=TRUE)
parDupli <- CreatePartition(unlist(lapply(whichdupli.metabol, length)), ngroup=5, uniform=TRUE)
no.assays <- (colnames(norm.metabol) %in% object.uplc_neg$metabolite_identification) +
  (colnames(norm.metabol) %in% object.uplc_pos$metabolite_identification) + 
  (colnames(norm.metabol) %in% object.gc$metabolite_identification) +
  (colnames(norm.metabol) %in% object.lc$metabolite_identification)
parNoassays <- CreatePartition(as.factor(no.assays))

# response data
resp <- object.bsi$Factor.Value.patient.group.[match(object.bsi$Sample.Name, rownames(norm.metabol))]
resp.nmfi <- as.numeric(resp=="non-malarial febrile illness")
resp.bsi <- as.numeric(resp=="bacterial bloodstream infection")

# final data dimensions
n <- length(resp.nmfi)
p <- ncol(norm.metabol)

### fitting models
fit.lasso <- cv.glmnet(norm.metabol, resp.nmfi, family="binomial", standardize=FALSE, alpha=1,
                       nfolds=n, grouped=FALSE)
fit.ridge <- cv.glmnet(norm.metabol, resp.nmfi, family="binomial", standardize=FALSE, alpha=0,
                       nfolds=n, grouped=FALSE)

fit1.grridge <- grridge(t(matrix(norm.metabol, ncol=p, dimnames=list(NULL, NULL))), resp.nmfi, 
                        list(assay=parNoassays), selectionEN=TRUE, maxsel=10)
fit2.grridge <- grridge(t(matrix(norm.metabol, ncol=p, dimnames=list(NULL, NULL))), resp.nmfi, 
                        list(duplicated=parDupli))
fit3.grridge <- grridge(t(matrix(norm.metabol, ncol=p, dimnames=list(NULL, NULL))), resp.nmfi, 
                        list(sds=parSds))
cv1.grridge <- grridgeCV(fit1.grridge, t(matrix(norm.metabol, ncol=p, dimnames=list(NULL, NULL))), 
                         resp.nmfi)
cv2.grridge <- grridgeCV(fit2.grridge, t(matrix(norm.metabol, ncol=p, dimnames=list(NULL, NULL))), 
                         resp.nmfi)
cv3.grridge <- grridgeCV(fit3.grridge, t(matrix(norm.metabol, ncol=p, dimnames=list(NULL, NULL))), 
                         resp.nmfi)

sellasso <- function(x, y, psel) {
  fit.lasso <- glmnet(x, y, family="binomial", alpha=1, standardize=FALSE, dfmax=psel)
  found <- any(fit.lasso$df==psel)
  while(!found) {
    lambdamin <- fit.lasso$lambda[head(which(fit.lasso$df > psel), 1L)]
    lambdamax <- fit.lasso$lambda[tail(which(fit.lasso$df < psel), 1L)]
    lambdaseq <- exp(seq(log(lambdamax), log(lambdamin), length.out=100))
    fit.lasso <- glmnet(x, y, family="binomial", alpha=1, standardize=FALSE, 
                        lambda=lambdaseq)
    found <- any(fit.lasso$df==psel)
  }
  fit.lasso$s <- tail(fit.lasso$lambda[fit.lasso$df==psel], 1L)
  return(fit.lasso)
}

sellassoCV <- function(x, y, psel, nfolds=nrow(x)) {
  n <- nrow(x)
  p <- ncol(x)
  rest <- n %% nfolds
  foldsize <- c(rep(n %/% nfolds + as.numeric(rest!=0), times=rest),
                rep(n %/% nfolds, times=nfolds - rest))
  foldid <- sample(rep(1:nfolds, times=foldsize))
  pred <- numeric(n) 
  for(k in 1:nfolds) {
    xtrain <- x[foldid!=k, ]
    ytrain <- y[foldid!=k]
    xtest <- x[foldid==k, ]
    
    fit.lasso <- sellasso(xtrain, ytrain, psel=psel)
    pred[k] <- predict(fit.lasso, matrix(xtest, ncol=p), s=fit.lasso$s, 
                       type="response")
  }
  return(pred)
}

pselseq <- c(5:10, 15, 20, 25, 30, 40, 50)
auc <- matrix(NA, nrow=pselseq, ncol=2)
for(csel in 1:length(pselseq)) {
  psel <- pselseq[csel]
  fit.grridge <- grridge(t(matrix(norm.metabol, ncol=p, dimnames=list(NULL, NULL))), resp.nmfi, 
                         list(assay=parNoassays), selectionEN=TRUE, maxsel=psel)
  cv.lasso <- sellassoCV(norm.metabol, resp.nmfi, psel=psel)
  cv.grridge <- grridgeCV(fit.grridge, t(matrix(norm.metabol, ncol=p, dimnames=list(NULL, NULL))), 
                          resp.nmfi)
  
  
}

pROC::roc(resp.nmfi, cv1.grridge$NoGroups)$auc
pROC::roc(resp.nmfi, cv1.grridge$GroupRegul)$auc
pROC::roc(resp.nmfi, cv1.grridge$EN)$auc
pROC::roc(resp.nmfi, cv2.grridge$GroupRegul)$auc
pROC::roc(resp.nmfi, cv3.grridge$GroupRegul)$auc

# cross-validate prediction under number of assays co-data
pred <- matrix(NA, nrow=n, ncol=4)
colnames(pred) <- c("lasso", "ridge.glmnet", "ridge.penalized", "GRridge")
for(k in c(1:n)) {
  
  cv.lasso <- glmnet(norm.metabol[-k, ], resp.nmfi[-k], family="binomial", standardize=FALSE, alpha=1,
                     lambda=fit.lasso$lambda.min)
  cv.ridge <- glmnet(norm.metabol[-k, ], resp.nmfi[-k], family="binomial", standardize=FALSE, alpha=0,
                     lambda=fit.ridge$lambda.min)
  cv.grridge <- grridge(t(matrix(norm.metabol[-k, ], ncol=p, nrow=n - 1, dimnames=list(NULL, NULL))),
                        resp.nmfi[-k], list(assay=parNoassays), optl=fit1.grridge$optl)
  
  pred[k, ] <- c(predict(cv.lasso, matrix(norm.metabol[k, ], nrow=1), type="response"),
                 predict(cv.ridge, matrix(norm.metabol[k, ], nrow=1), type="response"),
                 predict.grridge(cv.grridge, matrix(norm.metabol[k, ], ncol=1)))

}

pROC::roc(resp.nmfi, pred[, 1])$auc
pROC::roc(resp.nmfi, pred[, 2])$auc
pROC::roc(resp.nmfi, pred[, 3])$auc
pROC::roc(resp.nmfi, pred[, 4])$auc




bs <- c(8.39, 130.03, 15.48, 30, 5.71, 22.85, 23.87, 99, 32.8, 19.99, 102.2, 50, 34, 35, 14, 20, 31.10,
        17.63, 7.5, 9, 47.83, 8.84, 8.07, 15.52, 2.97, 111.12)
