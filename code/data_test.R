
path.data <- as.character("/Users/magnusmunch/Downloads/MTBLS315_20170922_111318/")
object.bsi <- read.table(paste(path.data, "s_NMFI and BSI diagnosis.txt", sep=""), 
                       header=TRUE, stringsAsFactors=FALSE, fill=TRUE)
object.gc <- read.table(paste(path.data, "m_GC_nmfi_and_bsi_diagnosis_v2_maf.tsv", sep=""), 
                      header=TRUE, stringsAsFactors=FALSE, fill=TRUE)
object.lc <- read.table(paste(path.data, "m_LC_nmfi_and_bsi_diagnosis_v2_maf.tsv", sep=""), 
                      header=TRUE, stringsAsFactors=FALSE, fill=TRUE, sep="\t", quote="", 
                      na.strings="")
object.uplc_neg <- read.table(paste(path.data, "m_UPLC_NEG_nmfi_and_bsi_diagnosis_v2_maf.tsv", sep=""), 
                            header=TRUE, stringsAsFactors=FALSE, fill=TRUE, sep="\t", quote="", 
                            na.strings="")
object.uplc_pos <- read.table(paste(path.data, "m_UPLC_POS_nmfi_and_bsi_diagnosis_v2_maf.tsv", sep=""), 
                            header=TRUE, stringsAsFactors=FALSE, fill=TRUE, sep="\t", quote="", 
                            na.strings="")

meta.uplc_neg <- read.table(paste(path.data, "a_UPLC_NEG_nmfi_and_bsi_diagnosis.txt", sep=""), 
                      header=TRUE, stringsAsFactors=FALSE, fill=TRUE)
meta.uplc_pos <- read.table(paste(path.data, "a_UPLC_POS_nmfi_and_bsi_diagnosis.txt", sep=""), 
                            header=TRUE, stringsAsFactors=FALSE, fill=TRUE)

samp.uplc_neg <- meta.uplc_neg$Sample.Name[match(sub("X", "", colnames(object.uplc_neg)[
  22:ncol(object.uplc_neg)]), meta.uplc_neg$MS.Assay.Name)]
samp.uplc_pos <- meta.uplc_pos$Sample.Name[match(sub("X", "", colnames(object.uplc_pos)[
  22:ncol(object.uplc_pos)]), meta.uplc_pos$MS.Assay.Name)]
samp.gc <- ifelse(nchar(sub("X", "", colnames(object.gc)[22:ncol(object.gc)]))==2, 
                  paste("MURB0", sub("X", "", colnames(object.gc)[22:ncol(object.gc)]), sep=""),
                  paste("MCMA", sub("X", "", colnames(object.gc)[22:ncol(object.gc)]), sep=""))
samp.lc <- sub("[^.]*.[^.]*.", "", sub("_([^_]*)$", "", colnames(object.lc)[22:ncol(object.lc)]))

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

data.metabol <- rawdata.metabol[, which(apply(rawdata.metabol, 2, function(metabol) {
  all(!is.na(metabol))}))]
select.metabol <- data.metabol[, apply(data.metabol, 2, sd)!=0]
trans.metabol <- sqrt(select.metabol)
norm.metabol <- apply(trans.metabol, 2, function(x) {(x - mean(x))/sd(x)})

assay.metabol <- rep(c(1:4), times=c(nrow(object.uplc_neg), nrow(object.uplc_pos), nrow(object.gc), 
                                     nrow(object.lc)))[which(apply(rawdata.metabol, 2, function(metabol) {
  all(!is.na(metabol))}))]

parAssay <- CreatePartition(as.factor(assay.metabol))

resp <- object.bsi$Factor.Value.patient.group.[match(object.bsi$Sample.Name, rownames(norm.metabol))]
resp.nmfi <- as.numeric(resp=="non-malarial febrile illness")
resp.bsi <- as.numeric(resp=="bacterial bloodstream infection")

n <- length(resp.nmfi)
p <- ncol(norm.metabol)

fit.lasso <- cv.glmnet(norm.metabol, resp.nmfi, family="binomial", standardize=FALSE)
fit.ridge <- cv.glmnet(norm.metabol, resp.nmfi, family="binomial", standardize=FALSE, alpha=0)

fit.grridge <- grridge(t(matrix(norm.metabol, ncol=p, dimnames=list(NULL, NULL))), resp.nmfi, 
                       list(assay=parAssay))
cv.grridge <- grridgeCV(fit.grridge, t(matrix(norm.metabol, ncol=p, dimnames=list(NULL, NULL))), 
                        resp.nmfi, list(assay=parAssay))

dupli.metabol <- vector(mode="list", length=sum(!duplicated(t(norm.metabol))))
um <- which(!duplicated(t(norm.metabol)))[1]
test <- sapply(which(!duplicated(t(norm.metabol))), function(um) {
  as.numeric(which(apply(norm.metabol, 2, function(metabol) {all(metabol==norm.metabol[, um])})))})


bs <- c(8.39, 130.03, 15.48, 30, 5.71, 22.85, 23.87, 99, 32.8, 19.99, 102.2, 50, 34, 35, 14, 20, 31.10,
        17.63, 7.5, 9, 47.83, 8.84, 8.07, 15.52, 2.97, 111.12)
