### paths
path.data <- "/Users/magnusmunch/Documents/PhD/EBEN/data/DataFelix/"
path.code <- as.character(ifelse(Sys.info()[1]=="Darwin","/Users/magnusmunch/Documents/PhD/EBEN/code/" ,"~/EBEN/code/"))
path.res <- as.character(ifelse(Sys.info()[1]=="Darwin","/Users/magnusmunch/Documents/PhD/EBEN/results/" ,"~/EBEN/results/"))
path.graph <- "/Users/magnusmunch/Documents/PhD/EBEN/graphs/"

### libraries
library(glmnet)
library(penalized)
library(GRridge)
library(pROC)

### functions
# source grENVB functions
source(paste(path.code, "grVBEM.R", sep=""))

### load TCGA data, PRAD is tissue data, PCBL is blood data
# $ctrlIndex - Control index (y = !ctrlIndex) 
# $isoRaw - The raw (but filtered) isomiR counts
# $isoNF - The isomiR normalization factors
# $isoDat - The filtered, normalized isomiR data 
# $mirRaw - Raw (filtered) miRNA counts
# $mirNF - miRNA normalization factors
# $mirDat - The filtered, normalized miRNA data 
# $isomirs - The isomiR variable names
# $mirs - The miRNA variable names
# $mirLookup - A lookup reference linking miRNA to isomiR (also seq length and chromosome)
# $iso_grro - The GRridge isomiR object
# $mir_grro - The GRridge miR object
# $barcodes - Patient reference
load(paste(path.data, "PRAD.Rdata", sep=""))
load(paste(path.data, "PCBL.Rdata", sep=""))

### blood data, 1=cancer, 0 = control
respblood <- 1 - PCBL$ctrlIndex

mirraw <- PCBL$isoRaw
mirblood <- PCBL$isoDat

rs <- rowSums(mirraw)
cs <- colSums(mirraw)

#filtering
mbfilt <- which(rs >= 100)
abund <- rs[mbfilt]
mirblst <- mirblood[mbfilt, ]
dim(mirblst)

parAbund <- CreatePartition(abund, ngroup=10, mingr=25, decreasing=TRUE)
mirblst <- mirblst[unlist(parAbund), ]
# rn <- rownames(mirblst)
# iso2mir <- sapply(rn, function(x) {strsplit(x[1],"\\|")[[1]][2]})
# pmir <- CreatePartition(factor(iso2mir))
# niso <- unlist(lapply(pmir,length))
# pmir2 <- pmir[which(niso>=20)]

groups <- rep(1:length(parAbund), unlist(lapply(parAbund, length)))
partitionRidge <- list(abundance=CreatePartition(as.factor(groups)))
partitionEnet <- list(abundance=groups)
# partitionFarkas <- list(cpg=firstPartition)

# source('C:/Synchr/Rscripts/Classification/GRridge/grridgeCV.R')
# source('C:/Synchr/RPackages/GRridge/temp.R')

fit1.grridge <- grridge(mirblst, respblood, partitionRidge, monotone=TRUE, savepredobj="all",
                        trace=TRUE, optl=10.69514)
fit1.grEBEN <- grVBEM2(t(mirblst), respblood, m=rep(1, length(respblood)), partitionEnet, 
                       lambda1=0.3875383, lambda2=0.9460494, intercept=TRUE, monotone=FALSE, 
                       iso=list(abundance=TRUE), posterior=TRUE, eps=0.001, maxiter=300, 
                       trace=TRUE)
fit1.2.grEBEN <- grVBEM2(t(mirblst), respblood, m=rep(1, length(respblood)), partitionEnet, 
                         lambda1=0.3875383, lambda2=0.9460494, intercept=TRUE, monotone=TRUE, 
                         iso=NULL, posterior=TRUE, eps=0.001, maxiter=300, 
                         trace=TRUE)
fit1.3.grEBEN <- grVBEM2(t(mirblst), respblood, m=rep(1, length(respblood)), partitionEnet, 
                         lambda1=0.3875383, lambda2=0.9460494, intercept=TRUE, monotone=FALSE, 
                         iso=NULL, posterior=TRUE, eps=0.001, maxiter=300, 
                         trace=TRUE)

### tissue data
resptissue <- 1 - PRAD$ctrlIndex
mirtissue <- PRAD$isoDat
mirrawtissue <- PRAD$isoRaw

rs <- rowSums(mirrawtissue)
cs <- colSums(mirrawtissue)

#filtering
mbfilt <- which(rs >= 100)
abund <- rs[mbfilt]
mirtist <- mirtissue[mbfilt, ]
dim(mirtist)

parAbund <- CreatePartition(abund, ngroup=10, mingr=25, decreasing=TRUE)
mirtist <- mirtist[unlist(parAbund), ]

groups <- rep(1:length(parAbund), unlist(lapply(parAbund, length)))
partitionRidge <- list(abundance=CreatePartition(as.factor(groups)))
partitionEnet <- list(abundance=groups)

# cv.ridge <- cv.glmnet(t(mirtist), resptissue, family="binomial", alpha=0, standardize=FALSE, 
#                       intercept=TRUE)
# optl <- cv.ridge$lambda.min*0.5*length(resptissue)

fit2.grridge <- grridge(mirtist, resptissue, partitionRidge, monotone=TRUE, savepredobj="all",
                        trace=TRUE, optl=602.8572, innfold=10)
fit2.grEBEN <- grVBEM2(t(mirtist), resptissue, m=rep(1, length(resptissue)), partitionEnet, 
                       lambda1=2.41, lambda2=38.98, intercept=TRUE, monotone=FALSE, 
                       iso=list(abundance=TRUE), posterior=TRUE, eps=0.001, maxiter=300, 
                       trace=TRUE)




grFcv2 <- grridgeCV(grF2,mirblst,respblood,outerfold=10)

save(grF,grFcv,grF2,grFcv2, file="Fresults.Rdata")
load("Fresults.Rdata")
grFcv
grFcv2

wl <- as.numeric(grF$reslasso$whichlasso)
wl2 <- as.numeric(grF2$resEN$whichEN)

odabund <- order(abund)
summary(match(wl,odabund)/length(abund))
summary(match(wl2,odabund)/length(abund))

dim(PCBL$mirLookup)

rn <- rownames(mirblst)
iso2mir <- sapply(rn, function(x) {strsplit(x[1],"\\|")[[1]][2]})
pmir <- CreatePartition(factor(iso2mir))
niso <- unlist(lapply(pmir,length))
pmir2 <- pmir[which(niso>=5)]


brier <- function(pr,tr) {return(mean((pr-tr)^2))}

cutoffs <- rev(seq(0,1,by=0.01))
rocridgeF <- roc(probs=grFcv[,2],true=grFcv[,1],cutoffs=cutoffs)
auc(rocridgeF)
brier(grFcv[,2],grFcv[,1])
rocgrridgeF <- roc(probs=grFcv[,3],true=grFcv[,1],cutoffs=cutoffs)
auc(rocgrridgeF)
brier(grFcv[,3],grFcv[,1])
roclasso <- roc(probs=grFcv[,4],true=grFcv[,1],cutoffs=cutoffs)
auc(roclasso)
brier(grFcv[,4],grFcv[,1])
rocENF <- roc(probs=grFcv[,5],true=grFcv[,1],cutoffs=cutoffs)
auc(rocENF)
brier(grFcv[,4],grFcv[,1])
rocENF2 <- roc(probs=grFcv2[,4],true=grFcv2[,1],cutoffs=cutoffs)
auc(rocENF2)
brier(grFcv2[,4],grFcv2[,1])
plot(rocridgeF[1,],rocridgeF[2,],type="l",lty=1,ann=F,col="black")
points(rocgrridgeF[1,],rocgrridgeF[2,],type="l",lty=1,col="grey")
points(roclasso[1,],roclasso[2,],type="l",lty=1,col="green")
points(rocENF[1,],rocENF[2,],type="l",lty=1,col="red")
points(rocENF2[1,],rocENF2[2,],type="l",lty=2,col="red")


mirtissue[1:5,1:5]
mirblood[1:5,1:5]

rntis <- rownames(mirtissue)
rnblood <- rownames(mirblood)
rnis <- intersect(rntis,rnblood)
length(rnis)

mtis <- match(rnis,rntis)
mbl <- match(rnis,rnblood)

mirtissue2 <- mirtissue[mtis,]
mirblood2 <- mirblood[mbl,]

rownames(mirtissue2) == rownames(mirblood2)

wilcoxp <- function(datarow,resp) {
  wilctest <- wilcox.test(mirx ~ resp,data=data.frame(mirx = datarow))
  pval <- wilctest$p.value
  return(pval)
}

sgn  <- function(datarow,resp) {
  #datarow <- mirblood2[1,];resp <- respblood
  return(sign(mean(datarow[resp==1]) - mean(datarow[resp==0])))
}

#Apply Wilcoxon to all data rows
pvalswilcbl <- apply(mirblood2,1,wilcoxp,resp=respblood)
pvalswilctis <- apply(mirtissue2,1,wilcoxp,resp=resptissue)

sgnbl <- apply(mirblood2,1,sgn,resp=respblood)
sgntis <- apply(mirtissue2,1,sgn,resp=resptissue)

plot(log(pvalswilcbl),log(pvalswilctis))
cor(log(pvalswilcbl),log(pvalswilctis),method="spearman")