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

### reading in data
object.mirna <- read.table(paste(path.data, "gastric_carcinoma_mirna/microarray20170618.txt", sep=""), 
                           header=TRUE, stringsAsFactors=FALSE, skip=1,
                           sep="\t", blank.lines.skip=TRUE)
sample.num <- read.table(paste(path.data, "gastric_carcinoma_mirna/microarray20170618.txt", sep=""), sep="\t", nrows=1)
sample.num <- sample.num[!is.na(sample.num)]
object.mirna <- object.mirna[which(object.mirna$miRBase_human_18th!=""), c(3:15)]
colnames(object.mirna)[c(4:ncol(object.mirna))] <- paste("Sample", sample.num)

object.gastric <- read.table(paste(path.data, "E-MTAB-5854.sdrf.txt", sep=""), header=TRUE, sep="\t",
                             stringsAsFactors=FALSE)

mirsfamily <- read.table(paste(path.data, "miR_Family_Info_edit.txt", sep=""), header=TRUE, 
                         fill=TRUE, as.is=c(1, 2, 4, 5, 7))

### data manipulation
# creating data table
rawdata.mirna <- as.data.frame(t(object.mirna[, c(4:ncol(object.mirna))]))
colnames(rawdata.mirna) <- object.mirna$miRBase_human_18th

# taking average of duplicate miRNAs
data.mirna <- as.data.frame(sapply(unique(names(rawdata.mirna)), function(col) {
  temp <- rowMeans(rawdata.mirna[colnames(rawdata.mirna)==col], na.rm=TRUE); 
  ifelse(is.nan(temp), NA, temp)}))

# selecting miRNAs with no missing values
select.mirna <- data.mirna[, apply(data.mirna, 2, function(var) {!any(is.na(var))})]

# standardize variables
norm.mirna <- apply(select.mirna, 2, function(var) {(var - mean(var))/sd(var)})

# response data, 0 is non-recurrent stage 3, 1 is recurrent stage 3
resp.stage <- as.numeric(object.gastric$Factor.Value.clinical.history.[
  match(object.gastric$Source.Name, rownames(norm.mirna))]=="Recurrent Stage II gastric carcinoma")

### co-data
# conservation status of miRNAs, 0: not conserved, 1: conserved accross vertebrates, 2: conserved
mirsfamily$Family.Conservation.2 <- (mirsfamily$Family.Conservation. > 0)*
  mirsfamily$Family.Conservation.
conservation <- mirsfamily$Family.Conservation.2[match(object.mirna$Target.Sequence[
  match(colnames(norm.mirna), object.mirna$miRBase_human_18th)],
  mirsfamily$Mature.sequence)] + 1

### select data with known co-data
cons.mirna <- norm.mirna[, !is.na(conservation)]
part.cons <- conservation[!is.na(conservation)]
partitions1 <- list(conservation=CreatePartition(as.factor(part.cons)))
partitions2 <- list(conservation=part.cons)

### fit models
enet.pen <- cv.pen(cons.mirna, resp.stage, unpenalized=NULL, intercept=TRUE, psel=NULL)
fit.ridge <- cv.glmnet(cons.mirna, resp.stage, family="binomial", alpha=0, grouped=FALSE,
                       standardize=FALSE)
fit.lasso <- cv.glmnet(cons.mirna, resp.stage, family="binomial", alpha=1, grouped=FALSE,
                       standardize=FALSE)
fit.enet <- glmnet(cons.mirna, resp.stage, family="binomial", standardize=FALSE,
                   alpha=enet.pen$alpha[which.min(enet.pen$cvll)],
                   lambda=enet.pen$lambda[which.min(enet.pen$cvll)])
fit.grridge <- grridge(t(cons.mirna), resp.stage, partitions1)
fit.greben <- grEBEN(cons.mirna, resp.stage, rep(1, length(resp.stage)), partitions=partitions2,
                     lambda1=enet.pen$lambda1bayes[which.min(enet.pen$cvll)],
                     lambda2=enet.pen$lambda2bayes[which.min(enet.pen$cvll)])


### cross-validate models
set.seed(1001)
p <- ncol(cons.mirna)
n <- nrow(cons.mirna)

nfolds <- n
rest <- n %% nfolds
foldid <- sample(rep(1:nfolds, times=round(c(rep(n %/% nfolds + as.numeric(rest!=0), times=rest),
                                             rep(n %/% nfolds, times=nfolds - rest)))))

methods <- c("ridge", "lasso", "enet")
pred <- matrix(NA, ncol=length(methods), nrow=n, dimnames=list(NULL, methods))
for(k in 1:nfolds) {
  print(paste("Fold", k, sep=" "))
  xtrain <- cons.mirna[foldid!=k, ]
  xtest <- cons.mirna[foldid==k, ]
  ytrain <- resp.stage[foldid!=k]
  ytest <- resp.stage[foldid==k]
  ntrain <- length(ytrain)
  ntest <- length(ytest)

  fitcv.ridge <- cv.glmnet(xtrain, ytrain, family="binomial", alpha=0, grouped=FALSE, 
                           standardize=FALSE)
  fitcv.lasso <- cv.glmnet(xtrain, ytrain, family="binomial", alpha=1, grouped=FALSE,
                           standardize=FALSE)
  fitcv.enet <- cv.glmnet(xtrain, ytrain, family="binomial", grouped=FALSE, standardize=FALSE,
                          alpha=enet.pen$alpha[which.min(enet.pen$cvll)])
  fit.grridge <- grridge(t(xtrain))
  
  pred[foldid==k, 1] <- predict(fitcv.ridge, matrix(xtest, ncol=p), s="lambda.min")
  pred[foldid==k, 2] <- predict(fitcv.lasso, matrix(xtest, ncol=p), s="lambda.min")
  pred[foldid==k, 3] <- predict(fitcv.enet, matrix(xtest, ncol=p), s="lambda.min")

}

auc <- setNames(sapply(1:ncol(pred), function(method) {
  return(pROC::roc(resp.stage, pred[, method])$auc)}), methods)





