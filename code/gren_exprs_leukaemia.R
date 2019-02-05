### paths
path.data <- "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/data/"
path.res <- ifelse(as.character(Sys.info()[1])!="Darwin", "~/EBEN/results/",
                   "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/results/")

### loading libraries
library(golubEsets)
library(hu6800.db)
library("KEGGREST")
library(gren) 
library(GRridge)
library(pROC)

### reading in data
data(Golub_Merge)
data(Golub_Test)
data(Golub_Train)

### extracting disease and pathway information from KEGG
# columns(hu6800.db) shows all available info
namemap <- as.data.frame(hu6800ALIAS2PROBE)
entrezmap <- as.data.frame(hu6800ENTREZID)
genenames <- namemap[match(rownames(fData(Golub_Merge)), namemap[, 1]), 2]
entrezids <- entrezmap[match(rownames(fData(Golub_Merge)), entrezmap[, 1]), 2]

# keeping only the known genes
entrezids2 <- entrezids[!is.na(entrezids)]
genenames2 <- genenames[!is.na(genenames)]

# extract number of cancer pathways and diseases in groups of ten
# cancerpathways <- numeric(length(entrezids2))
# cancers <- numeric(length(entrezids2))
# for(i in 1:ceiling(length(entrezids2)/10)) {
#   ind <- c(((i - 1)*10 + 1):(10*i))
#   ind <- ind[ind <= length(entrezids2)]
#   query <- keggGet(paste("hsa:", entrezids2[ind], sep=""))
#   query.names <- sapply(query, function(q) {q$ENTRY})
#   for(j in 1:length(ind)) {
#     if(any(query.names==entrezids2[ind][j])) {
#       paths <- query[[which(query.names==entrezids2[ind][j])[1]]]$PATHWAY
#       if(!is.null(query[[which(query.names==entrezids2[ind][j])[1]]]$DISEASE)) {
#         diseases <- query[[which(query.names==entrezids2[ind][j])[1]]]$DISEASE
#       } else {
#         diseases <- NA
#       }
#       cancers[ind][j] <- sum(grepl("cancer", diseases))
#       cancerpathways[ind][j] <- sum(grepl("cancer", paths))
#     } else {
#       cancerpathways[ind][j] <- NA
#       cancers[ind][j] <- NA
#     }
#   }
# }
# 
# # setting any genes without pathways and diseases to 0
# cancerpathways2 <- replace(cancerpathways, which(is.na(cancerpathways)), 0)
# cancers2 <- replace(cancers, which(is.na(cancers)), 0)
# save(cancerpathways2, cancers2,
#      file=paste(path.data, "gren_exprs_leukaemia_codata.Rdata", sep=""))
load(paste(path.data, "gren_exprs_leukaemia_codata.Rdata", sep=""))

### data preparation
# creating the partitioning
part1 <- as.numeric((cancerpathways2 > 1) + (cancerpathways2 > 0)) + 1
part2 <- as.numeric(cancers2 > 0) + 1

# keeping only the known genes
data <- exprs(Golub_Merge)
datatest <- exprs(Golub_Test)
datatrain <- exprs(Golub_Train)
data2 <- data[!is.na(entrezids), ]
datatest2 <- datatest[!is.na(entrezids), ]
datatrain2 <- datatrain[!is.na(entrezids), ]

# scaling data
x <- scale(t(data2))
xtest <- scale(t(datatest2))
xtrain <- scale(t(datatrain2))
y <- as.numeric(pData(Golub_Merge)$ALL.AML) - 1
ytest <- as.numeric(pData(Golub_Test)$ALL.AML) - 1
ytrain <- as.numeric(pData(Golub_Train)$ALL.AML) - 1

### fitting the models
psel <- c(seq(1, 10, 1), seq(12, 30, 2), seq(35, 50, 5))
train2.gren1 <- gren(xtrain, ytrain, partitions=list(pathways=part2), 
                     alpha=0.05, psel=psel)
train2.gren2 <- gren(xtrain, ytrain, partitions=list(pathways=part2), 
                     alpha=0.5, psel=psel)
train2.gren3 <- gren(xtrain, ytrain, partitions=list(pathways=part2), 
                     alpha=0.95, psel=psel)

pred2.gren1 <- predict(train2.gren1$freq.model$groupreg, xtest, 
                       type="response")
pred2.gren2 <- predict(train2.gren2$freq.model$groupreg, xtest, 
                       type="response")
pred2.gren3 <- predict(train2.gren3$freq.model$groupreg, xtest, 
                       type="response")

pred2.enet1 <- predict(train2.gren1$freq.model$regular, xtest, 
                       type="response")
pred2.enet2 <- predict(train2.gren2$freq.model$regular, xtest, 
                       type="response")
pred2.enet3 <- predict(train2.gren3$freq.model$regular, xtest, 
                       type="response")

auc2.gren1 <- apply(pred2.gren1, 2, function(s) {pROC::auc(ytest, s)})
auc2.gren2 <- apply(pred2.gren2, 2, function(s) {pROC::auc(ytest, s)})
auc2.gren3 <- apply(pred2.gren3, 2, function(s) {pROC::auc(ytest, s)})

auc2.enet1 <- apply(pred2.enet1, 2, function(s) {pROC::auc(ytest, s)})
auc2.enet2 <- apply(pred2.enet2, 2, function(s) {pROC::auc(ytest, s)})
auc2.enet3 <- apply(pred2.enet3, 2, function(s) {pROC::auc(ytest, s)})

psel2.gren1 <- colSums(coef(train2.gren1$freq.model$groupreg)!=0) - 1
psel2.gren2 <- colSums(coef(train2.gren2$freq.model$groupreg)!=0) - 1
psel2.gren3 <- colSums(coef(train2.gren3$freq.model$groupreg)!=0) - 1

psel2.enet1 <- colSums(coef(train2.gren1$freq.model$regular)!=0) - 1
psel2.enet2 <- colSums(coef(train2.gren2$freq.model$regular)!=0) - 1
psel2.enet3 <- colSums(coef(train2.gren3$freq.model$regular)!=0) - 1

ylim <- range(auc2.gren1, auc2.gren2, auc2.gren3, auc2.enet1, auc2.enet2, 
              auc2.enet3)
xlim <- range(psel2.gren1, psel2.gren2, psel2.gren3, psel2.enet1, psel2.enet2, 
              psel2.enet3)
plot(psel2.gren1, auc2.gren1, type="l", xlim=xlim, ylim=ylim)
lines(psel2.gren2, auc2.gren2, col=2)
lines(psel2.gren3, auc2.gren3, col=3)
lines(psel2.enet1, auc2.enet1, col=1, lty=2)
lines(psel2.enet2, auc2.enet2, col=2, lty=2)
lines(psel2.enet3, auc2.enet3, col=3, lty=2)

### fitting the models
psel <- c(seq(1, 10, 1), seq(12, 30, 2), seq(35, 50, 5))
fit2.gren1 <- gren(x, y, partitions=list(pathways=part2), alpha=0.05, psel=psel)
fit2.gren2 <- gren(x, y, partitions=list(pathways=part2), alpha=0.5, psel=psel)
fit2.gren3 <- gren(x, y, partitions=list(pathways=part2), alpha=0.95, psel=psel)

fit2.grridge <- vector("list", length(psel))
invisible(capture.output(
  fit2.grridge[[1]] <- grridge(t(x), y, partitions=list(
    pathways=CreatePartition(as.factor(part2))), selection=TRUE, maxsel=psel[1],
    trace=FALSE, standardizeX=FALSE)))
for(s in 2:length(psel)) {
  invisible(capture.output(
    fit2.grridge[[s]] <- grridge(t(x), y, partitions=list(
      pathways=CreatePartition(as.factor(part2))), selection=TRUE, 
      maxsel=psel[s], optl=fit2.grridge[[1]]$optl, trace=FALSE, 
      standardizeX=FALSE)))
}

### calculating performance by cross-validation
set.seed(2018)
n <- nrow(x)
p <- ncol(x)
nfolds <- n
rest <- n %% nfolds
foldid <- sample(rep(1:nfolds, times=round(c(rep(
  n %/% nfolds + as.numeric(rest!=0), times=rest),
  rep(n %/% nfolds, times=nfolds - rest)))))

pred2 <- list(ridge=numeric(n),
              grridge=matrix(NA, nrow=n, ncol=length(psel)),
              enet1=matrix(NA, nrow=n, ncol=length(psel)),
              enet2=matrix(NA, nrow=n, ncol=length(psel)),
              enet3=matrix(NA, nrow=n, ncol=length(psel)),
              gren1=matrix(NA, nrow=n, ncol=length(psel)),
              gren2=matrix(NA, nrow=n, ncol=length(psel)),
              gren3=matrix(NA, nrow=n, ncol=length(psel)))
psel2 <- list(grridge=matrix(NA, nrow=n, ncol=length(psel)),
              enet1=matrix(NA, nrow=n, ncol=length(psel)),
              enet2=matrix(NA, nrow=n, ncol=length(psel)),
              gren1=matrix(NA, nrow=n, ncol=length(psel)),
              enet3=matrix(NA, nrow=n, ncol=length(psel)),
              gren2=matrix(NA, nrow=n, ncol=length(psel)),
              gren3=matrix(NA, nrow=n, ncol=length(psel)))


for(k in 1:nfolds) {
  set.seed(2018 + k)
  cat(paste("Fold ", k, "\n"))
  
  xtrain <- scale(x[foldid!=k, ])
  xtest <- matrix(x[foldid==k, ], ncol=p, byrow=TRUE)
  ytrain <- y[foldid!=k]
  
  cv2.ridge <- cv.glmnet(xtrain, ytrain, alpha=0, standardize=FALSE)
  
  cv2.grridge <- vector("list", length(psel))
  invisible(capture.output(
    cv2.grridge[[1]] <- grridge(t(xtrain), ytrain, partitions=list(
      pathways=CreatePartition(as.factor(part2))), selection=TRUE, 
      maxsel=psel[1], trace=FALSE, standardizeX=FALSE)))
    for(s in 2:length(psel)) {
      invisible(capture.output(
        cv2.grridge[[s]] <- grridge(t(xtrain), ytrain, partitions=list(
          pathways=CreatePartition(as.factor(part2))), selection=TRUE, 
          maxsel=psel[s], optl=cv2.grridge[[1]]$optl, trace=FALSE, 
          standardizeX=FALSE)))
    }
  
  cv2.gren1 <- gren(xtrain, ytrain, partitions=list(pathways=part2), alpha=0.05, 
                    psel=psel, trace=FALSE)
  cv2.gren2 <- gren(xtrain, ytrain, partitions=list(pathways=part2), alpha=0.5, 
                    psel=psel, trace=FALSE)
  cv2.gren3 <- gren(xtrain, ytrain, partitions=list(pathways=part2), alpha=0.95, 
                    psel=psel, trace=FALSE)
  
  # predictions
  pred2$ridge[foldid==k] <- as.numeric(predict(cv2.ridge, xtest, "lambda.min"))
  
  pred2$grridge[foldid==k, ] <- sapply(cv2.grridge, function(s) {
    predict.grridge(s, t(xtest))[, 3]})
  
  pred2$gren1[foldid==k, ] <- predict(cv2.gren1, xtest, type="groupreg",
                                      s=cv2.gren1$freq.model$groupreg$lambda)
  pred2$gren2[foldid==k, ] <- predict(cv2.gren2, xtest, type="groupreg",
                                      s=cv2.gren2$freq.model$groupreg$lambda)
  pred2$gren3[foldid==k, ] <- predict(cv2.gren3, xtest, type="groupreg",
                                      s=cv2.gren3$freq.model$groupreg$lambda)
  
  pred2$enet1[foldid==k, ] <- predict(cv2.gren1, xtest, type="regular",
                                      s=cv2.gren1$freq.model$regular$lambda)
  pred2$enet2[foldid==k, ] <- predict(cv2.gren2, xtest, type="regular",
                                      s=cv2.gren2$freq.model$regular$lambda)
  pred2$enet3[foldid==k, ] <- predict(cv2.gren3, xtest, type="regular",
                                      s=cv2.gren3$freq.model$regular$lambda)
  
  # number of selected variables
  psel2$grridge[foldid==k, ] <- sapply(cv2.grridge, function(s) {
    length(s$resEN$whichEN)})
  
  psel2$gren1[foldid==k, ] <- cv2.gren1$freq.model$groupreg$df
  psel2$gren2[foldid==k, ] <- cv2.gren2$freq.model$groupreg$df
  psel2$gren3[foldid==k, ] <- cv2.gren3$freq.model$groupreg$df
  
  psel2$enet1[foldid==k, ] <- cv2.gren1$freq.model$regular$df
  psel2$enet2[foldid==k, ] <- cv2.gren2$freq.model$regular$df
  psel2$enet3[foldid==k, ] <- cv2.gren3$freq.model$regular$df
  
  results2 <- list(pred=pred2, psel=psel2)
  save(results2, file=paste(path.res, "gren_exprs_leukaemia_res2.Rdata", sep=""))
}

load(paste(path.res, "gren_exprs_leukaemia_res2.Rdata", sep=""))

str(results2)

plot(apply(results2$psel$enet1, 2, mean),
     apply(results2$pred$enet1, 2, function(s) {pROC::auc(y, s)}), type="l")
lines(apply(results2$psel$gren1, 2, mean),
     apply(results2$pred$gren1, 2, function(s) {pROC::auc(y, s)}), col=2)


# ### calculating predictive performance on test set
# cv2.gren1 <- gren(xtrain, ytrain, partitions=list(pathways=part2), alpha=0.5)
# cv2.grridge <- grridge(t(xtrain), as.factor(ytrain), 
#                        list(pathways=CreatePartition(as.factor(part2))))
# cv2.test1 <- glmnet(xtrain, ytrain, family="binomial", alpha=0.5)
# pred2.test1 <- predict(cv2.test1, xtest, type="response")
# 
# pred2.gren1 <- predict(cv2.gren1$freq.model$groupreg, xtest, type="response")
# pred2.enet1 <- predict(cv2.gren1$freq.model$regular, xtest, type="response")
# pred2.grridge <- predict.grridge(cv2.grridge, t(xtest))[, 2]
# 
# psel2.gren1 <- cv2.gren1$freq.model$groupreg$df
# psel2.enet1 <- cv2.gren1$freq.model$regular$df
# psel2.test1 <- cv2.test1$df
# 
# pred2.grridge <- predict.grridge(cv2.grridge, t(xtest))[, 2]
# pred2.ridge <- predict.grridge(cv2.grridge, t(xtest))[, 1]
# 
# auc2.gren1 <- apply(pred2.gren1, 2, function(s) {pROC::auc(ytest, s)})
# auc2.enet1 <- apply(pred2.enet1, 2, function(s) {pROC::auc(ytest, s)})
# auc2.grridge <- pROC::auc(ytest, pred2.grridge)
# auc2.test1 <- apply(pred2.test1, 2, function(s) {pROC::auc(ytest, s)})
# 
# plot(psel2.gren1, auc2.gren1, type="l")
# lines(psel2.enet1, auc2.enet1, col=2)
# lines(psel2.test1, auc2.test1, col=3)
