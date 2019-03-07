#!/usr/bin/env Rscript

### installation of gren if updated on github
if(substr(system('git log -n 1 --format="%h %aN %s %ad"', intern=TRUE), 1, 7)!=
   substr(packageDescription("gren")$GithubSHA1, 1, 7)) {
  if(!("devtools" %in% installed.packages())) {
    install.packages("devtools")
  }
  library(devtools)
  install_github("magnusmunch/gren/rpackage", local=FALSE,
                 auth_token=Sys.getenv("GITHUB_PAT"))
}

### libraries
library(gren)
library(GRridge)
library(grpreg)
library(SGL)

### load data
load("data/mirsData.RData")

### filter out the features with less than 2 unique counts
unique.counts <- apply(t(as.matrix(mirsData$countData)), 2, function(s) {
  length(unique(s))})
select <- which(unique.counts > 2)
filtered.micrornas <- t(as.matrix(mirsData$countData))[, select]

### transform the data
micrornas <- sqrt(filtered.micrornas + 3/8) - sqrt(3/8)
cervical <- as.numeric(as.character(mirsData$response)=="CIN3")

### conservation status of the filtered microRNAs
conservation <- (as.numeric(levels(mirsData$conservation))[
  mirsData$conservation] + 1)[select]
# 1 is not-conserved, 2 is conserved across mammals, 3 broadly conserved across
# most vertebrates

################################### model 1 ####################################
### fitting the models
set.seed(2019)
y <- cervical
x <- scale(micrornas)
part <- conservation
p <- ncol(x)
n <- nrow(x)

csel <- c(seq(1, 5, 1), seq(7, 15, 2), seq(18, 30, 3), seq(34, 50, 4))

fit.gren1 <- gren(x, y, partitions=list(part=part), alpha=0.05, 
                  standardize=TRUE, trace=FALSE, psel=csel)
fit.gren2 <- gren(x, y, partitions=list(part=part), alpha=0.5, 
                  standardize=TRUE, trace=FALSE, psel=csel)
fit.gren3 <- gren(x, y, partitions=list(part=part), alpha=0.95, 
                  standardize=TRUE, trace=FALSE, psel=csel)

fit.grridge <- grridge(t(x), y, list(part=split(1:p, part)))

save(fit.grridge, fit.gren1, fit.gren2, fit.gren3, 
     file="results/micrornaseq_cervical_cancer_fit1.Rdata")

nfolds <- 5
foldid <- sample(rep(1:nfolds, times=round(c(rep(
  n %/% nfolds + as.numeric((n %% nfolds)!=0), times=n %% nfolds),
  rep(n %/% nfolds, times=nfolds - n %% nfolds)))))

pred <- as.data.frame(matrix(NA, nrow=n, ncol=6*length(csel) + 2))
psel <- as.data.frame(matrix(NA, nrow=nfolds, ncol=6*length(csel) + 2))
colnames(pred) <- colnames(psel) <- 
  c("ridge", "grridge", paste0("gren1.psel", csel),
    paste0("gren2.psel", csel), paste0("gren3.psel", csel),
    paste0("enet1.psel", csel), paste0("enet2.psel", csel), 
    paste0("enet3.psel", csel))

for(k in sort(unique(foldid))) {
  cat(paste("fold ", k, "\n"))

  xtrain <- scale(matrix(x[foldid!=k, ], ncol=p))
  # find the constant features in the test data and set them to zero
  xtest <- scale(matrix(x[foldid==k, ], ncol=p))
  xtest[, apply(matrix(x[foldid==k, ], ncol=p), 2, sd)==0] <- 0
  ytrain <- y[foldid!=k]
  
  cv.gren1 <- gren(xtrain, ytrain, partitions=list(part=part), alpha=0.05, 
                   standardize=TRUE, trace=FALSE, psel=csel)
  cv.gren2 <- gren(xtrain, ytrain, partitions=list(part=part), alpha=0.5, 
                   standardize=TRUE, trace=FALSE, psel=csel)
  cv.gren3 <- gren(xtrain, ytrain, partitions=list(part=part), alpha=0.95, 
                   standardize=TRUE, trace=FALSE, psel=csel)
  
  cv.grridge <-  grridge(t(xtrain), ytrain, list(part=split(1:p, part)))

  pred$ridge[foldid==k] <- predict.grridge(cv.grridge, t(xtest))[, 1]
  pred$grridge[foldid==k] <- predict.grridge(cv.grridge, t(xtest))[, 2]
  
  pred[foldid==k, grep("gren1", colnames(pred))] <- 
    predict(cv.gren1, xtest, type="groupreg")
  pred[foldid==k, grep("gren2", colnames(pred))] <- 
    predict(cv.gren2, xtest, type="groupreg")
  pred[foldid==k, grep("gren3", colnames(pred))] <- 
    predict(cv.gren3, xtest, type="groupreg")
  
  pred[foldid==k, grep("enet1", colnames(pred))] <- 
    predict(cv.gren1, xtest, type="regular")
  pred[foldid==k, grep("enet2", colnames(pred))] <- 
    predict(cv.gren2, xtest, type="regular")
  pred[foldid==k, grep("enet3", colnames(pred))] <- 
    predict(cv.gren3, xtest, type="regular")

  psel$ridge <- psel$grridge <- p
  
  psel[k, grep("gren1", colnames(psel))] <- cv.gren1$freq.model$groupreg$df
  psel[k, grep("gren2", colnames(psel))] <- cv.gren2$freq.model$groupreg$df
  psel[k, grep("gren3", colnames(psel))] <- cv.gren3$freq.model$groupreg$df
  
  psel[k, grep("enet1", colnames(psel))] <- cv.gren1$freq.model$regular$df
  psel[k, grep("enet2", colnames(psel))] <- cv.gren2$freq.model$regular$df
  psel[k, grep("enet3", colnames(psel))] <- cv.gren3$freq.model$regular$df
  
}

auc <- apply(pred, 2, function(s) {pROC::auc(y, s)})
briers <- apply(pred, 2, function(m) {
  1 - mean((m - y)^2)/mean((mean(y) - y)^2)})
res <- rbind(pred, psel, auc, briers)
rownames(res) <- c(paste0("pred", c(1:length(y))), paste0("psel", c(1:nfolds)),
                   paste0("auc", 1), paste0("briers", 1))
write.table(res, file="results/micrornaseq_cervical_cancer_res1.csv")
