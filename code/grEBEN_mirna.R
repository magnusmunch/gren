##############################  preamble  #############################
# grMCEM on real data                                                 #
# version: 01                                                         #
# author: Magnus M?nch                                                #
# created: 21-03-2017                                                 #
# last edited: 21-03-2017                                             #
#######################################################################

###############################  notes  ###############################
# 21-03-2017: Results are used in presentation for channel meeting    #
#             2017                                                    #
#######################################################################

##########################  data description  #########################
# A deep sequencing analysis on small non-coding ribonucleic acid     #
# (miRNAseq) was performed on 56 samples (24 women with high-grade    #
# cervical intraepithelial neoplasia (CIN3) and 32 healthy women) for #
# the purpose of finding relevant screening markers for cervical      #
# cancer screening. The next generation sequencing analysis resulted  #
# in 2,576 transcripts. The data was normalized and pre-processed,    #
# rendering 772 transcripts. More detail description of the data sets #
# and the preprocessing step are available on the supplementary       #
# material of this following publication "Better diagnostic           #
# signatures from RNAseq data through use of auxiliary co-data",      #
# Bioinformatics (2017).                                              #
#######################################################################

### paths
path.code <- as.character(ifelse(Sys.info()[1]=="Darwin","/Users/magnusmunch/Documents/PhD/EBEN/code/" ,"~/EBEN/code/"))
path.res <- as.character(ifelse(Sys.info()[1]=="Darwin","/Users/magnusmunch/Documents/PhD/EBEN/results/" ,"~/EBEN/results/"))
path.data <- as.character(ifelse(Sys.info()[1]=="Darwin","/Users/magnusmunch/Documents/PhD/EBEN/data/" ,"~/EBEN/data/"))
path.graph <- "/Users/magnusmunch/Documents/PhD/EBEN/graphs/"

### libraries
library(glmnet)
library(penalized)
library(GRridge)
library(pROC)

### functions
# source grENVB functions
source(paste(path.code, "grVBEM.R", sep=""))

# functions to cross-validate elastic net
cv.pen <- function(x, y, intercept, psel) {
  n <- nrow(x)
  seq.alpha <- seq(0.01, 0.99, length.out=50)
  seq.lam <- seq.df <- seq.cvll <- numeric(length(seq.alpha))
  for(a in 1:length(seq.alpha)) {
    fit <- glmnet(x, y, family="binomial", alpha=seq.alpha[a], dfmax=psel, 
                  standardize=FALSE, intercept=intercept)
    cv.fit <- cv.glmnet(x, y, family="binomial", alpha=seq.alpha[a], 
                        lambda=tail(fit$lambda, n=2L), standardize=FALSE, 
                        intercept=intercept)
    seq.lam[a] <- tail(cv.fit$lambda.min, n=1L)
    seq.df[a] <- tail(cv.fit$nzero, n=1L)
    seq.cvll[a] <- tail(cv.fit$cvm, n=1L)
  }
  lambda1 <- 2*n*seq.alpha*seq.lam
  lambda2 <- n*(1 - seq.alpha)*seq.lam
  
  out <- list(alpha=seq.alpha, lambda=seq.lam, cvll=seq.cvll, df=seq.df)
  return(out)
}

### loading data
load(paste(path.data, "mirsData.RData", sep=""))

# ### creating partitions
parCons <- CreatePartition(mirsData$conservation) # using conservation status as grouping
parAbund <- CreatePartition(rowSums(mirsData$countData), mingr=25, ngroup=10, decreasing=TRUE)

### abundance on top of conservation status
conservation <- rep(1:length(parCons), times=unlist(lapply(parCons, length)))[order(unlist(parCons))]
abundance <- rep(1:length(parAbund), times=unlist(lapply(parAbund, length)))[order(unlist(parAbund))]
x <- apply(t(as.matrix(mirsData$transformedData)), 2, function(x) {
  (x - mean(x))/sd(x)})[, order(conservation, abundance)]
partitions <- list(conservation=conservation[order(conservation, abundance)],
                   abundance=abundance[order(conservation, abundance)])
y <- as.numeric(mirsData$response) - 1
n <- nrow(x)
p <- ncol(x)
m <- rep(1, n)

# number of selected variables
pselseq <- c(5, 10, 15, 20, 50, 100)

# estimating ridge model (fixed for number of selected variables)
cv.ridge <- cv.glmnet(x, y, family="binomial", alpha=0, standardize=FALSE, intercept=TRUE)
lambda2ridge <- 0.5*n*cv.ridge$lambda.min
lambdaridge <- cv.ridge$lambda.min

# determining the fixed folds
set.seed(1010)
nfolds <- n
rest <- n %% nfolds
foldsize <- c(rep(n %/% nfolds + as.numeric(rest!=0), times=rest),
              rep(n %/% nfolds, times=nfolds - rest))
foldid <- sample(rep(1:nfolds, times=foldsize))

# repeating the procedure for different numbers of selected variables
selmat <- matrix(NA, nrow=10*length(pselseq), ncol=3)
aucmat <- matrix(NA, nrow=length(pselseq), ncol=4)
for(cursel in 1:length(pselseq)) {

  psel <- pselseq[cursel]
  print(paste("Psel", psel, sep=" "))

  # estimating global penalty parameters by cross-validation
  cv.en <- cv.pen(x, y, intercept=TRUE, psel=psel)
  alphaglmnet <- cv.en$alpha[which.min(cv.en$cvll)]
  lambdaglmnet <- cv.en$lambda[which.min(cv.en$cvll)]
  lambda1gren <- 2*n*lambdaglmnet*alphaglmnet
  lambda2gren <- n*lambdaglmnet*(1 - alphaglmnet)

  # cross validating AUC
  pred.ridge <- pred.en <- pred.grEBEN <- pred.grridge <- numeric(n)
  for(k in 1:nfolds) {
    print(paste("Fold", k, sep=" "))
    xtrain <- x[foldid!=k, ]
    xtest <- x[foldid==k, ]
    ytrain <- y[foldid!=k]
    ytest <- y[foldid==k]
    mtrain <- m[foldid!=k]
    mtest <- m[foldid==k]
    ntrain <- length(ytrain)
    ntest <- length(ytest)

    # estimating the models
    fit.ridge <- glmnet(xtrain, ytrain, family="binomial", alpha=0, lambda=lambdaridge,
                        standardize=FALSE, intercept=TRUE)
    fit.en <- glmnet(xtrain, ytrain, family="binomial", alpha=alphaglmnet, lambda=lambdaglmnet,
                     standardize=FALSE, intercept=TRUE)
    fit.grEBEN <- grVBEM2(xtrain, ytrain, mtrain, partitions, lambda1=lambda1gren, lambda2=lambda2gren,
                          intercept=TRUE, posterior=FALSE, eps=0.001, maxiter=300, trace=TRUE)
    lambdag <- fit.grEBEN$lambdag$conservation[, fit.grEBEN$nouteriter + 1][conservation[order(conservation, abundance)]]*
      fit.grEBEN$lambdag$abundance[, fit.grEBEN$nouteriter + 1][abundance[order(conservation, abundance)]]
    fit.grEBEN2 <- penalized(ytrain, xtrain, unpenalized=~1, lambda1=0.5*fit.grEBEN$lambda1*sqrt(lambdag),
                             lambda2=0.5*fit.grEBEN$lambda2*lambdag, model="logistic")
    fit.grridge <- grridge(t(xtrain), ytrain, optl=lambda2ridge,
                           list(conservation=CreatePartition(as.factor(partitions$conservation)),
                                abundance=CreatePartition(as.factor(partitions$abundance))),
                           unpenal=~1, trace=FALSE, selectionEN=TRUE, maxsel=psel)

    # checking the number of selected variables
    selmat[(cursel - 1)*10 + k, ] <- c(sum(as.numeric(fit.en$beta)!=0),
                                       sum(fit.grEBEN2@penalized!=0),
                                       length(fit.grridge$resEN$whichEN))

    # getting cross-validated predictions
    pred.ridge[foldid==k] <- predict(fit.ridge, matrix(xtest, nrow=ntest), 
                                     type="response")
    pred.en[foldid==k] <- predict(fit.en, matrix(xtest, nrow=ntest), type="response")
    pred.grEBEN[foldid==k] <- predict(fit.grEBEN2, matrix(xtest, nrow=ntest))
    pred.grridge[foldid==k] <- predict(fit.grridge$predobj$EN, matrix(matrix(xtest, nrow=ntest)[, fit.grridge$resEN$whichEN],
                                                                      nrow=ntest))

  }

  auc.ridge <- pROC::roc(y, pred.ridge)$auc
  auc.en <- pROC::roc(y, pred.en)$auc
  auc.grEBEN <- pROC::roc(y, pred.grEBEN)$auc
  auc.grridge <- pROC::roc(y, pred.grridge)$auc

  aucmat[cursel, ] <- c(ridge=auc.ridge, en=auc.en, grEBEN=auc.grEBEN, grridge=auc.grridge)
  mirna1 <- list(auc=aucmat, sel=selmat)

}

save(mirna1, file=paste(path.res, "grVBEM_mirna_res1.Rdata", sep=""))

# load(paste(path.res, "grVBEM_mirna_res1.Rdata", sep=""))
# truesel <- t(sapply(1:length(pselseq), function(p) {colSums(mirna1$sel[(10*(p - 1) + 1):(10*p), ])/10}))
# colnames(mirna1$auc) <- c("ridge", "enet", "grEBEN", "GRridge")
# colnames(truesel) <- c("enet", "grEBEN", "GRridge")
# 
# plot(pselseq, mirna1$auc[, 1], col=4, type="l", ylim=range(mirna1$auc), xlim=range(truesel),
#      ylab="AUC", xlab="Average number of selected variables")
# lines(truesel[, 2], mirna1$auc[, 3], col=1)
# lines(truesel[, 3], mirna1$auc[, 4], col=2)
# lines(truesel[, 1], mirna1$auc[, 2], col=3)
# legend("topright", legend=c("grEBEN", "GRridge", "enet", "ridge"), col=c(1:4), lty=1)






