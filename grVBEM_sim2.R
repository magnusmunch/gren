##############################  preamble  #############################
# simulations for grVBEM                                              #
# version: 01                                                         #
# author: Magnus M?nch                                                #
# created: 15-03-2017                                                 #
# last edited: 15-03-2017                                             #
#######################################################################

###############################  notes  ###############################
#                                                                     #
#######################################################################

### paths
path.code <- as.character(ifelse(Sys.info()[1]=="Darwin","/Users/magnusmunch/Documents/PhD/EBEN/code/" ,"~/EBEN/code/"))
path.res <- as.character(ifelse(Sys.info()[1]=="Darwin","/Users/magnusmunch/Documents/PhD/EBEN/results/" ,"~/EBEN/results/"))
path.data <- as.character(ifelse(Sys.info()[1]=="Darwin","/Users/magnusmunch/Documents/PhD/EBEN/data/" ,"~/EBEN/data/"))
path.graph <- "/Users/magnusmunch/Documents/PhD/EBEN/graphs/"

### libraries
library(Rcpp)
library(glmnet)
library(penalized)
library(GRridge)
library(mvtnorm)
# library(grpreg)
library(pROC)
library(psych)

### functions
# source grENVB functions
source(paste(path.code, "grVBEM.R", sep=""))

# function to cross-validate elastic net
cv.en <- function(x, y, intercept=intercept) {
  fit.pen <- cv.pen(x, y, intercept=intercept)
  fit.en <- glmnet(x, y, family="binomial", alpha=fit.pen$alpha[which.min(fit.pen$cvll)],
                   lambda=fit.pen$lambda[which.min(fit.pen$cvll)])
  return(fit.en)
}

### the simulation
set.seed(6002)
settings <- cbind(expand.grid(f=c(1.3, 1.6, 2), q=c(0, 0.7, 0.9)), G=10, 
                  pg=100, meanbeta=0.05, nblock=20, rho=0.3, sigma2=1, m=1)
n <- 100
ntest <- 1000
reps <- 50
briermat <- aucmat <- matrix(NA, ncol=7, nrow=reps)
kappamat <- matrix(NA, ncol=3, nrow=reps)
for(s in 1:nrow(settings)) {
  
  # matrix with number of selected variables per group
  nselmat <- matrix(NA, ncol=3*settings$G[s], nrow=reps)
  
  # specifying parameters for current setting
  m <- rep(settings$m[s], n)
  p <- settings$G[s]*settings$pg[s]
  groups <- rep(1:settings$G[s], each=settings$pg[s])
  partition1 <- list(groups=groups)
  partition2 <- list(groups=CreatePartition(as.factor(groups)))
  pblock <- p/settings$nblock[s]
  ubeta <- rep(rev(sapply(0:(settings$G[s] - 1), function(g) {c(settings$f[s]^(-g), 0)})), 
               times=rep(c(settings$pg[s]*settings$q[s], settings$pg[s]*(1 - settings$q[s])), 
                         times=settings$G[s]))
  beta <- ubeta*settings$meanbeta[s]/mean(ubeta)
  bsigma <- matrix(settings$rho[s], ncol=pblock, nrow=pblock)
  diag(bsigma) <- settings$sigma2[s]
  
  for(r in 1:reps) {
    
    print(paste("rep ", r, sep=""))
    # creating the data
    x <- do.call(cbind, replicate(settings$nblock[s], 
                                  rmvnorm(n, mean=rep(0, pblock), 
                                          sigma=bsigma), simplify=FALSE))
    prob <- as.numeric(exp(x %*% beta)/(1 + exp(x %*% beta)))
    y <- rbinom(n, settings$m[s], prob)
    
    # estimating penalty parameters
    cv.ridge <- cv.glmnet(x, y, family="binomial", alpha=0, standardize=FALSE, intercept=TRUE)
    cv.en <- cv.pen(x, y, intercept=TRUE)
    
    lambda2ridge <- 0.5*n*cv.ridge$lambda.min
    lambdaridge <- cv.ridge$lambda.min
    alphaglmnet <- cv.en$alpha[which.min(cv.en$cvll)]
    lambdaglmnet <- cv.en$lambda[which.min(cv.en$cvll)]
    lambda1gren <- 2*n*lambdaglmnet*alphaglmnet
    lambda2gren <- n*lambdaglmnet*(1 - alphaglmnet)

    # fitting the models
    fit.ridge <- glmnet(x, y, family="binomial", alpha=0, lambda=lambdaridge, standardize=FALSE,
                        intercept=TRUE)
    fit.en <- glmnet(x, y, family="binomial", alpha=alphaglmnet, lambda=lambdaglmnet, standardize=FALSE,
                     intercept=TRUE)
    fit.grEBEN <- grVBEM(x, y, m=m, partitions=partition1, lambda1=lambda1gren, lambda2=lambda2gren,
                        intercept=TRUE, eps=0.001, maxiter=500, trace=FALSE)
    grlam1 <- 0.5*lambda1gren*rep(sqrt(fit.grEBEN$lambdag$groups[, fit.grEBEN$nouteriter + 1]), each=settings$pg[s])
    grlam2 <- 0.5*lambda2gren*rep(fit.grEBEN$lambdag$groups[, fit.grEBEN$nouteriter + 1], each=settings$pg[s])
    fit.grEBEN2 <- penalized(y, x, unpenalized=~1, lambda1=grlam1, grlam2, model="logistic")
    fit.grridge <- grridge(t(x), y, partition2, unpenal=~1, optl=lambda2ridge, trace=FALSE,
                           selectionEN=TRUE)
    # fit.grpreg <- cv.grpreg(x, y, group=groups, penalty="grLasso", family="binomial")
    
    # creating test data
    xtest <- do.call(cbind, replicate(settings$nblock[s], 
                                      rmvnorm(ntest, mean=rep(0, pblock), 
                                              sigma=bsigma), simplify=FALSE))
    probtest <- as.numeric(exp(xtest %*% beta)/(1 + exp(xtest %*% beta)))
    ytest <- rbinom(ntest, rep(settings$m[s], ntest), probtest)
    
    # computing predictions
    pred.truth <- as.numeric(exp(xtest %*% beta)/(1 + exp(xtest %*% beta)))
    pred.grridge <- predict.grridge(fit.grridge, t(xtest))[, 2]
    pred.grridgesel <- predict.grridge(fit.grridge, t(xtest))[, 3]
    pred.grEBEN <- as.numeric(exp(cbind(1, xtest) %*% fit.grEBEN$mu)/
                                (1 + exp(cbind(1, xtest) %*% fit.grEBEN$mu)))
    pred.grEBENsel <- as.numeric(exp(cbind(1, xtest) %*% coef(fit.grEBEN2, which="all"))/
                                   (1 + exp(cbind(1, xtest) %*% coef(fit.grEBEN2, which="all"))))
    # pred.grpreg <- predict(fit.grpreg, xtest, type="response")
    pred.ridge <- predict(fit.ridge, xtest, type="response")
    pred.en <- predict(fit.en, xtest, type="response")
    
    # number of selected variables
    nsel.grridge <- sapply(1:settings$G[s], function(g) {
      sum(as.numeric(cut(fit.grridge$resEN$whichEN, seq(0, 1000, 100), labels=1:10))==g)})
    nsel.grEBEN <- sapply(1:settings$G[s], function(g) {
      sum(as.numeric(cut(which(fit.grEBEN2@penalized!=0), seq(0, 1000, 100), labels=1:10))==g)})
    # nsel.grpreg <- sapply(1:settings$G[s], function(g) {
    #   sum(as.numeric(cut(which(coef(fit.grpreg)[-1]!=0), seq(0, 1000, 100), labels=1:10))==g)})
    nsel.en <- sapply(1:settings$G[s], function(g) {
      sum(as.numeric(cut(which(as.numeric(coef(fit.en))[-1]!=0), 
                         seq(0, 1000, 100), labels=1:10))==g)})
    
    # calculating Cohen's kappa
    kappa.grridge <- cohen.kappa(cbind(as.numeric(beta!=0), 
                                       as.numeric(c(1:p) %in% fit.grridge$resEN$whichEN)))$kappa
    kappa.grEBEN <- cohen.kappa(cbind(as.numeric(beta!=0), 
                                      as.numeric(fit.grEBEN2@penalized!=0)))$kappa
    # kappa.grpreg <- cohen.kappa(cbind(as.numeric(beta!=0), 
    #                                    as.numeric(coef(fit.grpreg)[-1]!=0)))$kappa
    kappa.en <- cohen.kappa(cbind(as.numeric(beta!=0), as.numeric(coef(fit.en)[-1, ]!=0)))$kappa
    
    # calculating AUCs on the predictions
    auc.truth <- pROC::roc(ytest, pred.truth)$auc
    auc.grridge <- pROC::roc(ytest, pred.grridge)$auc
    auc.grridgesel <- pROC::roc(ytest, pred.grridgesel)$auc
    auc.grEBEN <- pROC::roc(ytest, pred.grEBEN)$auc
    auc.grEBENsel <- pROC::roc(ytest, pred.grEBENsel)$auc
    # auc.grpreg <- pROC::roc(ytest, pred.grpreg)$auc
    auc.ridge <- pROC::roc(ytest, pred.ridge)$auc
    auc.en <- pROC::roc(ytest, pred.en)$auc
    
    # calculating Brier residuals on the predictions
    brier.truth <- mean((pred.truth - probtest)^2)
    brier.grridge <- mean((pred.grridge - probtest)^2)
    brier.grridgesel <- mean((pred.grridgesel - probtest)^2)
    brier.grEBEN <- mean((pred.grEBEN - probtest)^2)
    brier.grEBENsel <- mean((pred.grEBENsel - probtest)^2)
    # brier.grpreg <- mean((pred.grpreg - probtest)^2)
    brier.ridge <- mean((pred.ridge - probtest)^2)
    brier.en <- mean((pred.en - probtest)^2)
    
    # filling in the matrices with results
    aucmat[r, ] <- c(auc.truth, auc.grridge, auc.grridgesel, auc.grEBEN, auc.grEBENsel, 
                     auc.ridge, auc.en)
    briermat[r, ] <- c(brier.truth, brier.grridge, brier.grridgesel, brier.grEBEN, 
                       brier.grEBENsel, brier.ridge, brier.en)
    nselmat[r, ] <- c(nsel.grridge, nsel.grEBEN, nsel.en)
    kappamat[r, ] <- c(kappa.grridge, kappa.grEBEN, kappa.en)
    
  }
  
  # naming the columns and saving the objects
  colnames(aucmat) <- colnames(briermat) <- c("truth", "grridge", "grridgesel", "grEBEN", "grEBENsel", 
                                              "ridge", "enet")
  colnames(kappamat) <- c("grridge", "grEBEN", "enet")
  colnames(nselmat) <- paste(rep(c("grridge", "grEBEN", "enet"), each=settings$G[s]), 
                             rep(1:settings$G[s], 3), sep="")
  assign(paste("setting", s, sep=""), list(auc=aucmat, brier=briermat, kappa=kappamat,
                                           nsel=nselmat))
  save(list=paste("setting", s, sep=""), file=paste(path.res, "grVBEM_sim2_res2.2_setting", s, ".Rdata", sep=""))
  
}

# load("C:/Users/Magnus/Documents/phd/ENVB/results/grVBEM_res2_setting1.Rdata")
# load("C:/Users/Magnus/Documents/phd/ENVB/results/grVBEM_res2_setting2.Rdata")
# load("C:/Users/Magnus/Documents/phd/ENVB/results/grVBEM_res2_setting3.Rdata")
# load("C:/Users/Magnus/Documents/phd/ENVB/results/grVBEM_res2_setting4.Rdata")
# load("C:/Users/Magnus/Documents/phd/ENVB/results/grVBEM_res2_setting5.Rdata")
# load("C:/Users/Magnus/Documents/phd/ENVB/results/grVBEM_res2_setting6.Rdata")
# load("C:/Users/Magnus/Documents/phd/ENVB/results/grVBEM_res2_setting7.Rdata")
# load("C:/Users/Magnus/Documents/phd/ENVB/results/grVBEM_res2_setting8.Rdata")
# load("C:/Users/Magnus/Documents/phd/ENVB/results/grVBEM_res2_setting9.Rdata")
# 
# tab1 <- cbind(settings[order(settings$f), c(1:2)], 
#               t(sapply(1:nrow(settings), function(s) {
#                 apply(get(paste("setting", s, sep=""))$auc, 2, function(m) {
#                   paste(format(round(median(m), 2), nsmall=2), " (", 
#                         format(round(min(m), 2), nsmall=2), ")", sep="")})})))
# names(tab1)[3:7] <- c("True", "GRridge", "GRen", "Group-lasso", "Ridge")
# write.table(tab1, file=paste(path.tab, "grVBEM_sim2_tab1.csv", sep=""), quote=FALSE, sep=",", 
#             row.names=FALSE, col.names=TRUE)
# 
# boxplot(setting1$auc)
# boxplot(setting2$auc)
# boxplot(setting3$auc)
# boxplot(setting4$auc)
# boxplot(setting5$auc)
# boxplot(setting6$auc)
# boxplot(setting7$auc)
# boxplot(setting8$auc)
# boxplot(setting9$auc)
# 
# boxplot(setting1$brier)
# boxplot(setting2$brier)
# boxplot(setting3$brier)
# boxplot(setting4$brier)
# boxplot(setting5$brier)
# boxplot(setting6$brier)
# boxplot(setting7$brier)
# boxplot(setting8$brier)
# boxplot(setting9$brier)




load(paste(path.res, "grVBEM_sim2_res2.2_setting1.Rdata", sep=""))
load(paste(path.res, "grVBEM_sim2_res2.2_setting2.Rdata", sep=""))
load(paste(path.res, "grVBEM_sim2_res2.2_setting3.Rdata", sep=""))

tab2 <- cbind(settings[c(1:3), c(1, 2)], t(sapply(1:3, function(s) {
  apply(get(paste("setting", s, sep=""))$auc, 2, function(m) {
    paste(format(round(median(m), 2), nsmall=2), " (",
          format(round(min(m), 2), nsmall=2), ")", sep="")})})))[order(settings[c(1:3), ]$f), ]
names(tab2)[3:9] <- c("True", "GRridge", "GRridge+sel", "GReben", "GReben+sel", "Ridge", "Enet")

labels <- c("Truth", "GRridge", "GRridge+sel", "grEBEN", "grEBEN+sel", "ridge", "enet")
boxplot(setting1$auc, xaxt="n")
axis(1, at=1:7, labels=FALSE)
text(seq(1, 7, by=1), par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]), labels=labels, srt=45, pos=1, xpd=TRUE)

boxplot(setting2$auc, xaxt="n")
axis(1, at=1:7, labels=FALSE)
text(seq(1, 7, by=1), par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]), labels=labels, srt=45, pos=1, xpd=TRUE)

boxplot(setting3$auc, xaxt="n")
axis(1, at=1:7, labels=FALSE)
text(seq(1, 7, by=1), par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]), labels=labels, srt=45, pos=1, xpd=TRUE)

# number of selected variables
colSums(cbind(rowSums(setting1$nsel[, c(1:10)]), rowSums(setting1$nsel[, c(11:20)]), 
              rowSums(setting1$nsel[, c(21:30)])))/50

ncol(setting1$nsel)



