path.data <- as.character(ifelse(Sys.info()[1]=="Darwin", "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/data/" , "~/EBEN/data/"))
path.res <- ifelse(as.character(Sys.info()[1])!="Darwin", "~/EBEN/results/",
                   "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/results/")

library(gren)
library(pROC)
library(GRridge)
library(grpreg)
library(SGL)
library(Hmisc)
library(Biobase)

# colors <- bpy.colors(5)[-c(1, 5)]

load(paste(path.data, "ESetMbolCSFPR2.Rdata", sep=""))

x <- scale(t(exprs(ESetMbolCSFPR2)))
y <- as.numeric(pData(ESetMbolCSFPR2)$D_diag_name) - 1
colnames(pData(ESetMbolCSFPR2))
pData(ESetMbolCSFPR2)$APOE=="E4YES"


part1 <- replace(fData(ESetMbolCSFPR2)$PlatformCode, 
                 fData(ESetMbolCSFPR2)$PlatformCode==5, 4)
part2 <- fData(ESetMbolCSFPR2)$PlatformCode
part3 <- as.numeric(cut2(fData(ESetMbolCSFPR2)$RSDqc, g=3))
part4 <- as.numeric(cut2(fData(ESetMbolCSFPR2)$RSDqc, 0.06))

labels1 <- replace(levels(fData(ESetMbolCSFPR2)$Platform)[-5],
                   levels(fData(ESetMbolCSFPR2)$Platform)[-5]==
                     "Oxidative Stress high pH", "Oxidate Stress")
labels2 <- levels(fData(ESetMbolCSFPR2)$Platform)
labels3 <- levels(cut2(fData(ESetMbolCSFPR2)$RSDqc, g=3))
labels4 <- levels(cut2(fData(ESetMbolCSFPR2)$RSDqc, 0.06))

n <- nrow(x)
p <- ncol(x)
# set.seed(2016)
# idtrain <- sample(c(1:n), floor(n/2))
# xtrain <- x[idtrain, ]
# xtest <- x[-idtrain, ]
# ytrain <- y[idtrain]
# ytest <- y[-idtrain]
# 
# fit1.gren1 <- gren(xtrain, ytrain, partitions=list(platform=part1), alpha=0.05)
# fit1.gren2 <- gren(xtrain, ytrain, partitions=list(platform=part1), alpha=0.5)
# fit1.gren3 <- gren(xtrain, ytrain, partitions=list(platform=part1), alpha=0.95)
# 
# pred1.gren1 <- predict(fit1.gren1, xtest, type="groupreg",
#                        s=fit1.gren1$freq.model$groupreg$lambda)
# pred1.gren2 <- predict(fit1.gren2, xtest, type="groupreg",
#                        s=fit1.gren2$freq.model$groupreg$lambda)
# pred1.gren3 <- predict(fit1.gren3, xtest, type="groupreg",
#                        s=fit1.gren3$freq.model$groupreg$lambda)
# 
# pred1.enet1 <- predict(fit1.gren1, xtest, type="regular",
#                        s=fit1.gren1$freq.model$regular$lambda)
# pred1.enet2 <- predict(fit1.gren2, xtest, type="regular",
#                        s=fit1.gren2$freq.model$regular$lambda)
# pred1.enet3 <- predict(fit1.gren3, xtest, type="regular",
#                        s=fit1.gren3$freq.model$regular$lambda)
# 
# 
# 
# psel1.gren1 <- fit1.gren1$freq.model$groupreg$df
# psel1.gren2 <- fit1.gren2$freq.model$groupreg$df
# psel1.gren3 <- fit1.gren3$freq.model$groupreg$df
# 
# psel1.enet1 <- fit1.gren1$freq.model$regular$df
# psel1.enet2 <- fit1.gren2$freq.model$regular$df
# psel1.enet3 <- fit1.gren3$freq.model$regular$df
# 
# auc1.gren1 <- apply(pred1.gren1, 2, function(pred) {pROC::roc(ytest, pred)$auc})
# auc1.gren2 <- apply(pred1.gren2, 2, function(pred) {pROC::roc(ytest, pred)$auc})
# auc1.gren3 <- apply(pred1.gren3, 2, function(pred) {pROC::roc(ytest, pred)$auc})
# 
# auc1.enet1 <- apply(pred1.enet1, 2, function(pred) {pROC::roc(ytest, pred)$auc})
# auc1.enet2 <- apply(pred1.enet2, 2, function(pred) {pROC::roc(ytest, pred)$auc})
# auc1.enet3 <- apply(pred1.enet3, 2, function(pred) {pROC::roc(ytest, pred)$auc})
# 
# results1 <- list(idtrain=idtrain,
#                  pred=list(enet1=pred1.enet1, enet2=pred1.enet2, 
#                            enet3=pred1.enet3, gren1=pred1.gren1, 
#                            gren2=pred1.gren2, gren3=pred1.gren3),
#                  psel=list(enet1=psel1.enet1, enet2=psel1.enet2, 
#                            enet3=psel1.enet3, gren1=psel1.gren1, 
#                            gren2=psel1.gren2, gren3=psel1.gren3),
#                  auc=list(enet1=auc1.enet1, enet2=auc1.enet2, enet3=auc1.enet3, 
#                           gren1=auc1.gren1, gren2=auc1.gren2, gren3=auc1.gren3))
# 
# plot(sort(results1$psel$enet1), results1$auc$enet1[order(results1$psel$enet1)], 
#      type="l", xlim=range(results1$psel), ylim=range(results1$auc), lty=2,
#      col=colors[1])
# lines(sort(results1$psel$enet2), results1$auc$enet2[order(results1$psel$enet2)],
#       lty=2, col=colors[2])
# lines(sort(results1$psel$enet3), results1$auc$enet3[order(results1$psel$enet2)],
#       lty=2, col=colors[3])
# 
# lines(sort(results1$psel$gren1), results1$auc$gren1[order(results1$psel$gren1)],
#      lty=1, col=colors[1])
# lines(sort(results1$psel$gren2), results1$auc$gren2[order(results1$psel$gren2)],
#       lty=1, col=colors[2])
# lines(sort(results1$psel$gren3), results1$auc$gren3[order(results1$psel$gren3)],
#       lty=1, col=colors[3])
# 
# 
# #####
# fit2.gren1 <- gren(xtrain, ytrain, partitions=list(platform=part2), alpha=0.05)
# fit2.gren2 <- gren(xtrain, ytrain, partitions=list(platform=part2), alpha=0.5)
# fit2.gren3 <- gren(xtrain, ytrain, partitions=list(platform=part2), alpha=0.95)
# 
# pred2.gren1 <- predict(fit2.gren1, xtest, type="groupreg",
#                        s=fit2.gren1$freq.model$groupreg$lambda)
# pred2.gren2 <- predict(fit2.gren2, xtest, type="groupreg",
#                        s=fit2.gren2$freq.model$groupreg$lambda)
# pred2.gren3 <- predict(fit2.gren3, xtest, type="groupreg",
#                        s=fit2.gren3$freq.model$groupreg$lambda)
# 
# pred2.enet1 <- predict(fit2.gren1, xtest, type="regular",
#                        s=fit2.gren1$freq.model$regular$lambda)
# pred2.enet2 <- predict(fit2.gren2, xtest, type="regular",
#                        s=fit2.gren2$freq.model$regular$lambda)
# pred2.enet3 <- predict(fit2.gren3, xtest, type="regular",
#                        s=fit2.gren3$freq.model$regular$lambda)
# 
# psel2.gren1 <- fit2.gren1$freq.model$groupreg$df
# psel2.gren2 <- fit2.gren2$freq.model$groupreg$df
# psel2.gren3 <- fit2.gren3$freq.model$groupreg$df
# 
# psel2.enet1 <- fit2.gren1$freq.model$regular$df
# psel2.enet2 <- fit2.gren2$freq.model$regular$df
# psel2.enet3 <- fit2.gren3$freq.model$regular$df
# 
# auc2.gren1 <- apply(pred2.gren1, 2, function(pred) {pROC::roc(ytest, pred)$auc})
# auc2.gren2 <- apply(pred2.gren2, 2, function(pred) {pROC::roc(ytest, pred)$auc})
# auc2.gren3 <- apply(pred2.gren3, 2, function(pred) {pROC::roc(ytest, pred)$auc})
# 
# auc2.enet1 <- apply(pred2.enet1, 2, function(pred) {pROC::roc(ytest, pred)$auc})
# auc2.enet2 <- apply(pred2.enet2, 2, function(pred) {pROC::roc(ytest, pred)$auc})
# auc2.enet3 <- apply(pred2.enet3, 2, function(pred) {pROC::roc(ytest, pred)$auc})
# 
# results2 <- list(idtrain=idtrain,
#                  pred=list(enet1=pred2.enet1, enet2=pred2.enet2, 
#                            enet3=pred2.enet3, gren1=pred2.gren1, 
#                            gren2=pred2.gren2, gren3=pred2.gren3),
#                  psel=list(enet1=psel2.enet1, enet2=psel2.enet2, 
#                            enet3=psel2.enet3, gren1=psel2.gren1, 
#                            gren2=psel2.gren2, gren3=psel2.gren3),
#                  auc=list(enet1=auc2.enet1, enet2=auc2.enet2, enet3=auc2.enet3, 
#                           gren1=auc2.gren1, gren2=auc2.gren2, gren3=auc2.gren3))
# 
# plot(sort(results2$psel$enet1), results2$auc$enet1[order(results2$psel$enet1)], 
#      type="l", xlim=range(results2$psel), ylim=range(results2$auc), lty=2,
#      col=colors[1])
# lines(sort(results2$psel$enet2), results2$auc$enet2[order(results2$psel$enet2)],
#       lty=2, col=colors[2])
# lines(sort(results2$psel$enet3), results2$auc$enet3[order(results2$psel$enet2)],
#       lty=2, col=colors[3])
# 
# lines(sort(results2$psel$gren1), results2$auc$gren1[order(results2$psel$gren1)],
#       lty=1, col=colors[1])
# lines(sort(results2$psel$gren2), results2$auc$gren2[order(results2$psel$gren2)],
#       lty=1, col=colors[2])
# lines(sort(results2$psel$gren3), results2$auc$gren3[order(results2$psel$gren3)],
#       lty=1, col=colors[3])





# ####
# fit3.gren1 <- gren(xtrain, ytrain, partitions=list(RSDqc=part3), alpha=0.05)
# fit3.gren2 <- gren(xtrain, ytrain, partitions=list(RSDqc=part3), alpha=0.5)
# fit3.gren3 <- gren(xtrain, ytrain, partitions=list(RSDqc=part3), alpha=0.95)
# 


psel <- c(seq(1, 5, 1), seq(7, 15, 2), seq(20, 40, 5), seq(50, 90, 10), seq(110, 150, 20))


set.seed(2018)
fit3.ridge <- cv.glmnet(x, y, alpha=0, standardize=FALSE)
invisible(capture.output(fit3.grridge <- grridge(t(x), y, partitions=list(
  RSDqc=CreatePartition(as.factor(part3))), trace=FALSE, standardizeX=FALSE)))
fit3.gren1 <- gren(x, y, partitions=list(RSDqc=part3), alpha=0.05, psel=psel,
                   trace=FALSE)
fit3.gren2 <- gren(x, y, partitions=list(RSDqc=part3), alpha=0.5, psel=psel,
                   trace=FALSE)
fit3.gren3 <- gren(x, y, partitions=list(RSDqc=part3), alpha=0.95, psel=psel,
                   trace=FALSE)
fit3.sglasso1 <- SGL(list(x=x, y=y), part3, type="logit", alpha=0.05, 
                     standardize=FALSE, nlam=length(psel))
fit3.sglasso2 <- SGL(list(x=x, y=y), part3, type="logit", alpha=0.5, 
                     standardize=FALSE, nlam=length(psel))
fit3.sglasso3 <- SGL(list(x=x, y=y), part3, type="logit", alpha=0.95, 
                     standardize=FALSE, nlam=length(psel))
fit3.cmcp1 <- grpreg(x, y, part3, penalty="cMCP", alpha=0.05, 
                     nlambda=length(psel))
fit3.cmcp2 <- grpreg(x, y, part3, penalty="cMCP", alpha=0.5, 
                     nlambda=length(psel))
fit3.cmcp3 <- grpreg(x, y, part3, penalty="cMCP", alpha=0.95, 
                     nlambda=length(psel))

nfolds <- n
rest <- n %% nfolds
foldid <- sample(rep(1:nfolds, times=round(c(rep(
  n %/% nfolds + as.numeric(rest!=0), times=rest),
  rep(n %/% nfolds, times=nfolds - rest)))))

pred3 <- list(ridge=numeric(n),
              grridge=matrix(NA, nrow=n, ncol=length(psel)),
              enet1=matrix(NA, nrow=n, ncol=length(psel)),
              enet2=matrix(NA, nrow=n, ncol=length(psel)),
              enet3=matrix(NA, nrow=n, ncol=length(psel)),
              gren1=matrix(NA, nrow=n, ncol=length(psel)),
              gren2=matrix(NA, nrow=n, ncol=length(psel)),
              gren3=matrix(NA, nrow=n, ncol=length(psel)),
              sglasso1=vector("list", n),
              sglasso2=vector("list", n),
              sglasso3=vector("list", n),
              cmcp1=vector("list", n),
              cmcp2=vector("list", n),
              cmcp3=vector("list", n),
              gelasso2=vector("list", n),
              gelasso1=vector("list", n),
              gelasso3=vector("list", n))
psel3 <- list(grridge=matrix(NA, nrow=n, ncol=length(psel)),
              enet1=matrix(NA, nrow=n, ncol=length(psel)),
              enet2=matrix(NA, nrow=n, ncol=length(psel)),
              gren1=matrix(NA, nrow=n, ncol=length(psel)),
              enet3=matrix(NA, nrow=n, ncol=length(psel)),
              gren2=matrix(NA, nrow=n, ncol=length(psel)),
              gren3=matrix(NA, nrow=n, ncol=length(psel)),
              sglasso1=vector("list", n),
              sglasso2=vector("list", n),
              sglasso3=vector("list", n),
              cmcp1=vector("list", n),
              cmcp2=vector("list", n),
              cmcp3=vector("list", n),
              gelasso1=vector("list", n),
              gelasso2=vector("list", n),
              gelasso3=vector("list", n))


for(k in 1:nfolds) {
  cat(paste("Fold ", k, "\n"))
  
  xtrain <- scale(x[foldid!=k, ])
  xtest <- scale(matrix(x[foldid==k, ], ncol=p, byrow=TRUE))
  ytrain <- y[foldid!=k]
  
  cv3.ridge <- cv.glmnet(xtrain, ytrain, alpha=0, standardize=FALSE)
  
  cv3.grridge <- vector("list", length(psel))
  invisible(capture.output(
    cv3.grridge[[1]] <- grridge(t(xtrain), ytrain, partitions=list(
      RSDqc=CreatePartition(as.factor(part3))), selection=TRUE, maxsel=psel[1],
      trace=FALSE, standardizeX=FALSE)))
  for(s in 2:length(psel)) {
    invisible(capture.output(
      cv3.grridge[[s]] <- grridge(t(xtrain), ytrain, partitions=list(
        RSDqc=CreatePartition(as.factor(part3))), selection=TRUE,
        maxsel=psel[s], optl=cv3.grridge[[1]]$optl, trace=FALSE,
        standardizeX=FALSE)))
  }
  
  cv3.gren1 <- gren(xtrain, ytrain, partitions=list(RSDqc=part3), alpha=0.05,
                     psel=psel, trace=FALSE)
  cv3.gren2 <- gren(xtrain, ytrain, partitions=list(RSDqc=part3), alpha=0.5,
                     psel=psel, trace=FALSE)
  cv3.gren3 <- gren(xtrain, ytrain, partitions=list(RSDqc=part3), alpha=0.95,
                     psel=psel, trace=FALSE)
  
  cv3.sglasso1 <- SGL(list(x=xtrain, y=ytrain), part3, type="logit", 
                       alpha=0.05, standardize=FALSE, nlam=100)
  cv3.sglasso2 <- SGL(list(x=xtrain, y=ytrain), part3, type="logit", alpha=0.5,
                       standardize=FALSE, nlam=100)
  cv3.sglasso3 <- SGL(list(x=xtrain, y=ytrain), part3, type="logit", 
                       alpha=0.95, standardize=FALSE, nlam=100)
  
  cv3.cmcp1 <- grpreg(xtrain, ytrain, part3, penalty="cMCP", alpha=0.05)
  cv3.cmcp2 <- grpreg(xtrain, ytrain, part3, penalty="cMCP", alpha=0.5)
  cv3.cmcp3 <- grpreg(xtrain, ytrain, part3, penalty="cMCP", alpha=0.95)
  
  cv3.gelasso1 <- grpreg(xtrain, ytrain, part3, penalty="gel", alpha=0.05)
  cv3.gelasso2 <- grpreg(xtrain, ytrain, part3, penalty="gel", alpha=0.5)
  cv3.gelasso3 <- grpreg(xtrain, ytrain, part3, penalty="gel", alpha=0.95)
  
  # predictions
  pred3$ridge[foldid==k] <- as.numeric(predict(cv3.ridge, xtest, "lambda.min"))
  
  pred3$grridge[foldid==k, ] <- sapply(cv3.grridge, function(s) {
    predict.grridge(s, t(xtest))[, 3]})
  
  pred3$gren1[foldid==k, ] <- predict(cv3.gren1, xtest, type="groupreg",
                                      s=cv3.gren1$freq.model$groupreg$lambda)
  pred3$gren2[foldid==k, ] <- predict(cv3.gren2, xtest, type="groupreg",
                                      s=cv3.gren2$freq.model$groupreg$lambda)
  pred3$gren3[foldid==k, ] <- predict(cv3.gren3, xtest, type="groupreg",
                                      s=cv3.gren3$freq.model$groupreg$lambda)
  
  pred3$enet1[foldid==k, ] <- predict(cv3.gren1, xtest, type="regular",
                                      s=cv3.gren1$freq.model$regular$lambda)
  pred3$enet2[foldid==k, ] <- predict(cv3.gren2, xtest, type="regular",
                                      s=cv3.gren2$freq.model$regular$lambda)
  pred3$enet3[foldid==k, ] <- predict(cv3.gren3, xtest, type="regular",
                                      s=cv3.gren3$freq.model$regular$lambda)
  
  pred3$sglasso1[[which(foldid==k)]] <- predictSGL(cv3.sglasso1, xtest)
  pred3$sglasso2[[which(foldid==k)]] <- predictSGL(cv3.sglasso2, xtest)
  pred3$sglasso3[[which(foldid==k)]] <- predictSGL(cv3.sglasso3, xtest)
  
  pred3$cmcp1[[which(foldid==k)]] <- predict(cv3.cmcp1, xtest)
  pred3$cmcp2[[which(foldid==k)]] <- predict(cv3.cmcp2, xtest)
  pred3$cmcp3[[which(foldid==k)]] <- predict(cv3.cmcp3, xtest)
  
  pred3$gelasso1[[which(foldid==k)]] <- predict(cv3.gelasso1, xtest)
  pred3$gelasso2[[which(foldid==k)]] <- predict(cv3.gelasso2, xtest)
  pred3$gelasso3[[which(foldid==k)]] <- predict(cv3.gelasso3, xtest)
  
  # number of selected variables
  psel3$grridge[foldid==k, ] <- sapply(cv3.grridge, function(s) {
    length(s$resEN$whichEN)})
  
  psel3$gren1[foldid==k, ] <- cv3.gren1$freq.model$groupreg$df
  psel3$gren2[foldid==k, ] <- cv3.gren2$freq.model$groupreg$df
  psel3$gren3[foldid==k, ] <- cv3.gren3$freq.model$groupreg$df
  
  psel3$enet1[foldid==k, ] <- cv3.gren1$freq.model$regular$df
  psel3$enet2[foldid==k, ] <- cv3.gren2$freq.model$regular$df
  psel3$enet3[foldid==k, ] <- cv3.gren3$freq.model$regular$df
  
  psel3$sglasso1[[which(foldid==k)]] <- apply(
    cv3.sglasso1$beta, 2, function(b) {sum(b!=0)})
  psel3$sglasso2[[which(foldid==k)]] <- apply(
    cv3.sglasso2$beta, 2, function(b) {sum(b!=0)})
  psel3$sglasso3[[which(foldid==k)]] <- apply(
    cv3.sglasso3$beta, 2, function(b) {sum(b!=0)})
  
  psel3$cmcp1[[which(foldid==k)]] <- apply(
    cv3.cmcp1$beta, 2, function(b) {sum(b!=0)})
  psel3$cmcp2[[which(foldid==k)]] <- apply(
    cv3.cmcp2$beta, 2, function(b) {sum(b!=0)})
  psel3$cmcp3[[which(foldid==k)]] <- apply(
    cv3.cmcp3$beta, 2, function(b) {sum(b!=0)})
  
  psel3$gelasso1[[which(foldid==k)]] <- apply(
    cv3.gelasso1$beta, 2, function(b) {sum(b!=0)})
  psel3$gelasso2[[which(foldid==k)]] <- apply(
    cv3.gelasso2$beta, 2, function(b) {sum(b!=0)})
  psel3$gelasso3[[which(foldid==k)]] <- apply(
    cv3.gelasso3$beta, 2, function(b) {sum(b!=0)})
  
  results3 <- list(pred=pred3, psel=psel3)
  save(results3, file=paste(path.res, "gren_metabol_alzheimer_res1.Rdata", sep=""))
}



path.res <- ifelse(as.character(Sys.info()[1])!="Darwin", "~/EBEN/results/",
                   "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/results/")
path.graph <- "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/graphs/"
path.data <- as.character(ifelse(Sys.info()[1]=="Darwin", "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/data/" , "~/EBEN/data/"))

library(sp)

load(paste(path.res, "gren_metabol_alzheimer_res1.Rdata", sep=""))
load(paste(path.data, "ESetMbolCSFPR2.Rdata", sep=""))
y <- as.numeric(pData(ESetMbolCSFPR2)$D_diag_name) - 1

auc3.ridge <- pROC::auc(y, results3$pred$ridge)
auc3.grridge <- apply(results3$pred$grridge, 2, function(s) {pROC::auc(y, s)})
auc3.gren1 <- apply(results3$pred$gren1, 2, function(s) {pROC::auc(y, s)})
auc3.gren2 <- apply(results3$pred$gren2, 2, function(s) {pROC::auc(y, s)})
auc3.gren3 <- apply(results3$pred$gren3, 2, function(s) {pROC::auc(y, s)})
auc3.enet1 <- apply(results3$pred$enet1, 2, function(s) {pROC::auc(y, s)})
auc3.enet2 <- apply(results3$pred$enet2, 2, function(s) {pROC::auc(y, s)})
auc3.enet3 <- apply(results3$pred$enet3, 2, function(s) {pROC::auc(y, s)})

const <- sum((y - mean(y))^2)
briers3.ridge <- 1 - sum((y - results3$pred$ridge)^2)/const
briers3.grridge <- apply(results3$pred$grridge, 2, function(s) {
  1 - sum((y - s)^2)/const})
briers3.gren1 <- apply(results3$pred$gren1, 2, function(s) {
  1 - sum((y - s)^2)/const})
briers3.gren2 <- apply(results3$pred$gren2, 2, function(s) {
  1 - sum((y - s)^2)/const})
briers3.gren3 <- apply(results3$pred$gren3, 2, function(s) {
  1 - sum((y - s)^2)/const})
briers3.enet1 <- apply(results3$pred$enet1, 2, function(s) {
  1 - sum((y - s)^2)/const})
briers3.enet2 <- apply(results3$pred$enet2, 2, function(s) {
  1 - sum((y - s)^2)/const})
briers3.enet3 <- apply(results3$pred$enet3, 2, function(s) {
  1 - sum((y - s)^2)/const})

psel3.grridge <- colMeans(results3$psel$grridge)
psel3.gren1 <- colMeans(results3$psel$gren1)
psel3.gren2 <- colMeans(results3$psel$gren2)
psel3.gren3 <- colMeans(results3$psel$gren3)
psel3.enet1 <- colMeans(results3$psel$enet1)
psel3.enet2 <- colMeans(results3$psel$enet2)
psel3.enet3 <- colMeans(results3$psel$enet3)

colors <- bpy.colors(6)[-c(1, 6)]
png(paste(path.graph, "gren_metabol_alzheimer_res1_performance1.png", sep=""),
    units="in", width=14, height=6, res=120)
par(mfrow=c(1, 2), mar=c(5.1, 5.1, 4.1, 2.1))

# auc
ylim <- range(auc3.grridge, auc3.ridge, auc3.gren1, auc3.enet1, auc3.gren2, 
              auc3.enet2, auc3.gren3, auc3.enet3)
xlim <- range(psel3.grridge, psel3.gren1, psel3.enet1, psel3.gren2, psel3.enet2, 
              psel3.gren2, psel3.enet2)
plot(psel3.gren1, auc3.gren1, ylim=ylim, xlim=xlim, 
     main="a)", xlab="Number of selected features", ylab="AUC", 
     cex.axis=1.5, cex.lab=2, cex.main=2, lwd=1.5, col=colors[1], type="l")
lines(psel3.enet1, auc3.enet1, lwd=1.5, col=colors[1], lty=2)

lines(psel3.gren2, auc3.gren2, lwd=1.5, col=colors[2], lty=1)
lines(psel3.enet2, auc3.enet2, lwd=1.5, col=colors[2], lty=2)

lines(psel3.gren3, auc3.gren3, lwd=1.5, col=colors[3], lty=1)
lines(psel3.enet3, auc3.enet3, lwd=1.5, col=colors[3], lty=2)

lines(psel3.grridge, auc3.grridge, lwd=1.5, col=colors[4], lty=1)
abline(h=auc3.ridge, lwd=1.5, col=colors[4], lty=2)

# Briers
ylim <- range(briers3.grridge, briers3.ridge, briers3.gren1, briers3.enet1, 
              briers3.gren2, briers3.enet2, briers3.gren3, briers3.enet3)
xlim <- range(psel3.grridge, psel3.gren1, psel3.enet1, psel3.gren2, psel3.enet2, 
              psel3.gren2, psel3.enet2)
plot(psel3.gren1, briers3.gren1, ylim=ylim, xlim=xlim, 
     main="b)", xlab="Number of selected features", ylab="Brier skill score", 
     cex.axis=1.5, cex.lab=2, cex.main=2, lwd=1.5, col=colors[1], type="l")
lines(psel3.enet1, briers3.enet1, lwd=1.5, col=colors[1], lty=2)

lines(psel3.gren2, briers3.gren2, lwd=1.5, col=colors[2], lty=1)
lines(psel3.enet2, briers3.enet2, lwd=1.5, col=colors[2], lty=2)

lines(psel3.gren3, briers3.gren3, lwd=1.5, col=colors[3], lty=1)
lines(psel3.enet3, briers3.enet3, lwd=1.5, col=colors[3], lty=2)

lines(psel3.grridge, briers3.grridge, lwd=1.5, col=colors[4], lty=1)
abline(h=briers3.ridge, lwd=1.5, col=colors[4], lty=2)

# legend
leglabels <- c(expression(paste("enet, ", alpha==0.05)),
               expression(paste("enet, ", alpha==0.5)),
               expression(paste("enet, ", alpha==0.95)), "ridge",
               "group-regularized", "not group-regularized")
legend("bottomleft", legend=leglabels, fill=c(colors, 0, 0),
       lty=c(rep(NA, 4), 1, 2), lwd=c(rep(NA, 4), 1.5, 1.5),
       border=c(rep(1, 4), 0, 0), merge=TRUE, seg.len=1, cex=1.3)
dev.off()


# str(results3$pred$sglasso1)
# 
# 
# unlist(lapply(results3$pred$sglasso1, ncol))
# lapply(results3$pred$sglasso1)
# 
# method <- "sglasso1"
# measure <- "auc"
# data <- results3
# 
# data[["psel"]][[method]]
# pred.mean <- function(method, measure, data) {
#   ind <- unlist(data[["psel"]][[method]])
#   dep <- unlist(data[[measure]][[method]])
#   test <- sapply(sort(unique(ind)), function(s) {mean(dep[ind==s])})
#   
#   out <- tryCatch(list(x=sort(ind), y=predict(loess(dep[order(
#     ind)] ~ sort(ind)))), warning=function(w) {
#       list(x=unique(sort(ind)), y=sapply(unique(sort(ind)), function(psel) {
#         mean(dep[ind==psel], na.rm=TRUE)}))})
#   return(list(x=unique(out$x), y=sapply(unique(out$x), function(psel) {
#     mean(out$y[out$x==psel])})))
# }
# 
# 
# froot <- function(lambda, psel, x, y, part, alpha, type) {
#   if(type=="sglasso") {
#     fit <- SGL(list(x=x, y=y), part, type="logit", alpha=alpha, standardize=FALSE, 
#                lambdas=lambda, nlam=1, min.frac=0.001)
#     out <- psel - sum(fit$beta!=0)
#   } else if(type=="cmcp") {
#     fit <- grpreg(x, y, part, "cMCP", "binomial", lambda=lambda, alpha=alpha)
#     out <- psel - sum(fit$beta[-1, ]!=0)
#   } else if(type=="gelasso") {
#     fit <- grpreg(x, y, part, "gel", "binomial", lambda=lambda, alpha=alpha)
#     out <- psel - sum(fit$beta[-1, ]!=0)
#   }
#   return(out)
# }
# 
# 
# psel <- c(seq(1, 5, 1), seq(7, 15, 2), seq(20, 40, 5), seq(50, 90, 10), seq(110, 150, 20))
# 
# 
# test.init <- grpreg(x, y, part3, "gel", "binomial", nlambda=2, alpha=0.05)
# test.gelasso1 <- sapply(psel, function(csel) {
#   test.root <- uniroot(froot, range(test.init$lambda), psel=csel, x, y, part3, 
#                        alpha=0.05, type="gelasso", maxiter=100);
#   test.fit <- grpreg(x, y, part3, "gel", "binomial", lambda=test.root$root, 
#                      alpha=0.05);
#   return(sum(test.fit$beta[-1, ]!=0))})
# 
# test.init <- grpreg(x, y, part3, "gel", "binomial", nlambda=2, alpha=0.5)
# test.gelasso2 <- sapply(psel, function(csel) {
#   test.root <- uniroot(froot, range(test.init$lambda), psel=csel, x, y, part3, 
#                        alpha=0.5, type="gelasso", maxiter=100);
#   test.fit <- grpreg(x, y, part3, "gel", "binomial", lambda=test.root$root, 
#                      alpha=0.5);
#   return(sum(test.fit$beta[-1, ]!=0))})
# 
# test.init <- grpreg(x, y, part3, "gel", "binomial", nlambda=2, alpha=0.95)
# test.gelasso3 <- sapply(psel, function(csel) {
#   test.root <- uniroot(froot, range(test.init$lambda), psel=csel, x, y, part3, 
#                        alpha=0.95, type="gelasso", maxiter=100);
#   test.fit <- grpreg(x, y, part3, "gel", "binomial", lambda=test.root$root, 
#                      alpha=0.95);
#   return(sum(test.fit$beta[-1, ]!=0))})
# 
# 
# 
# 
# # SGL
# test.init <- SGL(list(x=x, y=y), part3, type="logit", alpha=0.05, 
#                  standardize=FALSE, nlam=2)
# test.sglasso1 <- sapply(psel, function(csel) {
#   test.root <- uniroot(froot, range(test.init$lambda), psel=csel, x, y, part3, 
#                        alpha=0.05, type="sglasso", maxiter=100);
#   test.fit <- SGL(list(x=x, y=y), part3, type="logit", alpha=0.05, 
#                   standardize=FALSE, lambdas=test.root$root, nlam=1);
#   return(sum(test.fit$beta!=0))})
# 
# test.init <- SGL(list(x=x, y=y), part3, type="logit", alpha=0.5, 
#                  standardize=FALSE, nlam=2)
# test.sglasso2 <- sapply(psel, function(csel) {
#   test.root <- uniroot(froot, range(test.init$lambda), psel=csel, x, y, part3, 
#                        alpha=0.5, type="sglasso", maxiter=100);
#   test.fit <- SGL(list(x=x, y=y), part3, type="logit", alpha=0.5, 
#                   standardize=FALSE, lambdas=test.root$root, nlam=1);
#   return(sum(test.fit$beta!=0))})
# 
# test.init <- SGL(list(x=x, y=y), part3, type="logit", alpha=0.95, 
#                  standardize=FALSE, nlam=2, min.frac=0.001)
# test.sglasso3 <- sapply(psel, function(csel) {
#   test.root <- uniroot(froot, range(test.init$lambda), psel=csel, x, y, part3, 
#                        alpha=0.95, type="sglasso", maxiter=100);
#   test.fit <- SGL(list(x=x, y=y), part3, type="logit", alpha=0.95, 
#                   standardize=FALSE, lambdas=test.root$root, nlam=1);
#   return(sum(test.fit$beta!=0))})
# 
# cbind(psel, test.gelasso1, test.gelasso2, test.gelasso3)
# cbind(psel, test.sglasso1, test.sglasso2, test.sglasso3)
# 
# 



# cv3.gren1 <- cv.gren(x, y, partitions=list(RSDqc=part3), alpha=0.05)
# cv3.gren2 <- cv.gren(x, y, partitions=list(RSDqc=part3), alpha=0.5)
# cv3.gren3 <- cv.gren(x, y, partitions=list(RSDqc=part3), alpha=0.95)
# 
# pred3.gren1 <- predict(fit3.gren1, xtest, type="groupreg",
#                        s=fit3.gren1$freq.model$groupreg$lambda)
# pred3.gren2 <- predict(fit3.gren2, xtest, type="groupreg",
#                        s=fit3.gren2$freq.model$groupreg$lambda)
# pred3.gren3 <- predict(fit3.gren3, xtest, type="groupreg",
#                        s=fit3.gren3$freq.model$groupreg$lambda)
# 
# pred3.enet1 <- predict(fit3.gren1, xtest, type="regular",
#                        s=fit3.gren1$freq.model$regular$lambda)
# pred3.enet2 <- predict(fit3.gren2, xtest, type="regular",
#                        s=fit3.gren2$freq.model$regular$lambda)
# pred3.enet3 <- predict(fit3.gren3, xtest, type="regular",
#                        s=fit3.gren3$freq.model$regular$lambda)
# 
# psel3.gren1 <- fit3.gren1$freq.model$groupreg$df
# psel3.gren2 <- fit3.gren2$freq.model$groupreg$df
# psel3.gren3 <- fit3.gren3$freq.model$groupreg$df
# 
# psel3.enet1 <- fit3.gren1$freq.model$regular$df
# psel3.enet2 <- fit3.gren2$freq.model$regular$df
# psel3.enet3 <- fit3.gren3$freq.model$regular$df
# 
# auc3.gren1 <- apply(pred3.gren1, 2, function(pred) {pROC::roc(ytest, pred)$auc})
# auc3.gren2 <- apply(pred3.gren2, 2, function(pred) {pROC::roc(ytest, pred)$auc})
# auc3.gren3 <- apply(pred3.gren3, 2, function(pred) {pROC::roc(ytest, pred)$auc})
# 
# auc3.enet1 <- apply(pred3.enet1, 2, function(pred) {pROC::roc(ytest, pred)$auc})
# auc3.enet2 <- apply(pred3.enet2, 2, function(pred) {pROC::roc(ytest, pred)$auc})
# auc3.enet3 <- apply(pred3.enet3, 2, function(pred) {pROC::roc(ytest, pred)$auc})
# 
# results3 <- list(idtrain=idtrain,
#                  pred=list(enet1=pred3.enet1, enet2=pred3.enet2, 
#                            enet3=pred3.enet3, gren1=pred3.gren1, 
#                            gren2=pred3.gren2, gren3=pred3.gren3),
#                  psel=list(enet1=psel3.enet1, enet2=psel3.enet2, 
#                            enet3=psel3.enet3, gren1=psel3.gren1, 
#                            gren2=psel3.gren2, gren3=psel3.gren3),
#                  auc=list(enet1=auc3.enet1, enet2=auc3.enet2, enet3=auc3.enet3, 
#                           gren1=auc3.gren1, gren2=auc3.gren2, gren3=auc3.gren3))
# 
# plot(sort(results3$psel$enet1), results3$auc$enet1[order(results3$psel$enet1)], 
#      type="l", xlim=range(results3$psel), ylim=range(results3$auc), lty=2,
#      col=colors[1])
# lines(sort(results3$psel$enet2), results3$auc$enet2[order(results3$psel$enet2)],
#       lty=2, col=colors[2])
# lines(sort(results3$psel$enet3), results3$auc$enet3[order(results3$psel$enet2)],
#       lty=2, col=colors[3])
# 
# lines(sort(results3$psel$gren1), results3$auc$gren1[order(results3$psel$gren1)],
#       lty=1, col=colors[1])
# lines(sort(results3$psel$gren2), results3$auc$gren2[order(results3$psel$gren2)],
#       lty=1, col=colors[2])
# lines(sort(results3$psel$gren3), results3$auc$gren3[order(results3$psel$gren3)],
#       lty=1, col=colors[3])
# 
# leglabels <- c(expression(alpha==0.05), expression(alpha==0.5),
#                expression(alpha==0.95), "gren", "enet")
# legend("bottomright", merge=TRUE, seg.len=1, cex=1.3, fill=c(colors, 0, 0), 
#        lty=c(rep(NA, length(colors)), 1, 2), 
#        lwd=c(rep(NA, length(colors)), 1.5, 1.5),
#        border=c(rep(1, length(colors)), 0, 0), legend=leglabels)
# 
# 
# 
# 
# 
# ####
# fit4.gren1 <- gren(xtrain, ytrain, partitions=list(RSDqc=part4), alpha=0.05)
# fit4.gren2 <- gren(xtrain, ytrain, partitions=list(RSDqc=part4), alpha=0.5)
# fit4.gren3 <- gren(xtrain, ytrain, partitions=list(RSDqc=part4), alpha=0.95)
# 
# pred4.gren1 <- predict(fit4.gren1, xtest, type="groupreg",
#                        s=fit4.gren1$freq.model$groupreg$lambda)
# pred4.gren2 <- predict(fit4.gren2, xtest, type="groupreg",
#                        s=fit4.gren2$freq.model$groupreg$lambda)
# pred4.gren3 <- predict(fit4.gren3, xtest, type="groupreg",
#                        s=fit4.gren3$freq.model$groupreg$lambda)
# 
# pred4.enet1 <- predict(fit4.gren1, xtest, type="regular",
#                        s=fit4.gren1$freq.model$regular$lambda)
# pred4.enet2 <- predict(fit4.gren2, xtest, type="regular",
#                        s=fit4.gren2$freq.model$regular$lambda)
# pred4.enet3 <- predict(fit4.gren3, xtest, type="regular",
#                        s=fit4.gren3$freq.model$regular$lambda)
# 
# psel4.gren1 <- fit4.gren1$freq.model$groupreg$df
# psel4.gren2 <- fit4.gren2$freq.model$groupreg$df
# psel4.gren3 <- fit4.gren3$freq.model$groupreg$df
# 
# psel4.enet1 <- fit4.gren1$freq.model$regular$df
# psel4.enet2 <- fit4.gren2$freq.model$regular$df
# psel4.enet3 <- fit4.gren3$freq.model$regular$df
# 
# auc4.gren1 <- apply(pred4.gren1, 2, function(pred) {pROC::roc(ytest, pred)$auc})
# auc4.gren2 <- apply(pred4.gren2, 2, function(pred) {pROC::roc(ytest, pred)$auc})
# auc4.gren3 <- apply(pred4.gren3, 2, function(pred) {pROC::roc(ytest, pred)$auc})
# 
# auc4.enet1 <- apply(pred4.enet1, 2, function(pred) {pROC::roc(ytest, pred)$auc})
# auc4.enet2 <- apply(pred4.enet2, 2, function(pred) {pROC::roc(ytest, pred)$auc})
# auc4.enet3 <- apply(pred4.enet3, 2, function(pred) {pROC::roc(ytest, pred)$auc})
# 
# results4 <- list(idtrain=idtrain,
#                  pred=list(enet1=pred4.enet1, enet2=pred4.enet2, 
#                            enet3=pred4.enet3, gren1=pred4.gren1, 
#                            gren2=pred4.gren2, gren3=pred4.gren3),
#                  psel=list(enet1=psel4.enet1, enet2=psel4.enet2, 
#                            enet3=psel4.enet3, gren1=psel4.gren1, 
#                            gren2=psel4.gren2, gren3=psel4.gren3),
#                  auc=list(enet1=auc4.enet1, enet2=auc4.enet2, enet3=auc4.enet3, 
#                           gren1=auc4.gren1, gren2=auc4.gren2, gren3=auc4.gren3))
# 
# plot(sort(results4$psel$enet1), results4$auc$enet1[order(results4$psel$enet1)], 
#      type="l", xlim=range(results4$psel), ylim=range(results4$auc), lty=2,
#      col=colors[1])
# lines(sort(results4$psel$enet2), results4$auc$enet2[order(results4$psel$enet2)],
#       lty=2, col=colors[2])
# lines(sort(results4$psel$enet3), results4$auc$enet3[order(results4$psel$enet2)],
#       lty=2, col=colors[3])
# 
# lines(sort(results4$psel$gren1), results4$auc$gren1[order(results4$psel$gren1)],
#       lty=1, col=colors[1])
# lines(sort(results4$psel$gren2), results4$auc$gren2[order(results4$psel$gren2)],
#       lty=1, col=colors[2])
# lines(sort(results4$psel$gren3), results4$auc$gren3[order(results4$psel$gren3)],
#       lty=1, col=colors[3])