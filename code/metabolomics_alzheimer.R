#!/usr/bin/env Rscript

### installation of packages
# if(!("devtools" %in% installed.packages())) {
#   install.packages("devtools")
# }
# library(devtools)
# install_github("magnusmunch/gren/rpackage", local=FALSE,
#                auth_token=Sys.getenv("GITHUB_PAT"))

### libraries
library(gren)
library(GRridge)
library(grpreg)
library(SGL)
library(Biobase)

### load data
load("data/ESetMbolCSFPR2.Rdata")

### select the patients with APOE E4 allel and controls without APOE E4
pheno <- pData(ESetMbolCSFPR2)
pheno.apoe <- pheno[(pheno$D_diag_name=="Probable AD" & pheno$APOE=="E4YES") |
                      (pheno$D_diag_name=="Subjectieve klachten" & 
                         pheno$APOE=="E4NO"), ]


metabol <- t(exprs(ESetMbolCSFPR2))
metabol.apoe <- metabol[(pheno$D_diag_name=="Probable AD" & 
                           pheno$APOE=="E4YES") |
                          (pheno$D_diag_name=="Subjectieve klachten" & 
                             pheno$APOE=="E4NO"), ]

### create co-data
feat <- fData(ESetMbolCSFPR2)

# platformcode gives the type of metabolites
platformcode <- feat$PlatformCode2

# quality gives groups based on a quality measure of the assayd features
part.RSDqc <- CreatePartition(feat$RSDqc, ngroup=3)
quality <- rep(c(1:3), sapply(part.RSDqc, length))[unlist(part.RSDqc)]

# sds gives groups based on the standard deviations of the features
part.sds.apoe <- CreatePartition(apply(metabol.apoe, 1, sd), ngroup=3)
sds.apoe <- rep(c(1:3), sapply(part.sds.apoe, length))[unlist(part.sds.apoe)]

### transform data
alzheim.apoe <- as.numeric(pheno.apoe$D_diag_name) - 1
metabol.apoe.scaled <- scale(metabol.apoe)

### fitting the models
y <- alzheim.apoe
x <- metabol.apoe.scaled
part <- platformcode
n <- nrow(x)
p <- ncol(x)
set.seed(2018)

fit1.gren1 <- gren(x, y, partitions=list(part=part), alpha=0.05)
fit1.gren2 <- gren(x, y, partitions=list(part=part), alpha=0.5)
fit1.gren3 <- gren(x, y, partitions=list(part=part), alpha=0.95)

fit1.grridge <- grridge(t(x), y, list(part=split(1:p, part)))
fit1.sgl1 <- cvSGL(list(x=x, y=y), part, type="logit", alpha=0.05)
fit1.sgl2 <- cvSGL(list(x=x, y=y), part, type="logit", alpha=0.5)
fit1.sgl3 <- cvSGL(list(x=x, y=y), part, type="logit", alpha=0.95)

fit1.cmcp1 <- cv.grpreg(x, y, part, penalty="cMCP", family="binomial", 
                        alpha=0.05)
fit1.cmcp2 <- cv.grpreg(x, y, part, penalty="cMCP", family="binomial", 
                        alpha=0.5)
fit1.cmcp3 <- cv.grpreg(x, y, part, penalty="cMCP", family="binomial", 
                        alpha=0.95)

fit1.gel1 <- cv.grpreg(x, y, part, penalty="gel", family="binomial", alpha=0.05)
fit1.gel2 <- cv.grpreg(x, y, part, penalty="gel", family="binomial", alpha=0.5)
fit1.gel3 <- cv.grpreg(x, y, part, penalty="gel", family="binomial", alpha=0.95)

save(fit1.gren1, file="results/metabolomics_alzheimer_fit1.Rdata")

### cross-validation of performance
nfolds <- 10
rest <- n %% nfolds
foldid <- sample(rep(1:nfolds, times=round(c(rep(
  n %/% nfolds + as.numeric(rest!=0), times=rest),
  rep(n %/% nfolds, times=nfolds - rest)))))
pred1 <- data.frame(ridge=rep(NA, n), grridge=rep(NA, n),
                    gren1=rep(NA, n), gren2=rep(NA, n), gren3=rep(NA, n),
                    enet1=rep(NA, n), enet2=rep(NA, n), enet3=rep(NA, n),
                    sgl1=rep(NA, n), sgl2=rep(NA, n), sgl3=rep(NA, n),
                    cmcp1=rep(NA, n), cmcp2=rep(NA, n), cmcp3=rep(NA, n),
                    gel1=rep(NA, n), gel2=rep(NA, n), gel3=rep(NA, n))
psel1 <- data.frame(ridge=rep(p, nfolds), grridge=rep(p, nfolds),
                    gren1=rep(NA, nfolds), gren2=rep(NA, nfolds), 
                    gren3=rep(NA, nfolds), enet1=rep(NA, nfolds), 
                    enet2=rep(NA, nfolds), enet3=rep(NA, nfolds),
                    sgl1=rep(NA, nfolds), sgl2=rep(NA, nfolds), 
                    sgl3=rep(NA, nfolds), cmcp1=rep(NA, nfolds), 
                    cmcp2=rep(NA, nfolds), cmcp3=rep(NA, nfolds),
                    gel1=rep(NA, nfolds), gel2=rep(NA, nfolds), 
                    gel3=rep(NA, nfolds))
for(k in 1:nfolds) {
  xtrain <- matrix(x[foldid!=k, ], ncol=p)
  xtest <- matrix(x[foldid==k, ], ncol=p)
  ytrain <- y[foldid!=k]
  ytest <- y[foldid==k]
  
  ### fit methods
  cv1.gren1 <- gren(xtrain, ytrain, partitions=list(part=part), alpha=0.05)
  cv1.gren2 <- gren(xtrain, ytrain, partitions=list(part=part), alpha=0.5)
  cv1.gren3 <- gren(xtrain, ytrain, partitions=list(part=part), alpha=0.95)
  
  cv1.grridge <- grridge(t(xtrain), ytrain, list(part=split(1:p, part)))
  cv1.sgl1 <- cvSGL(list(x=xtrain, y=ytrain), part, type="logit", alpha=0.05)
  cv1.sgl2 <- cvSGL(list(x=xtrain, y=ytrain), part, type="logit", alpha=0.5)
  cv1.sgl3 <- cvSGL(list(x=xtrain, y=ytrain), part, type="logit", alpha=0.95)
  
  cv1.cmcp1 <- cv.grpreg(xtrain, ytrain, part, penalty="cMCP", 
                         family="binomial", alpha=0.05)
  cv1.cmcp2 <- cv.grpreg(xtrain, ytrain, part, penalty="cMCP", 
                         family="binomial", alpha=0.5)
  cv1.cmcp3 <- cv.grpreg(xtrain, ytrain, part, penalty="cMCP", 
                         family="binomial", alpha=0.95)
  
  cv1.gel1 <- cv.grpreg(xtrain, ytrain, part, penalty="gel", family="binomial",
                        alpha=0.05)
  cv1.gel2 <- cv.grpreg(xtrain, ytrain, part, penalty="gel", family="binomial",
                        alpha=0.5)
  cv1.gel3 <- cv.grpreg(xtrain, ytrain, part, penalty="gel", family="binomial",
                        alpha=0.95)
  
  pred1$ridge[foldid==k] <- predict.grridge(cv1.grridge, t(xtest))[, 1]
  pred1$grridge[foldid==k] <- predict.grridge(cv1.grridge, t(xtest))[, 2]
  
  pred1$enet1[foldid==k] <- predict(cv1.gren1, xtest, type="regular")
  pred1$enet2[foldid==k] <- predict(cv1.gren2, xtest, type="regular")
  pred1$enet3[foldid==k] <- predict(cv1.gren3, xtest, type="regular")
  
  pred1$gren1[foldid==k] <- predict(cv1.gren1, xtest, type="groupreg")
  pred1$gren2[foldid==k] <- predict(cv1.gren2, xtest, type="groupreg")
  pred1$gren3[foldid==k] <- predict(cv1.gren3, xtest, type="groupreg")
  
  pred1$sgl1[foldid==k] <- 1/(1 + exp(-xtest %*% cv1.sgl1$fit$beta[, which.min(
    cv1.sgl1$lldiff)]))
  pred1$sgl2[foldid==k] <- 1/(1 + exp(-xtest %*% cv1.sgl2$fit$beta[, which.min(
    cv1.sgl2$lldiff)]))
  pred1$sgl3[foldid==k] <- 1/(1 + exp(-xtest %*% cv1.sgl3$fit$beta[, which.min(
    cv1.sgl3$lldiff)]))
  
  pred1$cmcp1[foldid==k] <- predict(cv1.cmcp1, xtest, type="response")
  pred1$cmcp2[foldid==k] <- predict(cv1.cmcp2, xtest, type="response")
  pred1$cmcp3[foldid==k] <- predict(cv1.cmcp3, xtest, type="response")
  
  pred1$gel1[foldid==k] <- predict(cv1.gel1, xtest, type="response")
  pred1$gel2[foldid==k] <- predict(cv1.gel2, xtest, type="response")
  pred1$gel3[foldid==k] <- predict(cv1.gel3, xtest, type="response")
  
  psel1$enet1[k] <- sum(coef(cv1.gren1, type="regular")!=0) - 1
  psel1$enet2[k] <- sum(coef(cv1.gren2, type="regular")!=0) - 1
  psel1$enet3[k] <- sum(coef(cv1.gren3, type="regular")!=0) - 1
  
  psel1$gren1[k] <- sum(coef(cv1.gren1, type="groupreg")!=0) - 1
  psel1$gren2[k] <- sum(coef(cv1.gren2, type="groupreg")!=0) - 1
  psel1$gren3[k] <- sum(coef(cv1.gren3, type="groupreg")!=0) - 1
  
  psel1$sgl1[k] <- sum(cv1.sgl1$fit$beta[, which.min(cv1.sgl1$lldiff)]!=0)
  psel1$sgl2[k] <- sum(cv1.sgl2$fit$beta[, which.min(cv1.sgl2$lldiff)]!=0)
  psel1$sgl3[k] <- sum(cv1.sgl3$fit$beta[, which.min(cv1.sgl3$lldiff)]!=0)
  
  psel1$cmcp1[k] <- sum(cv1.cmcp1$fit$beta[, cv1.cmcp1$min]!=0) - 1
  psel1$cmcp2[k] <- sum(cv1.cmcp2$fit$beta[, cv1.cmcp2$min]!=0) - 1
  psel1$cmcp3[k] <- sum(cv1.cmcp3$fit$beta[, cv1.cmcp3$min]!=0) - 1
  
  psel1$gel1[k] <- sum(cv1.gel1$fit$beta[, cv1.gel1$min]!=0) - 1
  psel1$gel2[k] <- sum(cv1.gel2$fit$beta[, cv1.gel2$min]!=0) - 1
  psel1$gel3[k] <- sum(cv1.gel3$fit$beta[, cv1.gel3$min]!=0) - 1
  
  res1 <- rbind(pred1, psel1)
  rownames(res1) <- c(paste0("pred", c(1:n)), paste0("psel", 1:nfolds))
  write.table(res1, file="results/metabolomics_alzheimer_res1.csv")
  
}

# 
# 
# 
# path.res <- ifelse(as.character(Sys.info()[1])!="Darwin", "~/EBEN/results/",
#                    "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/results/")
# path.graph <- "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/graphs/"
# path.data <- as.character(ifelse(Sys.info()[1]=="Darwin", "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/data/" , "~/EBEN/data/"))
# 
# library(sp)
# 
# load(paste(path.res, "gren_metabol_alzheimer_res1.Rdata", sep=""))
# load(paste(path.data, "ESetMbolCSFPR2.Rdata", sep=""))
# y <- as.numeric(pData(ESetMbolCSFPR2)$D_diag_name) - 1
# 
# auc3.ridge <- pROC::auc(y, results3$pred$ridge)
# auc3.grridge <- apply(results3$pred$grridge, 2, function(s) {pROC::auc(y, s)})
# auc3.gren1 <- apply(results3$pred$gren1, 2, function(s) {pROC::auc(y, s)})
# auc3.gren2 <- apply(results3$pred$gren2, 2, function(s) {pROC::auc(y, s)})
# auc3.gren3 <- apply(results3$pred$gren3, 2, function(s) {pROC::auc(y, s)})
# auc3.enet1 <- apply(results3$pred$enet1, 2, function(s) {pROC::auc(y, s)})
# auc3.enet2 <- apply(results3$pred$enet2, 2, function(s) {pROC::auc(y, s)})
# auc3.enet3 <- apply(results3$pred$enet3, 2, function(s) {pROC::auc(y, s)})
# 
# const <- sum((y - mean(y))^2)
# briers3.ridge <- 1 - sum((y - results3$pred$ridge)^2)/const
# briers3.grridge <- apply(results3$pred$grridge, 2, function(s) {
#   1 - sum((y - s)^2)/const})
# briers3.gren1 <- apply(results3$pred$gren1, 2, function(s) {
#   1 - sum((y - s)^2)/const})
# briers3.gren2 <- apply(results3$pred$gren2, 2, function(s) {
#   1 - sum((y - s)^2)/const})
# briers3.gren3 <- apply(results3$pred$gren3, 2, function(s) {
#   1 - sum((y - s)^2)/const})
# briers3.enet1 <- apply(results3$pred$enet1, 2, function(s) {
#   1 - sum((y - s)^2)/const})
# briers3.enet2 <- apply(results3$pred$enet2, 2, function(s) {
#   1 - sum((y - s)^2)/const})
# briers3.enet3 <- apply(results3$pred$enet3, 2, function(s) {
#   1 - sum((y - s)^2)/const})
# 
# psel3.grridge <- colMeans(results3$psel$grridge)
# psel3.gren1 <- colMeans(results3$psel$gren1)
# psel3.gren2 <- colMeans(results3$psel$gren2)
# psel3.gren3 <- colMeans(results3$psel$gren3)
# psel3.enet1 <- colMeans(results3$psel$enet1)
# psel3.enet2 <- colMeans(results3$psel$enet2)
# psel3.enet3 <- colMeans(results3$psel$enet3)
# 
# colors <- bpy.colors(6)[-c(1, 6)]
# png(paste(path.graph, "gren_metabol_alzheimer_res1_performance1.png", sep=""),
#     units="in", width=14, height=6, res=120)
# par(mfrow=c(1, 2), mar=c(5.1, 5.1, 4.1, 2.1))
# 
# # auc
# ylim <- range(auc3.grridge, auc3.ridge, auc3.gren1, auc3.enet1, auc3.gren2, 
#               auc3.enet2, auc3.gren3, auc3.enet3)
# xlim <- range(psel3.grridge, psel3.gren1, psel3.enet1, psel3.gren2, psel3.enet2, 
#               psel3.gren2, psel3.enet2)
# plot(psel3.gren1, auc3.gren1, ylim=ylim, xlim=xlim, 
#      main="a)", xlab="Number of selected features", ylab="AUC", 
#      cex.axis=1.5, cex.lab=2, cex.main=2, lwd=1.5, col=colors[1], type="l")
# lines(psel3.enet1, auc3.enet1, lwd=1.5, col=colors[1], lty=2)
# 
# lines(psel3.gren2, auc3.gren2, lwd=1.5, col=colors[2], lty=1)
# lines(psel3.enet2, auc3.enet2, lwd=1.5, col=colors[2], lty=2)
# 
# lines(psel3.gren3, auc3.gren3, lwd=1.5, col=colors[3], lty=1)
# lines(psel3.enet3, auc3.enet3, lwd=1.5, col=colors[3], lty=2)
# 
# lines(psel3.grridge, auc3.grridge, lwd=1.5, col=colors[4], lty=1)
# abline(h=auc3.ridge, lwd=1.5, col=colors[4], lty=2)
# 
# # Briers
# ylim <- range(briers3.grridge, briers3.ridge, briers3.gren1, briers3.enet1, 
#               briers3.gren2, briers3.enet2, briers3.gren3, briers3.enet3)
# xlim <- range(psel3.grridge, psel3.gren1, psel3.enet1, psel3.gren2, psel3.enet2, 
#               psel3.gren2, psel3.enet2)
# plot(psel3.gren1, briers3.gren1, ylim=ylim, xlim=xlim, 
#      main="b)", xlab="Number of selected features", ylab="Brier skill score", 
#      cex.axis=1.5, cex.lab=2, cex.main=2, lwd=1.5, col=colors[1], type="l")
# lines(psel3.enet1, briers3.enet1, lwd=1.5, col=colors[1], lty=2)
# 
# lines(psel3.gren2, briers3.gren2, lwd=1.5, col=colors[2], lty=1)
# lines(psel3.enet2, briers3.enet2, lwd=1.5, col=colors[2], lty=2)
# 
# lines(psel3.gren3, briers3.gren3, lwd=1.5, col=colors[3], lty=1)
# lines(psel3.enet3, briers3.enet3, lwd=1.5, col=colors[3], lty=2)
# 
# lines(psel3.grridge, briers3.grridge, lwd=1.5, col=colors[4], lty=1)
# abline(h=briers3.ridge, lwd=1.5, col=colors[4], lty=2)
# 
# # legend
# leglabels <- c(expression(paste("enet, ", alpha==0.05)),
#                expression(paste("enet, ", alpha==0.5)),
#                expression(paste("enet, ", alpha==0.95)), "ridge",
#                "group-regularized", "not group-regularized")
# legend("bottomleft", legend=leglabels, fill=c(colors, 0, 0),
#        lty=c(rep(NA, 4), 1, 2), lwd=c(rep(NA, 4), 1.5, 1.5),
#        border=c(rep(1, 4), 0, 0), merge=TRUE, seg.len=1, cex=1.3)
# dev.off()


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