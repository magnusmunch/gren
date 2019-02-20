################################## manuscript ##################################
# ---- figures ----
# ---- lines_metabolomics_alzheimer_res1_auc ----
library(pROC)
library(grpreg)
load("results/metabolomics_alzheimer_fit1.Rdata")
res1 <- read.table("results/metabolomics_alzheimer_res1.csv", 
                   stringsAsFactors=FALSE)

methods2 <- c("ridge", "grridge", "gren", "enet", "sgl", "cmcp", "gel")
pred1 <- as.matrix(res1[substr(rownames(res1), 1, 4)=="pred", ])
psel1 <- as.matrix(res1[substr(rownames(res1), 1, 4)=="psel", ])
auc1 <- as.matrix(res1[substr(rownames(res1), 1, 3)=="auc", ])
plot.data1 <- lapply(c(t(outer(methods2[-c(1, 2)], 1:3, paste, sep=""))), 
                     function(m) {aggregate(auc1[, grepl(m, colnames(auc1))], 
                                            list(psel=psel1[, grepl(m, colnames(
                                              psel1))]), mean)})
cv.psel1 <- c(sum(coef(fit1.gren1$freq.model$groupreg, 
                       s=fit1.gren1$lambda)[-1, ]!=0),
              sum(coef(fit1.gren2$freq.model$groupreg, 
                       s=fit1.gren2$lambda)[-1, ]!=0),
              sum(coef(fit1.gren3$freq.model$groupreg, 
                       s=fit1.gren3$lambda)[-1, ]!=0),
              sum(coef(fit1.gren1$freq.model$regular, 
                       s=fit1.gren1$lambda)[-1, ]!=0),
              sum(coef(fit1.gren2$freq.model$regular, 
                       s=fit1.gren2$lambda)[-1, ]!=0),
              sum(coef(fit1.gren3$freq.model$regular, 
                       s=fit1.gren3$lambda)[-1, ]!=0),
              psel1[, grepl("sgl1", colnames(psel1))][which.min(
                fit1.sgl1$lldiff)],
              psel1[, grepl("sgl2", colnames(psel1))][which.min(
                fit1.sgl2$lldiff)],
              psel1[, grepl("sgl3", colnames(psel1))][which.min(
                fit1.sgl3$lldiff)],
              psel1[, grepl("cmcp1", colnames(psel1))][fit1.cmcp1$min],
              psel1[, grepl("cmcp2", colnames(psel1))][fit1.cmcp2$min],
              psel1[, grepl("cmcp3", colnames(psel1))][fit1.cmcp3$min],
              psel1[, grepl("gel1", colnames(psel1))][fit1.gel1$min],
              psel1[, grepl("gel2", colnames(psel1))][fit1.gel2$min],
              psel1[, grepl("gel3", colnames(psel1))][fit1.gel3$min])
col1 <- c(1:length(methods2))
lty1 <- c(1:4)

par(mfrow=c(1, 2))
plot(plot.data1[[1]], type="l", xlab="Number of selected features", ylab="AUC", 
     main="a)", ylim=range(auc1), xlim=range(psel1), col=col1[3], lty=lty1[1])
lines(plot.data1[[2]], col=col1[3], lty=lty1[2])
lines(plot.data1[[3]], col=col1[3], lty=lty1[3])
lines(plot.data1[[4]], col=col1[4], lty=lty1[1])
lines(plot.data1[[5]], col=col1[4], lty=lty1[2])
lines(plot.data1[[6]], col=col1[4], lty=lty1[3])
lines(plot.data1[[7]], col=col1[5], lty=lty1[1])
lines(plot.data1[[8]], col=col1[5], lty=lty1[2])
lines(plot.data1[[9]], col=col1[5], lty=lty1[3])
lines(plot.data1[[10]], col=col1[6], lty=lty1[1])
lines(plot.data1[[11]], col=col1[6], lty=lty1[2])
lines(plot.data1[[12]], col=col1[6], lty=lty1[3])
lines(plot.data1[[13]], col=col1[7], lty=lty1[1])
lines(plot.data1[[14]], col=col1[8], lty=lty1[2])
lines(plot.data1[[15]], col=col1[9], lty=lty1[3])
abline(h=auc1[colnames(auc1)=="ridge"], col=col1[1], lty=lty1[4])
abline(h=auc1[colnames(auc1)=="grridge"], col=col1[2], lty=lty1[4])
points(t(sapply(1:length(plot.data1), function(m) {
  subset(plot.data1[[m]], psel==cv.psel1[m])})), 
  col=rep(col1[-c(1, 2)], each=3), pch=1)
legend("bottomright", 
       legend=c(methods2, expression(alpha==0.05, alpha==0.5, alpha==0.95),
                "CV model"), 
       fill=c(col1, rep(NA, 4)),
       border=c(rep(1, length(methods2)), rep(NA, 4)), 
       lty=c(rep(NA, length(methods2)), lty1[-4], NA),
       pch=c(rep(NA, length(methods2)), rep(NA, 3), 1),
       seg.len=1, merge=TRUE, cex=0.75)

plot(plot.data1[[1]], type="l", xlab="Number of selected features", ylab="AUC", 
     main="b)", ylim=c(0.5, max(auc1)), xlim=c(0, 50), col=col1[3], 
     lty=lty1[1])
lines(plot.data1[[2]], col=col1[3], lty=lty1[2])
lines(plot.data1[[3]], col=col1[3], lty=lty1[3])
lines(plot.data1[[4]], col=col1[4], lty=lty1[1])
lines(plot.data1[[5]], col=col1[4], lty=lty1[2])
lines(plot.data1[[6]], col=col1[4], lty=lty1[3])
lines(plot.data1[[7]], col=col1[5], lty=lty1[1])
lines(plot.data1[[8]], col=col1[5], lty=lty1[2])
lines(plot.data1[[9]], col=col1[5], lty=lty1[3])
lines(plot.data1[[10]], col=col1[6], lty=lty1[1])
lines(plot.data1[[11]], col=col1[6], lty=lty1[2])
lines(plot.data1[[12]], col=col1[6], lty=lty1[3])
lines(plot.data1[[13]], col=col1[7], lty=lty1[1])
lines(plot.data1[[14]], col=col1[8], lty=lty1[2])
lines(plot.data1[[15]], col=col1[9], lty=lty1[3])
abline(h=auc1[colnames(auc1)=="ridge"], col=col1[1], lty=lty1[4])
abline(h=auc1[colnames(auc1)=="grridge"], col=col1[2], lty=lty1[4])
points(t(sapply(1:length(plot.data1), function(m) {
  subset(plot.data1[[m]], psel==cv.psel1[m])})), 
  col=rep(col1[-c(1, 2)], each=3), pch=1)
par(mfrow=c(1, 1))

# ---- barplots_metabolomics_alzheimer_res1_auc ----
library(Biobase)
load("data/ESetMbolCSFPR2.Rdata")
load("results/metabolomics_alzheimer_fit1.Rdata")

methods1 <- c("grridge", paste0("gren", 1:3))
col1 <- grey.colors(length(methods1), start=0.3, end=0.9, gamma=2.2, alpha=NULL)
lty1 <- c(1:4)
labels1 <- list(platform=replace(unique(as.character(fData(
  ESetMbolCSFPR2)$Platform))[-5], 4, "Oxidative Stress"))

plot.data1 <- lapply(1:length(fit1.gren1$lambdag), function(s) {
  rbind(fit1.grridge$lambdamults[[s]], fit1.gren1$lambdag[[s]],
        fit1.gren2$lambdag[[s]], fit1.gren3$lambdag[[s]])})
opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
barplot(plot.data1[[1]], beside=TRUE, col=col1,
        legend.text=c(methods1[1],
                      expression("gren, "~alpha==0.05, "gren, "~alpha==0.5,
                                 "gren, "~alpha==0.95), "no co-data"),
        args.legend=list(x="topleft", fill=c(col1, NA),
                         border=c(rep(1, length(methods1)), NA),
                         lty=c(rep(NA, length(methods1)), 2),
                         seg.len=1, merge=TRUE),
        names.arg=labels1[[1]], ylab=expression(hat(lambda)~"'"[g]))
abline(h=1, lty=2)
par(opar)

# ---- lines_metabolomics_alzheimer_res2_auc ----
library(pROC)
library(grpreg)
load("results/metabolomics_alzheimer_fit2.Rdata")
res2 <- read.table("results/metabolomics_alzheimer_res2.csv", 
                   stringsAsFactors=FALSE)

methods2 <- c("ridge", "grridge", "gren", "enet", "sgl", "cmcp", "gel")
pred2 <- as.matrix(res2[substr(rownames(res2), 1, 4)=="pred", ])
psel2 <- as.matrix(res2[substr(rownames(res2), 1, 4)=="psel", ])
auc2 <- as.matrix(res2[substr(rownames(res2), 1, 3)=="auc", ])
plot.data2 <- lapply(c(t(outer(methods2[-c(1, 2)], 1:3, paste, sep=""))), 
                     function(m) {aggregate(auc2[, grepl(m, colnames(auc2))], 
                                            list(psel=psel2[, grepl(m, colnames(
                                              psel2))]), mean)})
cv.psel2 <- c(sum(coef(fit2.gren1$freq.model$groupreg, 
                       s=fit2.gren1$lambda)[-1, ]!=0),
              sum(coef(fit2.gren2$freq.model$groupreg, 
                       s=fit2.gren2$lambda)[-1, ]!=0),
              sum(coef(fit2.gren3$freq.model$groupreg, 
                       s=fit2.gren3$lambda)[-1, ]!=0),
              sum(coef(fit2.gren1$freq.model$regular, 
                       s=fit2.gren1$lambda)[-1, ]!=0),
              sum(coef(fit2.gren2$freq.model$regular, 
                       s=fit2.gren2$lambda)[-1, ]!=0),
              sum(coef(fit2.gren3$freq.model$regular, 
                       s=fit2.gren3$lambda)[-1, ]!=0),
              psel2[, grepl("sgl1", colnames(psel2))][which.min(
                fit2.sgl1$lldiff)],
              psel2[, grepl("sgl2", colnames(psel2))][which.min(
                fit2.sgl2$lldiff)],
              psel2[, grepl("sgl3", colnames(psel2))][which.min(
                fit2.sgl3$lldiff)],
              psel2[, grepl("cmcp1", colnames(psel2))][fit2.cmcp1$min],
              psel2[, grepl("cmcp2", colnames(psel2))][fit2.cmcp2$min],
              psel2[, grepl("cmcp3", colnames(psel2))][fit2.cmcp3$min],
              psel2[, grepl("gel1", colnames(psel2))][fit2.gel1$min],
              psel2[, grepl("gel2", colnames(psel2))][fit2.gel2$min],
              psel2[, grepl("gel3", colnames(psel2))][fit2.gel3$min])
col2 <- c(1:length(methods2))
lty2 <- c(1:4)

par(mfrow=c(1, 2))
plot(plot.data2[[1]], type="l", xlab="Number of selected features", ylab="AUC", 
     main="a)", ylim=range(auc2), xlim=range(psel2), col=col2[3], lty=lty2[1])
lines(plot.data2[[2]], col=col2[3], lty=lty2[2])
lines(plot.data2[[3]], col=col2[3], lty=lty2[3])
lines(plot.data2[[4]], col=col2[4], lty=lty2[1])
lines(plot.data2[[5]], col=col2[4], lty=lty2[2])
lines(plot.data2[[6]], col=col2[4], lty=lty2[3])
lines(plot.data2[[7]], col=col2[5], lty=lty2[1])
lines(plot.data2[[8]], col=col2[5], lty=lty2[2])
lines(plot.data2[[9]], col=col2[5], lty=lty2[3])
lines(plot.data2[[10]], col=col2[6], lty=lty2[1])
lines(plot.data2[[11]], col=col2[6], lty=lty2[2])
lines(plot.data2[[12]], col=col2[6], lty=lty2[3])
lines(plot.data2[[13]], col=col2[7], lty=lty2[1])
lines(plot.data2[[14]], col=col2[8], lty=lty2[2])
lines(plot.data2[[15]], col=col2[9], lty=lty2[3])
abline(h=auc2[colnames(auc2)=="ridge"], col=col2[1], lty=lty2[4])
abline(h=auc2[colnames(auc2)=="grridge"], col=col2[2], lty=lty2[4])
points(t(sapply(1:length(plot.data2), function(m) {
  subset(plot.data2[[m]], psel==cv.psel2[m])})), 
  col=rep(col2[-c(1, 2)], each=3), pch=1)
legend("bottomright", 
       legend=c(methods2, expression(alpha==0.05, alpha==0.5, alpha==0.95),
                "CV model"), 
       fill=c(col2, rep(NA, 4)),
       border=c(rep(1, length(methods2)), rep(NA, 4)), 
       lty=c(rep(NA, length(methods2)), lty2[-4], NA),
       pch=c(rep(NA, length(methods2)), rep(NA, 3), 1),
       seg.len=1, merge=TRUE, cex=0.75)

plot(plot.data2[[1]], type="l", xlab="Number of selected features", ylab="AUC", 
     main="b)", ylim=c(0.5, max(auc2)), xlim=c(0, 50), col=col2[3], 
     lty=lty2[1])
lines(plot.data2[[2]], col=col2[3], lty=lty2[2])
lines(plot.data2[[3]], col=col2[3], lty=lty2[3])
lines(plot.data2[[4]], col=col2[4], lty=lty2[1])
lines(plot.data2[[5]], col=col2[4], lty=lty2[2])
lines(plot.data2[[6]], col=col2[4], lty=lty2[3])
lines(plot.data2[[7]], col=col2[5], lty=lty2[1])
lines(plot.data2[[8]], col=col2[5], lty=lty2[2])
lines(plot.data2[[9]], col=col2[5], lty=lty2[3])
lines(plot.data2[[10]], col=col2[6], lty=lty2[1])
lines(plot.data2[[11]], col=col2[6], lty=lty2[2])
lines(plot.data2[[12]], col=col2[6], lty=lty2[3])
lines(plot.data2[[13]], col=col2[7], lty=lty2[1])
lines(plot.data2[[14]], col=col2[8], lty=lty2[2])
lines(plot.data2[[15]], col=col2[9], lty=lty2[3])
abline(h=auc2[colnames(auc2)=="ridge"], col=col2[1], lty=lty2[4])
abline(h=auc2[colnames(auc2)=="grridge"], col=col2[2], lty=lty2[4])
points(t(sapply(1:length(plot.data2), function(m) {
  subset(plot.data2[[m]], psel==cv.psel2[m])})), 
  col=rep(col2[-c(1, 2)], each=3), pch=1)
par(mfrow=c(1, 1))

# ---- barplots_metabolomics_alzheimer_res2_auc ----
library(Biobase)
load("data/ESetMbolCSFPR2.Rdata")
load("results/metabolomics_alzheimer_fit2.Rdata")

methods2 <- c("grridge", paste0("gren", 1:3))
col2 <- grey.colors(length(methods2), start=0.3, end=0.9, gamma=2.2, alpha=NULL)
lty2 <- c(1:4)
quality <- rep(c(1:length(fit2.grridge$arguments$partitions$part)), 
               times=sapply(fit2.grridge$arguments$partitions$part, 
                            length))[order(unlist(
                              fit2.grridge$arguments$partitions$part))]
feat <- fData(ESetMbolCSFPR2)
labels2 <- list(quality=c(
  paste("RSDqc >", round(min(feat$RSDqc[quality==1]), 3)), 
  paste(round(min(feat$RSDqc[quality==1]), 3), ">= RSDqc >", 
        round(min(feat$RSDqc[quality==2]), 3)),
  paste(round(min(feat$RSDqc[quality==2]), 3), ">= RSDqc")))
plot.data2 <- lapply(1:length(fit2.gren1$lambdag), function(s) {
  rbind(fit2.grridge$lambdamults[[s]], fit2.gren1$lambdag[[s]],
        fit2.gren2$lambdag[[s]], fit2.gren3$lambdag[[s]])})
opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
barplot(plot.data2[[1]], beside=TRUE, col=col2,
        legend.text=c(methods2[1],
                      expression("gren, "~alpha==0.05, "gren, "~alpha==0.5,
                                 "gren, "~alpha==0.95), "no co-data"),
        args.legend=list(x="topleft", fill=c(col2, NA),
                         border=c(rep(1, length(methods2)), NA),
                         lty=c(rep(NA, length(methods2)), 2),
                         seg.len=1, merge=TRUE),
        names.arg=labels2[[1]], ylab=expression(hat(lambda)~"'"[g]))
abline(h=1, lty=2)
par(opar)

# ---- lines_metabolomics_alzheimer_res3_auc ----
library(pROC)
library(grpreg)
load("results/metabolomics_alzheimer_fit3.Rdata")
res3 <- read.table("results/metabolomics_alzheimer_res3.csv", 
                   stringsAsFactors=FALSE)

methods3 <- c("ridge", "grridge", "gren", "enet", "sgl", "cmcp", "gel")
pred3 <- as.matrix(res3[substr(rownames(res3), 1, 4)=="pred", ])
psel3 <- as.matrix(res3[substr(rownames(res3), 1, 4)=="psel", ])
auc3 <- as.matrix(res3[substr(rownames(res3), 1, 3)=="auc", ])
plot.data3 <- lapply(c(t(outer(methods3[-c(1, 2)], 1:3, paste, sep=""))), 
                     function(m) {aggregate(auc3[, grepl(m, colnames(auc3))], 
                                            list(psel=psel3[, grepl(m, colnames(
                                              psel3))]), mean)})
cv.psel3 <- c(sum(coef(fit3.gren1$freq.model$groupreg, 
                       s=fit3.gren1$lambda)[-1, ]!=0),
              sum(coef(fit3.gren2$freq.model$groupreg, 
                       s=fit3.gren2$lambda)[-1, ]!=0),
              sum(coef(fit3.gren3$freq.model$groupreg, 
                       s=fit3.gren3$lambda)[-1, ]!=0),
              sum(coef(fit3.gren1$freq.model$regular, 
                       s=fit3.gren1$lambda)[-1, ]!=0),
              sum(coef(fit3.gren2$freq.model$regular, 
                       s=fit3.gren2$lambda)[-1, ]!=0),
              sum(coef(fit3.gren3$freq.model$regular, 
                       s=fit3.gren3$lambda)[-1, ]!=0),
              psel3[, grepl("sgl1", colnames(psel3))][which.min(
                fit3.sgl1$lldiff)],
              psel3[, grepl("sgl2", colnames(psel3))][which.min(
                fit3.sgl2$lldiff)],
              psel3[, grepl("sgl3", colnames(psel3))][which.min(
                fit3.sgl3$lldiff)],
              psel3[, grepl("cmcp1", colnames(psel3))][fit3.cmcp1$min],
              psel3[, grepl("cmcp2", colnames(psel3))][fit3.cmcp2$min],
              psel3[, grepl("cmcp3", colnames(psel3))][fit3.cmcp3$min],
              psel3[, grepl("gel1", colnames(psel3))][fit3.gel1$min],
              psel3[, grepl("gel2", colnames(psel3))][fit3.gel2$min],
              psel3[, grepl("gel3", colnames(psel3))][fit3.gel3$min])
col3 <- c(1:length(methods3))
lty3 <- c(1:4)

par(mfrow=c(1, 2))
plot(plot.data3[[1]], type="l", xlab="Number of selected features", ylab="AUC", 
     main="a)", ylim=range(auc3), xlim=range(psel3), col=col3[3], lty=lty3[1])
lines(plot.data3[[2]], col=col3[3], lty=lty3[2])
lines(plot.data3[[3]], col=col3[3], lty=lty3[3])
lines(plot.data3[[4]], col=col3[4], lty=lty3[1])
lines(plot.data3[[5]], col=col3[4], lty=lty3[2])
lines(plot.data3[[6]], col=col3[4], lty=lty3[3])
lines(plot.data3[[7]], col=col3[5], lty=lty3[1])
lines(plot.data3[[8]], col=col3[5], lty=lty3[2])
lines(plot.data3[[9]], col=col3[5], lty=lty3[3])
lines(plot.data3[[10]], col=col3[6], lty=lty3[1])
lines(plot.data3[[11]], col=col3[6], lty=lty3[2])
lines(plot.data3[[12]], col=col3[6], lty=lty3[3])
lines(plot.data3[[13]], col=col3[7], lty=lty3[1])
lines(plot.data3[[14]], col=col3[8], lty=lty3[2])
lines(plot.data3[[15]], col=col3[9], lty=lty3[3])
abline(h=auc3[colnames(auc3)=="ridge"], col=col3[1], lty=lty3[4])
abline(h=auc3[colnames(auc3)=="grridge"], col=col3[2], lty=lty3[4])
points(t(sapply(1:length(plot.data3), function(m) {
  subset(plot.data3[[m]], psel==cv.psel3[m])})), 
  col=rep(col3[-c(1, 2)], each=3), pch=1)
legend("bottomright", 
       legend=c(methods3, expression(alpha==0.05, alpha==0.5, alpha==0.95),
                "CV model"), 
       fill=c(col3, rep(NA, 4)),
       border=c(rep(1, length(methods3)), rep(NA, 4)), 
       lty=c(rep(NA, length(methods3)), lty3[-4], NA),
       pch=c(rep(NA, length(methods3)), rep(NA, 3), 1),
       seg.len=1, merge=TRUE, cex=0.75)

plot(plot.data3[[1]], type="l", xlab="Number of selected features", ylab="AUC", 
     main="b)", ylim=c(0.5, max(auc3)), xlim=c(0, 50), col=col3[3], 
     lty=lty3[1])
lines(plot.data3[[2]], col=col3[3], lty=lty3[2])
lines(plot.data3[[3]], col=col3[3], lty=lty3[3])
lines(plot.data3[[4]], col=col3[4], lty=lty3[1])
lines(plot.data3[[5]], col=col3[4], lty=lty3[2])
lines(plot.data3[[6]], col=col3[4], lty=lty3[3])
lines(plot.data3[[7]], col=col3[5], lty=lty3[1])
lines(plot.data3[[8]], col=col3[5], lty=lty3[2])
lines(plot.data3[[9]], col=col3[5], lty=lty3[3])
lines(plot.data3[[10]], col=col3[6], lty=lty3[1])
lines(plot.data3[[11]], col=col3[6], lty=lty3[2])
lines(plot.data3[[12]], col=col3[6], lty=lty3[3])
lines(plot.data3[[13]], col=col3[7], lty=lty3[1])
lines(plot.data3[[14]], col=col3[8], lty=lty3[2])
lines(plot.data3[[15]], col=col3[9], lty=lty3[3])
abline(h=auc3[colnames(auc3)=="ridge"], col=col3[1], lty=lty3[4])
abline(h=auc3[colnames(auc3)=="grridge"], col=col3[2], lty=lty3[4])
points(t(sapply(1:length(plot.data3), function(m) {
  subset(plot.data3[[m]], psel==cv.psel3[m])})), 
  col=rep(col3[-c(1, 2)], each=3), pch=1)
par(mfrow=c(1, 1))

# ---- barplots_metabolomics_alzheimer_res3_auc ----
library(Biobase)
load("data/ESetMbolCSFPR2.Rdata")
load("results/metabolomics_alzheimer_fit3.Rdata")

methods3 <- c("grridge", paste0("gren", 1:3))
col3 <- grey.colors(length(methods3), start=0.3, end=0.9, gamma=2.2, alpha=NULL)
lty3 <- c(1:4)
quality <- rep(c(1:length(fit3.grridge$arguments$partitions$part)), 
               times=sapply(fit3.grridge$arguments$partitions$part, 
                            length))[order(unlist(
                              fit3.grridge$arguments$partitions$part))]
feat <- fData(ESetMbolCSFPR2)
labels3 <- list(degree=c("degree 0", "0 <= degree < average", 
                         "average < degree"))
plot.data3 <- lapply(1:length(fit3.gren1$lambdag), function(s) {
  rbind(fit3.grridge$lambdamults[[s]], fit3.gren1$lambdag[[s]],
        fit3.gren2$lambdag[[s]], fit3.gren3$lambdag[[s]])})
opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
barplot(plot.data3[[1]], beside=TRUE, col=col3,
        legend.text=c(methods3[1],
                      expression("gren, "~alpha==0.05, "gren, "~alpha==0.5,
                                 "gren, "~alpha==0.95), "no co-data"),
        args.legend=list(x="topright", fill=c(col3, NA),
                         border=c(rep(1, length(methods3)), NA),
                         lty=c(rep(NA, length(methods3)), 2),
                         seg.len=1, merge=TRUE),
        names.arg=labels3[[1]], ylab=expression(hat(lambda)~"'"[g]))
abline(h=1, lty=2)
par(opar)

# ---- lines_metabolomics_alzheimer_res4_auc ----
library(pROC)
library(grpreg)
load("results/metabolomics_alzheimer_fit4.Rdata")
res4 <- read.table("results/metabolomics_alzheimer_res4.csv", 
                   stringsAsFactors=FALSE)

methods4 <- c("ridge", "grridge", "gren", "enet", "sgl", "cmcp", "gel", "ocmcp",
              "ogel")
pred4 <- as.matrix(res4[substr(rownames(res4), 1, 4)=="pred", ])
psel4 <- as.matrix(res4[substr(rownames(res4), 1, 4)=="psel", ])
auc4 <- as.matrix(res4[substr(rownames(res4), 1, 3)=="auc", ])
plot.data4 <- lapply(c(t(outer(methods4[-c(1, 2)], 1:3, paste, sep=""))), 
                     function(m) {aggregate(auc4[, grepl(m, colnames(auc4))], 
                                            list(psel=psel4[, grepl(m, colnames(
                                              psel4))]), mean)})
cv.psel4 <- c(sum(coef(fit4.gren1$freq.model$groupreg, 
                       s=fit4.gren1$lambda)[-1, ]!=0),
              sum(coef(fit4.gren2$freq.model$groupreg, 
                       s=fit4.gren2$lambda)[-1, ]!=0),
              sum(coef(fit4.gren3$freq.model$groupreg, 
                       s=fit4.gren3$lambda)[-1, ]!=0),
              sum(coef(fit4.gren1$freq.model$regular, 
                       s=fit4.gren1$lambda)[-1, ]!=0),
              sum(coef(fit4.gren2$freq.model$regular, 
                       s=fit4.gren2$lambda)[-1, ]!=0),
              sum(coef(fit4.gren3$freq.model$regular, 
                       s=fit4.gren3$lambda)[-1, ]!=0),
              psel4[, grepl("sgl1", colnames(psel4))][which.min(
                fit4.sgl1$lldiff)],
              psel4[, grepl("sgl2", colnames(psel4))][which.min(
                fit4.sgl2$lldiff)],
              psel4[, grepl("sgl3", colnames(psel4))][which.min(
                fit4.sgl3$lldiff)],
              psel4[, grepl("cmcp1", colnames(psel4))][fit4.cmcp1$min],
              psel4[, grepl("cmcp2", colnames(psel4))][fit4.cmcp2$min],
              psel4[, grepl("cmcp3", colnames(psel4))][fit4.cmcp3$min],
              psel4[, grepl("gel1", colnames(psel4))][fit4.gel1$min],
              psel4[, grepl("gel2", colnames(psel4))][fit4.gel2$min],
              psel4[, grepl("gel3", colnames(psel4))][fit4.gel3$min],
              psel4[, grepl("ocmcp1", colnames(psel4))][fit4.ocmcp1$min],
              psel4[, grepl("ocmcp2", colnames(psel4))][fit4.ocmcp2$min],
              psel4[, grepl("ocmcp3", colnames(psel4))][fit4.ocmcp3$min],
              psel4[, grepl("ogel1", colnames(psel4))][fit4.ogel1$min],
              psel4[, grepl("ogel2", colnames(psel4))][fit4.ogel2$min],
              psel4[, grepl("ogel3", colnames(psel4))][fit4.ogel3$min])
col4 <- c(1:length(methods4))
lty4 <- c(1:4)

par(mfrow=c(1, 2))
plot(plot.data4[[1]], type="l", xlab="Number of selected features", ylab="AUC", 
     main="a)", ylim=range(auc4), xlim=range(psel4), col=col4[3], lty=lty4[1])
lines(plot.data4[[2]], col=col4[3], lty=lty4[2])
lines(plot.data4[[3]], col=col4[3], lty=lty4[3])
lines(plot.data4[[4]], col=col4[4], lty=lty4[1])
lines(plot.data4[[5]], col=col4[4], lty=lty4[2])
lines(plot.data4[[6]], col=col4[4], lty=lty4[3])
lines(plot.data4[[7]], col=col4[5], lty=lty4[1])
lines(plot.data4[[8]], col=col4[5], lty=lty4[2])
lines(plot.data4[[9]], col=col4[5], lty=lty4[3])
lines(plot.data4[[10]], col=col4[6], lty=lty4[1])
lines(plot.data4[[11]], col=col4[6], lty=lty4[2])
lines(plot.data4[[12]], col=col4[6], lty=lty4[3])
lines(plot.data4[[13]], col=col4[7], lty=lty4[1])
lines(plot.data4[[14]], col=col4[7], lty=lty4[2])
lines(plot.data4[[15]], col=col4[7], lty=lty4[3])
lines(plot.data4[[16]], col=col4[8], lty=lty4[1])
lines(plot.data4[[17]], col=col4[8], lty=lty4[2])
lines(plot.data4[[18]], col=col4[8], lty=lty4[3])
lines(plot.data4[[19]], col=col4[9], lty=lty4[1])
lines(plot.data4[[20]], col=col4[9], lty=lty4[2])
lines(plot.data4[[21]], col=col4[9], lty=lty4[3])
abline(h=auc4[colnames(auc4)=="ridge"], col=col4[1], lty=lty4[4])
abline(h=auc4[colnames(auc4)=="grridge"], col=col4[2], lty=lty4[4])
points(t(sapply(1:length(plot.data4), function(m) {
  subset(plot.data4[[m]], psel==cv.psel4[m])})), 
  col=rep(col4[-c(1, 2)], each=3), pch=1)
legend("bottomright", 
       legend=c(methods4, expression(alpha==0.05, alpha==0.5, alpha==0.95),
                "CV model"), 
       fill=c(col4, rep(NA, 4)),
       border=c(rep(1, length(methods4)), rep(NA, 4)), 
       lty=c(rep(NA, length(methods4)), lty4[-4], NA),
       pch=c(rep(NA, length(methods4)), rep(NA, 3), 1),
       seg.len=1, merge=TRUE, cex=0.75)

plot(plot.data4[[1]], type="l", xlab="Number of selected features", ylab="AUC", 
     main="b)", ylim=c(0.5, max(auc4)), xlim=c(0, 50), col=col4[3], 
     lty=lty4[1])
lines(plot.data4[[2]], col=col4[3], lty=lty4[2])
lines(plot.data4[[3]], col=col4[3], lty=lty4[3])
lines(plot.data4[[4]], col=col4[4], lty=lty4[1])
lines(plot.data4[[5]], col=col4[4], lty=lty4[2])
lines(plot.data4[[6]], col=col4[4], lty=lty4[3])
lines(plot.data4[[7]], col=col4[5], lty=lty4[1])
lines(plot.data4[[8]], col=col4[5], lty=lty4[2])
lines(plot.data4[[9]], col=col4[5], lty=lty4[3])
lines(plot.data4[[10]], col=col4[6], lty=lty4[1])
lines(plot.data4[[11]], col=col4[6], lty=lty4[2])
lines(plot.data4[[12]], col=col4[6], lty=lty4[3])
lines(plot.data4[[13]], col=col4[7], lty=lty4[1])
lines(plot.data4[[14]], col=col4[7], lty=lty4[2])
lines(plot.data4[[15]], col=col4[7], lty=lty4[3])
lines(plot.data4[[16]], col=col4[8], lty=lty4[1])
lines(plot.data4[[17]], col=col4[8], lty=lty4[2])
lines(plot.data4[[18]], col=col4[8], lty=lty4[3])
lines(plot.data4[[19]], col=col4[9], lty=lty4[1])
lines(plot.data4[[20]], col=col4[9], lty=lty4[2])
lines(plot.data4[[21]], col=col4[9], lty=lty4[3])
abline(h=auc4[colnames(auc4)=="ridge"], col=col4[1], lty=lty4[4])
abline(h=auc4[colnames(auc4)=="grridge"], col=col4[2], lty=lty4[4])
points(t(sapply(1:length(plot.data4), function(m) {
  subset(plot.data4[[m]], psel==cv.psel4[m])})), 
  col=rep(col4[-c(1, 2)], each=3), pch=1)
par(mfrow=c(1, 1))

# ---- barplots_metabolomics_alzheimer_res4_auc ----
library(Biobase)
load("data/ESetMbolCSFPR2.Rdata")
load("results/metabolomics_alzheimer_fit4.Rdata")

methods4 <- c("grridge", paste0("gren", 1:3))
col4 <- grey.colors(length(methods4), start=0.3, end=0.9, gamma=2.2, alpha=NULL)
lty4 <- c(1:4)
quality <- rep(c(1:length(fit4.grridge$arguments$partitions$quality)), 
               times=sapply(fit4.grridge$arguments$partitions$quality, 
                            length))[order(unlist(
                              fit4.grridge$arguments$partitions$quality))]
feat <- fData(ESetMbolCSFPR2)
labels4 <- list(platform=replace(unique(as.character(fData(
  ESetMbolCSFPR2)$Platform))[-5], 4, "Oxidative Stress"),
  quality=c(paste("RSDqc >", round(min(feat$RSDqc[quality==1]), 3)), 
            paste(round(min(feat$RSDqc[quality==1]), 3), ">= RSDqc >", 
                  round(min(feat$RSDqc[quality==2]), 3)),
            paste(round(min(feat$RSDqc[quality==2]), 3), ">= RSDqc")),
  degree=c("degree 0", "0 <= degree < average", "average < degree"))

plot.data4 <- lapply(1:length(fit4.gren1$lambdag), function(s) {
  rbind(fit4.grridge$lambdamults[[s]], fit4.gren1$lambdag[[s]],
        fit4.gren2$lambdag[[s]], fit4.gren3$lambdag[[s]])})
opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
layout(matrix(c(rep(c(1, 1, 2, 2), 2), rep(c(0, 3, 3, 0), 2)), 4, 4, 
       byrow = TRUE))
barplot(plot.data4[[1]], beside=TRUE, col=col4,
        legend.text=c(methods4[1],
                      expression("gren, "~alpha==0.05, "gren, "~alpha==0.5,
                                 "gren, "~alpha==0.95), "no co-data"),
        args.legend=list(x="topleft", fill=c(col4, NA),
                         border=c(rep(1, length(methods4)), NA),
                         lty=c(rep(NA, length(methods4)), 2),
                         seg.len=1, merge=TRUE), main="a)",
        names.arg=labels4[[1]], ylab=expression(hat(lambda)~"'"[g]))
abline(h=1, lty=2)

barplot(plot.data4[[2]], beside=TRUE, col=col4, main="b)",
        names.arg=labels4[[2]], ylab=expression(hat(lambda)~"'"[g]))
abline(h=1, lty=2)

barplot(plot.data4[[3]], beside=TRUE, col=col4, main="c)",
        names.arg=labels4[[3]], ylab=expression(hat(lambda)~"'"[g]))
abline(h=1, lty=2)
par(opar)

# ---- lines_rnaseq_lymph_node_metastasis_res1_auc ----
library(CoRF)
library(pROC)
library(grpreg)
load("results/rnaseq_lymph_node_metastasis_fit1.Rdata")
res1 <- read.table("results/rnaseq_lymph_node_metastasis_res1.csv", 
                   stringsAsFactors=FALSE)

methods1 <- c("ridge", "grridge", "gren", "enet", "cmcp", "gel")
pred1 <- as.matrix(res1[substr(rownames(res1), 1, 4)=="pred", ])
psel1 <- as.matrix(res1[substr(rownames(res1), 1, 4)=="psel", ])
ytest <- as.numeric(levels(RespValidation))[RespValidation]
auc1 <- apply(pred1, 2, function(m) {pROC::auc(ytest, m)})
plot.data1 <- lapply(c(t(outer(methods1[-c(1, 2)], 1:3, paste, sep=""))), 
                     function(m) {aggregate(auc1[grepl(m, names(auc1))], list(
                       psel=psel1[, grepl(m, colnames(psel1))]), mean)})
cv.psel1 <- c(sum(coef(fit1.gren1$freq.model$groupreg, 
                       s=fit1.gren1$lambda)[-1, ]!=0),
              sum(coef(fit1.gren2$freq.model$groupreg, 
                       s=fit1.gren2$lambda)[-1, ]!=0),
              sum(coef(fit1.gren3$freq.model$groupreg, 
                       s=fit1.gren3$lambda)[-1, ]!=0),
              sum(coef(fit1.gren1$freq.model$regular, 
                       s=fit1.gren1$lambda)[-1, ]!=0),
              sum(coef(fit1.gren2$freq.model$regular, 
                       s=fit1.gren2$lambda)[-1, ]!=0),
              sum(coef(fit1.gren3$freq.model$regular, 
                       s=fit1.gren3$lambda)[-1, ]!=0),
              psel1[, grepl("cmcp1", colnames(psel1))][fit1.cmcp1$min],
              psel1[, grepl("cmcp2", colnames(psel1))][fit1.cmcp2$min],
              psel1[, grepl("cmcp3", colnames(psel1))][fit1.cmcp3$min],
              psel1[, grepl("gel1", colnames(psel1))][fit1.gel1$min],
              psel1[, grepl("gel2", colnames(psel1))][fit1.gel2$min],
              psel1[, grepl("gel3", colnames(psel1))][fit1.gel3$min])
col1 <- c(1:length(methods1))
lty1 <- c(1:4)

par(mfrow=c(1, 2))
plot(plot.data1[[1]], type="l", xlab="Number of selected features", ylab="AUC", 
     main="a)", ylim=range(auc1), xlim=range(psel1), col=col1[3], lty=lty1[1])
lines(plot.data1[[2]], col=col1[3], lty=lty1[2])
lines(plot.data1[[3]], col=col1[3], lty=lty1[3])
lines(plot.data1[[4]], col=col1[4], lty=lty1[1])
lines(plot.data1[[5]], col=col1[4], lty=lty1[2])
lines(plot.data1[[6]], col=col1[4], lty=lty1[3])
lines(plot.data1[[7]], col=col1[5], lty=lty1[1])
lines(plot.data1[[8]], col=col1[5], lty=lty1[2])
lines(plot.data1[[9]], col=col1[5], lty=lty1[3])
lines(plot.data1[[10]], col=col1[6], lty=lty1[1])
lines(plot.data1[[11]], col=col1[6], lty=lty1[2])
lines(plot.data1[[12]], col=col1[6], lty=lty1[3])
abline(h=auc1[names(auc1)=="ridge"], col=col1[1], lty=lty1[4])
abline(h=auc1[names(auc1)=="grridge"], col=col1[2], lty=lty1[4])
points(t(sapply(1:length(plot.data1), function(m) {
  subset(plot.data1[[m]], psel==cv.psel1[m])})), 
  col=rep(col1[-c(1, 2)], each=3), pch=1)
legend("bottomright", 
       legend=c(methods1, expression(alpha==0.05, alpha==0.5, alpha==0.95),
                "CV model"), 
       fill=c(col1, rep(NA, 4)),
       border=c(rep(1, length(methods1)), rep(NA, 4)), 
       lty=c(rep(NA, length(methods1)), lty1[-4], NA),
       pch=c(rep(NA, length(methods1)), rep(NA, 3), 1),
       seg.len=1, merge=TRUE, cex=0.75)

plot(plot.data1[[1]], type="l", xlab="Number of selected features", ylab="AUC", 
     main="b)", ylim=c(0.63, max(auc1)), xlim=c(0, 500), col=col1[3], 
     lty=lty1[1])
lines(plot.data1[[2]], col=col1[3], lty=lty1[2])
lines(plot.data1[[3]], col=col1[3], lty=lty1[3])
lines(plot.data1[[4]], col=col1[4], lty=lty1[1])
lines(plot.data1[[5]], col=col1[4], lty=lty1[2])
lines(plot.data1[[6]], col=col1[4], lty=lty1[3])
lines(plot.data1[[7]], col=col1[5], lty=lty1[1])
lines(plot.data1[[8]], col=col1[5], lty=lty1[2])
lines(plot.data1[[9]], col=col1[5], lty=lty1[3])
lines(plot.data1[[10]], col=col1[6], lty=lty1[1])
lines(plot.data1[[11]], col=col1[6], lty=lty1[2])
lines(plot.data1[[12]], col=col1[6], lty=lty1[3])
abline(h=auc1[names(auc1)=="ridge"], col=col1[1], lty=lty1[4])
abline(h=auc1[names(auc1)=="grridge"], col=col1[2], lty=lty1[4])
points(t(sapply(1:length(plot.data1), function(m) {
  subset(plot.data1[[m]], psel==cv.psel1[m])})), 
  col=rep(col1[-c(1, 2)], each=3), pch=1)
par(mfrow=c(1, 1))










