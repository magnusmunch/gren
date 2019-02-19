################################## manuscript ##################################
# ---- figures ----
# ---- boxplot_metabolomics_alzheimer_res1_auc ---- 
library(pROC)
library(Biobase)
res1 <- read.table("results/metabolomics_alzheimer_res1.csv", 
                   stringsAsFactors=FALSE)
load("data/ESetMbolCSFPR2.Rdata")
feat <- fData(ESetMbolCSFPR2)
pheno <- pData(ESetMbolCSFPR2)
pheno.apoe <- pheno[(pheno$D_diag_name=="Probable AD" & pheno$APOE=="E4YES") |
                      (pheno$D_diag_name=="Subjectieve klachten" & 
                         pheno$APOE=="E4NO"), ]
alzheim.apoe <- as.numeric(pheno.apoe$D_diag_name) - 1

auc1 <- sapply(res1[grepl("pred", rownames(res1)), ], function(s) {
  pROC::auc(alzheim.apoe, s)})
psel1 <- colMeans(res1[grepl("psel", rownames(res1)), ])

labels1 <- c(unique(as.character(feat$Platform))[-c(4, 5)], "Oxidative Stress")
plot(psel1, auc1)
boxplot(res1[grepl("psel", rownames(res1)), ])

barplot(rbind(fit1.gren1$lambdag$part, fit1.grridge$lambdamults$part), 
        beside=TRUE)
abline(h=1, lty=2)

# ---- lines_rnaseq_lymph_node_metastasis_res1_auc ----
library(CoRF)
library(pROC)
library(grpreg)
load("../results/rnaseq_lymph_node_metastasis_fit1.Rdata")
res1 <- read.table("../results/rnaseq_lymph_node_metastasis_res1.csv", 
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
       lty=c(rep(NA, length(methods1)), lty1, NA),
       pch=c(rep(NA, length(methods1)), rep(NA, 1), 4),
       seg.len=1, merge=TRUE)

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










