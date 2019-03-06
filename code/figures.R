################################## manuscript ##################################
# ---- barplot_lines_micrornaseq_colorectal_cancer_res2_auc ----
library(sp)
load("results/micrornaseq_colorectal_cancer_fit2.Rdata")
res <- read.table("results/micrornaseq_colorectal_cancer_res2.csv", 
                  stringsAsFactors=FALSE)
pred <- as.matrix(res[substr(rownames(res), 1, 4)=="pred", ])
psel <- as.matrix(res[substr(rownames(res), 1, 4)=="psel", ])
auc <- as.matrix(res[substr(rownames(res), 1, 3)=="auc", ])

methods1 <- c("grridge", expression("gren, "~alpha==0.05, "gren, "~alpha==0.5,
                                    "gren, "~alpha==0.95), "no co-data")
methods2 <- c("ridge", "grridge", "gren1", "gren2", "enet1", "cmcp1")
col1 <- grey.colors(length(methods1) - 1, start=0.3, end=0.9, gamma=2.2, 
                    alpha=NULL)
col2 <- bpy.colors(length(methods2), cutoff.tail=0.1)
lty1 <- 2
lty2 <- c(1, 2)
labels1 <- list(diff.threegroup=expression("FDR">0.05,
                                           atop(NA, atop(textstyle(
                                             0.05>=phantom(0)),
                                             textstyle("FDR" > 0.001))),
                                           "FDR"<=0.001))
labels2 <- expression("ridge", "GRridge", "gren,"~alpha==0.05, 
                      "gren,"~alpha==0.5, "enet"~alpha==0.05, 
                      "cMCP"~alpha==0.05)

plot.data1 <- lapply(1:length(fit2.gren1$lambdag), function(s) {
  rbind(fit2.grridge$lambdamults[[s]], fit2.gren1$lambdag[[s]],
        fit2.gren2$lambdag[[s]], fit2.gren3$lambdag[[s]])})
plot.data2 <- lapply(methods2, function(m) {
  aggregate(auc[, grepl(m, colnames(auc))], 
            list(psel=psel[, grepl(m, colnames(psel))]), mean)})
names(plot.data2) <- methods2

opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
layout(matrix(rep(c(1, 1, 2, 2), 2), 2, 4, byrow=TRUE))
barplot(plot.data1[[1]], beside=TRUE, col=col1,
        legend.text=methods1,
        args.legend=list(x="bottomleft", fill=c(col1, NA),
                         border=c(rep(1, length(methods1) - 1), NA),
                         lty=c(rep(NA, length(methods1) - 1), lty1[1]),
                         seg.len=1, merge=TRUE, bg="white"),
        names.arg=labels1[[1]], ylab=expression(hat(lambda)~"'"[g]),
        main="(a)")
abline(h=1, lty=2)

plot(plot.data2[[3]], type="l", xlab="Number of selected features", ylab="AUC", 
     main="(b)", ylim=range(auc), xlim=c(0, 600), col=col2[3], 
     lty=lty2[1])
lines(plot.data2[[4]], col=col2[4], lty=lty2[1])
lines(plot.data2[[5]], col=col2[5], lty=lty2[1])
lines(plot.data2[[6]], col=col2[6], lty=lty2[1])
abline(h=plot.data2[[1]][, 2], col=col2[1], lty=lty2[2])
abline(h=plot.data2[[2]][, 2], col=col2[2], lty=lty2[2])
legend("bottomright", legend=labels2, col=col2, 
       lty=c(rep(lty2[2], 2), rep(lty2[1], 4)), bg = "white")
par(opar)

# ---- barplot_lines_rnaseq_lymph_node_metastasis_res1_auc ----
library(sp)
load("results/rnaseq_lymph_node_metastasis_fit1.Rdata")
res <- read.table("results/rnaseq_lymph_node_metastasis_res1.csv", 
                  stringsAsFactors=FALSE)
pred <- as.matrix(res[substr(rownames(res), 1, 4)=="pred", ])
psel <- as.matrix(res[substr(rownames(res), 1, 4)=="psel", ])
auc <- as.matrix(res[substr(rownames(res), 1, 3)=="auc", ])

methods1 <- c("grridge", expression("gren, "~alpha==0.05, "gren, "~alpha==0.5,
                                    "gren, "~alpha==0.95), "no co-data")
methods2 <- c("ridge", "grridge", "gren1", "enet1", "ocmcp1")
col1 <- grey.colors(length(methods1) - 1, start=0.3, end=0.9, gamma=2.2, 
                    alpha=NULL)
col2 <- bpy.colors(length(methods2), cutoff.tail=0.1)
lty1 <- 2
lty2 <- c(1, 2)
labels1 <- list(corr=expression(tau > 0.54, 
                                atop(NA, atop(textstyle(0.45 <phantom(0)), 
                                              textstyle(tau <= 0.54))),
                                atop(NA, atop(textstyle(0.36 <phantom(0)), 
                                              textstyle(tau <= 0.45))),
                                atop(NA, atop(textstyle(0.21 <phantom(0)), 
                                              textstyle(tau <= 0.36))), 
                                tau <= 0.21),
                pv=expression(p <= 7.6e-5, 
                              atop(NA, atop(textstyle(7.6e-05 <phantom(0)), 
                                            textstyle(p <= 7.2e-04))),
                              atop(NA, atop(textstyle(7.2e-04 <phantom(0)), 
                                            textstyle(p <= 4.4e-03))),
                              atop(NA, atop(textstyle(4.4e-03 <phantom(0)), 
                                            textstyle(p <= 2.1e-02))), 
                              p > 2.1e-02))
labels2 <- expression("ridge", "GRridge", "gren,"~alpha==0.05, 
                      "enet"~alpha==0.05, "cMCP,"~alpha==0.05)

plot.data1 <- lapply(1:length(fit1.gren1$lambdag), function(s) {
  rbind(fit1.grridge$lambdamults[[s]], fit1.gren1$lambdag[[s]],
        fit1.gren2$lambdag[[s]], fit1.gren3$lambdag[[s]])})
plot.data2 <- lapply(methods2, function(m) {
  aggregate(auc[, grepl(m, colnames(auc))], 
            list(psel=psel[, grepl(m, colnames(psel))]), mean)})
names(plot.data2) <- methods2

opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
layout(matrix(c(rep(c(1, 1, 2, 2), 2), rep(c(0, 3, 3, 0), 2)), 4, 4, 
              byrow=TRUE))
barplot(plot.data1[[1]], beside=TRUE, col=col1,
        legend.text=methods1,
        args.legend=list(x="bottomright", fill=c(col1, NA),
                         border=c(rep(1, length(methods1) - 1), NA),
                         lty=c(rep(NA, length(methods1) - 1), lty1[1]),
                         seg.len=1, merge=TRUE, bg="white"),
        names.arg=labels1[[1]], ylab=expression(hat(lambda)~"'"[g]),
        main="(a)")
abline(h=1, lty=2)

barplot(plot.data1[[2]], beside=TRUE, col=col1, main="(b)", 
        names.arg=labels1[[2]], ylab=expression(hat(lambda)~"'"[g]), las=2)
abline(h=1, lty=2)

plot(plot.data2[[3]], type="l", xlab="Number of selected features", ylab="AUC", 
     main="(c)", ylim=c(0.65, max(auc)), xlim=c(0, 500), col=col2[3], 
     lty=lty2[1])
lines(plot.data2[[4]], col=col2[4], lty=lty2[1])
lines(plot.data2[[5]], col=col2[5], lty=lty2[1])
abline(h=plot.data2[[1]][, 2], col=col2[1], lty=lty2[2])
abline(h=plot.data2[[2]][, 2], col=col2[2], lty=lty2[2])
legend("bottomright", legend=labels2, col=col2, 
       lty=c(rep(lty2[2], 2), rep(lty2[1], 4)), bg = "white")
par(opar)

# ---- lines_colorectal_lymph_node_gren_auc ----
library(sp)
load("results/rnaseq_lymph_node_metastasis_fit1.Rdata")
res1 <- read.table("results/rnaseq_lymph_node_metastasis_res1.csv", 
                   stringsAsFactors=FALSE)
res2 <- read.table("results/micrornaseq_colorectal_cancer_res2.csv", 
                   stringsAsFactors=FALSE)
pred1 <- as.matrix(res1[substr(rownames(res1), 1, 4)=="pred", ])
pred2 <- as.matrix(res2[substr(rownames(res2), 1, 4)=="pred", ])
psel1 <- as.matrix(res1[substr(rownames(res1), 1, 4)=="psel", ])
psel2 <- as.matrix(res2[substr(rownames(res2), 1, 4)=="psel", ])
auc1 <- as.matrix(res1[substr(rownames(res1), 1, 3)=="auc", ])
auc2 <- as.matrix(res2[substr(rownames(res2), 1, 3)=="auc", ])

methods1 <- methods2 <- paste0("gren", 1:3)
col1 <- bpy.colors(length(methods1) + 2, cutoff.tail=0.1)[-c(1, 5)]
col2 <- bpy.colors(length(methods2) + 2, cutoff.tail=0.1)[-c(1, 5)]
lty1 <- lty2 <- 1
labels1 <- labels2 <- expression("gren,"~alpha==0.05, "gren,"~alpha==0.5, 
                                 "gren,"~alpha==0.95)

plot.data1 <- lapply(methods1, function(m) {
  aggregate(auc1[, grepl(m, colnames(auc1))], 
            list(psel=psel1[, grepl(m, colnames(psel1))]), mean)})
plot.data2 <- lapply(methods2, function(m) {
  aggregate(auc2[, grepl(m, colnames(auc2))], 
            list(psel=psel2[, grepl(m, colnames(psel2))]), mean)})
names(plot.data1) <- methods1
names(plot.data2) <- methods2

opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
layout(matrix(rep(c(1, 1, 2, 2), 2), 2, 4, byrow=TRUE))
plot(plot.data1[[1]], type="l", xlab="Number of selected features", ylab="AUC", 
     main="(a)", ylim=c(0.63, max(auc1)), xlim=c(0, 300), col=col1[1], 
     lty=lty1[1])
lines(plot.data1[[2]], col=col1[2], lty=lty1[1])
lines(plot.data1[[3]], col=col1[3], lty=lty1[1])
legend("bottomright", legend=labels1, col=col1, 
       lty=rep(lty1[1], 3), bg = "white")
plot(plot.data2[[1]], type="l", xlab="Number of selected features", ylab="AUC", 
     main="(b)", ylim=range(auc2), xlim=c(0, 400), col=col2[1], 
     lty=lty2[1])
lines(plot.data2[[2]], col=col2[2], lty=lty2[1])
lines(plot.data2[[3]], col=col2[3], lty=lty2[1])
par(opar)

################################## supplement ##################################
# ---- boxplots_micrornaseq_colorectal_cancer_res3_multipliers ----
library(sp)
res <- read.table("results/micrornaseq_colorectal_cancer_res3.csv", 
                  stringsAsFactors=FALSE)

methods <- c("grridge", "gren1", "gren2", "gren3")
labels <- c("GRridge", expression(paste("gren, ", alpha==0.05)),
            expression(paste("gren, ", alpha==0.5)),
            expression(paste("gren, ", alpha==0.95)))
col <- bpy.colors(length(methods), cutoff.tail=0.1)

plot.data <- data.frame(multipliers=unlist(res), 
                        group=rep(rep(1:3, each=100), 4),
                        method=rep(1:4, each=300))

opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
boxplot(multipliers ~ interaction(method, group), data=plot.data,
        at=c(c(1:4), c(6:9), c(11:14)), xaxt="n", col=rep(col, 3),
        outline=TRUE, cex.lab=2, cex.names=1.5,
        ylab=expression({lambda^{"'"}}[g]),
        cex.axis=1.5, boxlwd=0.5)
axis(1, c(2.5, 7.5, 12.5), c("Group 1", "Group 2", "Group 3"), tick=FALSE,
     cex.axis=1.5)
abline(h=1, lty=2, lwd=1.5)
legend("bottomleft", legend=labels, fill=col, border=rep(1, 4), seg.len=1, cex=1.3)
par(opar)

# ---- boxplots_micrornaseq_colorectal_cancer_res4_multipliers ----
res <- read.table("results/micrornaseq_colorectal_cancer_res4.csv", 
                  stringsAsFactors=FALSE)

methods <- c("grridge", "gren1", "gren2", "gren3")

labels <- expression("GRridge", 
                     atop(NA, atop(textstyle("gren"), 
                                   textstyle(alpha==0.05))),
                     atop(NA, atop(textstyle("gren"), 
                                   textstyle(alpha==0.5))),
                     atop(NA, atop(textstyle("gren"), 
                                   textstyle(alpha==0.95))))
col <- bpy.colors(length(methods), cutoff.tail=0.1)

plot.data <- data.frame(multipliers=unlist(res), 
                        method=rep(1:4, each=1000))

opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1/1.3))
layout(matrix(rep(c(1, 1, 2, 2), 2), 2, 4, byrow=TRUE))
boxplot(multipliers ~ method, data=plot.data, names=labels,
        col=col, outline=TRUE, main="(a)",
        ylab=expression({lambda^{"'"}}[g]), las=2,
        boxlwd=0.5, cex.lab=2, cex.names=1.5, cex.axis=1.3)
abline(h=1, lty=2, lwd=1.5)

boxplot(multipliers ~ method, data=plot.data, names=labels,
        col=col, outline=FALSE, main="(b)",
        ylab=expression({lambda^{"'"}}[g]), las=2,
        boxlwd=0.5, cex.lab=2, cex.names=1.5, cex.axis=1.3)
abline(h=1, lty=2, lwd=1.5)
par(opar)

# ---- hist_micrornaseq_colorectal_cancer_res5_overlap ----

# ---- lines_micrornaseq_colorectal_cancer_res2_auc ----
library(sp)
res <- read.table("results/micrornaseq_colorectal_cancer_res2.csv", 
                   stringsAsFactors=FALSE)

pred <- as.matrix(res[substr(rownames(res), 1, 4)=="pred", ])
psel <- as.matrix(res[substr(rownames(res), 1, 4)=="psel", ])
auc <- as.matrix(res[substr(rownames(res), 1, 3)=="auc", ])

methods <- c("ridge", "grridge", paste0("gren", c(1:3)), paste0("enet", c(1:3)),
             paste0("sgl", c(1:3)), paste0("cmcp", c(1:3)), 
             paste0("gel", c(1:3)))
col <- bpy.colors(length(methods), cutoff.tail=0.1)
lty <- 1:2
labels <- expression("ridge", "GRridge", "gren,"~alpha==0.05, 
                     "gren,"~alpha==0.5, "gren,"~alpha==0.95, 
                     "enet"~alpha==0.05, "enet"~alpha==0.5, "enet"~alpha==0.95, 
                     "SGL"~alpha==0.05, "SGL"~alpha==0.5, "SGL"~alpha==0.95,
                     "cMCP"~alpha==0.05, "cMCP"~alpha==0.5, "cMCP"~alpha==0.95,
                     "gel"~alpha==0.05, "gel"~alpha==0.5, "gel"~alpha==0.95)


plot.data <- lapply(methods, function(m) {
  aggregate(auc[, grepl(m, colnames(auc))], 
            list(psel=psel[, grepl(m, colnames(psel))]), mean)})
names(plot.data) <- methods

opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
plot(plot.data[[3]], type="l", xlab="Number of selected features", ylab="AUC", 
     main="(b)", ylim=range(auc), xlim=c(0, 500), col=col[3], 
     lty=lty[1])
lines(plot.data[[4]], col=col[4], lty=lty[1])
lines(plot.data[[5]], col=col[5], lty=lty[1])
lines(plot.data[[6]], col=col[6], lty=lty[1])
lines(plot.data[[7]], col=col[7], lty=lty[1])
lines(plot.data[[8]], col=col[8], lty=lty[1])
lines(plot.data[[9]], col=col[9], lty=lty[1])
lines(plot.data[[10]], col=col[10], lty=lty[1])
lines(plot.data[[11]], col=col[11], lty=lty[1])
lines(plot.data[[12]], col=col[12], lty=lty[1])
lines(plot.data[[13]], col=col[13], lty=lty[1])
lines(plot.data[[14]], col=col[14], lty=lty[1])
lines(plot.data[[15]], col=col[15], lty=lty[1])
lines(plot.data[[16]], col=col[16], lty=lty[1])
lines(plot.data[[17]], col=col[17], lty=lty[1])
abline(h=plot.data[[1]][, 2], col=col[1], lty=lty[2])
abline(h=plot.data[[2]][, 2], col=col[2], lty=lty[2])
legend("bottomright", legend=labels, col=col, 
       lty=c(rep(lty[2], 2), rep(lty[2], 15)), bg = "white")
par(opar)



str(plot.data)





# ---- lines_metabolomics_alzheimer_res1_auc ----
library(pROC)
library(grpreg)
library(sp)
load("results/metabolomics_alzheimer_fit1.Rdata")
res1 <- read.table("results/metabolomics_alzheimer_res1.csv", 
                   stringsAsFactors=FALSE)

methods1 <- c("ridge", "grridge", "gren", "enet", "sgl", "cmcp", "gel")
pred1 <- as.matrix(res1[substr(rownames(res1), 1, 4)=="pred", ])
psel1 <- as.matrix(res1[substr(rownames(res1), 1, 4)=="psel", ])
auc1 <- as.matrix(res1[substr(rownames(res1), 1, 3)=="auc", ])
plot.data1 <- lapply(c(t(outer(methods1[-c(1, 2)], 1:3, paste, sep=""))), 
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
col1 <- bpy.colors(length(methods1), cutoff.tail=0.1)
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
       legend=c(methods1, expression(alpha==0.05, alpha==0.5, alpha==0.95),
                "CV model"), 
       fill=c(col1, rep(NA, 4)),
       border=c(rep(1, length(methods1)), rep(NA, 4)), 
       lty=c(rep(NA, length(methods1)), lty1[-4], NA),
       pch=c(rep(NA, length(methods1)), rep(NA, 3), 1),
       seg.len=1, merge=TRUE, cex=0.75, bg="white")

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
        args.legend=list(x="topright", fill=c(col1, NA),
                         border=c(rep(1, length(methods1)), NA),
                         lty=c(rep(NA, length(methods1)), 2),
                         seg.len=1, merge=TRUE),
        names.arg=labels1[[1]], ylab=expression(hat(lambda)~"'"[g]))
abline(h=1, lty=2)
par(opar)

# ---- lines_metabolomics_alzheimer_res2_auc ----
library(pROC)
library(grpreg)
library(sp)
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
col2 <- bpy.colors(length(methods2), cutoff.tail=0.1)
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
        args.legend=list(x="topright", fill=c(col2, NA),
                         border=c(rep(1, length(methods2)), NA),
                         lty=c(rep(NA, length(methods2)), 2),
                         seg.len=1, merge=TRUE),
        names.arg=labels2[[1]], ylab=expression(hat(lambda)~"'"[g]))
abline(h=1, lty=2)
par(opar)

# ---- lines_metabolomics_alzheimer_res3_auc ----
library(pROC)
library(grpreg)
library(sp)
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
col3 <- bpy.colors(length(methods3), cutoff.tail=0.1)
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
library(sp)
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
col4 <- bpy.colors(length(methods4), cutoff.tail=0.1)
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

# ---- barplot_lines_metabolomics_alzheimer_res5_auc ----
library(sp)
load("results/metabolomics_alzheimer_fit5.Rdata")
res5 <- read.table("results/metabolomics_alzheimer_res5.csv", 
                   stringsAsFactors=FALSE)

methods1 <- c("grridge", expression("gren, "~alpha==0.05, "gren, "~alpha==0.5,
                                    "gren, "~alpha==0.95), "no co-data")
methods2 <- c("ridge", "grridge", "gren", "enet", "sgl", "cmcp", "gel", "ocmcp",
              "ogel")
pred <- as.matrix(res5[substr(rownames(res5), 1, 4)=="pred", ])
psel <- as.matrix(res5[substr(rownames(res5), 1, 4)=="psel", ])
auc <- as.matrix(res5[substr(rownames(res5), 1, 3)=="auc", ])
briers <- as.matrix(res5[substr(rownames(res5), 1, 6)=="briers", ])

plot.data1 <- lapply(1:length(fit5.gren1$lambdag), function(s) {
  rbind(fit5.grridge$lambdamults[[s]], fit5.gren1$lambdag[[s]],
        fit5.gren2$lambdag[[s]], fit5.gren3$lambdag[[s]])})
plot.data2 <- lapply(c(t(outer(methods2[-c(1, 2)], 1:3, paste, sep=""))), 
                     function(m) {aggregate(auc[, grepl(m, colnames(auc))], 
                                            list(psel=psel[, grepl(m, colnames(
                                              psel))]), mean)})
plot.data3 <- lapply(c(t(outer(methods2[-c(1, 2)], 1:3, paste, sep=""))), 
                     function(m) {aggregate(briers[, grepl(m, 
                                                           colnames(briers))], 
                                            list(psel=psel[, grepl(m, colnames(
                                              psel))]), mean)})
names(plot.data2) <- names(plot.data3) <- 
  c(t(outer(methods2[-c(1, 2)], 1:3, paste, sep="")))

col1 <- grey.colors(length(methods1) - 1, start=0.3, end=0.9, gamma=2.2, 
                    alpha=NULL)
col2 <- bpy.colors(length(methods2), cutoff.tail=0.1)
lty1 <- 2
lty2 <- c(1:4)

labels <- list(quality=expression("RSDqc" > 0.096, 
                                  atop(NA, atop(textstyle(0.096 >=phantom(0)), 
                                                textstyle("RSDqc" > 0.047))),
                                  "RSDqc" <= 0.047),
               degree=expression("degree" == 0, 
                                 atop(NA, atop(textstyle(0 < "degree"),
                                               textstyle(
                                                 phantom(0)<= "average"))),
                                 "degree" > "average"))

opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
layout(matrix(c(rep(c(1, 1, 2, 2), 2), rep(c(3, 3, 4, 4), 2)), 4, 4, 
              byrow=TRUE))
barplot(plot.data1[[1]], beside=TRUE, col=col1, main="(a)",
        names.arg=labels[[1]], ylab=expression(hat(lambda)~"'"[g]))
abline(h=1, lty=2)

barplot(plot.data1[[2]], beside=TRUE, col=col1, main="(b)", 
        legend.text=methods1, 
        args.legend=list(x="bottomright", fill=c(col1, NA),
                         border=c(rep(1, length(methods1)), NA),
                         lty=c(rep(NA, length(methods1)), 2),
                         seg.len=1, merge=TRUE, bg="white"), 
        names.arg=labels[[2]], ylab=expression(hat(lambda)~"'"[g]))
abline(h=1, lty=2)

plot(plot.data2[[1]], type="l", xlab="Number of selected features", ylab="AUC", 
     main="(c)", ylim=c(0.6, max(auc)), xlim=c(0, 100), col=col2[3], 
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
lines(plot.data2[[14]], col=col2[7], lty=lty2[2])
lines(plot.data2[[15]], col=col2[7], lty=lty2[3])
lines(plot.data2[[16]], col=col2[8], lty=lty2[1])
lines(plot.data2[[17]], col=col2[8], lty=lty2[2])
lines(plot.data2[[18]], col=col2[8], lty=lty2[3])
lines(plot.data2[[19]], col=col2[9], lty=lty2[1])
lines(plot.data2[[20]], col=col2[9], lty=lty2[2])
lines(plot.data2[[21]], col=col2[9], lty=lty2[3])
abline(h=auc[colnames(auc)=="ridge"], col=col2[1], lty=lty2[4])
abline(h=auc[colnames(auc)=="grridge"], col=col2[2], lty=lty2[4])
legend("bottomright", 
       legend=c(methods2, expression(alpha==0.05, alpha==0.5, 
                                     alpha==0.95, alpha==0)), 
       fill=c(col2, rep(NA, 4)),
       border=c(rep(1, length(methods2)), rep(NA, 4)), 
       lty=c(rep(NA, length(methods2)), lty2),
       seg.len=1, merge=TRUE, bg="white")

plot(plot.data3[[1]], type="l", xlab="Number of selected features", 
     ylab="Brier skill score", 
     main="(d)", ylim=range(briers), xlim=c(0, 100), col=col2[3], 
     lty=lty2[1])
lines(plot.data3[[2]], col=col2[3], lty=lty2[2])
lines(plot.data3[[3]], col=col2[3], lty=lty2[3])
lines(plot.data3[[4]], col=col2[4], lty=lty2[1])
lines(plot.data3[[5]], col=col2[4], lty=lty2[2])
lines(plot.data3[[6]], col=col2[4], lty=lty2[3])
lines(plot.data3[[7]], col=col2[5], lty=lty2[1])
lines(plot.data3[[8]], col=col2[5], lty=lty2[2])
lines(plot.data3[[9]], col=col2[5], lty=lty2[3])
lines(plot.data3[[10]], col=col2[6], lty=lty2[1])
lines(plot.data3[[11]], col=col2[6], lty=lty2[2])
lines(plot.data3[[12]], col=col2[6], lty=lty2[3])
lines(plot.data3[[13]], col=col2[7], lty=lty2[1])
lines(plot.data3[[14]], col=col2[7], lty=lty2[2])
lines(plot.data3[[15]], col=col2[7], lty=lty2[3])
lines(plot.data3[[16]], col=col2[8], lty=lty2[1])
lines(plot.data3[[17]], col=col2[8], lty=lty2[2])
lines(plot.data3[[18]], col=col2[8], lty=lty2[3])
lines(plot.data3[[19]], col=col2[9], lty=lty2[1])
lines(plot.data3[[20]], col=col2[9], lty=lty2[2])
lines(plot.data3[[21]], col=col2[9], lty=lty2[3])
abline(h=briers[colnames(briers)=="ridge"], col=col2[1], lty=lty2[4])
abline(h=briers[colnames(briers)=="grridge"], col=col2[2], lty=lty2[4])
par(opar)

# ---- lines_rnaseq_lymph_node_metastasis_res1_auc ----
library(CoRF)
library(pROC)
library(grpreg)
library(sp)
load("results/rnaseq_lymph_node_metastasis_fit1.Rdata")
res1 <- read.table("results/rnaseq_lymph_node_metastasis_res1.csv", 
                   stringsAsFactors=FALSE)

methods1 <- c("ridge", "grridge", "gren", "enet", "sgl", "cmcp", "gel")
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
col1 <- bpy.colors(length(methods1), cutoff.tail=0.1)
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

# ---- barplots_rnaseq_lymph_node_metastasis_res1_auc ----
library(sp)
library(Biobase)
load("results/rnaseq_lymph_node_metastasis_fit1.Rdata")

methods1 <- c("grridge", paste0("gren", 1:3))
col1 <- grey.colors(length(methods1), start=0.3, end=0.9, gamma=2.2, alpha=NULL)
lty1 <- c(1:4)
labels1 <- list(corr=paste0("corr", c(1:5)),
                pv=paste0("pv", c(1:5)))

plot.data1 <- lapply(1:length(fit1.gren1$lambdag), function(s) {
  rbind(fit1.grridge$lambdamults[[s]], fit1.gren1$lambdag[[s]],
        fit1.gren2$lambdag[[s]], fit1.gren3$lambdag[[s]])})
opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
layout(matrix(rep(c(1, 1, 2, 2), 2), 2, 4, byrow=TRUE))
barplot(plot.data1[[1]], beside=TRUE, col=col1, main="a)",
        names.arg=labels1[[1]], ylab=expression(hat(lambda)~"'"[g]))
abline(h=1, lty=2)

barplot(plot.data1[[2]], beside=TRUE, col=col1, 
        legend.text=c(methods1[1],
                      expression("gren, "~alpha==0.05, "gren, "~alpha==0.5,
                                 "gren, "~alpha==0.95), "no co-data"),
        args.legend=list(x="topleft", fill=c(col1, NA),
                         border=c(rep(1, length(methods1)), NA),
                         lty=c(rep(NA, length(methods1)), 2),
                         seg.len=1, merge=TRUE), main="b)",
        names.arg=labels1[[2]], ylab=expression(hat(lambda)~"'"[g]))
abline(h=1, lty=2)
par(opar)

# ---- lines_micrornaseq_colorectal_cancer_res1_auc ----
library(pROC)
library(grpreg)
library(sp)
load("results/micrornaseq_colorectal_cancer_fit1.Rdata")
res1 <- read.table("results/micrornaseq_colorectal_cancer_res1.csv", 
                   stringsAsFactors=FALSE)

methods1 <- c("ridge", "grridge", "gren", "enet", "sgl", "cmcp", "gel")
pred1 <- as.matrix(res1[substr(rownames(res1), 1, 4)=="pred", ])
psel1 <- as.matrix(res1[substr(rownames(res1), 1, 4)=="psel", ])
auc1 <- as.matrix(res1[substr(rownames(res1), 1, 3)=="auc", ])
plot.data1 <- lapply(c(t(outer(methods1[-c(1, 2)], 1:3, paste, sep=""))), 
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
col1 <- bpy.colors(length(methods1), cutoff.tail=0.1)
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
       legend=c(methods1, expression(alpha==0.05, alpha==0.5, alpha==0.95),
                "CV model"), 
       fill=c(col1, rep(NA, 4)),
       border=c(rep(1, length(methods1)), rep(NA, 4)), 
       lty=c(rep(NA, length(methods1)), lty1[-4], NA),
       pch=c(rep(NA, length(methods1)), rep(NA, 3), 1),
       seg.len=1, merge=TRUE, cex=0.75)

plot(plot.data1[[1]], type="l", xlab="Number of selected features", ylab="AUC", 
     main="b)", ylim=c(0.5, max(auc1)), xlim=c(0, 500), col=col1[3], 
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

# ---- barplots_micrornaseq_colorectal_cancer_res1_auc ----
library(Biobase)
load("results/micrornaseq_colorectal_cancer_fit1.Rdata")

methods1 <- c("grridge", paste0("gren", 1:3))
col1 <- grey.colors(length(methods1), start=0.3, end=0.9, gamma=2.2, alpha=NULL)
lty1 <- c(1:4)
labels1 <- list(diff.twogroup=c("FDR > 0.05", "FDR <= 0.05"))

plot.data1 <- lapply(1:length(fit1.gren1$lambdag), function(s) {
  rbind(fit1.grridge$lambdamults[[s]], fit1.gren1$lambdag[[s]],
        fit1.gren2$lambdag[[s]], fit1.gren3$lambdag[[s]])})
opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
barplot(plot.data1[[1]], beside=TRUE, col=col1,
        legend.text=c(methods1[1],
                      expression("gren, "~alpha==0.05, "gren, "~alpha==0.5,
                                 "gren, "~alpha==0.95), "no co-data"),
        args.legend=list(x="topright", fill=c(col1, NA),
                         border=c(rep(1, length(methods1)), NA),
                         lty=c(rep(NA, length(methods1)), 2),
                         seg.len=1, merge=TRUE),
        names.arg=labels1[[1]], ylab=expression(hat(lambda)~"'"[g]))
abline(h=1, lty=2)
par(opar)

# ---- lines_micrornaseq_colorectal_cancer_res2_auc ----
library(pROC)
library(grpreg)
library(sp)
load("results/micrornaseq_colorectal_cancer_fit2.Rdata")
res2 <- read.table("results/micrornaseq_colorectal_cancer_res2.csv", 
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
col2 <- bpy.colors(length(methods2), cutoff.tail=0.1)
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
     main="b)", ylim=c(0.5, max(auc2)), xlim=c(0, 500), col=col2[3], 
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

# ---- barplots_micrornaseq_colorectal_cancer_res2_auc ----
library(Biobase)
load("results/micrornaseq_colorectal_cancer_fit2.Rdata")

methods2 <- c("grridge", paste0("gren", 1:3))
col2 <- grey.colors(length(methods2), start=0.3, end=0.9, gamma=2.2, alpha=NULL)
lty2 <- c(1:4)
labels2 <- list(diff.threegroup=c("FDR > 0.05", "0.05 >= FDR > 0.001",
                                  "FDR <= 0.001"))

plot.data2 <- lapply(1:length(fit2.gren1$lambdag), function(s) {
  rbind(fit2.grridge$lambdamults[[s]], fit2.gren1$lambdag[[s]],
        fit2.gren2$lambdag[[s]], fit2.gren3$lambdag[[s]])})
opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
barplot(plot.data2[[1]], beside=TRUE, col=col2,
        legend.text=c(methods2[1],
                      expression("gren, "~alpha==0.05, "gren, "~alpha==0.5,
                                 "gren, "~alpha==0.95), "no co-data"),
        args.legend=list(x="topright", fill=c(col2, NA),
                         border=c(rep(1, length(methods2)), NA),
                         lty=c(rep(NA, length(methods2)), 2),
                         seg.len=1, merge=TRUE),
        names.arg=labels2[[1]], ylab=expression(hat(lambda)~"'"[g]))
abline(h=1, lty=2)
par(opar)

# ---- lines_micrornaseq_cervical_cancer_res1_auc ----
library(pROC)
library(grpreg)
library(sp)
load("results/micrornaseq_cervical_cancer_fit1.Rdata")
res1 <- read.table("results/micrornaseq_cervical_cancer_res1.csv", 
                   stringsAsFactors=FALSE)

methods1 <- c("ridge", "grridge", "gren", "enet")
pred1 <- as.matrix(res1[substr(rownames(res1), 1, 4)=="pred", ])
psel1 <- as.matrix(res1[substr(rownames(res1), 1, 4)=="psel", ])
auc1 <- as.matrix(res1[substr(rownames(res1), 1, 3)=="auc", ])

plot.data1 <- lapply(c(t(outer(methods1[-c(1, 2)], 1:3, paste, sep=""))), 
                     function(m) {cbind(colMeans(psel1[, grep(m, colnames(
                       psel1))]), auc1[, grep(m, colnames(auc1))])})
col1 <- bpy.colors(length(methods1) + 2, cutoff.tail=0.1)[
  -c(1, length(methods1) + 2)]
lty1 <- c(1:4)

plot(plot.data1[[1]], type="l", xlab="Number of selected features", ylab="AUC", 
     ylim=range(auc1), xlim=range(psel1[, -c(1, 2)]), col=col1[3], lty=lty1[1])
lines(plot.data1[[2]], col=col1[3], lty=lty1[2])
lines(plot.data1[[3]], col=col1[3], lty=lty1[3])
lines(plot.data1[[4]], col=col1[4], lty=lty1[1])
lines(plot.data1[[5]], col=col1[4], lty=lty1[2])
lines(plot.data1[[6]], col=col1[4], lty=lty1[3])
abline(h=auc1[colnames(auc1)=="ridge"], col=col1[1], lty=lty1[4])
abline(h=auc1[colnames(auc1)=="grridge"], col=col1[2], lty=lty1[4])
legend("bottomright", 
       legend=c(methods1, expression(alpha==0.05, alpha==0.5, alpha==0.95),
                "CV model"), 
       fill=c(col1, rep(NA, 4)),
       border=c(rep(1, length(methods1)), rep(NA, 4)), 
       lty=c(rep(NA, length(methods1)), lty1[-4], NA),
       pch=c(rep(NA, length(methods1)), rep(NA, 3), 1),
       seg.len=1, merge=TRUE, cex=0.75)

# ---- barplots_micrornaseq_cervical_cancer_res1_auc ----
load("results/micrornaseq_cervical_cancer_fit1.Rdata")

methods1 <- c("grridge", paste0("gren", 1:3))
col1 <- grey.colors(length(methods1), start=0.3, end=0.9, gamma=2.2, alpha=NULL)
lty1 <- c(1:4)
labels1 <- list(conservation=c("not-conserved", "mammals", "most vertebrates"))

plot.data1 <- lapply(1:length(fit1.gren1$lambdag), function(s) {
  rbind(fit1.grridge$lambdamults[[s]], fit1.gren1$lambdag[[s]],
        fit1.gren2$lambdag[[s]], fit1.gren3$lambdag[[s]])})
opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
barplot(plot.data1[[1]], beside=TRUE, col=col1,
        legend.text=c(methods1[1],
                      expression("gren, "~alpha==0.05, "gren, "~alpha==0.5,
                                 "gren, "~alpha==0.95), "no co-data"),
        args.legend=list(x="topright", fill=c(col1, NA),
                         border=c(rep(1, length(methods1)), NA),
                         lty=c(rep(NA, length(methods1)), 2),
                         seg.len=1, merge=TRUE),
        names.arg=labels1[[1]], ylab=expression(hat(lambda)~"'"[g]))
abline(h=1, lty=2)
par(opar)
