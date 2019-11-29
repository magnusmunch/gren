################################## manuscript ##################################
# ---- lines_simulations_res1_auc_briers ----
library(sp)
load(file="results/simulations_res1.Rdata")
plot.data <- c(sapply(c("grridge", "ridge"), function(m) {
  unname(cbind(NA, mean(sapply(res, function(s) {
    s$auc[substr(names(s$auc), 1, nchar(m))==m]}), na.rm=TRUE),
    mean(sapply(res, function(s) {
      s$auc[substr(names(s$briers), 1, nchar(m))==m]}), na.rm=TRUE)))}, 
  simplify=FALSE),
  sapply(c(paste0("gren", 1:3), paste0("enet", 1:3),
           paste0("cmcp", 1:3), paste0("gel", 1:3)), function(m) {
             cbind(rowMeans(sapply(res, function(s) {
               vec <- s$psel[substr(names(s$psel), 1, nchar(m))==m]
               replace(rep(NA, 8), 1:length(vec), vec)}), 
               na.rm=TRUE),
               rowMeans(sapply(res, function(s) {
                 vec <- s$auc[substr(names(s$auc), 1, nchar(m))==m]
                 replace(rep(NA, 8), 1:length(vec), vec)}), 
                 na.rm=TRUE),
               rowMeans(sapply(res, function(s) {
                 vec <- s$briers[substr(names(s$briers), 1, nchar(m))==m]
                 replace(rep(NA, 8), 1:length(vec), vec)}), 
                 na.rm=TRUE))}, simplify=FALSE))

methods <- c("GRridge", "ridge", "gren", "enet", "cmcp", "gel")
labels <- expression("GRridge", "ridge", "gren,"~alpha==0.05, 
                     "enet"~alpha==0.05, "cMCP,"~alpha==0.05, 
                     "gel,"~alpha==0.05)

col <- bpy.colors(length(methods), cutoff.tail=0.1)
lty <- c(1)

xlim <- range(sapply(plot.data[-c(1, 2, 4, 5, 7, 8, 10, 11, 13, 14)], 
                     function(s) {s[, 1]}), na.rm=TRUE)
ylim1 <- range(sapply(plot.data[-c(4, 5, 7, 8, 10, 11, 13, 14)], 
                      function(s) {s[, 2]}), na.rm=TRUE)
ylim2 <- range(sapply(plot.data[-c(4, 5, 7, 8, 10, 11, 13, 14, 9, 10, 11)], 
                      function(s) {s[, 3]}), na.rm=TRUE)

opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
layout(matrix(rep(c(1, 1, 2, 2), 2), 2, 4, byrow=TRUE))
plot(plot.data[[3]][, c(1, 2)], type="l", main="(a)", ylab="AUC", 
     xlab="Number of selected features", ylim=ylim1, xlim=xlim, col=col[3], 
     lty=lty[1])
lines(plot.data[[6]][, c(1, 2)], col=col[4], lty=lty[1])
lines(plot.data[[9]][, c(1, 2)], col=col[5], lty=lty[1])
lines(plot.data[[12]][, c(1, 2)], col=col[6], lty=lty[1])
abline(h=plot.data[[1]][, 2], col=col[1], lty=lty[1])
abline(h=plot.data[[2]][, 2], col=col[2], lty=lty[1])
legend("bottomright", legend=labels, col=col, lty=lty[1],
       seg.len=1, merge=TRUE, bg="white")
plot(plot.data[[3]][, c(1, 3)], type="l", main="(b)", ylab="Brier skill score", 
     xlab="Number of selected features", ylim=ylim2, xlim=xlim, col=col[3], 
     lty=lty[1])
lines(plot.data[[6]][, c(1, 3)], col=col[4], lty=lty[1])
lines(plot.data[[9]][, c(1, 3)], col=col[5], lty=lty[1])
lines(plot.data[[12]][, c(1, 3)], col=col[6], lty=lty[1])
abline(h=plot.data[[1]][, 3], col=col[1], lty=lty[1])
abline(h=plot.data[[2]][, 3], col=col[2], lty=lty[1])
par(opar)

# ---- lines_simulations_res2_auc_briers ----
library(sp)
load(file="results/simulations_res2.Rdata")
mults <- sapply(c("grridge", paste("gren", 1:3, sep="")), function(m) {
  sapply(res, function(r) {r$mults[grep(m, names(r$mults))]})}, simplify=FALSE)
plot.data <- c(sapply(c("grridge", "ridge"), function(m) {
  unname(cbind(NA, mean(sapply(res, function(s) {
    s$auc[grep(m, names(s$auc))[1]]}), na.rm=TRUE),
    mean(sapply(res, function(s) {
      s$auc[grep(m, names(s$briers))[1]]}), na.rm=TRUE)))}, simplify=FALSE),
  sapply(c(paste0("gren", 1:3), paste0("enet", 1:3),
           paste0("cmcp", 1:3), paste0("gel", 1:3)), function(m) {
             cbind(rowMeans(sapply(res, function(s) {
               vec <- s$psel[grep(m, names(s$psel))]
               replace(rep(NA, 8), 1:length(vec), vec)}), 
               na.rm=TRUE),
               rowMeans(sapply(res, function(s) {
                 vec <- s$auc[grep(m, names(s$auc))]
                 replace(rep(NA, 8), 1:length(vec), vec)}), 
                 na.rm=TRUE),
               rowMeans(sapply(res, function(s) {
                 vec <- s$briers[grep(m, names(s$briers))]
                 replace(rep(NA, 8), 1:length(vec), vec)}), 
                 na.rm=TRUE))}, simplify=FALSE))

methods <- c("GRridge", "ridge", "gren", "enet", "cmcp", "gel")
labels <- expression("GRridge", "ridge", "gren,"~alpha==0.05, 
                     "enet"~alpha==0.05, "cMCP,"~alpha==0.05, 
                     "gel,"~alpha==0.05)

col <- bpy.colors(length(methods), cutoff.tail=0.1)
lty <- c(1)

xlim <- range(sapply(plot.data[-c(1, 2, 4, 5, 7, 8, 10, 11, 13, 14)], 
                     function(s) {s[, 1]}), na.rm=TRUE)
ylim1 <- range(sapply(plot.data[-c(4, 5, 7, 8, 10, 11, 13, 14)], 
                      function(s) {s[, 2]}), na.rm=TRUE)
ylim2 <- range(sapply(plot.data[-c(4, 5, 7, 8, 10, 11, 13, 14, 9, 10, 11)], 
                      function(s) {s[, 3]}), na.rm=TRUE)

opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
layout(matrix(rep(c(1, 1, 2, 2), 2), 2, 4, byrow=TRUE))
plot(plot.data[[3]][, c(1, 2)], type="l", main="(a)", ylab="AUC", 
     xlab="Number of selected features", ylim=ylim1, xlim=xlim, col=col[3], 
     lty=lty[1])
lines(plot.data[[6]][, c(1, 2)], col=col[4], lty=lty[1])
lines(plot.data[[9]][, c(1, 2)], col=col[5], lty=lty[1])
lines(plot.data[[12]][, c(1, 2)], col=col[6], lty=lty[1])
abline(h=plot.data[[1]][, 2], col=col[1], lty=lty[1])
abline(h=plot.data[[2]][, 2], col=col[2], lty=lty[1])
legend("bottomright", legend=labels, col=col, lty=lty[1],
       seg.len=1, merge=TRUE, bg="white")
plot(plot.data[[3]][, c(1, 3)], type="l", main="(b)", ylab="Brier skill score", 
     xlab="Number of selected features", ylim=ylim2, xlim=xlim, col=col[3], 
     lty=lty[1])
lines(plot.data[[6]][, c(1, 3)], col=col[4], lty=lty[1])
lines(plot.data[[9]][, c(1, 3)], col=col[5], lty=lty[1])
lines(plot.data[[12]][, c(1, 3)], col=col[6], lty=lty[1])
abline(h=plot.data[[1]][, 3], col=col[1], lty=lty[1])
abline(h=plot.data[[2]][, 3], col=col[2], lty=lty[1])
par(opar)

# ---- lines_simulations_res3_auc_briers ----
library(sp)
load(file="results/simulations_res3.Rdata")
plot.data <- c(sapply(c("grridge", "ridge"), function(m) {
  unname(cbind(NA, mean(sapply(res, function(s) {
    s$auc[substr(names(s$auc), 1, nchar(m))==m]}), na.rm=TRUE),
    mean(sapply(res, function(s) {
      s$auc[substr(names(s$briers), 1, nchar(m))==m]}), na.rm=TRUE)))}, 
  simplify=FALSE),
  sapply(c(paste0("gren", 1:3), paste0("enet", 1:3),
           paste0("cmcp", 1:3), paste0("gel", 1:3)), function(m) {
             cbind(rowMeans(sapply(res, function(s) {
               vec <- s$psel[substr(names(s$psel), 1, nchar(m))==m]
               replace(rep(NA, 8), 1:length(vec), vec)}), 
               na.rm=TRUE),
               rowMeans(sapply(res, function(s) {
                 vec <- s$auc[substr(names(s$auc), 1, nchar(m))==m]
                 replace(rep(NA, 8), 1:length(vec), vec)}), 
                 na.rm=TRUE),
               rowMeans(sapply(res, function(s) {
                 vec <- s$briers[substr(names(s$briers), 1, nchar(m))==m]
                 replace(rep(NA, 8), 1:length(vec), vec)}), 
                 na.rm=TRUE))}, simplify=FALSE))

methods <- c("GRridge", "ridge", "gren", "enet", "cmcp", "gel")
labels <- expression("GRridge", "ridge", "gren,"~alpha==0.5, 
                     "enet"~alpha==0.5, "cMCP,"~alpha==0.05, 
                     "gel,"~alpha==0.05)

col <- bpy.colors(length(methods), cutoff.tail=0.1)
lty <- c(1)

xlim <- range(sapply(plot.data[-c(1, 2, 3, 5, 6, 8, 10, 11, 13, 14)], 
                     function(s) {s[, 1]}), na.rm=TRUE)
ylim1 <- range(sapply(plot.data[-c(3, 5, 6, 8, 10, 11, 13, 14)], 
                      function(s) {s[, 2]}), na.rm=TRUE)
ylim2 <- range(sapply(plot.data[-c(3, 5, 6, 8, 10, 11, 13, 14, 9, 10, 11)], 
                      function(s) {s[, 3]}), na.rm=TRUE)

opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
layout(matrix(rep(c(1, 1, 2, 2), 2), 2, 4, byrow=TRUE))
plot(plot.data[[4]][, c(1, 2)], type="l", main="(a)", ylab="AUC", 
     xlab="Number of selected features", ylim=ylim1, xlim=xlim, col=col[3], 
     lty=lty[1])
lines(plot.data[[7]][, c(1, 2)], col=col[4], lty=lty[1])
lines(plot.data[[9]][, c(1, 2)], col=col[5], lty=lty[1])
lines(plot.data[[12]][, c(1, 2)], col=col[6], lty=lty[1])
abline(h=plot.data[[1]][, 2], col=col[1], lty=lty[1])
abline(h=plot.data[[2]][, 2], col=col[2], lty=lty[1])
legend("bottomright", legend=labels, col=col, lty=lty[1],
       seg.len=1, merge=TRUE, bg="white")
plot(plot.data[[4]][, c(1, 3)], type="l", main="(b)", ylab="Brier skill score", 
     xlab="Number of selected features", ylim=ylim2, xlim=xlim, col=col[3], 
     lty=lty[1])
lines(plot.data[[7]][, c(1, 3)], col=col[4], lty=lty[1])
lines(plot.data[[9]][, c(1, 3)], col=col[5], lty=lty[1])
lines(plot.data[[12]][, c(1, 3)], col=col[6], lty=lty[1])
abline(h=plot.data[[1]][, 3], col=col[1], lty=lty[1])
abline(h=plot.data[[2]][, 3], col=col[2], lty=lty[1])
par(opar)

# ---- lines_simulations_res4_auc_briers ----
library(sp)
load(file="results/simulations_res4.Rdata")
plot.data <- c(sapply(c("grridge", "ridge"), function(m) {
  unname(cbind(NA, mean(sapply(res, function(s) {
    s$auc[substr(names(s$auc), 1, nchar(m))==m]}), na.rm=TRUE),
    mean(sapply(res, function(s) {
      s$auc[substr(names(s$briers), 1, nchar(m))==m]}), na.rm=TRUE)))}, 
  simplify=FALSE),
  sapply(c(paste0("gren", 1:3), paste0("enet", 1:3),
           paste0("cmcp", 1:3), paste0("gel", 1:3)), function(m) {
             cbind(rowMeans(sapply(res, function(s) {
               vec <- s$psel[substr(names(s$psel), 1, nchar(m))==m]
               replace(rep(NA, 8), 1:length(vec), vec)}), 
               na.rm=TRUE),
               rowMeans(sapply(res, function(s) {
                 vec <- s$auc[substr(names(s$auc), 1, nchar(m))==m]
                 replace(rep(NA, 8), 1:length(vec), vec)}), 
                 na.rm=TRUE),
               rowMeans(sapply(res, function(s) {
                 vec <- s$briers[substr(names(s$briers), 1, nchar(m))==m]
                 replace(rep(NA, 8), 1:length(vec), vec)}), 
                 na.rm=TRUE))}, simplify=FALSE))

methods <- c("GRridge", "ridge", "gren", "enet", "cmcp", "gel")
labels <- expression("GRridge", "ridge", "gren,"~alpha==0.5, 
                     "enet"~alpha==0.5, "cMCP,"~alpha==0.05, 
                     "gel,"~alpha==0.05)

col <- bpy.colors(length(methods), cutoff.tail=0.1)
lty <- c(1)

xlim <- range(sapply(plot.data[-c(1, 2, 3, 5, 6, 8, 10, 11, 13, 14)], 
                     function(s) {s[, 1]}), na.rm=TRUE)
ylim1 <- range(sapply(plot.data[-c(3, 5, 6, 8, 10, 11, 13, 14)], 
                      function(s) {s[, 2]}), na.rm=TRUE)
ylim2 <- range(sapply(plot.data[-c(3, 5, 6, 8, 10, 11, 13, 14, 9, 10, 11)], 
                      function(s) {s[, 3]}), na.rm=TRUE)

opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
layout(matrix(rep(c(1, 1, 2, 2), 2), 2, 4, byrow=TRUE))
plot(plot.data[[4]][, c(1, 2)], type="l", main="(a)", ylab="AUC", 
     xlab="Number of selected features", ylim=ylim1, xlim=xlim, col=col[3], 
     lty=lty[1])
lines(plot.data[[7]][, c(1, 2)], col=col[4], lty=lty[1])
lines(plot.data[[9]][, c(1, 2)], col=col[5], lty=lty[1])
lines(plot.data[[12]][, c(1, 2)], col=col[6], lty=lty[1])
abline(h=plot.data[[1]][, 2], col=col[1], lty=lty[1])
abline(h=plot.data[[2]][, 2], col=col[2], lty=lty[1])
legend("bottomright", legend=labels, col=col, lty=lty[1],
       seg.len=1, merge=TRUE, bg="white")
plot(plot.data[[4]][, c(1, 3)], type="l", main="(b)", ylab="Brier skill score", 
     xlab="Number of selected features", ylim=ylim2, xlim=xlim, col=col[3], 
     lty=lty[1])
lines(plot.data[[7]][, c(1, 3)], col=col[4], lty=lty[1])
lines(plot.data[[9]][, c(1, 3)], col=col[5], lty=lty[1])
lines(plot.data[[12]][, c(1, 3)], col=col[6], lty=lty[1])
abline(h=plot.data[[1]][, 3], col=col[1], lty=lty[1])
abline(h=plot.data[[2]][, 3], col=col[2], lty=lty[1])
par(opar)

# ---- lines_simulations_res5_auc_briers ----
library(sp)
load(file="results/simulations_res5.Rdata")
plot.data <- c(sapply(c("grridge", "ridge"), function(m) {
  unname(cbind(NA, mean(sapply(res, function(s) {
    s$auc[substr(names(s$auc), 1, nchar(m))==m]}), na.rm=TRUE),
    mean(sapply(res, function(s) {
      s$auc[substr(names(s$briers), 1, nchar(m))==m]}), na.rm=TRUE)))}, 
  simplify=FALSE),
  sapply(c(paste0("gren", 1:3), paste0("enet", 1:3),
           paste0("cmcp", 1:3), paste0("gel", 1:3)), function(m) {
             cbind(rowMeans(sapply(res, function(s) {
               vec <- s$psel[substr(names(s$psel), 1, nchar(m))==m]
               replace(rep(NA, 8), 1:length(vec), vec)}), 
               na.rm=TRUE),
               rowMeans(sapply(res, function(s) {
                 vec <- s$auc[substr(names(s$auc), 1, nchar(m))==m]
                 replace(rep(NA, 8), 1:length(vec), vec)}), 
                 na.rm=TRUE),
               rowMeans(sapply(res, function(s) {
                 vec <- s$briers[substr(names(s$briers), 1, nchar(m))==m]
                 replace(rep(NA, 8), 1:length(vec), vec)}), 
                 na.rm=TRUE))}, simplify=FALSE))

methods <- c("GRridge", "ridge", "gren", "enet", "cmcp", "gel")
labels <- expression("GRridge", "ridge", "gren,"~alpha==0.05, 
                     "enet"~alpha==0.05, "cMCP,"~alpha==0.05, 
                     "gel,"~alpha==0.05)

col <- bpy.colors(length(methods), cutoff.tail=0.1)
lty <- c(1)

xlim <- range(sapply(plot.data[-c(1, 2, 4, 5, 7, 8, 10, 11, 13, 14)], 
                     function(s) {s[, 1]}), na.rm=TRUE)
ylim1 <- range(sapply(plot.data[-c(4, 5, 7, 8, 10, 11, 13, 14)], 
                      function(s) {s[, 2]}), na.rm=TRUE)
ylim2 <- range(sapply(plot.data[-c(4, 5, 7, 8, 10, 11, 13, 14, 9, 10, 11)], 
                      function(s) {s[, 3]}), na.rm=TRUE)

opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
layout(matrix(rep(c(1, 1, 2, 2), 2), 2, 4, byrow=TRUE))
plot(plot.data[[3]][, c(1, 2)], type="l", main="(a)", ylab="AUC", 
     xlab="Number of selected features", ylim=ylim1, xlim=xlim, col=col[3], 
     lty=lty[1])
lines(plot.data[[6]][, c(1, 2)], col=col[4], lty=lty[1])
lines(plot.data[[9]][, c(1, 2)], col=col[5], lty=lty[1])
lines(plot.data[[12]][, c(1, 2)], col=col[6], lty=lty[1])
abline(h=plot.data[[1]][, 2], col=col[1], lty=lty[1])
abline(h=plot.data[[2]][, 2], col=col[2], lty=lty[1])
legend("bottomright", legend=labels, col=col, lty=lty[1],
       seg.len=1, merge=TRUE, bg="white")
plot(plot.data[[3]][, c(1, 3)], type="l", main="(b)", ylab="Brier skill score", 
     xlab="Number of selected features", ylim=ylim2, xlim=xlim, col=col[3], 
     lty=lty[1])
lines(plot.data[[6]][, c(1, 3)], col=col[4], lty=lty[1])
lines(plot.data[[9]][, c(1, 3)], col=col[5], lty=lty[1])
lines(plot.data[[12]][, c(1, 3)], col=col[6], lty=lty[1])
abline(h=plot.data[[1]][, 3], col=col[1], lty=lty[1])
abline(h=plot.data[[2]][, 3], col=col[2], lty=lty[1])
par(opar)

# ---- barplot_lines_micrornaseq_colorectal_cancer_res1_auc ----
library(sp)
load("results/micrornaseq_colorectal_cancer_fit1.Rdata")
res <- read.table("results/micrornaseq_colorectal_cancer_res1.csv", 
                  stringsAsFactors=FALSE)
pred <- as.matrix(res[substr(rownames(res), 1, 4)=="pred", ])
psel <- as.matrix(res[substr(rownames(res), 1, 4)=="psel", ])
auc <- as.matrix(res[substr(rownames(res), 1, 3)=="auc", ])

methods1 <- c("GRridge", expression("gren, "~alpha==0.05, "gren, "~alpha==0.5,
                                    "gren, "~alpha==0.95), "no co-data")
methods2 <- c("grridge2", "ridge", "rf", "gren5", "enet1", "cmcp1")
col1 <- grey.colors(length(methods1) - 1, start=0.3, end=0.9, gamma=2.2, 
                    alpha=NULL)
col2 <- bpy.colors(length(methods2) + 1, cutoff.tail=0.1)[
  -(length(methods2) + 1)]
lty1 <- 2
lty2 <- c(1, 2)
labels1 <- list(diff.expr=expression(atop(NA, atop(textstyle("FDR"), 
                                                   textstyle(""<=0.0001))),
                                     atop(NA, atop(textstyle(0.0001<phantom(0)),
                                                   textstyle("FDR"<=0.0186))),
                                     atop(NA, atop(textstyle("FDR"), 
                                                   textstyle("">0.0186))), 
                                     "no FDR"))
labels2 <- expression("GRridge", "ridge", "random forest",
                      "gren,"~alpha==0.5, "enet"~alpha==0.05, 
                      "cMCP"~alpha==0.05)

plot.data1 <- lapply(1:length(fit.gren4$lambdag), function(s) {
  rbind(fit.grridge2$lambdamults[[s]], fit.gren4$lambdag[[s]],
        fit.gren5$lambdag[[s]], fit.gren6$lambdag[[s]])})
plot.data2 <- lapply(methods2, function(m) {
  aggregate(auc[, substr(colnames(auc), 1, nchar(m))==m], 
            list(psel=psel[, substr(colnames(psel), 1, nchar(m))==m]), mean)})
names(plot.data2) <- methods2

opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
layout(matrix(rep(c(1, 1, 2, 2), 2), 2, 4, byrow=TRUE))
barplot(plot.data1[[1]], beside=TRUE, col=col1,
        legend.text=methods1,
        args.legend=list(x="bottomright", fill=c(col1, NA),
                         border=c(rep(1, length(methods1) - 1), NA),
                         lty=c(rep(NA, length(methods1) - 1), lty1[1]),
                         seg.len=1, merge=TRUE, bg="white"),
        names.arg=labels1[[1]], ylab=expression(hat(lambda)~"'"[g]),
        main="(a)")
abline(h=1, lty=2)

plot(plot.data2[[4]], type="l", xlab="Number of selected features", ylab="AUC", 
     main="(b)", ylim=range(sapply(plot.data2, "[", 2)), xlim=c(0, 600), 
     col=col2[4], lty=lty2[1])
lines(plot.data2[[5]], col=col2[5], lty=lty2[1])
lines(plot.data2[[6]], col=col2[6], lty=lty2[1])
abline(h=plot.data2[[1]][, 2], col=col2[1], lty=lty2[2])
abline(h=plot.data2[[2]][, 2], col=col2[2], lty=lty2[2])
abline(h=plot.data2[[3]][, 2], col=col2[3], lty=lty2[2])
legend("bottomright", legend=labels2, col=col2, 
       lty=c(rep(lty2[2], 3), rep(lty2[1], 3)), bg = "white")
par(opar)

# ---- barplot_lines_rnaseq_oral_cancer_metastasis_res1_auc ----
library(sp)
load("results/rnaseq_oral_cancer_metastasis_fit1.Rdata")
res <- read.table("results/rnaseq_oral_cancer_metastasis_res1.csv", 
                  stringsAsFactors=FALSE)
pred <- as.matrix(res[substr(rownames(res), 1, 4)=="pred", ])
psel <- as.matrix(res[substr(rownames(res), 1, 4)=="psel", ])
auc <- as.matrix(res[substr(rownames(res), 1, 3)=="auc", ])

methods1 <- c("grridge", expression("gren, "~alpha==0.05, "gren, "~alpha==0.5,
                                    "gren, "~alpha==0.95), "no co-data")
methods2 <- c("ridge", "rf", "grridge2", "gren4", "enet1", "ocmcp1")
col1 <- grey.colors(length(methods1) - 1, start=0.3, end=0.9, gamma=2.2, 
                    alpha=NULL)
col2 <- bpy.colors(length(methods2) + 1, cutoff.tail=0.1)[
  -(length(methods2) + 1)]
lty1 <- 2
lty2 <- c(1, 2)
labels1 <- list(corr=expression(tau > 0.460, 
                                atop(NA, atop(textstyle(0.316 <phantom(0)), 
                                              textstyle(tau <= 0.460))),
                                atop(NA, atop(textstyle(-0.010 <phantom(0)), 
                                              textstyle(tau <= 0.316))),
                                atop(NA, atop(textstyle(-0.087 <phantom(0)), 
                                              textstyle(tau <= -0.010))), 
                                tau <= -0.087),
                pv=expression(p <= 0.001, 
                              atop(NA, atop(textstyle(0.001 <phantom(0)), 
                                            textstyle(p <= 0.011))),
                              atop(NA, atop(textstyle(0.011 <phantom(0)), 
                                            textstyle(p <= 0.028))),
                              atop(NA, atop(textstyle(0.028 <phantom(0)), 
                                            textstyle(p <= 0.058))), 
                              p > 0.058))
labels2 <- expression("ridge", "GRridge", "random forest", "gren,"~alpha==0.05, 
                      "enet"~alpha==0.05, "cMCP,"~alpha==0.05)

plot.data1 <- lapply(1:length(fit.gren1$lambdag), function(s) {
  rbind(fit.grridge$lambdamults[[s]], fit.gren1$lambdag[[s]],
        fit.gren2$lambdag[[s]], fit.gren3$lambdag[[s]])})
plot.data2 <- lapply(methods2, function(m) {
  aggregate(auc[, substr(colnames(auc), 1, nchar(m))==m], 
            list(psel=psel[, substr(colnames(psel), 1, nchar(m))==m]), mean)})
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
                         seg.len=1, merge=TRUE, bg="white"), las=2,
        names.arg=labels1[[1]], ylab=expression(hat(lambda)~"'"[g]),
        main="(a)")
abline(h=1, lty=2)

barplot(plot.data1[[2]], beside=TRUE, col=col1, main="(b)", 
        names.arg=labels1[[2]], ylab=expression(hat(lambda)~"'"[g]), las=2)
abline(h=1, lty=2)

plot(plot.data2[[4]], type="l", xlab="Number of selected features", ylab="AUC", 
     main="(c)", ylim=c(0.65, max(auc)), xlim=c(0, 500), col=col2[4], 
     lty=lty2[1])
lines(plot.data2[[5]], col=col2[5], lty=lty2[1])
lines(plot.data2[[6]], col=col2[6], lty=lty2[1])
abline(h=plot.data2[[1]][, 2], col=col2[1], lty=lty2[2])
abline(h=plot.data2[[2]][, 2], col=col2[2], lty=lty2[2])
abline(h=plot.data2[[3]][, 2], col=col2[3], lty=lty2[2])
legend("bottomright", legend=labels2, col=col2, 
       lty=c(rep(lty2[2], 3), rep(lty2[1], 3)), bg = "white")
par(opar)

################################## supplement ##################################
# ---- lines2_simulations_res1_auc_briers ----
library(sp)
load(file="results/simulations_res1.Rdata")
plot.data <- c(sapply(c("grridge", "ridge"), function(m) {
  unname(cbind(NA, mean(sapply(res, function(s) {
    s$auc[substr(names(s$auc), 1, nchar(m))==m]}), na.rm=TRUE),
    mean(sapply(res, function(s) {
      s$auc[substr(names(s$briers), 1, nchar(m))==m]}), na.rm=TRUE)))}, 
  simplify=FALSE),
  sapply(c(paste0("gren", 1:3), paste0("enet", 1:3),
           paste0("cmcp", 1:3), paste0("gel", 1:3)), function(m) {
             cbind(rowMeans(sapply(res, function(s) {
               vec <- s$psel[substr(names(s$psel), 1, nchar(m))==m]
               replace(rep(NA, 8), 1:length(vec), vec)}), 
               na.rm=TRUE),
               rowMeans(sapply(res, function(s) {
                 vec <- s$auc[substr(names(s$auc), 1, nchar(m))==m]
                 replace(rep(NA, 8), 1:length(vec), vec)}), 
                 na.rm=TRUE),
               rowMeans(sapply(res, function(s) {
                 vec <- s$briers[substr(names(s$briers), 1, nchar(m))==m]
                 replace(rep(NA, 8), 1:length(vec), vec)}), 
                 na.rm=TRUE))}, simplify=FALSE))

methods <- c("GRridge", "ridge", "gren", "enet", "cmcp", "gel")
labels <- c(methods, expression(alpha==0.05), expression(alpha==0.5), 
            expression(alpha==0.95))

col <- bpy.colors(length(methods), cutoff.tail=0.1)
lty <- c(1:3)

xlim <- range(sapply(plot.data[-c(1, 2)], function(s) {s[, 1]}), na.rm=TRUE)
ylim1 <- range(sapply(plot.data, function(s) {s[, 2]}), na.rm=TRUE)
ylim2 <- range(sapply(plot.data[-c(9, 10, 11)], function(s) {s[, 3]}), 
               na.rm=TRUE)

opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
layout(matrix(rep(c(1, 1, 2, 2), 2), 2, 4, byrow=TRUE))
plot(plot.data[[3]][, c(1, 2)], type="l", main="(a)", ylab="AUC", 
     xlab="Number of selected features", ylim=ylim1, xlim=xlim, col=col[3], 
     lty=lty[1])
lines(plot.data[[4]][, c(1, 2)], col=col[3], lty=lty[2])
lines(plot.data[[5]][, c(1, 2)], col=col[3], lty=lty[3])
lines(plot.data[[6]][, c(1, 2)], col=col[4], lty=lty[1])
lines(plot.data[[7]][, c(1, 2)], col=col[4], lty=lty[2])
lines(plot.data[[8]][, c(1, 2)], col=col[4], lty=lty[3])
lines(plot.data[[9]][, c(1, 2)], col=col[5], lty=lty[1])
lines(plot.data[[10]][, c(1, 2)], col=col[5], lty=lty[2])
lines(plot.data[[11]][, c(1, 2)], col=col[5], lty=lty[3])
lines(plot.data[[12]][, c(1, 2)], col=col[6], lty=lty[1])
lines(plot.data[[13]][, c(1, 2)], col=col[6], lty=lty[2])
lines(plot.data[[14]][, c(1, 2)], col=col[6], lty=lty[3])
abline(h=plot.data[[1]][, 2], col=col[1], lty=lty[1])
abline(h=plot.data[[2]][, 2], col=col[2], lty=lty[1])
legend("bottomright", legend=labels, fill=c(col, rep(NA, 3)), 
       border=c(rep(1, length(methods)), rep(NA, 3)),
       lty=c(rep(NA, length(methods)), lty),
       seg.len=1, merge=TRUE, bg="white")
plot(plot.data[[3]][, c(1, 3)], type="l", main="(b)", ylab="Brier skill score", 
     xlab="Number of selected features", ylim=ylim2, xlim=xlim, col=col[3], 
     lty=lty[1])
lines(plot.data[[4]][, c(1, 3)], col=col[3], lty=lty[2])
lines(plot.data[[5]][, c(1, 3)], col=col[3], lty=lty[3])
lines(plot.data[[6]][, c(1, 3)], col=col[4], lty=lty[1])
lines(plot.data[[7]][, c(1, 3)], col=col[4], lty=lty[2])
lines(plot.data[[8]][, c(1, 3)], col=col[4], lty=lty[3])
lines(plot.data[[9]][, c(1, 3)], col=col[5], lty=lty[1])
lines(plot.data[[10]][, c(1, 3)], col=col[5], lty=lty[2])
lines(plot.data[[11]][, c(1, 3)], col=col[5], lty=lty[3])
lines(plot.data[[12]][, c(1, 3)], col=col[6], lty=lty[1])
lines(plot.data[[13]][, c(1, 3)], col=col[6], lty=lty[2])
lines(plot.data[[14]][, c(1, 3)], col=col[6], lty=lty[3])
abline(h=plot.data[[1]][, 3], col=col[1], lty=lty[1])
abline(h=plot.data[[2]][, 3], col=col[2], lty=lty[1])
par(opar)

# ---- lines2_simulations_res2_auc_briers ----
library(sp)
load(file="results/simulations_res2.Rdata")
mults <- sapply(c("grridge", paste("gren", 1:3, sep="")), function(m) {
  sapply(res, function(r) {r$mults[grep(m, names(r$mults))]})}, simplify=FALSE)
plot.data <- c(sapply(c("grridge", "ridge"), function(m) {
  unname(cbind(NA, mean(sapply(res, function(s) {
    s$auc[substr(names(s$auc), 1, nchar(m))==m]}), na.rm=TRUE),
    mean(sapply(res, function(s) {
      s$auc[substr(names(s$briers), 1, nchar(m))==m]}), na.rm=TRUE)))}, 
  simplify=FALSE),
  sapply(c(paste0("gren", 1:3), paste0("enet", 1:3),
           paste0("cmcp", 1:3), paste0("gel", 1:3)), function(m) {
             cbind(rowMeans(sapply(res, function(s) {
               vec <- s$psel[substr(names(s$psel), 1, nchar(m))==m]
               replace(rep(NA, 8), 1:length(vec), vec)}), 
               na.rm=TRUE),
               rowMeans(sapply(res, function(s) {
                 vec <- s$auc[substr(names(s$auc), 1, nchar(m))==m]
                 replace(rep(NA, 8), 1:length(vec), vec)}), 
                 na.rm=TRUE),
               rowMeans(sapply(res, function(s) {
                 vec <- s$briers[substr(names(s$briers), 1, nchar(m))==m]
                 replace(rep(NA, 8), 1:length(vec), vec)}), 
                 na.rm=TRUE))}, simplify=FALSE))

methods <- c("GRridge", "ridge", "gren", "enet", "cmcp", "gel")
labels <- c(methods, expression(alpha==0.05), expression(alpha==0.5), 
            expression(alpha==0.95))

col <- bpy.colors(length(methods), cutoff.tail=0.1)
lty <- c(1:3)

xlim <- range(sapply(plot.data[-c(1, 2)], function(s) {s[, 1]}), na.rm=TRUE)
ylim1 <- range(sapply(plot.data, function(s) {s[, 2]}), na.rm=TRUE)
ylim2 <- range(sapply(plot.data[-c(9, 10, 11)], function(s) {s[, 3]}), 
               na.rm=TRUE)

opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
layout(matrix(rep(c(1, 1, 2, 2), 2), 2, 4, byrow=TRUE))
plot(plot.data[[3]][, c(1, 2)], type="l", main="(a)", ylab="AUC", 
     xlab="Number of selected features", ylim=ylim1, xlim=xlim, col=col[3], 
     lty=lty[1])
lines(plot.data[[4]][, c(1, 2)], col=col[3], lty=lty[2])
lines(plot.data[[5]][, c(1, 2)], col=col[3], lty=lty[3])
lines(plot.data[[6]][, c(1, 2)], col=col[4], lty=lty[1])
lines(plot.data[[7]][, c(1, 2)], col=col[4], lty=lty[2])
lines(plot.data[[8]][, c(1, 2)], col=col[4], lty=lty[3])
lines(plot.data[[9]][, c(1, 2)], col=col[5], lty=lty[1])
lines(plot.data[[10]][, c(1, 2)], col=col[5], lty=lty[2])
lines(plot.data[[11]][, c(1, 2)], col=col[5], lty=lty[3])
lines(plot.data[[12]][, c(1, 2)], col=col[6], lty=lty[1])
lines(plot.data[[13]][, c(1, 2)], col=col[6], lty=lty[2])
lines(plot.data[[14]][, c(1, 2)], col=col[6], lty=lty[3])
abline(h=plot.data[[1]][, 2], col=col[1], lty=lty[1])
abline(h=plot.data[[2]][, 2], col=col[2], lty=lty[1])
legend("bottomright", legend=labels, fill=c(col, rep(NA, 3)), 
       border=c(rep(1, length(methods)), rep(NA, 3)),
       lty=c(rep(NA, length(methods)), lty),
       seg.len=1, merge=TRUE, bg="white")
plot(plot.data[[3]][, c(1, 3)], type="l", main="(b)", ylab="Brier skill score", 
     xlab="Number of selected features", ylim=ylim2, xlim=xlim, col=col[3], 
     lty=lty[1])
lines(plot.data[[4]][, c(1, 3)], col=col[3], lty=lty[2])
lines(plot.data[[5]][, c(1, 3)], col=col[3], lty=lty[3])
lines(plot.data[[6]][, c(1, 3)], col=col[4], lty=lty[1])
lines(plot.data[[7]][, c(1, 3)], col=col[4], lty=lty[2])
lines(plot.data[[8]][, c(1, 3)], col=col[4], lty=lty[3])
lines(plot.data[[9]][, c(1, 3)], col=col[5], lty=lty[1])
lines(plot.data[[10]][, c(1, 3)], col=col[5], lty=lty[2])
lines(plot.data[[11]][, c(1, 3)], col=col[5], lty=lty[3])
lines(plot.data[[12]][, c(1, 3)], col=col[6], lty=lty[1])
lines(plot.data[[13]][, c(1, 3)], col=col[6], lty=lty[2])
lines(plot.data[[14]][, c(1, 3)], col=col[6], lty=lty[3])
abline(h=plot.data[[1]][, 3], col=col[1], lty=lty[1])
abline(h=plot.data[[2]][, 3], col=col[2], lty=lty[1])
par(opar)

# ---- lines2_simulations_res3_auc_briers ----
library(sp)
load(file="results/simulations_res3.Rdata")
plot.data <- c(sapply(c("grridge", "ridge"), function(m) {
  unname(cbind(NA, mean(sapply(res, function(s) {
    s$auc[substr(names(s$auc), 1, nchar(m))==m]}), na.rm=TRUE),
    mean(sapply(res, function(s) {
      s$auc[substr(names(s$briers), 1, nchar(m))==m]}), na.rm=TRUE)))}, 
  simplify=FALSE),
  sapply(c(paste0("gren", 1:3), paste0("enet", 1:3),
           paste0("cmcp", 1:3), paste0("gel", 1:3)), function(m) {
             cbind(rowMeans(sapply(res, function(s) {
               vec <- s$psel[substr(names(s$psel), 1, nchar(m))==m]
               replace(rep(NA, 8), 1:length(vec), vec)}), 
               na.rm=TRUE),
               rowMeans(sapply(res, function(s) {
                 vec <- s$auc[substr(names(s$auc), 1, nchar(m))==m]
                 replace(rep(NA, 8), 1:length(vec), vec)}), 
                 na.rm=TRUE),
               rowMeans(sapply(res, function(s) {
                 vec <- s$briers[substr(names(s$briers), 1, nchar(m))==m]
                 replace(rep(NA, 8), 1:length(vec), vec)}), 
                 na.rm=TRUE))}, simplify=FALSE))

methods <- c("GRridge", "ridge", "gren", "enet", "cmcp", "gel")
labels <- c(methods, expression(alpha==0.05), expression(alpha==0.5), 
            expression(alpha==0.95))

col <- bpy.colors(length(methods), cutoff.tail=0.1)
lty <- c(1:3)

xlim <- range(sapply(plot.data[-c(1, 2)], function(s) {s[, 1]}), na.rm=TRUE)
ylim1 <- range(sapply(plot.data, function(s) {s[, 2]}), na.rm=TRUE)
ylim2 <- range(sapply(plot.data[-c(9, 10, 11)], function(s) {s[, 3]}), 
               na.rm=TRUE)

opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
layout(matrix(rep(c(1, 1, 2, 2), 2), 2, 4, byrow=TRUE))
plot(plot.data[[3]][, c(1, 2)], type="l", main="(a)", ylab="AUC", 
     xlab="Number of selected features", ylim=ylim1, xlim=xlim, col=col[3], 
     lty=lty[1])
lines(plot.data[[4]][, c(1, 2)], col=col[3], lty=lty[2])
lines(plot.data[[5]][, c(1, 2)], col=col[3], lty=lty[3])
lines(plot.data[[6]][, c(1, 2)], col=col[4], lty=lty[1])
lines(plot.data[[7]][, c(1, 2)], col=col[4], lty=lty[2])
lines(plot.data[[8]][, c(1, 2)], col=col[4], lty=lty[3])
lines(plot.data[[9]][, c(1, 2)], col=col[5], lty=lty[1])
lines(plot.data[[10]][, c(1, 2)], col=col[5], lty=lty[2])
lines(plot.data[[11]][, c(1, 2)], col=col[5], lty=lty[3])
lines(plot.data[[12]][, c(1, 2)], col=col[6], lty=lty[1])
lines(plot.data[[13]][, c(1, 2)], col=col[6], lty=lty[2])
lines(plot.data[[14]][, c(1, 2)], col=col[6], lty=lty[3])
abline(h=plot.data[[1]][, 2], col=col[1], lty=lty[1])
abline(h=plot.data[[2]][, 2], col=col[2], lty=lty[1])
legend("bottomright", legend=labels, fill=c(col, rep(NA, 3)), 
       border=c(rep(1, length(methods)), rep(NA, 3)),
       lty=c(rep(NA, length(methods)), lty),
       seg.len=1, merge=TRUE, bg="white")
plot(plot.data[[3]][, c(1, 3)], type="l", main="(b)", ylab="Brier skill score", 
     xlab="Number of selected features", ylim=ylim2, xlim=xlim, col=col[3], 
     lty=lty[1])
lines(plot.data[[4]][, c(1, 3)], col=col[3], lty=lty[2])
lines(plot.data[[5]][, c(1, 3)], col=col[3], lty=lty[3])
lines(plot.data[[6]][, c(1, 3)], col=col[4], lty=lty[1])
lines(plot.data[[7]][, c(1, 3)], col=col[4], lty=lty[2])
lines(plot.data[[8]][, c(1, 3)], col=col[4], lty=lty[3])
lines(plot.data[[9]][, c(1, 3)], col=col[5], lty=lty[1])
lines(plot.data[[10]][, c(1, 3)], col=col[5], lty=lty[2])
lines(plot.data[[11]][, c(1, 3)], col=col[5], lty=lty[3])
lines(plot.data[[12]][, c(1, 3)], col=col[6], lty=lty[1])
lines(plot.data[[13]][, c(1, 3)], col=col[6], lty=lty[2])
lines(plot.data[[14]][, c(1, 3)], col=col[6], lty=lty[3])
abline(h=plot.data[[1]][, 3], col=col[1], lty=lty[1])
abline(h=plot.data[[2]][, 3], col=col[2], lty=lty[1])
par(opar)

# ---- lines2_simulations_res4_auc_briers ----
library(sp)
load(file="results/simulations_res4.Rdata")
plot.data <- c(sapply(c("grridge", "ridge"), function(m) {
  unname(cbind(NA, mean(sapply(res, function(s) {
    s$auc[substr(names(s$auc), 1, nchar(m))==m]}), na.rm=TRUE),
    mean(sapply(res, function(s) {
      s$auc[substr(names(s$briers), 1, nchar(m))==m]}), na.rm=TRUE)))}, 
  simplify=FALSE),
  sapply(c(paste0("gren", 1:3), paste0("enet", 1:3),
           paste0("cmcp", 1:3), paste0("gel", 1:3)), function(m) {
             cbind(rowMeans(sapply(res, function(s) {
               vec <- s$psel[substr(names(s$psel), 1, nchar(m))==m]
               replace(rep(NA, 8), 1:length(vec), vec)}), 
               na.rm=TRUE),
               rowMeans(sapply(res, function(s) {
                 vec <- s$auc[substr(names(s$auc), 1, nchar(m))==m]
                 replace(rep(NA, 8), 1:length(vec), vec)}), 
                 na.rm=TRUE),
               rowMeans(sapply(res, function(s) {
                 vec <- s$briers[substr(names(s$briers), 1, nchar(m))==m]
                 replace(rep(NA, 8), 1:length(vec), vec)}), 
                 na.rm=TRUE))}, simplify=FALSE))

methods <- c("GRridge", "ridge", "gren", "enet", "cmcp", "gel")
labels <- c(methods, expression(alpha==0.05), expression(alpha==0.5), 
            expression(alpha==0.95))

col <- bpy.colors(length(methods), cutoff.tail=0.1)
lty <- c(1:3)

xlim <- range(sapply(plot.data[-c(1, 2)], function(s) {s[, 1]}), na.rm=TRUE)
ylim1 <- range(sapply(plot.data, function(s) {s[, 2]}), na.rm=TRUE)
ylim2 <- range(sapply(plot.data[-c(9, 10, 11)], function(s) {s[, 3]}), 
               na.rm=TRUE)

opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
layout(matrix(rep(c(1, 1, 2, 2), 2), 2, 4, byrow=TRUE))
plot(plot.data[[3]][, c(1, 2)], type="l", main="(a)", ylab="AUC", 
     xlab="Number of selected features", ylim=ylim1, xlim=xlim, col=col[3], 
     lty=lty[1])
lines(plot.data[[4]][, c(1, 2)], col=col[3], lty=lty[2])
lines(plot.data[[5]][, c(1, 2)], col=col[3], lty=lty[3])
lines(plot.data[[6]][, c(1, 2)], col=col[4], lty=lty[1])
lines(plot.data[[7]][, c(1, 2)], col=col[4], lty=lty[2])
lines(plot.data[[8]][, c(1, 2)], col=col[4], lty=lty[3])
lines(plot.data[[9]][, c(1, 2)], col=col[5], lty=lty[1])
lines(plot.data[[10]][, c(1, 2)], col=col[5], lty=lty[2])
lines(plot.data[[11]][, c(1, 2)], col=col[5], lty=lty[3])
lines(plot.data[[12]][, c(1, 2)], col=col[6], lty=lty[1])
lines(plot.data[[13]][, c(1, 2)], col=col[6], lty=lty[2])
lines(plot.data[[14]][, c(1, 2)], col=col[6], lty=lty[3])
abline(h=plot.data[[1]][, 2], col=col[1], lty=lty[1])
abline(h=plot.data[[2]][, 2], col=col[2], lty=lty[1])
legend("bottomright", legend=labels, fill=c(col, rep(NA, 3)), 
       border=c(rep(1, length(methods)), rep(NA, 3)),
       lty=c(rep(NA, length(methods)), lty),
       seg.len=1, merge=TRUE, bg="white")
plot(plot.data[[3]][, c(1, 3)], type="l", main="(b)", ylab="Brier skill score", 
     xlab="Number of selected features", ylim=ylim2, xlim=xlim, col=col[3], 
     lty=lty[1])
lines(plot.data[[4]][, c(1, 3)], col=col[3], lty=lty[2])
lines(plot.data[[5]][, c(1, 3)], col=col[3], lty=lty[3])
lines(plot.data[[6]][, c(1, 3)], col=col[4], lty=lty[1])
lines(plot.data[[7]][, c(1, 3)], col=col[4], lty=lty[2])
lines(plot.data[[8]][, c(1, 3)], col=col[4], lty=lty[3])
lines(plot.data[[9]][, c(1, 3)], col=col[5], lty=lty[1])
lines(plot.data[[10]][, c(1, 3)], col=col[5], lty=lty[2])
lines(plot.data[[11]][, c(1, 3)], col=col[5], lty=lty[3])
lines(plot.data[[12]][, c(1, 3)], col=col[6], lty=lty[1])
lines(plot.data[[13]][, c(1, 3)], col=col[6], lty=lty[2])
lines(plot.data[[14]][, c(1, 3)], col=col[6], lty=lty[3])
abline(h=plot.data[[1]][, 3], col=col[1], lty=lty[1])
abline(h=plot.data[[2]][, 3], col=col[2], lty=lty[1])
par(opar)

# ---- lines2_simulations_res5_auc_briers ----
library(sp)
load(file="results/simulations_res5.Rdata")
plot.data <- c(sapply(c("grridge", "ridge"), function(m) {
  unname(cbind(NA, mean(sapply(res, function(s) {
    s$auc[substr(names(s$auc), 1, nchar(m))==m]}), na.rm=TRUE),
    mean(sapply(res, function(s) {
      s$auc[substr(names(s$briers), 1, nchar(m))==m]}), na.rm=TRUE)))}, 
  simplify=FALSE),
  sapply(c(paste0("gren", 1:3), paste0("enet", 1:3),
           paste0("cmcp", 1:3), paste0("gel", 1:3)), function(m) {
             cbind(rowMeans(sapply(res, function(s) {
               vec <- s$psel[substr(names(s$psel), 1, nchar(m))==m]
               replace(rep(NA, 8), 1:length(vec), vec)}), 
               na.rm=TRUE),
               rowMeans(sapply(res, function(s) {
                 vec <- s$auc[substr(names(s$auc), 1, nchar(m))==m]
                 replace(rep(NA, 8), 1:length(vec), vec)}), 
                 na.rm=TRUE),
               rowMeans(sapply(res, function(s) {
                 vec <- s$briers[substr(names(s$briers), 1, nchar(m))==m]
                 replace(rep(NA, 8), 1:length(vec), vec)}), 
                 na.rm=TRUE))}, simplify=FALSE))

methods <- c("GRridge", "ridge", "gren", "enet", "cmcp", "gel")
labels <- c(methods, expression(alpha==0.05), expression(alpha==0.5), 
            expression(alpha==0.95))

col <- bpy.colors(length(methods), cutoff.tail=0.1)
lty <- c(1:3)

xlim <- range(sapply(plot.data[-c(1, 2)], function(s) {s[, 1]}), na.rm=TRUE)
ylim1 <- range(sapply(plot.data, function(s) {s[, 2]}), na.rm=TRUE)
ylim2 <- range(sapply(plot.data[-c(9, 10, 11)], function(s) {s[, 3]}), 
               na.rm=TRUE)

opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
layout(matrix(rep(c(1, 1, 2, 2), 2), 2, 4, byrow=TRUE))
plot(plot.data[[3]][, c(1, 2)], type="l", main="(a)", ylab="AUC", 
     xlab="Number of selected features", ylim=ylim1, xlim=xlim, col=col[3], 
     lty=lty[1])
lines(plot.data[[4]][, c(1, 2)], col=col[3], lty=lty[2])
lines(plot.data[[5]][, c(1, 2)], col=col[3], lty=lty[3])
lines(plot.data[[6]][, c(1, 2)], col=col[4], lty=lty[1])
lines(plot.data[[7]][, c(1, 2)], col=col[4], lty=lty[2])
lines(plot.data[[8]][, c(1, 2)], col=col[4], lty=lty[3])
lines(plot.data[[9]][, c(1, 2)], col=col[5], lty=lty[1])
lines(plot.data[[10]][, c(1, 2)], col=col[5], lty=lty[2])
lines(plot.data[[11]][, c(1, 2)], col=col[5], lty=lty[3])
lines(plot.data[[12]][, c(1, 2)], col=col[6], lty=lty[1])
lines(plot.data[[13]][, c(1, 2)], col=col[6], lty=lty[2])
lines(plot.data[[14]][, c(1, 2)], col=col[6], lty=lty[3])
abline(h=plot.data[[1]][, 2], col=col[1], lty=lty[1])
abline(h=plot.data[[2]][, 2], col=col[2], lty=lty[1])
legend("bottomright", legend=labels, fill=c(col, rep(NA, 3)), 
       border=c(rep(1, length(methods)), rep(NA, 3)),
       lty=c(rep(NA, length(methods)), lty),
       seg.len=1, merge=TRUE, bg="white")
plot(plot.data[[3]][, c(1, 3)], type="l", main="(b)", ylab="Brier skill score", 
     xlab="Number of selected features", ylim=ylim2, xlim=xlim, col=col[3], 
     lty=lty[1])
lines(plot.data[[4]][, c(1, 3)], col=col[3], lty=lty[2])
lines(plot.data[[5]][, c(1, 3)], col=col[3], lty=lty[3])
lines(plot.data[[6]][, c(1, 3)], col=col[4], lty=lty[1])
lines(plot.data[[7]][, c(1, 3)], col=col[4], lty=lty[2])
lines(plot.data[[8]][, c(1, 3)], col=col[4], lty=lty[3])
lines(plot.data[[9]][, c(1, 3)], col=col[5], lty=lty[1])
lines(plot.data[[10]][, c(1, 3)], col=col[5], lty=lty[2])
lines(plot.data[[11]][, c(1, 3)], col=col[5], lty=lty[3])
lines(plot.data[[12]][, c(1, 3)], col=col[6], lty=lty[1])
lines(plot.data[[13]][, c(1, 3)], col=col[6], lty=lty[2])
lines(plot.data[[14]][, c(1, 3)], col=col[6], lty=lty[3])
abline(h=plot.data[[1]][, 3], col=col[1], lty=lty[1])
abline(h=plot.data[[2]][, 3], col=col[2], lty=lty[1])
par(opar)

# ---- lines2_simulations_res6_auc_briers ----
library(sp)
load(file="results/simulations_res6.Rdata")
plot.data <- 
  sapply(c(paste0("gren", 1:3), paste0("enet", 1:3)), function(m) {
             cbind(rowMeans(sapply(res, function(s) {
               vec <- s$psel[substr(names(s$psel), 1, nchar(m))==m]
               replace(rep(NA, 8), 1:length(vec), vec)}), 
               na.rm=TRUE),
               rowMeans(sapply(res, function(s) {
                 vec <- s$auc[substr(names(s$auc), 1, nchar(m))==m]
                 replace(rep(NA, 8), 1:length(vec), vec)}), 
                 na.rm=TRUE),
               rowMeans(sapply(res, function(s) {
                 vec <- s$briers[substr(names(s$briers), 1, nchar(m))==m]
                 replace(rep(NA, 8), 1:length(vec), vec)}), 
                 na.rm=TRUE))}, simplify=FALSE)

methods <- c("gren", "enet")
labels <- c(methods, expression(alpha==0.05), expression(alpha==0.5), 
            expression(alpha==0.95))

col <- bpy.colors(length(methods), cutoff.tail=0.2)
lty <- c(1:3)

xlim <- range(sapply(plot.data, function(s) {s[, 1]}), na.rm=TRUE)
ylim1 <- range(sapply(plot.data, function(s) {s[, 2]}), na.rm=TRUE)
ylim2 <- range(sapply(plot.data, function(s) {s[, 3]}), 
               na.rm=TRUE)

opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
layout(matrix(rep(c(1, 1, 2, 2), 2), 2, 4, byrow=TRUE))
plot(plot.data[[1]][, c(1, 2)], type="l", main="(a)", ylab="AUC", 
     xlab="Number of selected features", ylim=ylim1, xlim=xlim, col=col[1], 
     lty=lty[1])
lines(plot.data[[2]][, c(1, 2)], col=col[1], lty=lty[2])
lines(plot.data[[3]][, c(1, 2)], col=col[1], lty=lty[3])
lines(plot.data[[4]][, c(1, 2)], col=col[2], lty=lty[1])
lines(plot.data[[5]][, c(1, 2)], col=col[2], lty=lty[2])
lines(plot.data[[6]][, c(1, 2)], col=col[2], lty=lty[3])

legend("bottomright", legend=labels, fill=c(col, rep(NA, 3)), 
       border=c(rep(1, length(methods)), rep(NA, 3)),
       lty=c(rep(NA, length(methods)), lty),
       seg.len=1, merge=TRUE, bg="white")
plot(plot.data[[1]][, c(1, 3)], type="l", main="(b)", ylab="Brier skill score", 
     xlab="Number of selected features", ylim=ylim2, xlim=xlim, col=col[1], 
     lty=lty[1])
lines(plot.data[[2]][, c(1, 3)], col=col[1], lty=lty[2])
lines(plot.data[[3]][, c(1, 3)], col=col[1], lty=lty[3])
lines(plot.data[[4]][, c(1, 3)], col=col[2], lty=lty[1])
lines(plot.data[[5]][, c(1, 3)], col=col[2], lty=lty[2])
lines(plot.data[[6]][, c(1, 3)], col=col[2], lty=lty[3])
par(opar)

# ---- lines2_simulations_res7_auc_briers ----
library(sp)
load(file="results/simulations_res7.Rdata")
plot.data <- 
  sapply(c(paste0("gren", 1:3), paste0("enet", 1:3)), function(m) {
    cbind(rowMeans(sapply(res, function(s) {
      vec <- s$psel[substr(names(s$psel), 1, nchar(m))==m]
      replace(rep(NA, 8), 1:length(vec), vec)}), 
      na.rm=TRUE),
      rowMeans(sapply(res, function(s) {
        vec <- s$auc[substr(names(s$auc), 1, nchar(m))==m]
        replace(rep(NA, 8), 1:length(vec), vec)}), 
        na.rm=TRUE),
      rowMeans(sapply(res, function(s) {
        vec <- s$briers[substr(names(s$briers), 1, nchar(m))==m]
        replace(rep(NA, 8), 1:length(vec), vec)}), 
        na.rm=TRUE))}, simplify=FALSE)

methods <- c("gren", "enet")
labels <- c(methods, expression(alpha==0.05), expression(alpha==0.5), 
            expression(alpha==0.95))

col <- bpy.colors(length(methods), cutoff.tail=0.2)
lty <- c(1:3)

xlim <- range(sapply(plot.data, function(s) {s[, 1]}), na.rm=TRUE)
ylim1 <- range(sapply(plot.data, function(s) {s[, 2]}), na.rm=TRUE)
ylim2 <- range(sapply(plot.data, function(s) {s[, 3]}), 
               na.rm=TRUE)

opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
layout(matrix(rep(c(1, 1, 2, 2), 2), 2, 4, byrow=TRUE))
plot(plot.data[[1]][, c(1, 2)], type="l", main="(a)", ylab="AUC", 
     xlab="Number of selected features", ylim=ylim1, xlim=xlim, col=col[1], 
     lty=lty[1])
lines(plot.data[[2]][, c(1, 2)], col=col[1], lty=lty[2])
lines(plot.data[[3]][, c(1, 2)], col=col[1], lty=lty[3])
lines(plot.data[[4]][, c(1, 2)], col=col[2], lty=lty[1])
lines(plot.data[[5]][, c(1, 2)], col=col[2], lty=lty[2])
lines(plot.data[[6]][, c(1, 2)], col=col[2], lty=lty[3])

legend("bottomright", legend=labels, fill=c(col, rep(NA, 3)), 
       border=c(rep(1, length(methods)), rep(NA, 3)),
       lty=c(rep(NA, length(methods)), lty),
       seg.len=1, merge=TRUE, bg="white")
plot(plot.data[[1]][, c(1, 3)], type="l", main="(b)", ylab="Brier skill score", 
     xlab="Number of selected features", ylim=ylim2, xlim=xlim, col=col[1], 
     lty=lty[1])
lines(plot.data[[2]][, c(1, 3)], col=col[1], lty=lty[2])
lines(plot.data[[3]][, c(1, 3)], col=col[1], lty=lty[3])
lines(plot.data[[4]][, c(1, 3)], col=col[2], lty=lty[1])
lines(plot.data[[5]][, c(1, 3)], col=col[2], lty=lty[2])
lines(plot.data[[6]][, c(1, 3)], col=col[2], lty=lty[3])
par(opar)

# ---- lines2_simulations_res8_auc_briers ----
library(sp)
load(file="results/simulations_res8.Rdata")
plot.data <- 
  sapply(c(paste0("gren", 1:3), paste0("enet", 1:3)), function(m) {
    cbind(rowMeans(sapply(res, function(s) {
      vec <- s$psel[substr(names(s$psel), 1, nchar(m))==m]
      replace(rep(NA, 8), 1:length(vec), vec)}), 
      na.rm=TRUE),
      rowMeans(sapply(res, function(s) {
        vec <- s$auc[substr(names(s$auc), 1, nchar(m))==m]
        replace(rep(NA, 8), 1:length(vec), vec)}), 
        na.rm=TRUE),
      rowMeans(sapply(res, function(s) {
        vec <- s$briers[substr(names(s$briers), 1, nchar(m))==m]
        replace(rep(NA, 8), 1:length(vec), vec)}), 
        na.rm=TRUE))}, simplify=FALSE)

methods <- c("gren", "enet")
labels <- c(methods, expression(alpha==0.05), expression(alpha==0.5), 
            expression(alpha==0.95))

col <- bpy.colors(length(methods), cutoff.tail=0.2)
lty <- c(1:3)

xlim <- range(sapply(plot.data, function(s) {s[, 1]}), na.rm=TRUE)
ylim1 <- range(sapply(plot.data, function(s) {s[, 2]}), na.rm=TRUE)
ylim2 <- range(sapply(plot.data, function(s) {s[, 3]}), 
               na.rm=TRUE)

opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
layout(matrix(rep(c(1, 1, 2, 2), 2), 2, 4, byrow=TRUE))
plot(plot.data[[1]][, c(1, 2)], type="l", main="(a)", ylab="AUC", 
     xlab="Number of selected features", ylim=ylim1, xlim=xlim, col=col[1], 
     lty=lty[1])
lines(plot.data[[2]][, c(1, 2)], col=col[1], lty=lty[2])
lines(plot.data[[3]][, c(1, 2)], col=col[1], lty=lty[3])
lines(plot.data[[4]][, c(1, 2)], col=col[2], lty=lty[1])
lines(plot.data[[5]][, c(1, 2)], col=col[2], lty=lty[2])
lines(plot.data[[6]][, c(1, 2)], col=col[2], lty=lty[3])

legend("bottomright", legend=labels, fill=c(col, rep(NA, 3)), 
       border=c(rep(1, length(methods)), rep(NA, 3)),
       lty=c(rep(NA, length(methods)), lty),
       seg.len=1, merge=TRUE, bg="white")
plot(plot.data[[1]][, c(1, 3)], type="l", main="(b)", ylab="Brier skill score", 
     xlab="Number of selected features", ylim=ylim2, xlim=xlim, col=col[1], 
     lty=lty[1])
lines(plot.data[[2]][, c(1, 3)], col=col[1], lty=lty[2])
lines(plot.data[[3]][, c(1, 3)], col=col[1], lty=lty[3])
lines(plot.data[[4]][, c(1, 3)], col=col[2], lty=lty[1])
lines(plot.data[[5]][, c(1, 3)], col=col[2], lty=lty[2])
lines(plot.data[[6]][, c(1, 3)], col=col[2], lty=lty[3])
par(opar)

# ---- lines2_simulations_res9_auc_briers ----
library(sp)
load(file="results/simulations_res9.Rdata")
plot.data <- 
  sapply(c(paste0("gren", 1:3), paste0("enet", 1:3)), function(m) {
    cbind(rowMeans(sapply(res, function(s) {
      vec <- s$psel[substr(names(s$psel), 1, nchar(m))==m]
      replace(rep(NA, 8), 1:length(vec), vec)}), 
      na.rm=TRUE),
      rowMeans(sapply(res, function(s) {
        vec <- s$auc[substr(names(s$auc), 1, nchar(m))==m]
        replace(rep(NA, 8), 1:length(vec), vec)}), 
        na.rm=TRUE),
      rowMeans(sapply(res, function(s) {
        vec <- s$briers[substr(names(s$briers), 1, nchar(m))==m]
        replace(rep(NA, 8), 1:length(vec), vec)}), 
        na.rm=TRUE))}, simplify=FALSE)

methods <- c("gren", "enet")
labels <- c(methods, expression(alpha==0.05), expression(alpha==0.5), 
            expression(alpha==0.95))

col <- bpy.colors(length(methods), cutoff.tail=0.2)
lty <- c(1:3)

xlim <- range(sapply(plot.data, function(s) {s[, 1]}), na.rm=TRUE)
ylim1 <- range(sapply(plot.data, function(s) {s[, 2]}), na.rm=TRUE)
ylim2 <- range(sapply(plot.data, function(s) {s[, 3]}), 
               na.rm=TRUE)

opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
layout(matrix(rep(c(1, 1, 2, 2), 2), 2, 4, byrow=TRUE))
plot(plot.data[[1]][, c(1, 2)], type="l", main="(a)", ylab="AUC", 
     xlab="Number of selected features", ylim=ylim1, xlim=xlim, col=col[1], 
     lty=lty[1])
lines(plot.data[[2]][, c(1, 2)], col=col[1], lty=lty[2])
lines(plot.data[[3]][, c(1, 2)], col=col[1], lty=lty[3])
lines(plot.data[[4]][, c(1, 2)], col=col[2], lty=lty[1])
lines(plot.data[[5]][, c(1, 2)], col=col[2], lty=lty[2])
lines(plot.data[[6]][, c(1, 2)], col=col[2], lty=lty[3])

legend("bottomright", legend=labels, fill=c(col, rep(NA, 3)), 
       border=c(rep(1, length(methods)), rep(NA, 3)),
       lty=c(rep(NA, length(methods)), lty),
       seg.len=1, merge=TRUE, bg="white")
plot(plot.data[[1]][, c(1, 3)], type="l", main="(b)", ylab="Brier skill score", 
     xlab="Number of selected features", ylim=ylim2, xlim=xlim, col=col[1], 
     lty=lty[1])
lines(plot.data[[2]][, c(1, 3)], col=col[1], lty=lty[2])
lines(plot.data[[3]][, c(1, 3)], col=col[1], lty=lty[3])
lines(plot.data[[4]][, c(1, 3)], col=col[2], lty=lty[1])
lines(plot.data[[5]][, c(1, 3)], col=col[2], lty=lty[2])
lines(plot.data[[6]][, c(1, 3)], col=col[2], lty=lty[3])
par(opar)

# ---- lines2_simulations_res10_auc_briers ----
library(sp)
load(file="results/simulations_res10.Rdata")
plot.data <- 
  sapply(c(paste0("gren", 1:3), paste0("enet", 1:3)), function(m) {
    cbind(rowMeans(sapply(res, function(s) {
      vec <- s$psel[substr(names(s$psel), 1, nchar(m))==m]
      replace(rep(NA, 8), 1:length(vec), vec)}), 
      na.rm=TRUE),
      rowMeans(sapply(res, function(s) {
        vec <- s$auc[substr(names(s$auc), 1, nchar(m))==m]
        replace(rep(NA, 8), 1:length(vec), vec)}), 
        na.rm=TRUE),
      rowMeans(sapply(res, function(s) {
        vec <- s$briers[substr(names(s$briers), 1, nchar(m))==m]
        replace(rep(NA, 8), 1:length(vec), vec)}), 
        na.rm=TRUE))}, simplify=FALSE)

methods <- c("gren", "enet")
labels <- c(methods, expression(alpha==0.05), expression(alpha==0.5), 
            expression(alpha==0.95))

col <- bpy.colors(length(methods), cutoff.tail=0.2)
lty <- c(1:3)

xlim <- range(sapply(plot.data, function(s) {s[, 1]}), na.rm=TRUE)
ylim1 <- range(sapply(plot.data, function(s) {s[, 2]}), na.rm=TRUE)
ylim2 <- range(sapply(plot.data, function(s) {s[, 3]}), 
               na.rm=TRUE)

opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
layout(matrix(rep(c(1, 1, 2, 2), 2), 2, 4, byrow=TRUE))
plot(plot.data[[1]][, c(1, 2)], type="l", main="(a)", ylab="AUC", 
     xlab="Number of selected features", ylim=ylim1, xlim=xlim, col=col[1], 
     lty=lty[1])
lines(plot.data[[2]][, c(1, 2)], col=col[1], lty=lty[2])
lines(plot.data[[3]][, c(1, 2)], col=col[1], lty=lty[3])
lines(plot.data[[4]][, c(1, 2)], col=col[2], lty=lty[1])
lines(plot.data[[5]][, c(1, 2)], col=col[2], lty=lty[2])
lines(plot.data[[6]][, c(1, 2)], col=col[2], lty=lty[3])

legend("bottomright", legend=labels, fill=c(col, rep(NA, 3)), 
       border=c(rep(1, length(methods)), rep(NA, 3)),
       lty=c(rep(NA, length(methods)), lty),
       seg.len=1, merge=TRUE, bg="white")
plot(plot.data[[1]][, c(1, 3)], type="l", main="(b)", ylab="Brier skill score", 
     xlab="Number of selected features", ylim=ylim2, xlim=xlim, col=col[1], 
     lty=lty[1])
lines(plot.data[[2]][, c(1, 3)], col=col[1], lty=lty[2])
lines(plot.data[[3]][, c(1, 3)], col=col[1], lty=lty[3])
lines(plot.data[[4]][, c(1, 3)], col=col[2], lty=lty[1])
lines(plot.data[[5]][, c(1, 3)], col=col[2], lty=lty[2])
lines(plot.data[[6]][, c(1, 3)], col=col[2], lty=lty[3])
par(opar)

# ---- lines2_simulations_res11_auc_briers ----
library(sp)
load(file="results/simulations_res11.Rdata")
plot.data <- 
  sapply(c(paste0("gren", 1:3), paste0("enet", 1:3)), function(m) {
    cbind(rowMeans(sapply(res, function(s) {
      vec <- s$psel[substr(names(s$psel), 1, nchar(m))==m]
      replace(rep(NA, 8), 1:length(vec), vec)}), 
      na.rm=TRUE),
      rowMeans(sapply(res, function(s) {
        vec <- s$auc[substr(names(s$auc), 1, nchar(m))==m]
        replace(rep(NA, 8), 1:length(vec), vec)}), 
        na.rm=TRUE),
      rowMeans(sapply(res, function(s) {
        vec <- s$briers[substr(names(s$briers), 1, nchar(m))==m]
        replace(rep(NA, 8), 1:length(vec), vec)}), 
        na.rm=TRUE))}, simplify=FALSE)

methods <- c("gren", "enet")
labels <- c(methods, expression(alpha==0.05), expression(alpha==0.5), 
            expression(alpha==0.95))

col <- bpy.colors(length(methods), cutoff.tail=0.2)
lty <- c(1:3)

xlim <- range(sapply(plot.data, function(s) {s[, 1]}), na.rm=TRUE)
ylim1 <- range(sapply(plot.data, function(s) {s[, 2]}), na.rm=TRUE)
ylim2 <- range(sapply(plot.data, function(s) {s[, 3]}), 
               na.rm=TRUE)

opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
layout(matrix(rep(c(1, 1, 2, 2), 2), 2, 4, byrow=TRUE))
plot(plot.data[[1]][, c(1, 2)], type="l", main="(a)", ylab="AUC", 
     xlab="Number of selected features", ylim=ylim1, xlim=xlim, col=col[1], 
     lty=lty[1])
lines(plot.data[[2]][, c(1, 2)], col=col[1], lty=lty[2])
lines(plot.data[[3]][, c(1, 2)], col=col[1], lty=lty[3])
lines(plot.data[[4]][, c(1, 2)], col=col[2], lty=lty[1])
lines(plot.data[[5]][, c(1, 2)], col=col[2], lty=lty[2])
lines(plot.data[[6]][, c(1, 2)], col=col[2], lty=lty[3])

legend("bottomright", legend=labels, fill=c(col, rep(NA, 3)), 
       border=c(rep(1, length(methods)), rep(NA, 3)),
       lty=c(rep(NA, length(methods)), lty),
       seg.len=1, merge=TRUE, bg="white")
plot(plot.data[[1]][, c(1, 3)], type="l", main="(b)", ylab="Brier skill score", 
     xlab="Number of selected features", ylim=ylim2, xlim=xlim, col=col[1], 
     lty=lty[1])
lines(plot.data[[2]][, c(1, 3)], col=col[1], lty=lty[2])
lines(plot.data[[3]][, c(1, 3)], col=col[1], lty=lty[3])
lines(plot.data[[4]][, c(1, 3)], col=col[2], lty=lty[1])
lines(plot.data[[5]][, c(1, 3)], col=col[2], lty=lty[2])
lines(plot.data[[6]][, c(1, 3)], col=col[2], lty=lty[3])
par(opar)

# ---- lines2_simulations_res12_auc_briers ----
library(sp)
load(file="results/simulations_res12.Rdata")
plot.data <- 
  sapply(c(paste0("gren", 1:3), paste0("enet", 1:3)), function(m) {
    cbind(rowMeans(sapply(res, function(s) {
      vec <- s$psel[substr(names(s$psel), 1, nchar(m))==m]
      replace(rep(NA, 8), 1:length(vec), vec)}), 
      na.rm=TRUE),
      rowMeans(sapply(res, function(s) {
        vec <- s$auc[substr(names(s$auc), 1, nchar(m))==m]
        replace(rep(NA, 8), 1:length(vec), vec)}), 
        na.rm=TRUE),
      rowMeans(sapply(res, function(s) {
        vec <- s$briers[substr(names(s$briers), 1, nchar(m))==m]
        replace(rep(NA, 8), 1:length(vec), vec)}), 
        na.rm=TRUE))}, simplify=FALSE)

methods <- c("gren", "enet")
labels <- c(methods, expression(alpha==0.05), expression(alpha==0.5), 
            expression(alpha==0.95))

col <- bpy.colors(length(methods), cutoff.tail=0.2)
lty <- c(1:3)

xlim <- range(sapply(plot.data, function(s) {s[, 1]}), na.rm=TRUE)
ylim1 <- range(sapply(plot.data, function(s) {s[, 2]}), na.rm=TRUE)
ylim2 <- range(sapply(plot.data, function(s) {s[, 3]}), 
               na.rm=TRUE)

opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
layout(matrix(rep(c(1, 1, 2, 2), 2), 2, 4, byrow=TRUE))
plot(plot.data[[1]][, c(1, 2)], type="l", main="(a)", ylab="AUC", 
     xlab="Number of selected features", ylim=ylim1, xlim=xlim, col=col[1], 
     lty=lty[1])
lines(plot.data[[2]][, c(1, 2)], col=col[1], lty=lty[2])
lines(plot.data[[3]][, c(1, 2)], col=col[1], lty=lty[3])
lines(plot.data[[4]][, c(1, 2)], col=col[2], lty=lty[1])
lines(plot.data[[5]][, c(1, 2)], col=col[2], lty=lty[2])
lines(plot.data[[6]][, c(1, 2)], col=col[2], lty=lty[3])

legend("bottomright", legend=labels, fill=c(col, rep(NA, 3)), 
       border=c(rep(1, length(methods)), rep(NA, 3)),
       lty=c(rep(NA, length(methods)), lty),
       seg.len=1, merge=TRUE, bg="white")
plot(plot.data[[1]][, c(1, 3)], type="l", main="(b)", ylab="Brier skill score", 
     xlab="Number of selected features", ylim=ylim2, xlim=xlim, col=col[1], 
     lty=lty[1])
lines(plot.data[[2]][, c(1, 3)], col=col[1], lty=lty[2])
lines(plot.data[[3]][, c(1, 3)], col=col[1], lty=lty[3])
lines(plot.data[[4]][, c(1, 3)], col=col[2], lty=lty[1])
lines(plot.data[[5]][, c(1, 3)], col=col[2], lty=lty[2])
lines(plot.data[[6]][, c(1, 3)], col=col[2], lty=lty[3])
par(opar)

# ---- lines2_simulations_res13_auc_briers ----
library(sp)
load(file="results/simulations_res13.Rdata")
plot.data <- 
  sapply(c(paste0("gren", 1:3), paste0("enet", 1:3)), function(m) {
    cbind(rowMeans(sapply(res, function(s) {
      vec <- s$psel[substr(names(s$psel), 1, nchar(m))==m]
      replace(rep(NA, 8), 1:length(vec), vec)}), 
      na.rm=TRUE),
      rowMeans(sapply(res, function(s) {
        vec <- s$auc[substr(names(s$auc), 1, nchar(m))==m]
        replace(rep(NA, 8), 1:length(vec), vec)}), 
        na.rm=TRUE),
      rowMeans(sapply(res, function(s) {
        vec <- s$briers[substr(names(s$briers), 1, nchar(m))==m]
        replace(rep(NA, 8), 1:length(vec), vec)}), 
        na.rm=TRUE))}, simplify=FALSE)

methods <- c("gren", "enet")
labels <- c(methods, expression(alpha==0.05), expression(alpha==0.5), 
            expression(alpha==0.95))

col <- bpy.colors(length(methods), cutoff.tail=0.2)
lty <- c(1:3)

xlim <- range(sapply(plot.data, function(s) {s[, 1]}), na.rm=TRUE)
ylim1 <- range(sapply(plot.data, function(s) {s[, 2]}), na.rm=TRUE)
ylim2 <- range(sapply(plot.data, function(s) {s[, 3]}), 
               na.rm=TRUE)

opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
layout(matrix(rep(c(1, 1, 2, 2), 2), 2, 4, byrow=TRUE))
plot(plot.data[[1]][, c(1, 2)], type="l", main="(a)", ylab="AUC", 
     xlab="Number of selected features", ylim=ylim1, xlim=xlim, col=col[1], 
     lty=lty[1])
lines(plot.data[[2]][, c(1, 2)], col=col[1], lty=lty[2])
lines(plot.data[[3]][, c(1, 2)], col=col[1], lty=lty[3])
lines(plot.data[[4]][, c(1, 2)], col=col[2], lty=lty[1])
lines(plot.data[[5]][, c(1, 2)], col=col[2], lty=lty[2])
lines(plot.data[[6]][, c(1, 2)], col=col[2], lty=lty[3])

legend("bottomright", legend=labels, fill=c(col, rep(NA, 3)), 
       border=c(rep(1, length(methods)), rep(NA, 3)),
       lty=c(rep(NA, length(methods)), lty),
       seg.len=1, merge=TRUE, bg="white")
plot(plot.data[[1]][, c(1, 3)], type="l", main="(b)", ylab="Brier skill score", 
     xlab="Number of selected features", ylim=ylim2, xlim=xlim, col=col[1], 
     lty=lty[1])
lines(plot.data[[2]][, c(1, 3)], col=col[1], lty=lty[2])
lines(plot.data[[3]][, c(1, 3)], col=col[1], lty=lty[3])
lines(plot.data[[4]][, c(1, 3)], col=col[2], lty=lty[1])
lines(plot.data[[5]][, c(1, 3)], col=col[2], lty=lty[2])
lines(plot.data[[6]][, c(1, 3)], col=col[2], lty=lty[3])
par(opar)

# ---- lines2_simulations_res14_auc_briers ----
library(sp)
load(file="results/simulations_res14.Rdata")
plot.data <- 
  sapply(c(paste0("gren", 1:3), paste0("enet", 1:3)), function(m) {
    cbind(rowMeans(sapply(res, function(s) {
      vec <- s$psel[substr(names(s$psel), 1, nchar(m))==m]
      replace(rep(NA, 8), 1:length(vec), vec)}), 
      na.rm=TRUE),
      rowMeans(sapply(res, function(s) {
        vec <- s$auc[substr(names(s$auc), 1, nchar(m))==m]
        replace(rep(NA, 8), 1:length(vec), vec)}), 
        na.rm=TRUE),
      rowMeans(sapply(res, function(s) {
        vec <- s$briers[substr(names(s$briers), 1, nchar(m))==m]
        replace(rep(NA, 8), 1:length(vec), vec)}), 
        na.rm=TRUE))}, simplify=FALSE)

methods <- c("gren", "enet")
labels <- c(methods, expression(alpha==0.05), expression(alpha==0.5), 
            expression(alpha==0.95))

col <- bpy.colors(length(methods), cutoff.tail=0.2)
lty <- c(1:3)

xlim <- range(sapply(plot.data, function(s) {s[, 1]}), na.rm=TRUE)
ylim1 <- range(sapply(plot.data, function(s) {s[, 2]}), na.rm=TRUE)
ylim2 <- range(sapply(plot.data, function(s) {s[, 3]}), 
               na.rm=TRUE)

opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
layout(matrix(rep(c(1, 1, 2, 2), 2), 2, 4, byrow=TRUE))
plot(plot.data[[1]][, c(1, 2)], type="l", main="(a)", ylab="AUC", 
     xlab="Number of selected features", ylim=ylim1, xlim=xlim, col=col[1], 
     lty=lty[1])
lines(plot.data[[2]][, c(1, 2)], col=col[1], lty=lty[2])
lines(plot.data[[3]][, c(1, 2)], col=col[1], lty=lty[3])
lines(plot.data[[4]][, c(1, 2)], col=col[2], lty=lty[1])
lines(plot.data[[5]][, c(1, 2)], col=col[2], lty=lty[2])
lines(plot.data[[6]][, c(1, 2)], col=col[2], lty=lty[3])

legend("bottomright", legend=labels, fill=c(col, rep(NA, 3)), 
       border=c(rep(1, length(methods)), rep(NA, 3)),
       lty=c(rep(NA, length(methods)), lty),
       seg.len=1, merge=TRUE, bg="white")
plot(plot.data[[1]][, c(1, 3)], type="l", main="(b)", ylab="Brier skill score", 
     xlab="Number of selected features", ylim=ylim2, xlim=xlim, col=col[1], 
     lty=lty[1])
lines(plot.data[[2]][, c(1, 3)], col=col[1], lty=lty[2])
lines(plot.data[[3]][, c(1, 3)], col=col[1], lty=lty[3])
lines(plot.data[[4]][, c(1, 3)], col=col[2], lty=lty[1])
lines(plot.data[[5]][, c(1, 3)], col=col[2], lty=lty[2])
lines(plot.data[[6]][, c(1, 3)], col=col[2], lty=lty[3])
par(opar)

# ---- lines2_simulations_res15_auc_briers ----
library(sp)
load(file="results/simulations_res15.Rdata")
plot.data <- 
  sapply(c(paste0("gren", 1:3), paste0("enet", 1:3)), function(m) {
    cbind(rowMeans(sapply(res, function(s) {
      vec <- s$psel[substr(names(s$psel), 1, nchar(m))==m]
      replace(rep(NA, 8), 1:length(vec), vec)}), 
      na.rm=TRUE),
      rowMeans(sapply(res, function(s) {
        vec <- s$auc[substr(names(s$auc), 1, nchar(m))==m]
        replace(rep(NA, 8), 1:length(vec), vec)}), 
        na.rm=TRUE),
      rowMeans(sapply(res, function(s) {
        vec <- s$briers[substr(names(s$briers), 1, nchar(m))==m]
        replace(rep(NA, 8), 1:length(vec), vec)}), 
        na.rm=TRUE))}, simplify=FALSE)

methods <- c("gren", "enet")
labels <- c(methods, expression(alpha==0.05), expression(alpha==0.5), 
            expression(alpha==0.95))

col <- bpy.colors(length(methods), cutoff.tail=0.2)
lty <- c(1:3)

xlim <- range(sapply(plot.data, function(s) {s[, 1]}), na.rm=TRUE)
ylim1 <- range(sapply(plot.data, function(s) {s[, 2]}), na.rm=TRUE)
ylim2 <- range(sapply(plot.data, function(s) {s[, 3]}), 
               na.rm=TRUE)

opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
layout(matrix(rep(c(1, 1, 2, 2), 2), 2, 4, byrow=TRUE))
plot(plot.data[[1]][, c(1, 2)], type="l", main="(a)", ylab="AUC", 
     xlab="Number of selected features", ylim=ylim1, xlim=xlim, col=col[1], 
     lty=lty[1])
lines(plot.data[[2]][, c(1, 2)], col=col[1], lty=lty[2])
lines(plot.data[[3]][, c(1, 2)], col=col[1], lty=lty[3])
lines(plot.data[[4]][, c(1, 2)], col=col[2], lty=lty[1])
lines(plot.data[[5]][, c(1, 2)], col=col[2], lty=lty[2])
lines(plot.data[[6]][, c(1, 2)], col=col[2], lty=lty[3])

legend("bottomright", legend=labels, fill=c(col, rep(NA, 3)), 
       border=c(rep(1, length(methods)), rep(NA, 3)),
       lty=c(rep(NA, length(methods)), lty),
       seg.len=1, merge=TRUE, bg="white")
plot(plot.data[[1]][, c(1, 3)], type="l", main="(b)", ylab="Brier skill score", 
     xlab="Number of selected features", ylim=ylim2, xlim=xlim, col=col[1], 
     lty=lty[1])
lines(plot.data[[2]][, c(1, 3)], col=col[1], lty=lty[2])
lines(plot.data[[3]][, c(1, 3)], col=col[1], lty=lty[3])
lines(plot.data[[4]][, c(1, 3)], col=col[2], lty=lty[1])
lines(plot.data[[5]][, c(1, 3)], col=col[2], lty=lty[2])
lines(plot.data[[6]][, c(1, 3)], col=col[2], lty=lty[3])
par(opar)

# ---- lines2_simulations_res16_auc_briers ----
library(sp)
load(file="results/simulations_res16.Rdata")
plot.data <- 
  sapply(c(paste0("gren", 1:3), paste0("enet", 1:3)), function(m) {
    cbind(rowMeans(sapply(res, function(s) {
      vec <- s$psel[substr(names(s$psel), 1, nchar(m))==m]
      replace(rep(NA, 8), 1:length(vec), vec)}), 
      na.rm=TRUE),
      rowMeans(sapply(res, function(s) {
        vec <- s$auc[substr(names(s$auc), 1, nchar(m))==m]
        replace(rep(NA, 8), 1:length(vec), vec)}), 
        na.rm=TRUE),
      rowMeans(sapply(res, function(s) {
        vec <- s$briers[substr(names(s$briers), 1, nchar(m))==m]
        replace(rep(NA, 8), 1:length(vec), vec)}), 
        na.rm=TRUE))}, simplify=FALSE)

methods <- c("gren", "enet")
labels <- c(methods, expression(alpha==0.05), expression(alpha==0.5), 
            expression(alpha==0.95))

col <- bpy.colors(length(methods), cutoff.tail=0.2)
lty <- c(1:3)

xlim <- range(sapply(plot.data, function(s) {s[, 1]}), na.rm=TRUE)
ylim1 <- range(sapply(plot.data, function(s) {s[, 2]}), na.rm=TRUE)
ylim2 <- range(sapply(plot.data, function(s) {s[, 3]}), 
               na.rm=TRUE)

opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
layout(matrix(rep(c(1, 1, 2, 2), 2), 2, 4, byrow=TRUE))
plot(plot.data[[1]][, c(1, 2)], type="l", main="(a)", ylab="AUC", 
     xlab="Number of selected features", ylim=ylim1, xlim=xlim, col=col[1], 
     lty=lty[1])
lines(plot.data[[2]][, c(1, 2)], col=col[1], lty=lty[2])
lines(plot.data[[3]][, c(1, 2)], col=col[1], lty=lty[3])
lines(plot.data[[4]][, c(1, 2)], col=col[2], lty=lty[1])
lines(plot.data[[5]][, c(1, 2)], col=col[2], lty=lty[2])
lines(plot.data[[6]][, c(1, 2)], col=col[2], lty=lty[3])

legend("bottomright", legend=labels, fill=c(col, rep(NA, 3)), 
       border=c(rep(1, length(methods)), rep(NA, 3)),
       lty=c(rep(NA, length(methods)), lty),
       seg.len=1, merge=TRUE, bg="white")
plot(plot.data[[1]][, c(1, 3)], type="l", main="(b)", ylab="Brier skill score", 
     xlab="Number of selected features", ylim=ylim2, xlim=xlim, col=col[1], 
     lty=lty[1])
lines(plot.data[[2]][, c(1, 3)], col=col[1], lty=lty[2])
lines(plot.data[[3]][, c(1, 3)], col=col[1], lty=lty[3])
lines(plot.data[[4]][, c(1, 3)], col=col[2], lty=lty[1])
lines(plot.data[[5]][, c(1, 3)], col=col[2], lty=lty[2])
lines(plot.data[[6]][, c(1, 3)], col=col[2], lty=lty[3])
par(opar)

# ---- lines2_simulations_res17_auc_briers ----
library(sp)
load(file="results/simulations_res17.Rdata")
plot.data <- 
  sapply(c(paste0("gren", 1:3), paste0("enet", 1:3)), function(m) {
    cbind(rowMeans(sapply(res, function(s) {
      vec <- s$psel[substr(names(s$psel), 1, nchar(m))==m]
      replace(rep(NA, 8), 1:length(vec), vec)}), 
      na.rm=TRUE),
      rowMeans(sapply(res, function(s) {
        vec <- s$auc[substr(names(s$auc), 1, nchar(m))==m]
        replace(rep(NA, 8), 1:length(vec), vec)}), 
        na.rm=TRUE),
      rowMeans(sapply(res, function(s) {
        vec <- s$briers[substr(names(s$briers), 1, nchar(m))==m]
        replace(rep(NA, 8), 1:length(vec), vec)}), 
        na.rm=TRUE))}, simplify=FALSE)

methods <- c("gren", "enet")
labels <- c(methods, expression(alpha==0.05), expression(alpha==0.5), 
            expression(alpha==0.95))

col <- bpy.colors(length(methods), cutoff.tail=0.2)
lty <- c(1:3)

xlim <- range(sapply(plot.data, function(s) {s[, 1]}), na.rm=TRUE)
ylim1 <- range(sapply(plot.data, function(s) {s[, 2]}), na.rm=TRUE)
ylim2 <- range(sapply(plot.data, function(s) {s[, 3]}), 
               na.rm=TRUE)

opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
layout(matrix(rep(c(1, 1, 2, 2), 2), 2, 4, byrow=TRUE))
plot(plot.data[[1]][, c(1, 2)], type="l", main="(a)", ylab="AUC", 
     xlab="Number of selected features", ylim=ylim1, xlim=xlim, col=col[1], 
     lty=lty[1])
lines(plot.data[[2]][, c(1, 2)], col=col[1], lty=lty[2])
lines(plot.data[[3]][, c(1, 2)], col=col[1], lty=lty[3])
lines(plot.data[[4]][, c(1, 2)], col=col[2], lty=lty[1])
lines(plot.data[[5]][, c(1, 2)], col=col[2], lty=lty[2])
lines(plot.data[[6]][, c(1, 2)], col=col[2], lty=lty[3])

legend("bottomright", legend=labels, fill=c(col, rep(NA, 3)), 
       border=c(rep(1, length(methods)), rep(NA, 3)),
       lty=c(rep(NA, length(methods)), lty),
       seg.len=1, merge=TRUE, bg="white")
plot(plot.data[[1]][, c(1, 3)], type="l", main="(b)", ylab="Brier skill score", 
     xlab="Number of selected features", ylim=ylim2, xlim=xlim, col=col[1], 
     lty=lty[1])
lines(plot.data[[2]][, c(1, 3)], col=col[1], lty=lty[2])
lines(plot.data[[3]][, c(1, 3)], col=col[1], lty=lty[3])
lines(plot.data[[4]][, c(1, 3)], col=col[2], lty=lty[1])
lines(plot.data[[5]][, c(1, 3)], col=col[2], lty=lty[2])
lines(plot.data[[6]][, c(1, 3)], col=col[2], lty=lty[3])
par(opar)


# ---- lines_micrornaseq_colorectal_cancer_res1_auc ----
library(sp)
res <- read.table("results/micrornaseq_colorectal_cancer_res1.csv", 
                  stringsAsFactors=FALSE)

pred <- as.matrix(res[substr(rownames(res), 1, 4)=="pred", ])
psel <- as.matrix(res[substr(rownames(res), 1, 4)=="psel", ])
auc <- as.matrix(res[substr(rownames(res), 1, 3)=="auc", ])

methods <- c("grridge2", "ridge", "rf", paste0("gren", c(4:6)), 
             paste0("enet", c(1:3)), paste0("sgl", c(1:3)), 
             paste0("cmcp", c(1:3)), paste0("gel", c(1:3)))
labels <- list("ridge", "GRridge", "random forest", "gren", "enet", "SGL", 
               "cMCP", "gel")
col <- bpy.colors(length(labels), cutoff.tail=0.1)
lty <- 1:length(methods)

plot.data <- lapply(methods, function(m) {
  aggregate(auc[, substr(colnames(auc), 1, nchar(m))==m], 
            list(psel=psel[, substr(colnames(psel), 1, nchar(m))==m]), mean)})
names(plot.data) <- methods

opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
layout(matrix(c(rep(c(1, 1, 2, 2), 2), rep(c(0, 3, 3, 0), 2)), 
              4, 4, byrow=TRUE))
plot(plot.data[[4]], type="l", xlab="Number of selected features", 
     ylab="AUC", main="(a)", ylim=range(auc), xlim=c(0, 500), 
     col=col[4], lty=lty[4])
lines(plot.data[[7]], col=col[5], lty=lty[5])
lines(plot.data[[10]], col=col[6], lty=lty[6])
lines(plot.data[[13]], col=col[7], lty=lty[7])
lines(plot.data[[16]], col=col[8], lty=lty[8])
abline(h=plot.data[[1]][, 2], col=col[1], lty=lty[1])
abline(h=plot.data[[2]][, 2], col=col[2], lty=lty[2])
abline(h=plot.data[[3]][, 2], col=col[3], lty=lty[3])
legend("bottomright", legend=labels, col=col, 
       lty=lty, bg = "white")

plot(plot.data[[5]], type="l", xlab="Number of selected features", 
     ylab="AUC", main="(b)", ylim=range(auc), xlim=c(0, 500), 
     col=col[4], lty=lty[4])
lines(plot.data[[8]], col=col[5], lty=lty[5])
lines(plot.data[[11]], col=col[6], lty=lty[6])
lines(plot.data[[14]], col=col[7], lty=lty[7])
lines(plot.data[[17]], col=col[8], lty=lty[8])
abline(h=plot.data[[1]][, 2], col=col[1], lty=lty[1])
abline(h=plot.data[[2]][, 2], col=col[2], lty=lty[2])
abline(h=plot.data[[3]][, 2], col=col[3], lty=lty[3])

plot(plot.data[[6]], type="l", xlab="Number of selected features", ylab="AUC", 
     main="(c)", ylim=range(auc), xlim=c(0, 500), col=col[4], 
     lty=lty[4])
lines(plot.data[[9]], col=col[5], lty=lty[5])
lines(plot.data[[12]], col=col[6], lty=lty[6])
lines(plot.data[[15]], col=col[7], lty=lty[7])
lines(plot.data[[18]], col=col[8], lty=lty[8])
abline(h=plot.data[[1]][, 2], col=col[1], lty=lty[1])
abline(h=plot.data[[2]][, 2], col=col[2], lty=lty[2])
abline(h=plot.data[[3]][, 2], col=col[3], lty=lty[3])
par(opar)

# ---- lines_micrornaseq_colorectal_cancer_res1_briers ----
library(sp)
res <- read.table("results/micrornaseq_colorectal_cancer_res1.csv", 
                  stringsAsFactors=FALSE)

pred <- as.matrix(res[substr(rownames(res), 1, 4)=="pred", ])
psel <- as.matrix(res[substr(rownames(res), 1, 4)=="psel", ])
briers <- as.matrix(res[substr(rownames(res), 1, 6)=="briers", ])

methods <- c("grridge2", "rf", "ridge", paste0("gren", c(4:6)), 
             paste0("enet", c(1:3)), paste0("sgl", c(1:3)), 
             paste0("cmcp", c(1:3)), paste0("gel", c(1:3)))
labels <- list("ridge", "GRridge", "random forest", "gren", "enet", "SGL", 
               "cMCP", "gel")
col <- bpy.colors(length(labels), cutoff.tail=0.1)
lty <- 1:length(methods)

plot.data <- lapply(methods, function(m) {
  aggregate(briers[, substr(colnames(briers), 1, nchar(m))==m], 
            list(psel=psel[, substr(colnames(psel), 1, nchar(m))==m]), mean)})
names(plot.data) <- methods

opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
layout(matrix(c(rep(c(1, 1, 2, 2), 2), rep(c(0, 3, 3, 0), 2)), 
              4, 4, byrow=TRUE))
plot(plot.data[[4]], type="l", xlab="Number of selected features", 
     ylab="Brier skill score", main="(a)", ylim=range(briers), xlim=c(0, 500), 
     col=col[4], lty=lty[4])
lines(plot.data[[7]], col=col[5], lty=lty[5])
lines(plot.data[[10]], col=col[6], lty=lty[6])
lines(plot.data[[13]], col=col[7], lty=lty[7])
lines(plot.data[[16]], col=col[8], lty=lty[8])
abline(h=plot.data[[1]][, 2], col=col[1], lty=lty[1])
abline(h=plot.data[[2]][, 2], col=col[2], lty=lty[2])
abline(h=plot.data[[3]][, 2], col=col[3], lty=lty[3])
legend("bottomright", legend=labels, col=col, 
       lty=lty, bg = "white")

plot(plot.data[[5]], type="l", xlab="Number of selected features", 
     ylab="Brier skill score", main="(b)", ylim=range(briers), xlim=c(0, 500), 
     col=col[4], lty=lty[4])
lines(plot.data[[8]], col=col[5], lty=lty[5])
lines(plot.data[[11]], col=col[6], lty=lty[6])
lines(plot.data[[14]], col=col[7], lty=lty[7])
lines(plot.data[[17]], col=col[8], lty=lty[8])
abline(h=plot.data[[1]][, 2], col=col[1], lty=lty[1])
abline(h=plot.data[[2]][, 2], col=col[2], lty=lty[2])
abline(h=plot.data[[3]][, 2], col=col[3], lty=lty[3])

plot(plot.data[[6]], type="l", xlab="Number of selected features", 
     ylab="Brier skill score", 
     main="(c)", ylim=range(briers), xlim=c(0, 500), col=col[4], 
     lty=lty[4])
lines(plot.data[[9]], col=col[5], lty=lty[5])
lines(plot.data[[12]], col=col[6], lty=lty[6])
lines(plot.data[[15]], col=col[7], lty=lty[7])
lines(plot.data[[18]], col=col[8], lty=lty[8])
abline(h=plot.data[[1]][, 2], col=col[1], lty=lty[1])
abline(h=plot.data[[2]][, 2], col=col[2], lty=lty[2])
abline(h=plot.data[[3]][, 2], col=col[3], lty=lty[3])
par(opar)

# ---- boxplots_micrornaseq_colorectal_cancer_res2_multipliers ----
library(sp)
res <- read.table("results/micrornaseq_colorectal_cancer_res2.csv", 
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

# ---- boxplots_micrornaseq_colorectal_cancer_res3_multipliers ----
res <- read.table("results/micrornaseq_colorectal_cancer_res3.csv", 
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

# ---- hist_micrornaseq_colorectal_cancer_res4_overlap ----
library(sp)
res <- read.table("results/micrornaseq_colorectal_cancer_res4.csv", 
                  stringsAsFactors=FALSE)

methods <- colnames(res)
labels <- c("gren", "enet")
col <- bpy.colors(length(labels), cutoff.tail=0.1)

plot.data <- sapply(c(1:3), function(m) {
  sapply(1:max(res), function(s) {
    c(sum(res[[paste0("gren", m)]]==s), sum(res[[paste0("enet", m)]]==s))/
      nrow(res)})}, simplify=FALSE)

opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
layout(matrix(c(rep(c(1, 1, 2, 2), 2), rep(c(0, 3, 3, 0), 2)), 
              nrow=4, ncol=4, byrow=TRUE))
barplot(plot.data[[1]], beside=TRUE, xlab="Number of overlapping features",
        ylab="Density", main="(a)", col=col, names.arg=1:max(res), 
        legend.text=labels, args.legend=list(x="topright", fill=col, 
                                             border=col))
barplot(plot.data[[2]], beside=TRUE, xlab="Number of overlapping features",
        ylab="Density", main="(b)", col=col, names.arg=1:max(res))
barplot(plot.data[[3]], beside=TRUE, xlab="Number of overlapping features",
        ylab="Density", main="(c)", col=col, names.arg=1:max(res))
par(opar)

# ---- lines_rnaseq_oral_cancer_metastasis_res1_auc ----
library(sp)
res <- read.table("results/rnaseq_oral_cancer_metastasis_res1.csv", 
                  stringsAsFactors=FALSE)

pred <- as.matrix(res[substr(rownames(res), 1, 4)=="pred", ])
psel <- as.matrix(res[substr(rownames(res), 1, 4)=="psel", ])
auc <- as.matrix(res[substr(rownames(res), 1, 3)=="auc", ])

methods <- c("ridge", "grridge", "rf", paste0("gren", c(1:3)), 
             paste0("enet", c(1:3)), paste0("sgl", c(1:3)), 
             paste0("cmcp", c(1:3)), paste0("gel", c(1:3)), 
             paste0("ocmcp", c(1:3)), paste0("ogel", c(1:3)))
labels <- list("ridge", "GRridge", "random forest", "gren", "enet", "SGL", 
               "cMCP", "gel", "OcMCP", "Ogel")
col <- bpy.colors(length(labels), cutoff.tail=0.1)
lty <- 1:length(labels)

plot.data <- lapply(methods, function(m) {
  aggregate(auc[, substr(colnames(auc), 1, nchar(m))==m], 
            list(psel=psel[, substr(colnames(psel), 1, nchar(m))==m]), mean)})
names(plot.data) <- methods

opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
layout(matrix(c(rep(c(1, 1, 2, 2), 2), rep(c(0, 3, 3, 0), 2)), 
       4, 4, byrow=TRUE))
plot(plot.data[[4]], type="l", xlab="Number of selected features", 
     ylab="AUC", main="(a)", ylim=range(auc), xlim=c(0, 500), 
     col=col[4], lty=lty[4])
lines(plot.data[[7]], col=col[5], lty=lty[5])
lines(plot.data[[10]], col=col[6], lty=lty[6])
lines(plot.data[[13]], col=col[7], lty=lty[7])
lines(plot.data[[16]], col=col[8], lty=lty[8])
lines(plot.data[[19]], col=col[9], lty=lty[9])
lines(plot.data[[22]], col=col[10], lty=lty[10])
abline(h=plot.data[[1]][, 2], col=col[1], lty=lty[1])
abline(h=plot.data[[2]][, 2], col=col[2], lty=lty[2])
abline(h=plot.data[[3]][, 2], col=col[3], lty=lty[3])
legend("bottomright", legend=labels, col=col, 
       lty=lty, bg = "white")

plot(plot.data[[5]], type="l", xlab="Number of selected features", 
     ylab="AUC", main="(b)", ylim=range(auc), xlim=c(0, 500), 
     col=col[4], lty=lty[4])
lines(plot.data[[8]], col=col[5], lty=lty[5])
lines(plot.data[[11]], col=col[6], lty=lty[6])
lines(plot.data[[14]], col=col[7], lty=lty[7])
lines(plot.data[[17]], col=col[8], lty=lty[8])
lines(plot.data[[20]], col=col[9], lty=lty[9])
lines(plot.data[[23]], col=col[10], lty=lty[10])
abline(h=plot.data[[1]][, 2], col=col[1], lty=lty[1])
abline(h=plot.data[[2]][, 2], col=col[2], lty=lty[2])
abline(h=plot.data[[3]][, 2], col=col[3], lty=lty[3])

plot(plot.data[[6]], type="l", xlab="Number of selected features", ylab="AUC", 
     main="(c)", ylim=range(auc), xlim=c(0, 500), col=col[4], 
     lty=lty[4])
lines(plot.data[[9]], col=col[5], lty=lty[5])
lines(plot.data[[10]], col=col[6], lty=lty[6])
lines(plot.data[[15]], col=col[7], lty=lty[7])
lines(plot.data[[18]], col=col[8], lty=lty[8])
lines(plot.data[[21]], col=col[9], lty=lty[9])
lines(plot.data[[24]], col=col[10], lty=lty[10])
abline(h=plot.data[[1]][, 2], col=col[1], lty=lty[1])
abline(h=plot.data[[2]][, 2], col=col[2], lty=lty[2])
abline(h=plot.data[[3]][, 2], col=col[3], lty=lty[3])

# ---- lines_rnaseq_oral_cancer_metastasis_res1_briers ----
library(sp)
res <- read.table("results/rnaseq_oral_cancer_metastasis_res1.csv", 
                  stringsAsFactors=FALSE)

pred <- as.matrix(res[substr(rownames(res), 1, 4)=="pred", ])
psel <- as.matrix(res[substr(rownames(res), 1, 4)=="psel", ])
briers <- as.matrix(res[substr(rownames(res), 1, 6)=="briers", ])

methods <- c("ridge", "grridge", "rf", paste0("gren", c(1:3)), 
             paste0("enet", c(1:3)), paste0("sgl", c(1:3)), 
             paste0("cmcp", c(1:3)), paste0("gel", c(1:3)), 
             paste0("ocmcp", c(1:3)), paste0("ogel", c(1:3)))
labels <- list("ridge", "GRridge", "random forest", "gren", "enet", "SGL", 
               "cMCP", "gel", "OcMCP", "Ogel")
col <- bpy.colors(length(labels), cutoff.tail=0.1)
lty <- 1:length(labels)

plot.data <- lapply(methods, function(m) {
  aggregate(briers[, substr(colnames(briers), 1, nchar(m))==m], 
            list(psel=psel[, substr(colnames(psel), 1, nchar(m))==m]), mean)})
names(plot.data) <- methods

opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
layout(matrix(c(rep(c(1, 1, 2, 2), 2), rep(c(0, 3, 3, 0), 2)), 
              4, 4, byrow=TRUE))
plot(plot.data[[4]], type="l", xlab="Number of selected features", 
     ylab="Brier skill score", main="(a)", ylim=range(briers), xlim=c(0, 500), 
     col=col[4], lty=lty[4])
lines(plot.data[[7]], col=col[5], lty=lty[5])
lines(plot.data[[10]], col=col[6], lty=lty[6])
lines(plot.data[[13]], col=col[7], lty=lty[7])
lines(plot.data[[16]], col=col[8], lty=lty[8])
lines(plot.data[[19]], col=col[9], lty=lty[9])
lines(plot.data[[22]], col=col[10], lty=lty[10])
abline(h=plot.data[[1]][, 2], col=col[1], lty=lty[1])
abline(h=plot.data[[2]][, 2], col=col[2], lty=lty[2])
abline(h=plot.data[[3]][, 2], col=col[3], lty=lty[3])
legend("bottomright", legend=labels, col=col, 
       lty=lty, bg = "white")

plot(plot.data[[5]], type="l", xlab="Number of selected features", 
     ylab="Brier skill score", main="(b)", ylim=range(briers), xlim=c(0, 500), 
     col=col[4], lty=lty[4])
lines(plot.data[[8]], col=col[5], lty=lty[5])
lines(plot.data[[11]], col=col[6], lty=lty[6])
lines(plot.data[[14]], col=col[7], lty=lty[7])
lines(plot.data[[17]], col=col[8], lty=lty[8])
lines(plot.data[[20]], col=col[9], lty=lty[9])
lines(plot.data[[23]], col=col[10], lty=lty[10])
abline(h=plot.data[[1]][, 2], col=col[1], lty=lty[1])
abline(h=plot.data[[2]][, 2], col=col[2], lty=lty[2])
abline(h=plot.data[[3]][, 2], col=col[3], lty=lty[3])

plot(plot.data[[6]], type="l", xlab="Number of selected features", 
     ylab="Brier skill score", main="(c)", ylim=range(briers), xlim=c(0, 500), 
     col=col[4], lty=lty[4])
lines(plot.data[[9]], col=col[5], lty=lty[5])
lines(plot.data[[10]], col=col[6], lty=lty[6])
lines(plot.data[[15]], col=col[7], lty=lty[7])
lines(plot.data[[18]], col=col[8], lty=lty[8])
lines(plot.data[[21]], col=col[9], lty=lty[9])
lines(plot.data[[24]], col=col[10], lty=lty[10])
abline(h=plot.data[[1]][, 2], col=col[1], lty=lty[1])
abline(h=plot.data[[2]][, 2], col=col[2], lty=lty[2])
abline(h=plot.data[[3]][, 2], col=col[3], lty=lty[3])
par(opar)

# ---- barplot_metabolomics_alzheimer_res1 ----
library(sp)
load("results/metabolomics_alzheimer_fit1.Rdata")

methods <- c("grridge", expression("gren, "~alpha==0.05, "gren, "~alpha==0.5,
                                   "gren, "~alpha==0.95), "no co-data")
col <- grey.colors(length(methods) - 1, start=0.3, end=0.9, gamma=2.2, 
                   alpha=NULL)
lty <- 2

labels <- list(quality=expression("RSDqc" > 0.113, 
                                  atop(NA, atop(textstyle(0.113 >=phantom(0)), 
                                                textstyle("RSDqc" > 0.055))),
                                  "RSDqc" <= 0.055),
               degree=expression("degree" == 0, 
                                 atop(NA, atop(textstyle(0 < "degree"),
                                               textstyle(
                                                 phantom(0)<= "average"))),
                                 atop(NA, atop(textstyle("degree" > phantom(0)),
                                               textstyle("average")))))

plot.data <- lapply(1:length(fit.gren4$lambdag), function(s) {
  rbind(fit.grridge$lambdamults[[s]], fit.gren4$lambdag[[s]],
        fit.gren5$lambdag[[s]], fit.gren6$lambdag[[s]])})

opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
layout(matrix(rep(c(1, 1, 2, 2), 2), 2, 4, byrow=TRUE))
barplot(plot.data[[1]], beside=TRUE, col=col, legend.text=methods,
        args.legend=list(x="bottomright", fill=c(col, NA),
                         border=c(rep(1, length(methods) - 1), NA),
                         lty=c(rep(NA, length(methods) - 1), lty[1]),
                         seg.len=1, merge=TRUE, bg="white"),
        names.arg=labels[[1]], ylab=expression(hat(lambda)~"'"[g]),
        main="")
abline(h=1, lty=2)

barplot(plot.data[[2]], beside=TRUE, col=col, main="(b)", 
        names.arg=labels[[2]], ylab=expression(hat(lambda)~"'"[g]), las=2)
abline(h=1, lty=2)
par(opar)

# ---- lines_metabolomics_alzheimer_res1_auc ----
library(sp)
res <- read.table("results/metabolomics_alzheimer_res1.csv", 
                  stringsAsFactors=FALSE)

pred <- as.matrix(res[substr(rownames(res), 1, 4)=="pred", ])
psel <- as.matrix(res[substr(rownames(res), 1, 4)=="psel", ])
auc <- as.matrix(res[substr(rownames(res), 1, 3)=="auc", ])

methods <- c("ridge", "grridge", "rf", paste0("gren", c(4:6)), 
             paste0("enet", c(1:3)),
             paste0("sgl", c(1:3)), paste0("cmcp", c(1:3)), 
             paste0("gel", c(1:3)), paste0("ocmcp", c(1:3)), 
             paste0("ogel", c(1:3)))
labels <- list("ridge", "GRridge", "random forest", "gren", "enet", "SGL", 
               "cMCP", "gel", "OcMCP", "Ogel")
col <- bpy.colors(length(labels), cutoff.tail=0.1)
lty <- 1:length(labels)

plot.data <- lapply(methods, function(m) {
  aggregate(auc[, substr(colnames(auc), 1, nchar(m))==m], 
            list(psel=psel[, substr(colnames(psel), 1, nchar(m))==m]), mean)})
names(plot.data) <- methods

opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
layout(matrix(c(rep(c(1, 1, 2, 2), 2), rep(c(0, 3, 3, 0), 2)), 
              4, 4, byrow=TRUE))
plot(plot.data[[4]], type="l", xlab="Number of selected features", 
     ylab="AUC", main="(a)", ylim=range(auc), xlim=range(psel), 
     col=col[4], lty=lty[4])
lines(plot.data[[7]], col=col[5], lty=lty[5])
lines(plot.data[[10]], col=col[6], lty=lty[6])
lines(plot.data[[13]], col=col[7], lty=lty[7])
lines(plot.data[[16]], col=col[8], lty=lty[8])
lines(plot.data[[19]], col=col[9], lty=lty[9])
lines(plot.data[[22]], col=col[10], lty=lty[10])
abline(h=plot.data[[1]][, 2], col=col[1], lty=lty[1])
abline(h=plot.data[[2]][, 2], col=col[2], lty=lty[2])
abline(h=plot.data[[3]][, 2], col=col[3], lty=lty[3])
legend("bottomright", legend=labels, col=col, 
       lty=lty, bg = "white")

plot(plot.data[[5]], type="l", xlab="Number of selected features", 
     ylab="AUC", main="(b)", ylim=range(auc), xlim=range(psel), 
     col=col[4], lty=lty[4])
lines(plot.data[[8]], col=col[5], lty=lty[5])
lines(plot.data[[11]], col=col[6], lty=lty[6])
lines(plot.data[[14]], col=col[7], lty=lty[7])
lines(plot.data[[17]], col=col[8], lty=lty[8])
lines(plot.data[[20]], col=col[9], lty=lty[9])
lines(plot.data[[23]], col=col[10], lty=lty[10])
abline(h=plot.data[[1]][, 2], col=col[1], lty=lty[1])
abline(h=plot.data[[2]][, 2], col=col[2], lty=lty[2])
abline(h=plot.data[[3]][, 2], col=col[3], lty=lty[3])

plot(plot.data[[6]], type="l", xlab="Number of selected features", ylab="AUC", 
     main="(c)", ylim=range(auc), xlim=range(psel), col=col[4], 
     lty=lty[4])
lines(plot.data[[9]], col=col[5], lty=lty[5])
lines(plot.data[[12]], col=col[6], lty=lty[6])
lines(plot.data[[15]], col=col[7], lty=lty[7])
lines(plot.data[[18]], col=col[8], lty=lty[8])
lines(plot.data[[21]], col=col[9], lty=lty[9])
lines(plot.data[[24]], col=col[10], lty=lty[10])
abline(h=plot.data[[1]][, 2], col=col[1], lty=lty[1])
abline(h=plot.data[[2]][, 2], col=col[2], lty=lty[2])
abline(h=plot.data[[3]][, 2], col=col[3], lty=lty[3])
par(opar)

# ---- lines_metabolomics_alzheimer_res1_briers ----
library(sp)
res <- read.table("results/metabolomics_alzheimer_res1.csv", 
                  stringsAsFactors=FALSE)

pred <- as.matrix(res[substr(rownames(res), 1, 4)=="pred", ])
psel <- as.matrix(res[substr(rownames(res), 1, 4)=="psel", ])
briers <- as.matrix(res[substr(rownames(res), 1, 6)=="briers", ])

methods <- c("ridge", "grridge", "rf", paste0("gren", c(4:6)), 
             paste0("enet", c(1:3)),
             paste0("sgl", c(1:3)), paste0("cmcp", c(1:3)), 
             paste0("gel", c(1:3)), paste0("ocmcp", c(1:3)), 
             paste0("ogel", c(1:3)))
labels <- list("ridge", "GRridge", "random forest", "gren", "enet", "SGL", 
               "cMCP", 
               "gel", "OcMCP", "Ogel")
col <- bpy.colors(length(labels), cutoff.tail=0.1)
lty <- 1:length(labels)

plot.data <- lapply(methods, function(m) {
  aggregate(briers[, substr(colnames(briers), 1, nchar(m))==m], 
            list(psel=psel[, substr(colnames(psel), 1, nchar(m))==m]), mean)})
names(plot.data) <- methods

opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
layout(matrix(c(rep(c(1, 1, 2, 2), 2), rep(c(0, 3, 3, 0), 2)), 
              4, 4, byrow=TRUE))
plot(plot.data[[4]], type="l", xlab="Number of selected features", 
     ylab="Brier skill score", main="(a)", ylim=range(briers), xlim=range(psel), 
     col=col[4], lty=lty[4])
lines(plot.data[[7]], col=col[5], lty=lty[5])
lines(plot.data[[10]], col=col[6], lty=lty[6])
lines(plot.data[[13]], col=col[7], lty=lty[7])
lines(plot.data[[16]], col=col[8], lty=lty[8])
lines(plot.data[[19]], col=col[9], lty=lty[9])
lines(plot.data[[22]], col=col[10], lty=lty[10])
abline(h=plot.data[[1]][, 2], col=col[1], lty=lty[1])
abline(h=plot.data[[2]][, 2], col=col[2], lty=lty[2])
abline(h=plot.data[[3]][, 2], col=col[3], lty=lty[3])
legend("bottomright", legend=labels, col=col, 
       lty=lty, bg = "white")

plot(plot.data[[5]], type="l", xlab="Number of selected features", 
     ylab="Brier skill score", main="(b)", ylim=range(briers), xlim=range(psel), 
     col=col[4], lty=lty[4])
lines(plot.data[[8]], col=col[5], lty=lty[5])
lines(plot.data[[11]], col=col[6], lty=lty[6])
lines(plot.data[[14]], col=col[7], lty=lty[7])
lines(plot.data[[17]], col=col[8], lty=lty[8])
lines(plot.data[[20]], col=col[9], lty=lty[9])
lines(plot.data[[23]], col=col[10], lty=lty[10])
abline(h=plot.data[[1]][, 2], col=col[1], lty=lty[1])
abline(h=plot.data[[2]][, 2], col=col[2], lty=lty[2])
abline(h=plot.data[[3]][, 2], col=col[3], lty=lty[3])

plot(plot.data[[6]], type="l", xlab="Number of selected features", 
     ylab="Brier skill score", main="(c)", ylim=range(briers), xlim=range(psel), 
     col=col[4], lty=lty[4])
lines(plot.data[[9]], col=col[5], lty=lty[5])
lines(plot.data[[12]], col=col[6], lty=lty[6])
lines(plot.data[[15]], col=col[7], lty=lty[7])
lines(plot.data[[18]], col=col[8], lty=lty[8])
lines(plot.data[[21]], col=col[9], lty=lty[9])
lines(plot.data[[24]], col=col[10], lty=lty[10])
abline(h=plot.data[[1]][, 2], col=col[1], lty=lty[1])
abline(h=plot.data[[2]][, 2], col=col[2], lty=lty[2])
abline(h=plot.data[[3]][, 2], col=col[3], lty=lty[3])
par(opar)

# ---- barplot_micrornaseq_cervical_cancer_res1 ----
library(sp)
load("results/micrornaseq_cervical_cancer_fit1.Rdata")

methods <- c("grridge", expression("gren, "~alpha==0.05, "gren, "~alpha==0.5,
                                   "gren, "~alpha==0.95), "no co-data")
col <- grey.colors(length(methods) - 1, start=0.3, end=0.9, gamma=2.2, 
                   alpha=NULL)
lty <- 2
labels <- list(conservation=c("not-conserved", "mammals", "most vertebrates"))

plot.data <- lapply(1:length(fit.gren1$lambdag), function(s) {
  rbind(fit.grridge$lambdamults[[s]], fit.gren1$lambdag[[s]],
        fit.gren2$lambdag[[s]], fit.gren3$lambdag[[s]])})

opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
barplot(plot.data[[1]], beside=TRUE, col=col, legend.text=methods,
        args.legend=list(x="topright", fill=c(col, NA),
                         border=c(rep(1, length(methods) - 1), NA),
                         lty=c(rep(NA, length(methods) - 1), lty[1]),
                         seg.len=1, merge=TRUE, bg="white"),
        names.arg=labels[[1]], ylab=expression(hat(lambda)~"'"[g]),
        main="")
abline(h=1, lty=2)
par(opar)

# ---- lines_micrornaseq_cervical_cancer_res1_auc ----
library(sp)
res <- read.table("results/micrornaseq_cervical_cancer_res1.csv", 
                  stringsAsFactors=FALSE)

pred <- as.matrix(res[substr(rownames(res), 1, 4)=="pred", ])
psel <- as.matrix(res[substr(rownames(res), 1, 4)=="psel", ])
auc <- as.matrix(res[substr(rownames(res), 1, 3)=="auc", ])

methods <- c("ridge", "rf", "grridge", paste0("gren", c(1:3)), 
             paste0("enet", c(1:3)))
col <- bpy.colors(length(methods), cutoff.tail=0.1)
lty <- 1:length(methods)
labels <- expression("ridge", "GRridge", "random forest", "gren,"~alpha==0.05, 
                     "gren,"~alpha==0.5, "gren,"~alpha==0.95, 
                     "enet,"~alpha==0.05, "enet,"~alpha==0.5, 
                     "enet,"~alpha==0.95)

plot.data <- lapply(methods, function(m) {
  aggregate(auc[, grepl(m, colnames(auc))], 
            list(psel=colMeans(psel)[grepl(m, colnames(psel))]), mean)})
names(plot.data) <- methods

opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
plot(plot.data[[4]], type="l", xlab="Number of selected features", ylab="AUC", 
     main="", ylim=range(auc), xlim=c(0, max(psel[, -c(1, 2, ncol(psel))])), 
     col=col[4], lty=lty[4])
lines(plot.data[[5]], col=col[5], lty=lty[5])
lines(plot.data[[6]], col=col[6], lty=lty[6])
lines(plot.data[[7]], col=col[7], lty=lty[7])
lines(plot.data[[8]], col=col[8], lty=lty[8])
lines(plot.data[[9]], col=col[9], lty=lty[9])
abline(h=plot.data[[1]][, 2], col=col[1], lty=lty[1])
abline(h=plot.data[[2]][, 2], col=col[2], lty=lty[2])
abline(h=plot.data[[3]][, 2], col=col[3], lty=lty[3])
legend("bottomright", legend=labels, col=col, 
       lty=lty, bg = "white")
par(opar)

# ---- lines_micrornaseq_cervical_cancer_res1_briers ----
library(sp)
res <- read.table("results/micrornaseq_cervical_cancer_res1.csv", 
                  stringsAsFactors=FALSE)

pred <- as.matrix(res[substr(rownames(res), 1, 4)=="pred", ])
psel <- as.matrix(res[substr(rownames(res), 1, 4)=="psel", ])
briers <- as.matrix(res[substr(rownames(res), 1, 6)=="briers", ])

methods <- c("ridge", "rf", "grridge", paste0("gren", c(1:3)), 
             paste0("enet", c(1:3)))
col <- bpy.colors(length(methods), cutoff.tail=0.1)
lty <- 1:length(methods)
labels <- expression("ridge", "GRridge", "random forest", "gren,"~alpha==0.05, 
                     "gren,"~alpha==0.5, "gren,"~alpha==0.95, 
                     "enet,"~alpha==0.05, "enet,"~alpha==0.5, 
                     "enet,"~alpha==0.95)

plot.data <- lapply(methods, function(m) {
  aggregate(briers[, grepl(m, colnames(briers))], 
            list(psel=colMeans(psel)[grepl(m, colnames(psel))]), mean)})
names(plot.data) <- methods

opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
plot(plot.data[[4]], type="l", xlab="Number of selected features", 
     ylab="Brier skill score", 
     main="", ylim=range(briers), xlim=c(0, max(psel[, -c(1, 2, ncol(psel))])), 
     col=col[4], lty=lty[4])
lines(plot.data[[5]], col=col[5], lty=lty[5])
lines(plot.data[[6]], col=col[6], lty=lty[6])
lines(plot.data[[7]], col=col[7], lty=lty[7])
lines(plot.data[[8]], col=col[8], lty=lty[8])
lines(plot.data[[9]], col=col[9], lty=lty[9])
abline(h=plot.data[[1]][, 2], col=col[1], lty=lty[1])
abline(h=plot.data[[2]][, 2], col=col[2], lty=lty[2])
abline(h=plot.data[[3]][, 2], col=col[3], lty=lty[3])
legend("bottomright", legend=labels, col=col, 
       lty=lty, bg = "white")
par(opar)

################################## response 2 ##################################
# ---- comparison_micrornaseq_colorectal_cancer_res1 ----
library(sp)
load("results/micrornaseq_colorectal_cancer_fit1.Rdata")
res <- read.table("results/micrornaseq_colorectal_cancer_res1.csv", 
                  stringsAsFactors=FALSE)
pred <- as.matrix(res[substr(rownames(res), 1, 4)=="pred", ])
psel <- as.matrix(res[substr(rownames(res), 1, 4)=="psel", ])
auc <- as.matrix(res[substr(rownames(res), 1, 3)=="auc", ])

methods1 <- c("automatic groups", "manual groups", "no co-data")
methods2 <- c(paste0("gren", 1:6))
col1 <- grey.colors(length(methods1) - 1, start=0.3, end=0.9, gamma=2.2, 
                    alpha=NULL)
col2 <- bpy.colors(3, cutoff.tail=0.1)
lty1 <- 2
lty2 <- c(1, 2)
labels1 <- list(diff.expr=c("group 1", "group 2", "group 3", "no FDR"))
labels2 <- expression(alpha==0.05, alpha==0.5, alpha==0.95, "automatic groups",
                      "manual groups")

plot.data1 <- lapply(1:length(fit.gren5$lambdag), function(s) {
  rbind(fit.gren5$lambdag[[s]], fit.gren2$lambdag[[s]])})
plot.data2 <- lapply(methods2, function(m) {
  aggregate(auc[, substr(colnames(auc), 1, nchar(m))==m], 
            list(psel=psel[, substr(colnames(psel), 1, nchar(m))==m]), mean)})
names(plot.data2) <- methods2

opar <- par(no.readonly=TRUE)
par(mar=opar$mar*c(1, 1.3, 1, 1))
layout(matrix(rep(c(1, 1, 2, 2), 2), 2, 4, byrow=TRUE))
barplot(plot.data1[[1]], beside=TRUE, col=col1,
        legend.text=methods1,
        args.legend=list(x="bottomright", fill=c(col1, NA),
                         border=c(rep(1, length(methods1) - 1), NA),
                         lty=c(rep(NA, length(methods1) - 1), lty1[1]),
                         seg.len=1, merge=TRUE, bg="white"),
        names.arg=labels1[[1]], ylab=expression(hat(lambda)~"'"[g]),
        main="(a)")
abline(h=1, lty=2)

plot(plot.data2[[1]], type="l", xlab="Number of selected features", ylab="AUC", 
     main="(b)", ylim=range(sapply(plot.data2, "[", 2)), xlim=c(0, 600), 
     col=col2[1], lty=lty2[2])
lines(plot.data2[[2]], col=col2[2], lty=lty2[2])
lines(plot.data2[[3]], col=col2[3], lty=lty2[2])
lines(plot.data2[[4]], col=col2[1], lty=lty2[1])
lines(plot.data2[[5]], col=col2[2], lty=lty2[1])
lines(plot.data2[[6]], col=col2[3], lty=lty2[1])
legend("bottomright", legend=labels2, fill=c(col2, NA, NA), 
       border=c(rep(1, 3), NA, NA), lty=c(rep(NA, 3), lty2),
       seg.len=1, merge=TRUE, bg="white")
par(opar)