### paths
path.code <- "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/code/"
path.graph <- "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/graphs/"
path.res <- "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/results/"
path.data <- "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/data/"

### libraries
library(pROC)
library(gren)
library(GRridge)
library(sp)
# library(myglmnet)
# install.packages("/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/rpackages/glmnet_2.0-13.tar.gz", 
#                  repos=NULL)

### load data
load(paste(path.data, "mirsData.RData", sep=""))

### data manipulation
x <- apply(t(as.matrix(mirsData$transformedData)), 2,function(v) {
  (v - mean(v))/sd(v)})[, order(as.numeric(mirsData$conservation))]
y <- as.numeric(mirsData$response) - 1 # 1=CIN3, 0=normal
n <- nrow(x)
p <- ncol(x)
m <- rep(1, n)

### partitioning
conservation <- sort(as.numeric(mirsData$conservation))
parCons <- CreatePartition(as.factor(conservation))

# ## fit full models
# set.seed(2018)
# fit.gren1 <- gren(x, y, m, NULL, list(conservation=conservation), 0.05)
# fit.gren2 <- gren(x, y, m, NULL, list(conservation=conservation), 0.5)
# fit.gren3 <- gren(x, y, m, NULL, list(conservation=conservation), 0.95)
# fit.grridge <- grridge(t(x), y, parCons)
# fit.enet1 <- cv.glmnet(x, y, family="binomial", alpha=0.05, standardize=FALSE)
# fit.enet2 <- cv.glmnet(x, y, family="binomial", alpha=0.5, standardize=FALSE)
# fit.enet3 <- cv.glmnet(x, y, family="binomial", alpha=0.95, standardize=FALSE)
# fit.ridge <- cv.glmnet(x, y, family="binomial", alpha=0, standardize=FALSE)

# save(fit.ridge, fit.grridge, fit.enet1, fit.enet2, fit.enet3, fit.gren1, 
#      fit.gren2, fit.gren3,
#      file=paste(path.res, "gren_mirna_prostate_fitted1.Rdata", sep=""))

# ### cross-validate predictions
# set.seed(2018)
# nfolds <- n
# rest <- n %% nfolds
# foldid <- sample(rep(1:nfolds, times=round(c(rep(
#   n %/% nfolds + as.numeric(rest!=0), times=rest),
#   rep(n %/% nfolds, times=nfolds - rest)))))
# 
# psel <- c(seq(1, 5, 1), seq(7, 15, 2), seq(18, 30, 3), seq(34, 50, 4))
# methods <- c("ridge", "grridge", "enet1", "enet2", "enet3", "gren1",
#              "gren2", "gren3")
# pred.ridge <- numeric(n)
# for(s in 2:length(methods)) {
#   assign(paste("pred.", methods[s], sep=""), matrix(NA, nrow=n,
#                                                     ncol=length(psel)))
#   assign(paste("psel.", methods[s], sep=""), matrix(NA, nrow=n,
#                                                     ncol=length(psel)))
# }
# 
# for(k in sort(unique(foldid))) {
#   cat(paste("Fold ", k, "\n"))
# 
#   xtrain <- x[foldid!=k, ]
#   xtest <- matrix(x[foldid==k, ], ncol=p, byrow=TRUE)
#   ytrain <- y[foldid!=k]
#   mtrain <- m[foldid!=k]
# 
#   cv.ridge <- cv.glmnet(xtrain, ytrain, family="binomial", alpha=0, 
#                         standardize=FALSE, lambda=fit.ridge$lambda.min)
#   cv.grridge <- vector("list", length(psel))
#   for(s in 1:length(psel)) {
#     cv.grridge[[s]] <- grridge(t(xtrain), ytrain, parCons, 
#                                optl=fit.grridge$optl, selection=TRUE,
#                                maxsel=psel[s])
#   }
#   cv.gren1 <- gren(xtrain, ytrain, mtrain, NULL, 
#                    list(conservation=conservation), alpha=0.05, psel=psel,
#                    compare=TRUE, lambda=fit.enet1$lambda.min)
#   cv.gren2 <- gren(xtrain, ytrain, mtrain, NULL, 
#                    list(conservation=conservation), alpha=0.5, psel=psel,
#                    compare=TRUE, lambda=fit.enet2$lambda.min)
#   cv.gren3 <- gren(xtrain, ytrain, mtrain, NULL, 
#                    list(conservation=conservation), alpha=0.95, psel=psel,
#                    compare=TRUE, lambda=fit.enet3$lambda.min)
#   
#   pred.ridge[foldid==k] <- predict(cv.ridge, xtest, s="lambda.min", 
#                                    type="response")
#   pred.grridge[foldid==k, ] <- sapply(cv.grridge, function(s) {
#     predict.grridge(s, t(xtest))[, 3]})
#   
#   pred.enet1[foldid==k, ] <- predict(cv.gren1$freq.model$regular, xtest, 
#                                      type="response")
#   pred.enet2[foldid==k, ] <- predict(cv.gren2$freq.model$regular, xtest, 
#                                      type="response")
#   pred.enet3[foldid==k, ] <- predict(cv.gren3$freq.model$regular, xtest, 
#                                      type="response")
#   pred.gren1[foldid==k, ] <- predict(cv.gren1$freq.model$groupreg, xtest, 
#                                      type="response")
#   pred.gren2[foldid==k, ] <- predict(cv.gren2$freq.model$groupreg, xtest, 
#                                      type="response")
#   pred.gren3[foldid==k, ] <- predict(cv.gren3$freq.model$groupreg, xtest, 
#                                      type="response")
#   
#   psel.grridge[foldid==k, ] <- sapply(cv.grridge, function(s) {
#     length(s$resEN$whichEN)})
#   psel.enet1[foldid==k, ] <- cv.gren1$freq.model$regular$df
#   psel.enet2[foldid==k, ] <- cv.gren2$freq.model$regular$df
#   psel.enet3[foldid==k, ] <- cv.gren3$freq.model$regular$df
#   psel.gren1[foldid==k, ] <- cv.gren1$freq.model$groupreg$df
#   psel.gren2[foldid==k, ] <- cv.gren2$freq.model$groupreg$df
#   psel.gren3[foldid==k, ] <- cv.gren3$freq.model$groupreg$df
# 
# }
# 
# results1 <- list(pred=list(pred.ridge, pred.grridge, pred.enet1, pred.enet2,
#                            pred.enet3, pred.gren1, pred.gren2, pred.gren3),
#                  psel=list(psel.grridge, psel.enet1, psel.enet2, psel.enet3,
#                            psel.gren1, psel.gren2, psel.gren3))
# save(results1, file=paste(path.res, "gren_mirna_prostate_res1.Rdata", sep=""))

### stability selection
set.seed(2018)
K <- 50
psel <- 15

methods <- c("grridge", "enet1", "enet2", "enet3", "gren1",
             "gren2", "gren3")
for(s in 1:length(methods)) {
  assign(paste("sel.", methods[s], sep=""), vector("list", K))
}

for(k in 1:K) {
  cat(paste("Sample ", k, "\n"))
  
  id <- sample(c(1:n), n, replace=TRUE)
  xtrain <- x[id, ]
  id.const <- which(apply(xtrain, 2, sd)==0)
  if(length(id.const)!=0) {
    xtrain <- xtrain[, -id.const]
  }
  ytrain <- y[id]
  mtrain <- m[id]
  
  boot.conservation <- sort(as.numeric(mirsData$conservation))[-id.const]
  boot.parCons <- CreatePartition(as.factor(boot.conservation))
  
  boot.grridge <- grridge(t(xtrain), ytrain, list(conservation=boot.parCons), 
                          selection=TRUE, maxsel=psel)
  boot.gren1 <- gren(xtrain, ytrain, mtrain, NULL, 
                     list(conservation=boot.conservation), alpha=0.05, 
                     psel=psel, compare=TRUE)
  boot.gren2 <- gren(xtrain, ytrain, mtrain, NULL, 
                     list(conservation=boot.conservation), alpha=0.5, psel=psel,
                     compare=TRUE)
  boot.gren3 <- gren(xtrain, ytrain, mtrain, NULL, 
                     list(conservation=boot.conservation), alpha=0.95, psel=psel,
                     compare=TRUE)
  
  sel.grridge[[k]] <- names(boot.grridge$resEN$whichEN)
  sel.gren1[[k]] <- rownames(coef(boot.gren1$freq.model$groupreg))[
    as.numeric(coef(boot.gren1$freq.model$groupreg))!=0][-1]
  sel.gren2[[k]] <- rownames(coef(boot.gren2$freq.model$groupreg))[
    as.numeric(coef(boot.gren2$freq.model$groupreg))!=0][-1]
  sel.gren3[[k]] <- rownames(coef(boot.gren3$freq.model$groupreg))[
    as.numeric(coef(boot.gren3$freq.model$groupreg))!=0][-1]
  sel.enet1[[k]] <- rownames(coef(boot.gren1$freq.model$regular))[
    as.numeric(coef(boot.gren1$freq.model$regular))!=0][-1]
  sel.enet2[[k]] <- rownames(coef(boot.gren2$freq.model$regular))[
    as.numeric(coef(boot.gren2$freq.model$regular))!=0][-1]
  sel.enet3[[k]] <- rownames(coef(boot.gren3$freq.model$regular))[
    as.numeric(coef(boot.gren3$freq.model$regular))!=0][-1]
  
}

results2 <- list(sel=list(sel.grridge, sel.enet1, sel.enet2, sel.enet3,
                          sel.gren1, sel.gren2, sel.gren3))
save(results2, file=paste(path.res, "gren_mirna_prostate_res2.Rdata", sep=""))

### penalty multiplier plot
load(paste(path.res, "gren_mirna_prostate_fitted1.Rdata", sep=""))
leg.lab <- c("GRridge", expression(paste("gren, ", alpha==0.05)),
             expression(paste("gren, ", alpha==0.5)),
             expression(paste("gren, ", alpha==0.95)),
             "not group-regularized")
png(paste(path.graph, "gren_mirna_prostate_bar1.png", sep=""), units="in",
    width=7, height=5, res=200)
par(mar=c(5.1, 5.6, 4.1, 2.1))
barplot(rbind(fit.grridge$lambdamults$group,
              fit.gren1$lambdag$conservation,
              fit.gren2$lambdag$conservation,
              fit.gren3$lambdag$conservation), beside=TRUE,
        names.arg=rep("", 3), col=colors, 
        ylab=expression({lambda^{"'"}}[g]), xlab="", 
        main="", legend.text=leg.lab, cex.axis=1.5, cex.names=1.5, cex.lab=2,
        args.legend=list(x="topright", fill=c(colors, 0),
                         lty=c(rep(NA, 4), 2), lwd=c(rep(NA, 4), 1.5),
                         border=c(rep(1, 4), 0), merge=TRUE, seg.len=1, 
                         cex=1.3), cex.main=2)
axis(1, c(3, 8, 13), c("Not conserved", "\n Conserved \n in mammals", 
                       "Broadly \n conserved"), pos=-0.15, tck=0, lty=0, 
     cex.axis=1.5)
abline(h=1, lty=2, lwd=1.5)
dev.off()

### performance plots
leglabels <- c("ridge", expression(paste("enet, ", alpha==0.05)),
               expression(paste("enet, ", alpha==0.5)),
               expression(paste("enet, ", alpha==0.95)),
               "group-regularized", "not group-regularized")
colors <- bpy.colors(6)[-c(1, 6)]

load(paste(path.res, "gren_mirna_prostate_res1.Rdata", sep=""))

psel1 <- lapply(results1$psel, function(l) {colMeans(l)})

auc1 <- lapply(results1$pred[-1], function(l) {
  apply(l, 2, function(preds) {pROC::roc(y, preds, direction="<")$auc})})
auc1 <- lapply(1:length(psel1), function(l) {auc1[[l]][psel1[[l]] >= 2]})
auc1 <- c(pROC::roc(y, results1$pred[[1]], direction="<")$auc, auc1)

briers1 <- lapply(results1$pred[-1], function(l) {
  apply(l, 2, function(preds) {1 - sum((y - preds)^2)/sum((y - mean(y))^2)})})
briers1 <- lapply(1:length(psel1), function(l) {briers1[[l]][psel1[[l]] >= 2]})
briers1 <- c(1 - sum((y - results1$pred[[1]])^2)/sum((y - mean(y))^2), briers1)

psel1 <- lapply(psel1, function(l) {l[l >= 2]})

png(paste(path.graph, "gren_mirna_prostate_performance.png", sep=""),
    units="in", width=14, height=6, res=120)
par(mfrow=c(1, 2), mar=c(5.1, 5.1, 4.1, 2.1))
plot(sort(psel1[[1]]), auc1[[2]][order(psel1[[1]])], type="l", 
     xlim=range(psel1), ylim=range(auc1), col=colors[1], 
     xlab="Number of selected features", ylab="AUC", main="a)", 
     lwd=1.5, cex.axis=1.5, cex.lab=2, cex.main=2)
lines(range(psel1), rep(auc1[[1]], 2), col=colors[1], lty=2, lwd=1.5)
lines(psel1[[2]], auc1[[3]], col=colors[2], lty=2, lwd=1.5)
lines(psel1[[3]], auc1[[4]], col=colors[3], lty=2, lwd=1.5)
lines(psel1[[4]], auc1[[5]], col=colors[4], lty=2, lwd=1.5)
lines(psel1[[5]], auc1[[6]], col=colors[2], lwd=1.5)
lines(psel1[[6]], auc1[[7]], col=colors[3], lwd=1.5)
lines(psel1[[7]], auc1[[8]], col=colors[4], lwd=1.5)
legend("bottomright", legend=leglabels, fill=c(colors, 0, 0),
       lty=c(rep(NA, 4), 1, 2), lwd=c(rep(NA, 4), 1.5, 1.5), 
       border=c(rep(1, 4), 0, 0), merge=TRUE, seg.len=1, cex=1.3)

plot(sort(psel1[[1]]), briers1[[2]][order(psel1[[1]])], type="l", 
     xlim=range(psel1), ylim=range(briers1), col=colors[1], 
     xlab="Number of selected features", ylab="Brier skill score", main="b)", 
     cex.axis=1.5, cex.lab=2, lwd=1.5, cex.main=2)
lines(range(psel1), rep(briers1[[1]], 2), col=colors[1], lty=2, lwd=1.5)
lines(psel1[[2]], briers1[[3]], col=colors[2], lty=2, lwd=1.5)
lines(psel1[[3]], briers1[[4]], col=colors[3], lty=2, lwd=1.5)
lines(psel1[[4]], briers1[[5]], col=colors[4], lty=2, lwd=1.5)
lines(psel1[[5]], briers1[[6]], col=colors[2], lwd=1.5)
lines(psel1[[6]], briers1[[7]], col=colors[3], lwd=1.5)
lines(psel1[[7]], briers1[[8]], col=colors[4], lwd=1.5)
dev.off()


# plaatjes bayescomp 2018
leglabels <- c(expression(paste("enet, ", alpha==0.05)),
               expression(paste("enet, ", alpha==0.5)),
               expression(paste("enet, ", alpha==0.95)),
               "group-regularized", "not group-regularized")

pdf(paste(path.graph, "gren_mirna_prostate_performance.pdf", sep=""))
par(mar=c(5.1, 5.1, 4.1, 2.1))
plot(psel1[[2]], auc1[[3]], type="l", 
     xlim=range(psel1), ylim=range(auc1), col=colors[2], 
     xlab="Number of selected features", ylab="AUC", main="", 
     lwd=2.5, cex.axis=1.5, cex.lab=2, cex.main=2, lty=2)
lines(psel1[[3]], auc1[[4]], col=colors[3], lty=2, lwd=2.5)
lines(psel1[[4]], auc1[[5]], col=colors[4], lty=2, lwd=2.5)
lines(psel1[[5]], auc1[[6]], col=colors[2], lwd=2.5)
lines(psel1[[6]], auc1[[7]], col=colors[3], lwd=2.5)
lines(psel1[[7]], auc1[[8]], col=colors[4], lwd=2.5)
legend("bottomright", legend=leglabels, fill=c(colors[c(2:4)], 0, 0),
       lty=c(rep(NA, 3), 1, 2), lwd=c(rep(NA, 3), 1.5, 1.5), 
       border=c(rep(1, 3), 0, 0), merge=TRUE, seg.len=1, cex=1.3)
dev.off()

load(paste(path.res, "gren_mirna_prostate_fitted1.Rdata", sep=""))
leg.lab <- c(expression(paste("gren, ", alpha==0.05)),
             expression(paste("gren, ", alpha==0.5)),
             expression(paste("gren, ", alpha==0.95)),
             "not group-regularized")
pdf(paste(path.graph, "gren_mirna_prostate_bar1.pdf", sep=""))
par(mar=c(5.1, 5.6, 4.1, 2.1))
barplot(rbind(fit.gren1$lambdag$conservation,
              fit.gren2$lambdag$conservation,
              fit.gren3$lambdag$conservation), beside=TRUE,
        names.arg=rep("", 3), col=colors[c(2:4)], 
        ylab=expression({lambda^{"'"}}[g]), xlab="", 
        main="", legend.text=leg.lab, cex.axis=1.5, cex.names=1.5, cex.lab=2,
        args.legend=list(x="topright", fill=c(colors[c(2:4)], 0),
                         lty=c(rep(NA, 3), 2), lwd=c(rep(NA, 3), 1.5),
                         border=c(rep(1, 3), 0), merge=TRUE, seg.len=1, 
                         cex=1.3), cex.main=2)
axis(1, c(2.5, 6.5, 10.5), c("Not conserved", "\n Conserved \n in mammals", 
                       "Broadly \n conserved"), pos=-0.15, tck=0, lty=0, 
     cex.axis=1.5)
abline(h=1, lty=2, lwd=2.5)
dev.off()

### histogram with bootstrap variable selection
colors <- bpy.colors(6)[-c(1, 6)]

leg.lab <- c("not group-regularized", "group-regularized")

int2 <- lapply(results2$sel, function(l) {
  sapply(combn(l, 2, simplify=FALSE), function(x) {
    length(intersect(x[[1]], x[[2]]))}, simplify=TRUE)})

png(paste(path.graph, "gren_mirna_prostate_overlap.png", sep=""),
    units="in", width=12, height=12, res=120)
par(mfrow=c(2, 2), mar=c(5.1, 5.1, 4.1, 2.1))
barplot(rbind(rep(0, 15),
              sapply(c(1:max(unlist(int2))), function(s) {sum(int2[[1]]==s)}))/
          choose(50, 2), xlab="Number of overlapping features", 
        ylab="Density", beside=TRUE, col=colors[c(2:3)], names.arg=c(1:15), 
        main="a)", cex.axis=1.5, cex.names=1.5, cex.lab=2, cex.main=2,
        legend.text=leg.lab, 
        args.legend=list(x="topright", fill=colors[c(2:3)], 
                         border=colors[c(2:3)], cex=1.3))
barplot(rbind(sapply(c(1:max(unlist(int2))), function(s) {sum(int2[[2]]==s)}),
              sapply(c(1:max(unlist(int2))), function(s) {sum(int2[[5]]==s)}))/
          choose(50, 2), 
        xlab="Number of overlapping features", ylab="Density", beside=TRUE, 
        col=colors[c(2:3)], names.arg=c(1:15), main="b)", cex.axis=1.5, 
        cex.names=1.5, cex.lab=2, cex.main=2)
barplot(rbind(sapply(c(1:max(unlist(int2))), function(s) {sum(int2[[3]]==s)}),
              sapply(c(1:max(unlist(int2))), function(s) {sum(int2[[6]]==s)}))/
          choose(50, 2), 
        xlab="Number of overlapping features", ylab="Density", beside=TRUE, 
        col=colors[c(2:3)], names.arg=c(1:15), main="c)", cex.axis=1.5, 
        cex.names=1.5, cex.lab=2, cex.main=2)
barplot(rbind(sapply(c(1:max(unlist(int2))), function(s) {sum(int2[[4]]==s)}),
              sapply(c(1:max(unlist(int2))), function(s) {sum(int2[[7]]==s)}))/
          choose(50, 2), 
        xlab="Number of overlapping features", ylab="Density", beside=TRUE, 
        col=colors[c(2:3)], names.arg=c(1:15), main="d)", cex.axis=1.5, 
        cex.names=1.5, cex.lab=2, cex.main=2)
dev.off()

# plaatjes short course milaan 2018
leglabels <- c("ridge", expression(paste("enet, ", alpha==0.05)),
               expression(paste("enet, ", alpha==0.5)),
               expression(paste("enet, ", alpha==0.95)),
               "group-regularized", "not group-regularized")
colors <- bpy.colors(6)[-c(1, 6)]

pdf(paste(path.graph, "gren_mirna_prostate_performance2.pdf", sep=""))
par(mar=c(5.1, 5.1, 4.1, 2.1))
plot(psel1[[2]], auc1[[3]], type="l", 
     xlim=range(psel1), ylim=range(auc1), col=colors[2], 
     xlab="Number of selected features", ylab="AUC", main="", 
     lwd=2.5, cex.axis=1.5, cex.lab=2, cex.main=2, lty=2)
lines(sort(psel1[[1]]), auc1[[2]][order(psel1[[1]])], col=colors[1], lwd=2.5)
lines(range(psel1), rep(auc1[[1]], 2), col=colors[1], lty=2, lwd=2.5)
lines(psel1[[3]], auc1[[4]], col=colors[3], lty=2, lwd=2.5)
lines(psel1[[4]], auc1[[5]], col=colors[4], lty=2, lwd=2.5)
lines(psel1[[5]], auc1[[6]], col=colors[2], lwd=2.5)
lines(psel1[[6]], auc1[[7]], col=colors[3], lwd=2.5)
lines(psel1[[7]], auc1[[8]], col=colors[4], lwd=2.5)
legend("bottomright", legend=leglabels, fill=c(colors[c(1:4)], 0, 0),
       lty=c(rep(NA, 4), 1, 2), lwd=c(rep(NA, 4), 1.5, 1.5), 
       border=c(rep(1, 4), 0, 0), merge=TRUE, seg.len=1, cex=1.3)
dev.off()

load(paste(path.res, "gren_mirna_prostate_fitted1.Rdata", sep=""))
leg.lab <- c("GRridge", expression(paste("gren, ", alpha==0.05)),
             expression(paste("gren, ", alpha==0.5)),
             expression(paste("gren, ", alpha==0.95)),
             "not group-regularized")
pdf(paste(path.graph, "gren_mirna_prostate_bar2.pdf", sep=""))
par(mar=c(5.1, 5.6, 4.1, 2.1))
barplot(rbind(fit.grridge$lambdamults$group, 
              fit.gren1$lambdag$conservation,
              fit.gren2$lambdag$conservation,
              fit.gren3$lambdag$conservation), beside=TRUE,
        names.arg=rep("", 3), col=colors[c(1:4)], 
        ylab=expression({lambda^{"'"}}[g]), xlab="", 
        main="", legend.text=leg.lab, cex.axis=1.5, cex.names=1.5, cex.lab=2,
        args.legend=list(x="topright", fill=c(colors[c(1:4)], 0),
                         lty=c(rep(NA, 4), 2), lwd=c(rep(NA, 4), 1.5),
                         border=c(rep(1, 4), 0), merge=TRUE, seg.len=1, 
                         cex=1.3), cex.main=2)
axis(1, c(2.5, 6.5, 10.5), c("Not conserved", "\n Conserved \n in mammals", 
                             "Broadly \n conserved"), pos=-0.15, tck=0, lty=0, 
     cex.axis=1.5)
abline(h=1, lty=2, lwd=2.5)
dev.off()





# plaatjes poster warschau 2018
library(sp)
library(pROC)

path.data <- "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/data/"
path.graph <- "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/graphs/"
path.res <- "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/results/"

load(paste(path.res, "gren_mirna_prostate_res1.Rdata", sep=""))
load(paste(path.res, "gren_mirna_prostate_fitted1.Rdata", sep=""))
load(paste(path.data, "mirsData.RData", sep=""))

x <- apply(t(as.matrix(mirsData$transformedData)), 2,function(v) {
  (v - mean(v))/sd(v)})[, order(as.numeric(mirsData$conservation))]
y <- as.numeric(mirsData$response) - 1 # 1=CIN3, 0=normal
n <- nrow(x)
p <- ncol(x)
m <- rep(1, n)

psel1 <- lapply(results1$psel, function(l) {colMeans(l)})

auc1 <- lapply(results1$pred[-1], function(l) {
  apply(l, 2, function(preds) {pROC::roc(y, preds, direction="<")$auc})})
auc1 <- lapply(1:length(psel1), function(l) {auc1[[l]][psel1[[l]] >= 2]})
auc1 <- c(pROC::roc(y, results1$pred[[1]], direction="<")$auc, auc1)

briers1 <- lapply(results1$pred[-1], function(l) {
  apply(l, 2, function(preds) {1 - sum((y - preds)^2)/sum((y - mean(y))^2)})})
briers1 <- lapply(1:length(psel1), function(l) {briers1[[l]][psel1[[l]] >= 2]})
briers1 <- c(1 - sum((y - results1$pred[[1]])^2)/sum((y - mean(y))^2), briers1)

psel1 <- lapply(psel1, function(l) {l[l >= 2]})

leglabels <- c("ridge", expression(paste("enet, ", alpha==0.05)),
               expression(paste("enet, ", alpha==0.5)),
               expression(paste("enet, ", alpha==0.95)),
               "group-regularized", "not group-regularized")
colors <- bpy.colors(6)[-c(1, 6)]

pdf(paste(path.graph, "gren_mirna_cervical_performance2.pdf", sep=""))
par(mar=c(5.1, 5.1, 4.1, 2.1))
plot(psel1[[2]], auc1[[3]], type="l", 
     xlim=range(psel1), ylim=range(auc1), col=colors[2], 
     xlab="Number of selected features", ylab="AUC", main="", 
     lwd=2.5, cex.axis=1.5, cex.lab=2, cex.main=2, lty=2)
lines(sort(psel1[[1]]), auc1[[2]][order(psel1[[1]])], col=colors[1], lwd=2.5)
lines(range(psel1), rep(auc1[[1]], 2), col=colors[1], lty=2, lwd=2.5)
lines(psel1[[3]], auc1[[4]], col=colors[3], lty=2, lwd=2.5)
lines(psel1[[4]], auc1[[5]], col=colors[4], lty=2, lwd=2.5)
lines(psel1[[5]], auc1[[6]], col=colors[2], lwd=2.5)
lines(psel1[[6]], auc1[[7]], col=colors[3], lwd=2.5)
lines(psel1[[7]], auc1[[8]], col=colors[4], lwd=2.5)
legend("bottomright", legend=leglabels, fill=c(colors[c(1:4)], 0, 0),
       lty=c(rep(NA, 4), 1, 2), lwd=c(rep(NA, 4), 1.5, 1.5), 
       border=c(rep(1, 4), 0, 0), merge=TRUE, seg.len=1, cex=1.3)
dev.off()

leg.lab <- c("GRridge", expression(paste("gren, ", alpha==0.05)),
             expression(paste("gren, ", alpha==0.5)),
             expression(paste("gren, ", alpha==0.95)),
             "not group-regularized")
pdf(paste(path.graph, "gren_mirna_cervical_bar2.pdf", sep=""))
par(mar=c(5.1, 5.6, 4.1, 2.1))
barplot(rbind(fit.grridge$lambdamults$group, 
              fit.gren1$lambdag$conservation,
              fit.gren2$lambdag$conservation,
              fit.gren3$lambdag$conservation), beside=TRUE,
        names.arg=rep("", 3), col=colors[c(1:4)], 
        ylab=expression({lambda^{"'"}}[g]), xlab="", 
        main="", legend.text=leg.lab, cex.axis=1.5, cex.names=1.5, cex.lab=2,
        args.legend=list(x="topright", fill=c(colors[c(1:4)], 0),
                         lty=c(rep(NA, 4), 2), lwd=c(rep(NA, 4), 1.5),
                         border=c(rep(1, 4), 0), merge=TRUE, seg.len=1, 
                         cex=1.3), cex.main=2)
axis(1, c(3, 8, 13), c("Not conserved", "Conserved \n in mammals", 
                             "Broadly \n conserved"), pos=-0.15, tck=0, lty=0, 
     cex.axis=1.5)
abline(h=1, lty=2, lwd=2.5)
dev.off()
