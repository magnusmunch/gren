##############################  preamble  #############################
# testing grEBEN                                                      #
# version: 01                                                         #
# author: Magnus M?nch                                                #
# created: 27-10-2017                                                 #
# last edited: 27-10-2017                                             #
#######################################################################

###############################  notes  ###############################
#######################################################################

# paths
path.code <- as.character(ifelse(Sys.info()[1]=="Darwin", 
                                 "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/code/" ,
                                 "~/EBEN/code/"))
path.graph <- "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/graphs/"
path.res <- as.character(ifelse(Sys.info()[1]=="Darwin", 
                                "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/results/" ,
                                "~/EBEN/results/"))

### libraries
library(mvtnorm)
library(GRridge)
library(pROC)

### source grEBEN
source(paste(path.code, "grVBEM.R", sep=""))
source(paste(path.code, "mygrridge.R", sep=""))

## simulations
# ## simulation 1
# # create data
# n <- 100
# p <- 1000
# G <- 10
# rhow <- 0.7
# rhob <- 0.5
# blockw <- blockb <- diag(p/G)
# for(i in 1:ncol(blockw)) {
#   for(j in 1:ncol(blockw)) {
#     blockw[i, j] <- rhow^abs(i - j)
#     blockb[i, j] <- rhob^(abs(i - j) + 1)
#   }
# }
# Sigma <- rbind(cbind(blockw, blockb), cbind(blockb, blockw))
# q <- 0.7 # proportion of non-zero coefficients
# beta.mean <- 0.05
# f <- 1.4
# m <- rep(1, n)
# betag <- beta.mean*G*f^(c(1:G) - 1)/(q*sum(f^(c(1:G) - 1)))
# beta <- as.numeric(sapply(betag, function(b) {
#   rep(c(0, b), times=c((1 - q)*p/G, q*p/G))}))
# part.greben <- list(groups=rep(1:G, each=p/G))
# part.grridge <- list(groups=CreatePartition(as.factor(part.greben$groups)))
# ntest <- 1000
# mtest <- rep(1, ntest)
# 
# methods <- c("ridge", "GRridge", "GRridge+sel",
#              "enet+a=0.05", "enet+a=0.5", "enet+a=0.95",
#              "grEBEN+a=0.05", "grEBEN+a=0.5","grEBEN+a=0.95")
# nreps <- 50
# auc1 <- briers1 <- mse1 <- rep(list(vector(mode="list", length=nreps)),
#                                length(methods))
# names(auc1) <- names(briers1) <- names(mse1) <- methods
# lambdag1 <- rep(list(matrix(NA, ncol=G, nrow=nreps)), 4)
# names(lambdag1) <- methods[c(2, 7, 8, 9)]
# # the simulations
# for(r in 1:nreps) {
# 
#   set.seed(400 + r)
#   x <- do.call(cbind, replicate(G/2, rmvnorm(n, rep(0, 2*p/G), Sigma),
#                                 simplify=FALSE))
#   prob <- as.numeric(exp(x %*% beta)/(1 + exp(x %*% beta)))
#   y <- rbinom(n, m, prob)
# 
#   xtest <- do.call(cbind, replicate(G/2, rmvnorm(ntest, rep(0, 2*p/G), Sigma),
#                                     simplify=FALSE))
#   probtest <- as.numeric(exp(xtest %*% beta)/(1 + exp(xtest %*% beta)))
#   ytest <- rbinom(ntest, mtest, probtest)
#   
#   fit1.greben <- grEBEN3(x, y, m, partitions=part.greben, alpha=0.05, psel=TRUE)
#   fit2.greben <- grEBEN3(x, y, m, partitions=part.greben, alpha=0.5, psel=TRUE)
#   fit3.greben <- grEBEN3(x, y, m, partitions=part.greben, alpha=0.95, psel=TRUE)
# 
#   # number of selected variables
#   psel1.enet <- apply(fit1.greben$beta.nogroups, 2, function(b) {sum(b!=0) - 1})
#   psel2.enet <- apply(fit2.greben$beta.nogroups, 2, function(b) {sum(b!=0) - 1})
#   psel3.enet <- apply(fit3.greben$beta.nogroups, 2, function(b) {sum(b!=0) - 1})
# 
#   psel1.greben <- apply(fit1.greben$beta, 2, function(b) {sum(b!=0) - 1})
#   psel2.greben <- apply(fit2.greben$beta, 2, function(b) {sum(b!=0) - 1})
#   psel3.greben <- apply(fit3.greben$beta, 2, function(b) {sum(b!=0) - 1})
# 
#   psel.all <- unique(c(psel1.enet, psel2.enet, psel3.enet,
#     psel1.greben, psel2.greben, psel3.greben))
#   psel.in <- floor(quantile(psel.all[psel.all!=0],
#     prob=c(seq(0.01, 0.05, 0.01), seq(0.07, 0.25, 0.03), seq(0.3, 1, 0.1))))
# 
#   fit1.grridge <- vector(mode="list", length=length(psel.in))
#   fit1.grridge[[1]] <- grridge(t(x), y, part.grridge, selectionEN=TRUE,
#                                maxsel=psel.in[1])
#   for(s in 2:length(psel.in)) {
#     fit1.grridge[[s]] <- grridge(t(x), y, part.grridge, selectionEN=TRUE,
#                                  optl=fit1.grridge[[1]]$optl,
#                                  maxsel=psel.in[s])
#   }
# 
#   psel1.grridge <- sapply(1:length(fit1.grridge), function(s) {
#     return(length(fit1.grridge[[s]]$resEN$whichEN))})
# 
#   # estimates
#   est1.ridge <- coef(fit1.grridge[[1]]$predobj$NoGroups, "all")
# 
#   est1.grridge <- coef(fit1.grridge[[1]]$predobj$GroupRegul, "all")
# 
#   est2.grridge <- sapply(fit1.grridge, function(s) {
#     replace(rep(0, p + 1), c(1, s$resEN$whichEN + 1), 
#             coef(s$predobj$EN, "all"))})
#   
#   est1.enet <- fit1.greben$beta.nogroups
#   est2.enet <- fit2.greben$beta.nogroups
#   est3.enet <- fit3.greben$beta.nogroups
# 
#   est1.greben <- fit1.greben$beta
#   est2.greben <- fit2.greben$beta
#   est3.greben <- fit3.greben$beta
# 
#   # predictions on fit data
#   pred1.ridge <- predict.grridge(fit1.grridge[[1]], t(xtest))[, 1]
# 
#   pred1.grridge <- predict.grridge(fit1.grridge[[1]], t(xtest))[, 2]
#   pred2.grridge <- sapply(fit1.grridge, function(s) {
#     predict.grridge(s, t(xtest))[, 3]})
# 
#   pred1.enet <- 1/(1 + exp(-xtest %*% est1.enet[-1, ]))
#   pred2.enet <- 1/(1 + exp(-xtest %*% est2.enet[-1, ]))
#   pred3.enet <- 1/(1 + exp(-xtest %*% est3.enet[-1, ]))
# 
#   pred1.greben <- 1/(1 + exp(-xtest %*% est1.greben[-1, ]))
#   pred2.greben <- 1/(1 + exp(-xtest %*% est2.greben[-1, ]))
#   pred3.greben <- 1/(1 + exp(-xtest %*% est3.greben[-1, ]))
# 
#    # AUCs
#   auc.true <- pROC::roc(ytest, probtest)$auc
# 
#   auc1.ridge <- pROC::roc(ytest, pred1.ridge)$auc
# 
#   auc1.grridge <- pROC::roc(ytest, pred1.grridge)$auc
#   auc2.grridge <- apply(pred2.grridge, 2, function(r) {pROC::roc(ytest, r)$auc})
# 
#   auc1.enet <- apply(pred1.enet, 2, function(r) {pROC::roc(ytest, r)$auc})
#   auc2.enet <- apply(pred2.enet, 2, function(r) {pROC::roc(ytest, r)$auc})
#   auc3.enet <- apply(pred3.enet, 2, function(r) {pROC::roc(ytest, r)$auc})
# 
#   auc1.greben <- apply(pred1.greben, 2, function(r) {pROC::roc(ytest, r)$auc})
#   auc2.greben <- apply(pred2.greben, 2, function(r) {pROC::roc(ytest, r)$auc})
#   auc3.greben <- apply(pred3.greben, 2, function(r) {pROC::roc(ytest, r)$auc})
# 
#   # Brier scores
#   brier.null <- sum((ytest - mean(ytest))^2)
#   briers.true <- 1 - sum((ytest - probtest)^2)/brier.null
# 
#   briers1.ridge <- 1 - sum((ytest - pred1.ridge)^2)/brier.null
# 
#   briers1.grridge <- 1 - sum((ytest - pred1.grridge)^2)/brier.null
#   briers2.grridge <- apply(pred2.grridge, 2, function(pred) {
#     1 - sum((ytest - pred)^2)/brier.null})
# 
#   briers1.enet <- apply(pred1.enet, 2, function(pred) {
#     1 - sum((ytest - pred)^2)/brier.null})
#   briers2.enet <- apply(pred2.enet, 2, function(pred) {
#     1 - sum((ytest - pred)^2)/brier.null})
#   briers3.enet <- apply(pred3.enet, 2, function(pred) {
#     1 - sum((ytest - pred)^2)/brier.null})
#   
#   briers1.greben <- apply(pred1.greben, 2, function(pred) {
#     1 - sum((ytest - pred)^2)/brier.null})
#   briers2.greben <- apply(pred2.greben, 2, function(pred) {
#     1 - sum((ytest - pred)^2)/brier.null})
#   briers3.greben <- apply(pred3.greben, 2, function(pred) {
#     1 - sum((ytest - pred)^2)/brier.null})
#   
#   # MSE
#   mse.true <- 0
# 
#   mse1.ridge <- mean((c(0, beta) - est1.ridge)^2)
# 
#   mse1.grridge <- mean((c(0, beta) - est1.grridge)^2)
#   mse2.grridge <- apply(est2.grridge, 2, function(b) {mean((c(0, beta) - b)^2)})
# 
#   mse1.enet <- apply(est1.enet, 2, function(b) {mean((c(0, beta) - b)^2)})
#   mse2.enet <- apply(est2.enet, 2, function(b) {mean((c(0, beta) - b)^2)})
#   mse3.enet <- apply(est3.enet, 2, function(b) {mean((c(0, beta) - b)^2)})
#   
#   mse1.greben <- apply(est1.greben, 2, function(b) {mean((c(0, beta) - b)^2)})
#   mse2.greben <- apply(est2.greben, 2, function(b) {mean((c(0, beta) - b)^2)})
#   mse3.greben <- apply(est3.greben, 2, function(b) {mean((c(0, beta) - b)^2)})
#   
#   auc1[[1]][[r]] <- cbind(psel=p, auc=auc1.ridge)
#   auc1[[2]][[r]] <- cbind(psel=p, auc=auc1.grridge)
#   auc1[[3]][[r]] <- cbind(psel=psel1.grridge, auc=auc2.grridge)
#   auc1[[4]][[r]] <- cbind(psel=psel1.enet, auc=auc1.enet)
#   auc1[[5]][[r]] <- cbind(psel=psel2.enet, auc=auc2.enet)
#   auc1[[6]][[r]] <- cbind(psel=psel3.enet, auc=auc3.enet)
#   auc1[[7]][[r]] <- cbind(psel=psel1.greben, auc=auc1.greben)
#   auc1[[8]][[r]] <- cbind(psel=psel2.greben, auc=auc2.greben)
#   auc1[[9]][[r]] <- cbind(psel=psel3.greben, auc=auc3.greben)
#   
#   briers1[[1]][[r]] <- cbind(psel=p, briers=briers1.ridge)
#   briers1[[2]][[r]] <- cbind(psel=p, briers=briers1.grridge)
#   briers1[[3]][[r]] <- cbind(psel=psel1.grridge, briers=briers2.grridge)
#   briers1[[4]][[r]] <- cbind(psel=psel1.enet, briers=briers1.enet)
#   briers1[[5]][[r]] <- cbind(psel=psel2.enet, briers=briers2.enet)
#   briers1[[6]][[r]] <- cbind(psel=psel3.enet, briers=briers3.enet)
#   briers1[[7]][[r]] <- cbind(psel=psel1.greben, briers=briers1.greben)
#   briers1[[8]][[r]] <- cbind(psel=psel2.greben, briers=briers2.greben)
#   briers1[[9]][[r]] <- cbind(psel=psel3.greben, briers=briers3.greben)
#   
#   mse1[[1]][[r]] <- cbind(psel=p, mse=mse1.ridge)
#   mse1[[2]][[r]] <- cbind(psel=p, mse=mse1.grridge)
#   mse1[[3]][[r]] <- cbind(psel=psel1.grridge, mse=mse2.grridge)
#   mse1[[4]][[r]] <- cbind(psel=psel1.enet, mse=mse1.enet)
#   mse1[[5]][[r]] <- cbind(psel=psel2.enet, mse=mse2.enet)
#   mse1[[6]][[r]] <- cbind(psel=psel3.enet, mse=mse3.enet)
#   mse1[[7]][[r]] <- cbind(psel=psel1.greben, mse=mse1.greben)
#   mse1[[8]][[r]] <- cbind(psel=psel2.greben, mse=mse2.greben)
#   mse1[[9]][[r]] <- cbind(psel=psel3.greben, mse=mse3.greben)
# 
# 
#   lambdag1[[1]][r, ] <- fit1.grridge[[1]]$lambdamults$groups
#   lambdag1[[2]][r, ] <- fit1.greben$lambdag$groups[, fit1.greben$nouteriter + 1]
#   lambdag1[[3]][r, ] <- fit2.greben$lambdag$groups[, fit2.greben$nouteriter + 1]
#   lambdag1[[4]][r, ] <- fit3.greben$lambdag$groups[, fit3.greben$nouteriter + 1]
# 
#   results1 <- list(auc=auc1, briers=briers1, mse=mse1, lambdag=lambdag1)
#   save(results1, file=paste(path.res, "grEBEN_test2_res1.Rdata", sep=""))
# 
# }
# 
### plots
load(paste(path.res, "grEBEN_test2_res1.Rdata", sep=""))

# preparing data
auc1.ridge <- do.call(rbind, results1$auc[[1]])[
  order(do.call(rbind, results1$auc[[1]])[, 1]), ]

auc1.grridge <- do.call(rbind, results1$auc[[2]])[
  order(do.call(rbind, results1$auc[[2]])[, 1]), ]
auc2.grridge <- do.call(rbind, results1$auc[[3]])[
  order(do.call(rbind, results1$auc[[3]])[, 1]), ]

auc1.enet <- do.call(rbind, results1$auc[[4]])[
  order(do.call(rbind, results1$auc[[4]])[, 1]), ]
auc2.enet <- do.call(rbind, results1$auc[[5]])[
  order(do.call(rbind, results1$auc[[5]])[, 1]), ]
auc3.enet <- do.call(rbind, results1$auc[[6]])[
  order(do.call(rbind, results1$auc[[6]])[, 1]), ]

auc1.greben <- do.call(rbind, results1$auc[[7]])[
  order(do.call(rbind, results1$auc[[7]])[, 1]), ]
auc2.greben <- do.call(rbind, results1$auc[[8]])[
  order(do.call(rbind, results1$auc[[8]])[, 1]), ]
auc3.greben <- do.call(rbind, results1$auc[[9]])[
  order(do.call(rbind, results1$auc[[9]])[, 1]), ]

lauc1.enet <- lowess(auc1.enet[, 1], auc1.enet[, 2])
lauc2.enet <- lowess(auc2.enet[, 1], auc2.enet[, 2])
lauc3.enet <- lowess(auc3.enet[, 1], auc3.enet[, 2])

lauc1.greben <- lowess(auc1.greben[, 1], auc1.greben[, 2])
lauc2.greben <- lowess(auc2.greben[, 1], auc2.greben[, 2])
lauc3.greben <- lowess(auc3.greben[, 1], auc3.greben[, 2])

lauc1.grridge <- lowess(auc1.grridge[, 1], auc1.grridge[, 2])
lauc2.grridge <- lowess(auc2.grridge[, 1], auc2.grridge[, 2])

lauc1.ridge <- lowess(auc1.ridge[, 1], auc1.ridge[, 2])

briers1.ridge <- do.call(rbind, results1$briers[[1]])[
  order(do.call(rbind, results1$briers[[1]])[, 1]), ]

briers1.grridge <- do.call(rbind, results1$briers[[2]])[
  order(do.call(rbind, results1$briers[[2]])[, 1]), ]
briers2.grridge <- do.call(rbind, results1$briers[[3]])[
  order(do.call(rbind, results1$briers[[3]])[, 1]), ]

briers1.enet <- do.call(rbind, results1$briers[[4]])[
  order(do.call(rbind, results1$briers[[4]])[, 1]), ]
briers2.enet <- do.call(rbind, results1$briers[[5]])[
  order(do.call(rbind, results1$briers[[5]])[, 1]), ]
briers3.enet <- do.call(rbind, results1$briers[[6]])[
  order(do.call(rbind, results1$briers[[6]])[, 1]), ]

briers1.greben <- do.call(rbind, results1$briers[[7]])[
  order(do.call(rbind, results1$briers[[7]])[, 1]), ]
briers2.greben <- do.call(rbind, results1$briers[[8]])[
  order(do.call(rbind, results1$briers[[8]])[, 1]), ]
briers3.greben <- do.call(rbind, results1$briers[[9]])[
  order(do.call(rbind, results1$briers[[9]])[, 1]), ]

lbriers1.enet <- lowess(briers1.enet[, 1], briers1.enet[, 2])
lbriers2.enet <- lowess(briers2.enet[, 1], briers2.enet[, 2])
lbriers3.enet <- lowess(briers3.enet[, 1], briers3.enet[, 2])

lbriers1.greben <- lowess(briers1.greben[, 1], briers1.greben[, 2])
lbriers2.greben <- lowess(briers2.greben[, 1], briers2.greben[, 2])
lbriers3.greben <- lowess(briers3.greben[, 1], briers3.greben[, 2])

lbriers1.grridge <- lowess(briers1.grridge[, 1], briers1.grridge[, 2])
lbriers2.grridge <- lowess(briers2.grridge[, 1], briers2.grridge[, 2])

lbriers1.ridge <- lowess(briers1.ridge[, 1], briers1.ridge[, 2])

mse1.ridge <- do.call(rbind, results1$mse[[1]])[
  order(do.call(rbind, results1$mse[[1]])[, 1]), ]

mse1.grridge <- do.call(rbind, results1$mse[[2]])[
  order(do.call(rbind, results1$mse[[2]])[, 1]), ]
mse2.grridge <- do.call(rbind, results1$mse[[3]])[
  order(do.call(rbind, results1$mse[[3]])[, 1]), ]

mse1.enet <- do.call(rbind, results1$mse[[4]])[
  order(do.call(rbind, results1$mse[[4]])[, 1]), ]
mse2.enet <- do.call(rbind, results1$mse[[5]])[
  order(do.call(rbind, results1$mse[[5]])[, 1]), ]
mse3.enet <- do.call(rbind, results1$mse[[6]])[
  order(do.call(rbind, results1$mse[[6]])[, 1]), ]

mse1.greben <- do.call(rbind, results1$mse[[7]])[
  order(do.call(rbind, results1$mse[[7]])[, 1]), ]
mse2.greben <- do.call(rbind, results1$mse[[8]])[
  order(do.call(rbind, results1$mse[[8]])[, 1]), ]
mse3.greben <- do.call(rbind, results1$mse[[9]])[
  order(do.call(rbind, results1$mse[[9]])[, 1]), ]

lmse1.enet <- lowess(mse1.enet[, 1], mse1.enet[, 2])
lmse2.enet <- lowess(mse2.enet[, 1], mse2.enet[, 2])
lmse3.enet <- lowess(mse3.enet[, 1], mse3.enet[, 2])

lmse1.greben <- lowess(mse1.greben[, 1], mse1.greben[, 2])
lmse2.greben <- lowess(mse2.greben[, 1], mse2.greben[, 2])
lmse3.greben <- lowess(mse3.greben[, 1], mse3.greben[, 2])

lmse1.grridge <- lowess(mse1.grridge[, 1], mse1.grridge[, 2])
lmse2.grridge <- lowess(mse2.grridge[, 1], mse2.grridge[, 2])

lmse1.ridge <- lowess(mse1.ridge[, 1], mse1.ridge[, 2])

### diagnostics: checking lowess fits
# AUC
xlim2.auc <- range(auc1.enet[, 1], auc2.enet[, 1], auc3.enet[, 1],
                   auc1.greben[, 1], auc2.greben[, 1], auc3.greben[, 1],
                   auc2.grridge[, 1])
ylim2.auc <- range(auc1.enet[, 2], auc2.enet[, 2], auc3.enet[, 2],
                   auc1.greben[, 2], auc2.greben[, 2], auc3.greben[, 2],
                   auc2.grridge[, 2])

png(paste(path.graph, "grEBEN_test2_res1_lowess_auc.png", sep=""),
    units="in", width=10, height=6, res=120)
par(mfrow=c(2, 4))

plot(auc1.greben[, 1], auc1.greben[, 2], col=1, ylab="AUC", main="a)",
     xlab="Number of selected variables", xlim=xlim2.auc, ylim=ylim2.auc)
lines(lauc1.greben, col=3)

plot(auc2.greben[, 1], auc2.greben[, 2], col=1, ylab="AUC", main="b)",
     xlab="Number of selected variables", xlim=xlim2.auc, ylim=ylim2.auc)
lines(lauc2.greben, col=4)

plot(auc3.greben[, 1], auc3.greben[, 2], col=1, ylab="AUC", main="c)",
     xlab="Number of selected variables", xlim=xlim2.auc, ylim=ylim2.auc)
lines(lauc3.greben, col=5)

plot(auc2.grridge[, 1], auc2.grridge[, 2], col=1, ylab="AUC", main="d)",
     xlab="Number of selected variables", xlim=xlim2.auc, ylim=ylim2.auc)
lines(lauc2.grridge, col=2)

plot(auc1.enet[, 1], auc1.enet[, 2], col=1, ylab="AUC", main="e)",
     xlab="Number of selected variables", xlim=xlim2.auc, ylim=ylim2.auc)
lines(lauc1.enet, col=3)

plot(auc2.enet[, 1], auc2.enet[, 2], col=1, ylab="AUC", main="f)",
     xlab="Number of selected variables", xlim=xlim2.auc, ylim=ylim2.auc)
lines(lauc2.enet, col=4)

plot(auc3.enet[, 1], auc3.enet[, 2], col=1, ylab="AUC", main="g)",
     xlab="Number of selected variables", xlim=xlim2.auc, ylim=ylim2.auc)
lines(lauc3.enet, col=5)
dev.off()

# briers
xlim2.briers <- range(briers1.enet[, 1], briers2.enet[, 1], briers3.enet[, 1],
                   briers1.greben[, 1], briers2.greben[, 1], briers3.greben[, 1],
                   briers2.grridge[, 1])
ylim2.briers <- range(briers1.enet[, 2], briers2.enet[, 2], briers3.enet[, 2],
                   briers1.greben[, 2], briers2.greben[, 2], briers3.greben[, 2],
                   briers2.grridge[, 2])

png(paste(path.graph, "grEBEN_test2_res1_lowess_briers.png", sep=""),
    units="in", width=10, height=6, res=120)
par(mfrow=c(2, 4))

plot(briers1.greben[, 1], briers1.greben[, 2], col=1, ylab="briers", main="a)",
     xlab="Number of selected variables", xlim=xlim2.briers, ylim=ylim2.briers)
lines(lbriers1.greben, col=3)

plot(briers2.greben[, 1], briers2.greben[, 2], col=1, ylab="briers", main="b)",
     xlab="Number of selected variables", xlim=xlim2.briers, ylim=ylim2.briers)
lines(lbriers2.greben, col=4)

plot(briers3.greben[, 1], briers3.greben[, 2], col=1, ylab="briers", main="c)",
     xlab="Number of selected variables", xlim=xlim2.briers, ylim=ylim2.briers)
lines(lbriers3.greben, col=5)

plot(briers2.grridge[, 1], briers2.grridge[, 2], col=1, ylab="briers", main="d)",
     xlab="Number of selected variables", xlim=xlim2.briers, ylim=ylim2.briers)
lines(lbriers2.grridge, col=2)

plot(briers1.enet[, 1], briers1.enet[, 2], col=1, ylab="briers", main="e)",
     xlab="Number of selected variables", xlim=xlim2.briers, ylim=ylim2.briers)
lines(lbriers1.enet, col=4)

plot(briers2.enet[, 1], briers2.enet[, 2], col=1, ylab="briers", main="f)",
     xlab="Number of selected variables", xlim=xlim2.briers, ylim=ylim2.briers)
lines(lbriers2.enet, col=5)

plot(briers3.enet[, 1], briers3.enet[, 2], col=1, ylab="briers", main="g)",
     xlab="Number of selected variables", xlim=xlim2.briers, ylim=ylim2.briers)
lines(lbriers3.enet, col=6)
dev.off()

# mse
xlim2.mse <- range(mse1.enet[, 1], mse2.enet[, 1], mse3.enet[, 1],
                   mse1.greben[, 1], mse2.greben[, 1], mse3.greben[, 1],
                   mse2.grridge[, 1])
ylim2.mse <- range(mse1.enet[, 2], mse2.enet[, 2], mse3.enet[, 2],
                   mse1.greben[, 2], mse2.greben[, 2], mse3.greben[, 2],
                   mse2.grridge[, 2])

png(paste(path.graph, "grEBEN_test2_res1_lowess_mse.png", sep=""),
    units="in", width=10, height=6, res=120)
par(mfrow=c(2, 4))

plot(mse1.greben[, 1], mse1.greben[, 2], col=1, ylab="mse", main="a)",
     xlab="Number of selected variables", xlim=xlim2.mse, ylim=ylim2.mse)
lines(lmse1.greben, col=3)

plot(mse2.greben[, 1], mse2.greben[, 2], col=1, ylab="mse", main="b)",
     xlab="Number of selected variables", xlim=xlim2.mse, ylim=ylim2.mse)
lines(lmse2.greben, col=4)

plot(mse3.greben[, 1], mse3.greben[, 2], col=1, ylab="mse", main="c)",
     xlab="Number of selected variables", xlim=xlim2.mse, ylim=ylim2.mse)
lines(lmse3.greben, col=5)

plot(mse2.grridge[, 1], mse2.grridge[, 2], col=1, ylab="mse", main="d)",
     xlab="Number of selected variables", xlim=xlim2.mse, ylim=ylim2.mse)
lines(lmse2.grridge, col=2)

plot(mse1.enet[, 1], mse1.enet[, 2], col=1, ylab="mse", main="e)",
     xlab="Number of selected variables", xlim=xlim2.mse, ylim=ylim2.mse)
lines(lmse1.enet, col=3)

plot(mse2.enet[, 1], mse2.enet[, 2], col=1, ylab="mse", main="f)",
     xlab="Number of selected variables", xlim=xlim2.mse, ylim=ylim2.mse)
lines(lmse2.enet, col=4)

plot(mse3.enet[, 1], mse3.enet[, 2], col=1, ylab="mse", main="g)",
     xlab="Number of selected variables", xlim=xlim2.mse, ylim=ylim2.mse)
lines(lmse3.enet, col=5)
dev.off()



### Performance measures
# AUC
png(paste(path.graph, "grEBEN_test2_res1_performance.png", sep=""),
    units="in", width=12, height=4, res=120)
par(mfrow=c(1, 3))
xlim1.auc <- range(lauc1.enet$x, lauc2.enet$x, lauc3.enet$x, lauc1.greben$x,
                   lauc2.greben$x, lauc3.greben$x, lauc2.grridge$x)
ylim1.auc <- range(lauc1.enet$y, lauc2.enet$y, lauc3.enet$y, lauc1.greben$y,
                   lauc2.greben$y, lauc3.greben$y, lauc1.grridge$y,
                   lauc2.grridge$y, lauc1.ridge$y)

plot(0, 0, col=2, type="n", ylab="AUC", main="a)",
     xlab="Number of selected variables", xlim=xlim1.auc, ylim=ylim1.auc)

lines(lauc2.grridge, col=2, lty=1)
abline(h=lauc1.grridge$y, col=2, lty=2)
abline(h=lauc1.ridge$y, col=2, lty=3)

lines(lauc1.greben, col=3, lty=1)
lines(lauc1.enet, col=3, lty=3)

lines(lauc2.greben, col=4, lty=1)
lines(lauc2.enet, col=4, lty=3)

lines(lauc3.greben, col=5, lty=1)
lines(lauc3.enet, col=5, lty=3)

# Brier skilll
xlim1.briers <- range(lbriers1.enet$x, lbriers2.enet$x, lbriers3.enet$x, lbriers1.greben$x,
                   lbriers2.greben$x, lbriers3.greben$x, lbriers2.grridge$x)
ylim1.briers <- range(lbriers1.enet$y, lbriers2.enet$y, lbriers3.enet$y, lbriers1.greben$y,
                   lbriers2.greben$y, lbriers3.greben$y, lbriers1.grridge$y,
                   lbriers2.grridge$y, lbriers1.ridge$y)

plot(0, 0, col=2, type="n", ylab="Brier skill score", main="a)",
     xlab="Number of selected variables", xlim=xlim1.briers, ylim=ylim1.briers)

lines(lbriers2.grridge, col=2, lty=1)
abline(h=lbriers1.grridge$y, col=2, lty=2)
abline(h=lbriers1.ridge$y, col=2, lty=3)

lines(lbriers1.greben, col=3, lty=1)
lines(lbriers1.enet, col=3, lty=3)

lines(lbriers2.greben, col=4, lty=1)
lines(lbriers2.enet, col=4, lty=3)

lines(lbriers3.greben, col=5, lty=1)
lines(lbriers3.enet, col=5, lty=3)

leglabels <- c("ridge",
               expression(paste("enet, ", alpha==0.05)),
               expression(paste("enet, ", alpha==0.5)),
               expression(paste("enet, ", alpha==0.95)),
               "group-regularized + selection", "group-regularized",
               "not group-regularized")
legend("bottomright", legend=leglabels, fill=c(2:5, 0, 0, 0),
       lty=c(rep(NA, 4), 1, 2, 3), border=c(rep(1, 4), 0 ,0, 0), merge=TRUE,
       seg.len=1)

# MSE
xlim1.mse <- range(lmse1.enet$x, lmse2.enet$x, lmse3.enet$x, lmse1.greben$x,
                      lmse2.greben$x, lmse3.greben$x, lmse2.grridge$x)
ylim1.mse <- range(lmse1.enet$y, lmse2.enet$y, lmse3.enet$y, lmse1.greben$y,
                      lmse2.greben$y, lmse3.greben$y, lmse1.grridge$y,
                      lmse2.grridge$y, lmse1.ridge$y)

plot(0, 0, col=2, type="n", ylab="MSE", main="a)",
     xlab="Number of selected variables", xlim=xlim1.mse, ylim=ylim1.mse)

lines(lmse2.grridge, col=2, lty=1)
abline(h=lmse1.grridge$y, col=2, lty=2)
abline(h=lmse1.ridge$y, col=2, lty=3)

lines(lmse1.greben, col=3, lty=1)
lines(lmse1.enet, col=3, lty=3)

lines(lmse2.greben, col=4, lty=1)
lines(lmse2.enet, col=4, lty=3)

lines(lmse3.greben, col=5, lty=1)
lines(lmse3.enet, col=5, lty=3)
dev.off()

### estimated penalty parameters
png(paste(path.graph, "grEBEN_test2_res1_penalties.png", sep=""),
    units="in", width=8, height=6, res=120)
par(mfrow=c(2, 2))
boxplot(results1$lambdag[[1]], main="a)",  xlab="", ylab="")
title(ylab=expression(paste(lambda[g], "'")), 
      xlab="Groups in increasing effect size", line=2.5)
boxplot(results1$lambdag[[2]], main="b)",  xlab="", ylab="")
title(ylab=expression(paste(lambda[g], "'")), 
      xlab="Groups in increasing effect size", line=2.5)
boxplot(results1$lambdag[[3]], main="c)",  xlab="", ylab="")
title(ylab=expression(paste(lambda[g], "'")), 
      xlab="Groups in increasing effect size", line=2.5)
boxplot(results1$lambdag[[4]], main="d)",  xlab="", ylab="")
title(ylab=expression(paste(lambda[g], "'")), 
      xlab="Groups in increasing effect size", line=2.5)
par(fig=c(0.15, 0.5, 0.6, 1), new=TRUE)
boxplot(results1$lambdag[[1]], outline=FALSE, xaxt="n")
dev.off()




## simulation 2
# create data
n <- 100
p <- 1000
G <- 10
rhow <- 0.7
rhob <- 0.5
blockw <- blockb <- diag(p/G)
for(i in 1:ncol(blockw)) {
  for(j in 1:ncol(blockw)) {
    blockw[i, j] <- rhow^abs(i - j)
    blockb[i, j] <- rhob^(abs(i - j) + 1)
  }
}
Sigma <- rbind(cbind(blockw, blockb), cbind(blockb, blockw))
q <- 0.3 # proportion of non-zero coefficients
beta.mean <- 0.05
f <- 1.4
m <- rep(1, n)
betag <- beta.mean*G*f^(c(1:G) - 1)/(q*sum(f^(c(1:G) - 1)))
beta <- as.numeric(sapply(betag, function(b) {
  rep(c(0, b), times=c((1 - q)*p/G, q*p/G))}))
part.greben <- list(groups=rep(1:G, each=p/G))
part.grridge <- list(groups=CreatePartition(as.factor(part.greben$groups)))
ntest <- 1000
mtest <- rep(1, ntest)

methods <- c("ridge", "GRridge", "GRridge+sel",
             "enet+a=0.05", "enet+a=0.5", "enet+a=0.95",
             "grEBEN+a=0.05", "grEBEN+a=0.5","grEBEN+a=0.95")
nreps <- 50
auc2 <- briers2 <- mse2 <- rep(list(vector(mode="list", length=nreps)),
                               length(methods))
names(auc2) <- names(briers2) <- names(mse2) <- methods
lambdag2 <- rep(list(matrix(NA, ncol=G, nrow=nreps)), 4)
names(lambdag2) <- methods[c(2, 7, 8, 9)]
# the simulations
for(r in 1:nreps) {

  set.seed(400 + r)
  x <- do.call(cbind, replicate(G/2, rmvnorm(n, rep(0, 2*p/G), Sigma),
                                simplify=FALSE))
  prob <- as.numeric(exp(x %*% beta)/(1 + exp(x %*% beta)))
  y <- rbinom(n, m, prob)
  
  xtest <- do.call(cbind, replicate(G/2, rmvnorm(ntest, rep(0, 2*p/G), Sigma),
                                    simplify=FALSE))
  probtest <- as.numeric(exp(xtest %*% beta)/(1 + exp(xtest %*% beta)))
  ytest <- rbinom(ntest, mtest, probtest)
  
  fit1.greben <- grEBEN3(x, y, m, partitions=part.greben, alpha=0.05, psel=TRUE)
  fit2.greben <- grEBEN3(x, y, m, partitions=part.greben, alpha=0.5, psel=TRUE)
  fit3.greben <- grEBEN3(x, y, m, partitions=part.greben, alpha=0.95, psel=TRUE)
  
  # number of selected variables
  psel1.enet <- apply(fit1.greben$beta.nogroups, 2, function(b) {sum(b!=0) - 1})
  psel2.enet <- apply(fit2.greben$beta.nogroups, 2, function(b) {sum(b!=0) - 1})
  psel3.enet <- apply(fit3.greben$beta.nogroups, 2, function(b) {sum(b!=0) - 1})
  
  psel1.greben <- apply(fit1.greben$beta, 2, function(b) {sum(b!=0) - 1})
  psel2.greben <- apply(fit2.greben$beta, 2, function(b) {sum(b!=0) - 1})
  psel3.greben <- apply(fit3.greben$beta, 2, function(b) {sum(b!=0) - 1})
  
  psel.all <- unique(c(psel1.enet, psel2.enet, psel3.enet,
                       psel1.greben, psel2.greben, psel3.greben))
  psel.in <- floor(quantile(psel.all[psel.all!=0],
                            prob=c(seq(0.01, 0.05, 0.01), seq(0.07, 0.25, 0.03), seq(0.3, 1, 0.1))))
  
  fit1.grridge <- vector(mode="list", length=length(psel.in))
  fit1.grridge[[1]] <- grridge(t(x), y, part.grridge, selectionEN=TRUE,
                               maxsel=psel.in[1])
  for(s in 2:length(psel.in)) {
    fit1.grridge[[s]] <- grridge(t(x), y, part.grridge, selectionEN=TRUE,
                                 optl=fit1.grridge[[1]]$optl,
                                 maxsel=psel.in[s])
  }
  
  psel1.grridge <- sapply(1:length(fit1.grridge), function(s) {
    return(length(fit1.grridge[[s]]$resEN$whichEN))})
  
  # estimates
  est1.ridge <- coef(fit1.grridge[[1]]$predobj$NoGroups, "all")
  
  est1.grridge <- coef(fit1.grridge[[1]]$predobj$GroupRegul, "all")
  
  est2.grridge <- sapply(fit1.grridge, function(s) {
    replace(rep(0, p + 1), c(1, s$resEN$whichEN + 1), 
            coef(s$predobj$EN, "all"))})
  
  est1.enet <- fit1.greben$beta.nogroups
  est2.enet <- fit2.greben$beta.nogroups
  est3.enet <- fit3.greben$beta.nogroups
  
  est1.greben <- fit1.greben$beta
  est2.greben <- fit2.greben$beta
  est3.greben <- fit3.greben$beta
  
  # predictions on fit data
  pred1.ridge <- predict.grridge(fit1.grridge[[1]], t(xtest))[, 1]
  
  pred1.grridge <- predict.grridge(fit1.grridge[[1]], t(xtest))[, 2]
  pred2.grridge <- sapply(fit1.grridge, function(s) {
    predict.grridge(s, t(xtest))[, 3]})
  
  pred1.enet <- 1/(1 + exp(-xtest %*% est1.enet[-1, ]))
  pred2.enet <- 1/(1 + exp(-xtest %*% est2.enet[-1, ]))
  pred3.enet <- 1/(1 + exp(-xtest %*% est3.enet[-1, ]))
  
  pred1.greben <- 1/(1 + exp(-xtest %*% est1.greben[-1, ]))
  pred2.greben <- 1/(1 + exp(-xtest %*% est2.greben[-1, ]))
  pred3.greben <- 1/(1 + exp(-xtest %*% est3.greben[-1, ]))
  
  # AUCs
  auc.true <- pROC::roc(ytest, probtest)$auc
  
  auc1.ridge <- pROC::roc(ytest, pred1.ridge)$auc
  
  auc1.grridge <- pROC::roc(ytest, pred1.grridge)$auc
  auc2.grridge <- apply(pred2.grridge, 2, function(r) {pROC::roc(ytest, r)$auc})
  
  auc1.enet <- apply(pred1.enet, 2, function(r) {pROC::roc(ytest, r)$auc})
  auc2.enet <- apply(pred2.enet, 2, function(r) {pROC::roc(ytest, r)$auc})
  auc3.enet <- apply(pred3.enet, 2, function(r) {pROC::roc(ytest, r)$auc})
  
  auc1.greben <- apply(pred1.greben, 2, function(r) {pROC::roc(ytest, r)$auc})
  auc2.greben <- apply(pred2.greben, 2, function(r) {pROC::roc(ytest, r)$auc})
  auc3.greben <- apply(pred3.greben, 2, function(r) {pROC::roc(ytest, r)$auc})
  
  # Brier scores
  brier.null <- sum((ytest - mean(ytest))^2)
  briers.true <- 1 - sum((ytest - probtest)^2)/brier.null
  
  briers1.ridge <- 1 - sum((ytest - pred1.ridge)^2)/brier.null
  
  briers1.grridge <- 1 - sum((ytest - pred1.grridge)^2)/brier.null
  briers2.grridge <- apply(pred2.grridge, 2, function(pred) {
    1 - sum((ytest - pred)^2)/brier.null})
  
  briers1.enet <- apply(pred1.enet, 2, function(pred) {
    1 - sum((ytest - pred)^2)/brier.null})
  briers2.enet <- apply(pred2.enet, 2, function(pred) {
    1 - sum((ytest - pred)^2)/brier.null})
  briers3.enet <- apply(pred3.enet, 2, function(pred) {
    1 - sum((ytest - pred)^2)/brier.null})
  
  briers1.greben <- apply(pred1.greben, 2, function(pred) {
    1 - sum((ytest - pred)^2)/brier.null})
  briers2.greben <- apply(pred2.greben, 2, function(pred) {
    1 - sum((ytest - pred)^2)/brier.null})
  briers3.greben <- apply(pred3.greben, 2, function(pred) {
    1 - sum((ytest - pred)^2)/brier.null})
  
  # MSE
  mse.true <- 0
  
  mse1.ridge <- mean((c(0, beta) - est1.ridge)^2)
  
  mse1.grridge <- mean((c(0, beta) - est1.grridge)^2)
  mse2.grridge <- apply(est2.grridge, 2, function(b) {mean((c(0, beta) - b)^2)})
  
  mse1.enet <- apply(est1.enet, 2, function(b) {mean((c(0, beta) - b)^2)})
  mse2.enet <- apply(est2.enet, 2, function(b) {mean((c(0, beta) - b)^2)})
  mse3.enet <- apply(est3.enet, 2, function(b) {mean((c(0, beta) - b)^2)})
  
  mse1.greben <- apply(est1.greben, 2, function(b) {mean((c(0, beta) - b)^2)})
  mse2.greben <- apply(est2.greben, 2, function(b) {mean((c(0, beta) - b)^2)})
  mse3.greben <- apply(est3.greben, 2, function(b) {mean((c(0, beta) - b)^2)})
  
  auc2[[1]][[r]] <- cbind(psel=p, auc=auc1.ridge)
  auc2[[2]][[r]] <- cbind(psel=p, auc=auc1.grridge)
  auc2[[3]][[r]] <- cbind(psel=psel1.grridge, auc=auc2.grridge)
  auc2[[4]][[r]] <- cbind(psel=psel1.enet, auc=auc1.enet)
  auc2[[5]][[r]] <- cbind(psel=psel2.enet, auc=auc2.enet)
  auc2[[6]][[r]] <- cbind(psel=psel3.enet, auc=auc3.enet)
  auc2[[7]][[r]] <- cbind(psel=psel1.greben, auc=auc1.greben)
  auc2[[8]][[r]] <- cbind(psel=psel2.greben, auc=auc2.greben)
  auc2[[9]][[r]] <- cbind(psel=psel3.greben, auc=auc3.greben)
  
  briers2[[1]][[r]] <- cbind(psel=p, briers=briers1.ridge)
  briers2[[2]][[r]] <- cbind(psel=p, briers=briers1.grridge)
  briers2[[3]][[r]] <- cbind(psel=psel1.grridge, briers=briers2.grridge)
  briers2[[4]][[r]] <- cbind(psel=psel1.enet, briers=briers1.enet)
  briers2[[5]][[r]] <- cbind(psel=psel2.enet, briers=briers2.enet)
  briers2[[6]][[r]] <- cbind(psel=psel3.enet, briers=briers3.enet)
  briers2[[7]][[r]] <- cbind(psel=psel1.greben, briers=briers1.greben)
  briers2[[8]][[r]] <- cbind(psel=psel2.greben, briers=briers2.greben)
  briers2[[9]][[r]] <- cbind(psel=psel3.greben, briers=briers3.greben)
  
  mse2[[1]][[r]] <- cbind(psel=p, mse=mse1.ridge)
  mse2[[2]][[r]] <- cbind(psel=p, mse=mse1.grridge)
  mse2[[3]][[r]] <- cbind(psel=psel1.grridge, mse=mse2.grridge)
  mse2[[4]][[r]] <- cbind(psel=psel1.enet, mse=mse1.enet)
  mse2[[5]][[r]] <- cbind(psel=psel2.enet, mse=mse2.enet)
  mse2[[6]][[r]] <- cbind(psel=psel3.enet, mse=mse3.enet)
  mse2[[7]][[r]] <- cbind(psel=psel1.greben, mse=mse1.greben)
  mse2[[8]][[r]] <- cbind(psel=psel2.greben, mse=mse2.greben)
  mse2[[9]][[r]] <- cbind(psel=psel3.greben, mse=mse3.greben)
  
  
  lambdag2[[1]][r, ] <- fit1.grridge[[1]]$lambdamults$groups
  lambdag2[[2]][r, ] <- fit1.greben$lambdag$groups[, fit1.greben$nouteriter + 1]
  lambdag2[[3]][r, ] <- fit2.greben$lambdag$groups[, fit2.greben$nouteriter + 1]
  lambdag2[[4]][r, ] <- fit3.greben$lambdag$groups[, fit3.greben$nouteriter + 1]
  
  results2 <- list(auc=auc2, briers=briers2, mse=mse2, lambdag=lambdag2)
  save(results2, file=paste(path.res, "grEBEN_test2_res2.Rdata", sep=""))
  
}

### plots
load(paste(path.res, "grEBEN_test2_res2.Rdata", sep=""))

# preparing data
auc1.ridge <- do.call(rbind, results2$auc[[1]])[
  order(do.call(rbind, results2$auc[[1]])[, 1]), ]

auc1.grridge <- do.call(rbind, results2$auc[[2]])[
  order(do.call(rbind, results2$auc[[2]])[, 1]), ]
auc2.grridge <- do.call(rbind, results2$auc[[3]])[
  order(do.call(rbind, results2$auc[[3]])[, 1]), ]

auc1.enet <- do.call(rbind, results2$auc[[4]])[
  order(do.call(rbind, results2$auc[[4]])[, 1]), ]
auc2.enet <- do.call(rbind, results2$auc[[5]])[
  order(do.call(rbind, results2$auc[[5]])[, 1]), ]
auc3.enet <- do.call(rbind, results2$auc[[6]])[
  order(do.call(rbind, results2$auc[[6]])[, 1]), ]

auc1.greben <- do.call(rbind, results2$auc[[7]])[
  order(do.call(rbind, results2$auc[[7]])[, 1]), ]
auc2.greben <- do.call(rbind, results2$auc[[8]])[
  order(do.call(rbind, results2$auc[[8]])[, 1]), ]
auc3.greben <- do.call(rbind, results2$auc[[9]])[
  order(do.call(rbind, results2$auc[[9]])[, 1]), ]

lauc1.enet <- lowess(auc1.enet[, 1], auc1.enet[, 2])
lauc2.enet <- lowess(auc2.enet[, 1], auc2.enet[, 2])
lauc3.enet <- lowess(auc3.enet[, 1], auc3.enet[, 2])

lauc1.greben <- lowess(auc1.greben[, 1], auc1.greben[, 2])
lauc2.greben <- lowess(auc2.greben[, 1], auc2.greben[, 2])
lauc3.greben <- lowess(auc3.greben[, 1], auc3.greben[, 2])

lauc1.grridge <- lowess(auc1.grridge[, 1], auc1.grridge[, 2])
lauc2.grridge <- lowess(auc2.grridge[, 1], auc2.grridge[, 2])

lauc1.ridge <- lowess(auc1.ridge[, 1], auc1.ridge[, 2])

briers1.ridge <- do.call(rbind, results2$briers[[1]])[
  order(do.call(rbind, results2$briers[[1]])[, 1]), ]

briers1.grridge <- do.call(rbind, results2$briers[[2]])[
  order(do.call(rbind, results2$briers[[2]])[, 1]), ]
briers2.grridge <- do.call(rbind, results2$briers[[3]])[
  order(do.call(rbind, results2$briers[[3]])[, 1]), ]

briers1.enet <- do.call(rbind, results2$briers[[4]])[
  order(do.call(rbind, results2$briers[[4]])[, 1]), ]
briers2.enet <- do.call(rbind, results2$briers[[5]])[
  order(do.call(rbind, results2$briers[[5]])[, 1]), ]
briers3.enet <- do.call(rbind, results2$briers[[6]])[
  order(do.call(rbind, results2$briers[[6]])[, 1]), ]

briers1.greben <- do.call(rbind, results2$briers[[7]])[
  order(do.call(rbind, results2$briers[[7]])[, 1]), ]
briers2.greben <- do.call(rbind, results2$briers[[8]])[
  order(do.call(rbind, results2$briers[[8]])[, 1]), ]
briers3.greben <- do.call(rbind, results2$briers[[9]])[
  order(do.call(rbind, results2$briers[[9]])[, 1]), ]

lbriers1.enet <- lowess(briers1.enet[, 1], briers1.enet[, 2])
lbriers2.enet <- lowess(briers2.enet[, 1], briers2.enet[, 2])
lbriers3.enet <- lowess(briers3.enet[, 1], briers3.enet[, 2])

lbriers1.greben <- lowess(briers1.greben[, 1], briers1.greben[, 2])
lbriers2.greben <- lowess(briers2.greben[, 1], briers2.greben[, 2])
lbriers3.greben <- lowess(briers3.greben[, 1], briers3.greben[, 2])

lbriers1.grridge <- lowess(briers1.grridge[, 1], briers1.grridge[, 2])
lbriers2.grridge <- lowess(briers2.grridge[, 1], briers2.grridge[, 2])

lbriers1.ridge <- lowess(briers1.ridge[, 1], briers1.ridge[, 2])

mse1.ridge <- do.call(rbind, results2$mse[[1]])[
  order(do.call(rbind, results2$mse[[1]])[, 1]), ]

mse1.grridge <- do.call(rbind, results2$mse[[2]])[
  order(do.call(rbind, results2$mse[[2]])[, 1]), ]
mse2.grridge <- do.call(rbind, results2$mse[[3]])[
  order(do.call(rbind, results2$mse[[3]])[, 1]), ]

mse1.enet <- do.call(rbind, results2$mse[[4]])[
  order(do.call(rbind, results2$mse[[4]])[, 1]), ]
mse2.enet <- do.call(rbind, results2$mse[[5]])[
  order(do.call(rbind, results2$mse[[5]])[, 1]), ]
mse3.enet <- do.call(rbind, results2$mse[[6]])[
  order(do.call(rbind, results2$mse[[6]])[, 1]), ]

mse1.greben <- do.call(rbind, results2$mse[[7]])[
  order(do.call(rbind, results2$mse[[7]])[, 1]), ]
mse2.greben <- do.call(rbind, results2$mse[[8]])[
  order(do.call(rbind, results2$mse[[8]])[, 1]), ]
mse3.greben <- do.call(rbind, results2$mse[[9]])[
  order(do.call(rbind, results2$mse[[9]])[, 1]), ]

lmse1.enet <- lowess(mse1.enet[, 1], mse1.enet[, 2])
lmse2.enet <- lowess(mse2.enet[, 1], mse2.enet[, 2])
lmse3.enet <- lowess(mse3.enet[, 1], mse3.enet[, 2])

lmse1.greben <- lowess(mse1.greben[, 1], mse1.greben[, 2])
lmse2.greben <- lowess(mse2.greben[, 1], mse2.greben[, 2])
lmse3.greben <- lowess(mse3.greben[, 1], mse3.greben[, 2])

lmse1.grridge <- lowess(mse1.grridge[, 1], mse1.grridge[, 2])
lmse2.grridge <- lowess(mse2.grridge[, 1], mse2.grridge[, 2])

lmse1.ridge <- lowess(mse1.ridge[, 1], mse1.ridge[, 2])

### diagnostics: checking lowess fits
# AUC
xlim2.auc <- range(auc1.enet[, 1], auc2.enet[, 1], auc3.enet[, 1],
                   auc1.greben[, 1], auc2.greben[, 1], auc3.greben[, 1],
                   auc2.grridge[, 1])
ylim2.auc <- range(auc1.enet[, 2], auc2.enet[, 2], auc3.enet[, 2],
                   auc1.greben[, 2], auc2.greben[, 2], auc3.greben[, 2],
                   auc2.grridge[, 2])

png(paste(path.graph, "grEBEN_test2_res2_lowess_auc.png", sep=""),
    units="in", width=10, height=6, res=120)
par(mfrow=c(2, 4))

plot(auc1.greben[, 1], auc1.greben[, 2], col=1, ylab="AUC", main="a)",
     xlab="Number of selected variables", xlim=xlim2.auc, ylim=ylim2.auc)
lines(lauc1.greben, col=3)

plot(auc2.greben[, 1], auc2.greben[, 2], col=1, ylab="AUC", main="b)",
     xlab="Number of selected variables", xlim=xlim2.auc, ylim=ylim2.auc)
lines(lauc2.greben, col=4)

plot(auc3.greben[, 1], auc3.greben[, 2], col=1, ylab="AUC", main="c)",
     xlab="Number of selected variables", xlim=xlim2.auc, ylim=ylim2.auc)
lines(lauc3.greben, col=5)

plot(auc2.grridge[, 1], auc2.grridge[, 2], col=1, ylab="AUC", main="d)",
     xlab="Number of selected variables", xlim=xlim2.auc, ylim=ylim2.auc)
lines(lauc2.grridge, col=2)

plot(auc1.enet[, 1], auc1.enet[, 2], col=1, ylab="AUC", main="e)",
     xlab="Number of selected variables", xlim=xlim2.auc, ylim=ylim2.auc)
lines(lauc1.enet, col=3)

plot(auc2.enet[, 1], auc2.enet[, 2], col=1, ylab="AUC", main="f)",
     xlab="Number of selected variables", xlim=xlim2.auc, ylim=ylim2.auc)
lines(lauc2.enet, col=4)

plot(auc3.enet[, 1], auc3.enet[, 2], col=1, ylab="AUC", main="g)",
     xlab="Number of selected variables", xlim=xlim2.auc, ylim=ylim2.auc)
lines(lauc3.enet, col=5)
dev.off()

# briers
xlim2.briers <- range(briers1.enet[, 1], briers2.enet[, 1], briers3.enet[, 1],
                      briers1.greben[, 1], briers2.greben[, 1], briers3.greben[, 1],
                      briers2.grridge[, 1])
ylim2.briers <- range(briers1.enet[, 2], briers2.enet[, 2], briers3.enet[, 2],
                      briers1.greben[, 2], briers2.greben[, 2], briers3.greben[, 2],
                      briers2.grridge[, 2])

png(paste(path.graph, "grEBEN_test2_res2_lowess_briers.png", sep=""),
    units="in", width=10, height=6, res=120)
par(mfrow=c(2, 4))

plot(briers1.greben[, 1], briers1.greben[, 2], col=1, ylab="briers", main="a)",
     xlab="Number of selected variables", xlim=xlim2.briers, ylim=ylim2.briers)
lines(lbriers1.greben, col=3)

plot(briers2.greben[, 1], briers2.greben[, 2], col=1, ylab="briers", main="b)",
     xlab="Number of selected variables", xlim=xlim2.briers, ylim=ylim2.briers)
lines(lbriers2.greben, col=4)

plot(briers3.greben[, 1], briers3.greben[, 2], col=1, ylab="briers", main="c)",
     xlab="Number of selected variables", xlim=xlim2.briers, ylim=ylim2.briers)
lines(lbriers3.greben, col=5)

plot(briers2.grridge[, 1], briers2.grridge[, 2], col=1, ylab="briers", main="d)",
     xlab="Number of selected variables", xlim=xlim2.briers, ylim=ylim2.briers)
lines(lbriers2.grridge, col=2)

plot(briers1.enet[, 1], briers1.enet[, 2], col=1, ylab="briers", main="e)",
     xlab="Number of selected variables", xlim=xlim2.briers, ylim=ylim2.briers)
lines(lbriers1.enet, col=4)

plot(briers2.enet[, 1], briers2.enet[, 2], col=1, ylab="briers", main="f)",
     xlab="Number of selected variables", xlim=xlim2.briers, ylim=ylim2.briers)
lines(lbriers2.enet, col=5)

plot(briers3.enet[, 1], briers3.enet[, 2], col=1, ylab="briers", main="g)",
     xlab="Number of selected variables", xlim=xlim2.briers, ylim=ylim2.briers)
lines(lbriers3.enet, col=6)
dev.off()

# mse
xlim2.mse <- range(mse1.enet[, 1], mse2.enet[, 1], mse3.enet[, 1],
                   mse1.greben[, 1], mse2.greben[, 1], mse3.greben[, 1],
                   mse2.grridge[, 1])
ylim2.mse <- range(mse1.enet[, 2], mse2.enet[, 2], mse3.enet[, 2],
                   mse1.greben[, 2], mse2.greben[, 2], mse3.greben[, 2],
                   mse2.grridge[, 2])

png(paste(path.graph, "grEBEN_test2_res2_lowess_mse.png", sep=""),
    units="in", width=10, height=6, res=120)
par(mfrow=c(2, 4))

plot(mse1.greben[, 1], mse1.greben[, 2], col=1, ylab="mse", main="a)",
     xlab="Number of selected variables", xlim=xlim2.mse, ylim=ylim2.mse)
lines(lmse1.greben, col=3)

plot(mse2.greben[, 1], mse2.greben[, 2], col=1, ylab="mse", main="b)",
     xlab="Number of selected variables", xlim=xlim2.mse, ylim=ylim2.mse)
lines(lmse2.greben, col=4)

plot(mse3.greben[, 1], mse3.greben[, 2], col=1, ylab="mse", main="c)",
     xlab="Number of selected variables", xlim=xlim2.mse, ylim=ylim2.mse)
lines(lmse3.greben, col=5)

plot(mse2.grridge[, 1], mse2.grridge[, 2], col=1, ylab="mse", main="d)",
     xlab="Number of selected variables", xlim=xlim2.mse, ylim=ylim2.mse)
lines(lmse2.grridge, col=2)

plot(mse1.enet[, 1], mse1.enet[, 2], col=1, ylab="mse", main="e)",
     xlab="Number of selected variables", xlim=xlim2.mse, ylim=ylim2.mse)
lines(lmse1.enet, col=3)

plot(mse2.enet[, 1], mse2.enet[, 2], col=1, ylab="mse", main="f)",
     xlab="Number of selected variables", xlim=xlim2.mse, ylim=ylim2.mse)
lines(lmse2.enet, col=4)

plot(mse3.enet[, 1], mse3.enet[, 2], col=1, ylab="mse", main="g)",
     xlab="Number of selected variables", xlim=xlim2.mse, ylim=ylim2.mse)
lines(lmse3.enet, col=5)
dev.off()



### Performance measures
# AUC
png(paste(path.graph, "grEBEN_test2_res2_performance.png", sep=""),
    units="in", width=12, height=4, res=120)
par(mfrow=c(1, 3))
xlim1.auc <- range(lauc1.enet$x, lauc2.enet$x, lauc3.enet$x, lauc1.greben$x,
                   lauc2.greben$x, lauc3.greben$x, lauc2.grridge$x)
ylim1.auc <- range(lauc1.enet$y, lauc2.enet$y, lauc3.enet$y, lauc1.greben$y,
                   lauc2.greben$y, lauc3.greben$y, lauc1.grridge$y,
                   lauc2.grridge$y, lauc1.ridge$y)

plot(0, 0, col=2, type="n", ylab="AUC", main="a)",
     xlab="Number of selected variables", xlim=xlim1.auc, ylim=ylim1.auc)

lines(lauc2.grridge, col=2, lty=1)
abline(h=lauc1.grridge$y, col=2, lty=2)
abline(h=lauc1.ridge$y, col=2, lty=3)

lines(lauc1.greben, col=3, lty=1)
lines(lauc1.enet, col=3, lty=3)

lines(lauc2.greben, col=4, lty=1)
lines(lauc2.enet, col=4, lty=3)

lines(lauc3.greben, col=5, lty=1)
lines(lauc3.enet, col=5, lty=3)

# Brier skilll
xlim1.briers <- range(lbriers1.enet$x, lbriers2.enet$x, lbriers3.enet$x, lbriers1.greben$x,
                      lbriers2.greben$x, lbriers3.greben$x, lbriers2.grridge$x)
ylim1.briers <- range(lbriers1.enet$y, lbriers2.enet$y, lbriers3.enet$y, lbriers1.greben$y,
                      lbriers2.greben$y, lbriers3.greben$y, lbriers1.grridge$y,
                      lbriers2.grridge$y, lbriers1.ridge$y)

plot(0, 0, col=2, type="n", ylab="Brier skill score", main="a)",
     xlab="Number of selected variables", xlim=xlim1.briers, ylim=ylim1.briers)

lines(lbriers2.grridge, col=2, lty=1)
abline(h=lbriers1.grridge$y, col=2, lty=2)
abline(h=lbriers1.ridge$y, col=2, lty=3)

lines(lbriers1.greben, col=3, lty=1)
lines(lbriers1.enet, col=3, lty=3)

lines(lbriers2.greben, col=4, lty=1)
lines(lbriers2.enet, col=4, lty=3)

lines(lbriers3.greben, col=5, lty=1)
lines(lbriers3.enet, col=5, lty=3)

leglabels <- c("ridge",
               expression(paste("enet, ", alpha==0.05)),
               expression(paste("enet, ", alpha==0.5)),
               expression(paste("enet, ", alpha==0.95)),
               "group-regularized + selection", "group-regularized",
               "not group-regularized")
legend("bottomright", legend=leglabels, fill=c(2:5, 0, 0, 0),
       lty=c(rep(NA, 4), 1, 2, 3), border=c(rep(1, 4), 0 ,0, 0), merge=TRUE,
       seg.len=1)

# MSE
xlim1.mse <- range(lmse1.enet$x, lmse2.enet$x, lmse3.enet$x, lmse1.greben$x,
                   lmse2.greben$x, lmse3.greben$x, lmse2.grridge$x)
ylim1.mse <- range(lmse1.enet$y, lmse2.enet$y, lmse3.enet$y, lmse1.greben$y,
                   lmse2.greben$y, lmse3.greben$y, lmse1.grridge$y,
                   lmse2.grridge$y, lmse1.ridge$y)

plot(0, 0, col=2, type="n", ylab="MSE", main="a)",
     xlab="Number of selected variables", xlim=xlim1.mse, ylim=ylim1.mse)

lines(lmse2.grridge, col=2, lty=1)
abline(h=lmse1.grridge$y, col=2, lty=2)
abline(h=lmse1.ridge$y, col=2, lty=3)

lines(lmse1.greben, col=3, lty=1)
lines(lmse1.enet, col=3, lty=3)

lines(lmse2.greben, col=4, lty=1)
lines(lmse2.enet, col=4, lty=3)

lines(lmse3.greben, col=5, lty=1)
lines(lmse3.enet, col=5, lty=3)
dev.off()

### estimated penalty parameters
png(paste(path.graph, "grEBEN_test2_res2_penalties.png", sep=""),
    units="in", width=8, height=6, res=120)
par(mfrow=c(2, 2))
boxplot(results2$lambdag[[1]], main="a)",  xlab="", ylab="")
title(ylab=expression(paste(lambda[g], "'")), 
      xlab="Groups in increasing effect size", line=2.5)
boxplot(results2$lambdag[[2]], main="b)",  xlab="", ylab="")
title(ylab=expression(paste(lambda[g], "'")), 
      xlab="Groups in increasing effect size", line=2.5)
boxplot(results2$lambdag[[3]], main="c)",  xlab="", ylab="")
title(ylab=expression(paste(lambda[g], "'")), 
      xlab="Groups in increasing effect size", line=2.5)
boxplot(results2$lambdag[[4]], main="d)",  xlab="", ylab="")
title(ylab=expression(paste(lambda[g], "'")), 
      xlab="Groups in increasing effect size", line=2.5)
par(fig=c(0.15, 0.5, 0.6, 1), new=TRUE)
boxplot(results2$lambdag[[1]], outline=FALSE, xaxt="n")
dev.off()
