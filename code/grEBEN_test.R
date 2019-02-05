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

# ## simulations
# ## simulation 1
# # create data
# n <- 200
# p <- 1000
# G <- 5
# pblock <- 20
# rho <- 0.7
# sigma <- matrix(rho, ncol=pblock, nrow=pblock); diag(sigma) <- 1
# lambda <- 0.02
# alpha <- 0.05
# lambdag <- exp(seq(-1, 1, length.out=G))
# m <- rep(1, n)
# part.greben <- list(groups=rep(1:G, each=p/G))
# part.grridge <- list(groups=CreatePartition(as.factor(part.greben$groups)))
# ntest <- 1000
# 
# methods <- c("ridge+truelambda", "ridge", "GRridge+truelambda", "GRridge", 
#              "GRridge+truelambda+sel", "GRridge+sel",
#              "enet+truelambda", "enet+a=0.05", "enet+a=0.5", "enet+a=0.95", 
#              "grEBEN+truelambda", "grEBEN+a=0.05", "grEBEN+a=0.5","grEBEN+a=0.95")
# nreps <- 50
# auc1 <- briers1 <- mse1 <- rep(list(vector(mode="list", length=nreps)), 
#                                length(methods)) 
# names(auc1) <- names(briers1) <- names(mse1) <- methods
# lambdag1 <- rep(list(matrix(NA, ncol=G, nrow=nreps)), 6)
# names(lambdag1) <- methods[c(3, 4, 11, 12, 13, 14)]
# varbeta1 <- matrix(nrow=nreps, ncol=G)
# # the simulations
# for(r in 1:nreps) {
#   
#   set.seed(200 + r)
#   x <- do.call(cbind, replicate(p/pblock, rmvnorm(n, mean=rep(0, pblock), sigma=sigma), 
#                                 simplify=FALSE))
#   beta <- as.numeric(sapply(1:G, function(g) {
#     renbeta(p/G, 2*n*lambda*alpha*sqrt(lambdag[g]), n*lambda*(1 - alpha)*lambdag[g])}))
#   prob <- as.numeric(exp(x %*% beta)/(1 + exp(x %*% beta)))
#   y <- rbinom(n, 1, prob)
#   optl <- n*lambda*alpha*sum(abs(beta))/sum(beta^2) + 0.5*n*lambda*(1 - alpha)
#   
#   xtest <- do.call(cbind, replicate(p/pblock, rmvnorm(ntest, mean=rep(0, pblock),
#                                                       sigma=sigma), simplify=FALSE))
#   probtest <- as.numeric(exp(xtest %*% beta)/(1 + exp(xtest %*% beta)))
#   ytest <- rbinom(ntest, 1, probtest)
#   
#   fit1.greben <- grEBEN3(x, y, m, partitions=part.greben, alpha=alpha, 
#                          lambda=lambda, psel=TRUE)
#   fit2.greben <- grEBEN3(x, y, m, partitions=part.greben, alpha=0.05, psel=TRUE)
#   fit3.greben <- grEBEN3(x, y, m, partitions=part.greben, alpha=0.5, psel=TRUE)
#   fit4.greben <- grEBEN3(x, y, m, partitions=part.greben, alpha=0.95, psel=TRUE)
#   
#   # number of selected variables
#   psel1.enet <- apply(fit1.greben$beta.nogroups, 2, function(b) {sum(b!=0) - 1})
#   psel2.enet <- apply(fit2.greben$beta.nogroups, 2, function(b) {sum(b!=0) - 1})
#   psel3.enet <- apply(fit3.greben$beta.nogroups, 2, function(b) {sum(b!=0) - 1})
#   psel4.enet <- apply(fit4.greben$beta.nogroups, 2, function(b) {sum(b!=0) - 1})
#   
#   psel1.greben <- apply(fit1.greben$beta, 2, function(b) {sum(b!=0) - 1})
#   psel2.greben <- apply(fit2.greben$beta, 2, function(b) {sum(b!=0) - 1})
#   psel3.greben <- apply(fit3.greben$beta, 2, function(b) {sum(b!=0) - 1})
#   psel4.greben <- apply(fit4.greben$beta, 2, function(b) {sum(b!=0) - 1})
#   
#   psel.all <- unique(c(psel1.enet, psel2.enet, psel3.enet, psel4.enet,
#     psel1.greben, psel2.greben, psel3.greben, psel4.greben))
#   psel.in <- floor(quantile(psel.all[psel.all!=0], 
#     prob=c(seq(0.01, 0.05, 0.01), seq(0.07, 0.25, 0.03), seq(0.3, 1, 0.1))))
#   
#   fit1.grridge <- vector(mode="list", length=length(psel.in))
#   for(s in 1:length(psel.in)) {
#     fit1.grridge[[s]] <- grridge(t(x), y, part.grridge, selectionEN=TRUE, 
#                                  optl=optl, maxsel=psel.in[s])
#   }
#   
#   fit2.grridge <- vector(mode="list", length=length(psel.in))
#   fit2.grridge[[1]] <- grridge(t(x), y, part.grridge, selectionEN=TRUE,
#                                maxsel=psel.in[1])
#   for(s in 2:length(psel.in)) {
#     fit2.grridge[[s]] <- grridge(t(x), y, part.grridge, selectionEN=TRUE,
#                                  optl=fit2.grridge[[1]]$optl, 
#                                  maxsel=psel.in[s])
#   }
#   
#   psel1.grridge <- sapply(1:length(fit1.grridge), function(s) {
#     return(length(fit1.grridge[[s]]$resEN$whichEN))})
#   psel2.grridge <- sapply(1:length(fit2.grridge), function(s) {
#     return(length(fit2.grridge[[s]]$resEN$whichEN))})
#   
#   # estimates
#   est1.ridge <- coef(fit1.grridge[[1]]$predobj$NoGroups, "all")
#   est2.ridge <- coef(fit2.grridge[[1]]$predobj$NoGroups, "all")
#   
#   est1.grridge <- coef(fit1.grridge[[1]]$predobj$GroupRegul, "all")
#   est2.grridge <- coef(fit2.grridge[[1]]$predobj$GroupRegul, "all")
#   
#   est3.grridge <- sapply(fit1.grridge, function(s) {
#     replace(rep(0, p + 1), c(1, s$resEN$whichEN + 1), coef(s$predobj$EN, "all"))})
#   est4.grridge <- sapply(fit2.grridge, function(s) {
#     replace(rep(0, p + 1), c(1, s$resEN$whichEN + 1), coef(s$predobj$EN, "all"))})
#   
#   est1.enet <- fit1.greben$beta.nogroups
#   est2.enet <- fit2.greben$beta.nogroups
#   est3.enet <- fit3.greben$beta.nogroups
#   est4.enet <- fit4.greben$beta.nogroups
#   
#   est1.greben <- fit1.greben$beta
#   est2.greben <- fit2.greben$beta
#   est3.greben <- fit3.greben$beta
#   est4.greben <- fit4.greben$beta
#   
#   # predictions on fit data
#   pred1.ridge <- predict.grridge(fit1.grridge[[1]], t(xtest))[, 1]
#   pred2.ridge <- predict.grridge(fit2.grridge[[1]], t(xtest))[, 1]
#   
#   pred1.grridge <- predict.grridge(fit1.grridge[[1]], t(xtest))[, 2]
#   pred2.grridge <- predict.grridge(fit2.grridge[[1]], t(xtest))[, 2]
#   pred3.grridge <- sapply(fit1.grridge, function(s) {
#     predict.grridge(s, t(xtest))[, 3]})
#   pred4.grridge <- sapply(fit2.grridge, function(s) {
#     predict.grridge(s, t(xtest))[, 3]})
#   
#   pred1.enet <- 1/(1 + exp(-xtest %*% est1.enet[-1, ]))
#   pred2.enet <- 1/(1 + exp(-xtest %*% est2.enet[-1, ]))
#   pred3.enet <- 1/(1 + exp(-xtest %*% est3.enet[-1, ]))
#   pred4.enet <- 1/(1 + exp(-xtest %*% est4.enet[-1, ]))
#   
#   pred1.greben <- 1/(1 + exp(-xtest %*% est1.greben[-1, ]))
#   pred2.greben <- 1/(1 + exp(-xtest %*% est2.greben[-1, ]))
#   pred3.greben <- 1/(1 + exp(-xtest %*% est3.greben[-1, ]))
#   pred4.greben <- 1/(1 + exp(-xtest %*% est4.greben[-1, ]))
#   
#    # AUCs
#   auc.true <- pROC::roc(ytest, probtest)$auc
#   
#   auc1.ridge <- pROC::roc(ytest, pred1.ridge)$auc
#   auc2.ridge <- pROC::roc(ytest, pred2.ridge)$auc
#   
#   auc1.grridge <- pROC::roc(ytest, pred1.grridge)$auc
#   auc2.grridge <- pROC::roc(ytest, pred2.grridge)$auc
#   auc3.grridge <- apply(pred3.grridge, 2, function(r) {pROC::roc(ytest, r)$auc})
#   auc4.grridge <- apply(pred4.grridge, 2, function(r) {pROC::roc(ytest, r)$auc})
#   
#   auc1.enet <- apply(pred1.enet, 2, function(r) {pROC::roc(ytest, r)$auc})
#   auc2.enet <- apply(pred2.enet, 2, function(r) {pROC::roc(ytest, r)$auc})
#   auc3.enet <- apply(pred3.enet, 2, function(r) {pROC::roc(ytest, r)$auc})
#   auc4.enet <- apply(pred4.enet, 2, function(r) {pROC::roc(ytest, r)$auc})
#   
#   auc1.greben <- apply(pred1.greben, 2, function(r) {pROC::roc(ytest, r)$auc})
#   auc2.greben <- apply(pred2.greben, 2, function(r) {pROC::roc(ytest, r)$auc})
#   auc3.greben <- apply(pred3.greben, 2, function(r) {pROC::roc(ytest, r)$auc})
#   auc4.greben <- apply(pred4.greben, 2, function(r) {pROC::roc(ytest, r)$auc})
#   
#   # Brier scores
#   brier.null <- sum((ytest - mean(ytest))^2)
#   briers.true <- 1 - sum((ytest - probtest)^2)/brier.null
#   
#   briers1.ridge <- 1 - sum((ytest - pred1.ridge)^2)/brier.null
#   briers2.ridge <- 1 - sum((ytest - pred2.ridge)^2)/brier.null
# 
#   briers1.grridge <- 1 - sum((ytest - pred1.grridge)^2)/brier.null
#   briers2.grridge <- 1 - sum((ytest - pred2.grridge)^2)/brier.null
#   briers3.grridge <- apply(pred3.grridge, 2, function(pred) {
#     1 - sum((ytest - pred)^2)/brier.null})
#   briers4.grridge <- apply(pred4.grridge, 2, function(pred) {
#     1 - sum((ytest - pred)^2)/brier.null})
#   
#   briers1.enet <- apply(pred1.enet, 2, function(pred) {
#     1 - sum((ytest - pred)^2)/brier.null})
#   briers2.enet <- apply(pred2.enet, 2, function(pred) {
#     1 - sum((ytest - pred)^2)/brier.null})
#   briers3.enet <- apply(pred3.enet, 2, function(pred) {
#     1 - sum((ytest - pred)^2)/brier.null})
#   briers4.enet <- apply(pred4.enet, 2, function(pred) {
#     1 - sum((ytest - pred)^2)/brier.null})
#   
#   briers1.greben <- apply(pred1.greben, 2, function(pred) {
#     1 - sum((ytest - pred)^2)/brier.null})
#   briers2.greben <- apply(pred2.greben, 2, function(pred) {
#     1 - sum((ytest - pred)^2)/brier.null})
#   briers3.greben <- apply(pred3.greben, 2, function(pred) {
#     1 - sum((ytest - pred)^2)/brier.null})
#   briers4.greben <- apply(pred4.greben, 2, function(pred) {
#     1 - sum((ytest - pred)^2)/brier.null})
#   
#   # MSE
#   mse.true <- 0
#   
#   mse1.ridge <- mean((c(0, beta) - est1.ridge)^2)
#   mse2.ridge <- mean((c(0, beta) - est2.ridge)^2)
#   
#   mse1.grridge <- mean((c(0, beta) - est1.grridge)^2)
#   mse2.grridge <- mean((c(0, beta) - est2.grridge)^2)
#   mse3.grridge <- apply(est3.grridge, 2, function(b) {mean((c(0, beta) - b)^2)})
#   mse4.grridge <- apply(est4.grridge, 2, function(b) {mean((c(0, beta) - b)^2)})
#   
#   mse1.enet <- apply(est1.enet, 2, function(b) {mean((c(0, beta) - b)^2)})
#   mse2.enet <- apply(est2.enet, 2, function(b) {mean((c(0, beta) - b)^2)})
#   mse3.enet <- apply(est3.enet, 2, function(b) {mean((c(0, beta) - b)^2)})
#   mse4.enet <- apply(est4.enet, 2, function(b) {mean((c(0, beta) - b)^2)})
#   
#   mse1.greben <- apply(est1.greben, 2, function(b) {mean((c(0, beta) - b)^2)})
#   mse2.greben <- apply(est2.greben, 2, function(b) {mean((c(0, beta) - b)^2)})
#   mse3.greben <- apply(est3.greben, 2, function(b) {mean((c(0, beta) - b)^2)})
#   mse4.greben <- apply(est4.greben, 2, function(b) {mean((c(0, beta) - b)^2)})
#   
#   auc1[[1]][[r]] <- cbind(psel=p, auc=auc1.ridge)
#   auc1[[2]][[r]] <- cbind(psel=p, auc=auc2.ridge)
#   auc1[[3]][[r]] <- cbind(psel=p, auc=auc1.grridge)
#   auc1[[4]][[r]] <- cbind(psel=p, auc=auc2.grridge)
#   auc1[[5]][[r]] <- cbind(psel=psel1.grridge, auc=auc3.grridge)
#   auc1[[6]][[r]] <- cbind(psel=psel2.grridge, auc=auc4.grridge)
#   auc1[[7]][[r]] <- cbind(psel=psel1.enet, auc=auc1.enet)
#   auc1[[8]][[r]] <- cbind(psel=psel2.enet, auc=auc2.enet)
#   auc1[[9]][[r]] <- cbind(psel=psel3.enet, auc=auc3.enet)
#   auc1[[10]][[r]] <- cbind(psel=psel4.enet, auc=auc4.enet)
#   auc1[[11]][[r]] <- cbind(psel=psel1.greben, auc=auc1.greben)
#   auc1[[12]][[r]] <- cbind(psel=psel2.greben, auc=auc2.greben)
#   auc1[[13]][[r]] <- cbind(psel=psel3.greben, auc=auc3.greben)
#   auc1[[14]][[r]] <- cbind(psel=psel4.greben, auc=auc4.greben)
#   
#   briers1[[1]][[r]] <- cbind(psel=p, briers=briers1.ridge)
#   briers1[[2]][[r]] <- cbind(psel=p, briers=briers2.ridge)
#   briers1[[3]][[r]] <- cbind(psel=p, briers=briers1.grridge)
#   briers1[[4]][[r]] <- cbind(psel=p, briers=briers2.grridge)
#   briers1[[5]][[r]] <- cbind(psel=psel1.grridge, briers=briers3.grridge)
#   briers1[[6]][[r]] <- cbind(psel=psel2.grridge, briers=briers4.grridge)
#   briers1[[7]][[r]] <- cbind(psel=psel1.enet, briers=briers1.enet)
#   briers1[[8]][[r]] <- cbind(psel=psel2.enet, briers=briers2.enet)
#   briers1[[9]][[r]] <- cbind(psel=psel3.enet, briers=briers3.enet)
#   briers1[[10]][[r]] <- cbind(psel=psel4.enet, briers=briers4.enet)
#   briers1[[11]][[r]] <- cbind(psel=psel1.greben, briers=briers1.greben)
#   briers1[[12]][[r]] <- cbind(psel=psel2.greben, briers=briers2.greben)
#   briers1[[13]][[r]] <- cbind(psel=psel3.greben, briers=briers3.greben)
#   briers1[[14]][[r]] <- cbind(psel=psel4.greben, briers=briers4.greben)
#   
#   mse1[[1]][[r]] <- cbind(psel=p, mse=mse1.ridge)
#   mse1[[2]][[r]] <- cbind(psel=p, mse=mse2.ridge)
#   mse1[[3]][[r]] <- cbind(psel=p, mse=mse1.grridge)
#   mse1[[4]][[r]] <- cbind(psel=p, mse=mse2.grridge)
#   mse1[[5]][[r]] <- cbind(psel=psel1.grridge, mse=mse3.grridge)
#   mse1[[6]][[r]] <- cbind(psel=psel2.grridge, mse=mse4.grridge)
#   mse1[[7]][[r]] <- cbind(psel=psel1.enet, mse=mse1.enet)
#   mse1[[8]][[r]] <- cbind(psel=psel2.enet, mse=mse2.enet)
#   mse1[[9]][[r]] <- cbind(psel=psel3.enet, mse=mse3.enet)
#   mse1[[10]][[r]] <- cbind(psel=psel4.enet, mse=mse4.enet)
#   mse1[[11]][[r]] <- cbind(psel=psel1.greben, mse=mse1.greben)
#   mse1[[12]][[r]] <- cbind(psel=psel2.greben, mse=mse2.greben)
#   mse1[[13]][[r]] <- cbind(psel=psel3.greben, mse=mse3.greben)
#   mse1[[14]][[r]] <- cbind(psel=psel4.greben, mse=mse4.greben)
#   
#   
#   lambdag1[[1]][r, ] <- fit1.grridge[[1]]$lambdamults$groups
#   lambdag1[[2]][r, ] <- fit2.grridge[[1]]$lambdamults$groups
#   lambdag1[[3]][r, ] <- fit1.greben$lambdag$groups[, fit1.greben$nouteriter + 1]
#   lambdag1[[4]][r, ] <- fit2.greben$lambdag$groups[, fit2.greben$nouteriter + 1]
#   lambdag1[[5]][r, ] <- fit3.greben$lambdag$groups[, fit3.greben$nouteriter + 1]
#   lambdag1[[6]][r, ] <- fit4.greben$lambdag$groups[, fit4.greben$nouteriter + 1]
#   
#   varbeta1[r, ] <- sapply(1:G, function(g) {var(beta[(p*(g - 1)/G + 1):(p*g/G)])})
#   
#   results1 <- list(auc=auc1, briers=briers1, mse=mse1, lambdag=lambdag1,
#                    varbeta=varbeta1)
#   save(results1, file=paste(path.res, "grEBEN_test_res1.Rdata", sep=""))
#   
# }


### plots
load(paste(path.res, "grEBEN_test_res1.Rdata", sep=""))

# preparing data
auc1.enet <- do.call(rbind, results1$auc[[7]])[
  order(do.call(rbind, results1$auc[[7]])[, 1]), ]
auc2.enet <- do.call(rbind, results1$auc[[8]])[
  order(do.call(rbind, results1$auc[[8]])[, 1]), ]
auc3.enet <- do.call(rbind, results1$auc[[9]])[
  order(do.call(rbind, results1$auc[[9]])[, 1]), ]
auc4.enet <- do.call(rbind, results1$auc[[10]])[
  order(do.call(rbind, results1$auc[[10]])[, 1]), ]

auc1.greben <- do.call(rbind, results1$auc[[11]])[
  order(do.call(rbind, results1$auc[[11]])[, 1]), ]
auc2.greben <- do.call(rbind, results1$auc[[12]])[
  order(do.call(rbind, results1$auc[[12]])[, 1]), ]
auc3.greben <- do.call(rbind, results1$auc[[13]])[
  order(do.call(rbind, results1$auc[[13]])[, 1]), ]
auc4.greben <- do.call(rbind, results1$auc[[14]])[
  order(do.call(rbind, results1$auc[[14]])[, 1]), ]

auc1.grridge <- do.call(rbind, results1$auc[[3]])[
  order(do.call(rbind, results1$auc[[3]])[, 1]), ]
auc2.grridge <- do.call(rbind, results1$auc[[4]])[
  order(do.call(rbind, results1$auc[[4]])[, 1]), ]
auc3.grridge <- do.call(rbind, results1$auc[[5]])[
  order(do.call(rbind, results1$auc[[5]])[, 1]), ]
auc4.grridge <- do.call(rbind, results1$auc[[6]])[
  order(do.call(rbind, results1$auc[[6]])[, 1]), ]

auc1.ridge <- do.call(rbind, results1$auc[[1]])[
  order(do.call(rbind, results1$auc[[1]])[, 1]), ]
auc2.ridge <- do.call(rbind, results1$auc[[2]])[
  order(do.call(rbind, results1$auc[[2]])[, 1]), ]

lauc1.enet <- lowess(auc1.enet[, 1], auc1.enet[, 2])
lauc2.enet <- lowess(auc2.enet[, 1], auc2.enet[, 2])
lauc3.enet <- lowess(auc3.enet[, 1], auc3.enet[, 2])
lauc4.enet <- lowess(auc4.enet[, 1], auc4.enet[, 2])

lauc1.greben <- lowess(auc1.greben[, 1], auc1.greben[, 2])
lauc2.greben <- lowess(auc2.greben[, 1], auc2.greben[, 2])
lauc3.greben <- lowess(auc3.greben[, 1], auc3.greben[, 2])
lauc4.greben <- lowess(auc4.greben[, 1], auc4.greben[, 2])

lauc1.grridge <- lowess(auc1.grridge[, 1], auc1.grridge[, 2])
lauc2.grridge <- lowess(auc2.grridge[, 1], auc2.grridge[, 2])
lauc3.grridge <- lowess(auc3.grridge[, 1], auc3.grridge[, 2])
lauc4.grridge <- lowess(auc4.grridge[, 1], auc4.grridge[, 2])

lauc1.ridge <- lowess(auc1.ridge[, 1], auc1.ridge[, 2])
lauc2.ridge <- lowess(auc2.ridge[, 1], auc2.ridge[, 2])

briers1.enet <- do.call(rbind, results1$briers[[7]])[
  order(do.call(rbind, results1$briers[[7]])[, 1]), ]
briers2.enet <- do.call(rbind, results1$briers[[8]])[
  order(do.call(rbind, results1$briers[[8]])[, 1]), ]
briers3.enet <- do.call(rbind, results1$briers[[9]])[
  order(do.call(rbind, results1$briers[[9]])[, 1]), ]
briers4.enet <- do.call(rbind, results1$briers[[10]])[
  order(do.call(rbind, results1$briers[[10]])[, 1]), ]

briers1.greben <- do.call(rbind, results1$briers[[11]])[
  order(do.call(rbind, results1$briers[[11]])[, 1]), ]
briers2.greben <- do.call(rbind, results1$briers[[12]])[
  order(do.call(rbind, results1$briers[[12]])[, 1]), ]
briers3.greben <- do.call(rbind, results1$briers[[13]])[
  order(do.call(rbind, results1$briers[[13]])[, 1]), ]
briers4.greben <- do.call(rbind, results1$briers[[14]])[
  order(do.call(rbind, results1$briers[[14]])[, 1]), ]

briers1.grridge <- do.call(rbind, results1$briers[[3]])[
  order(do.call(rbind, results1$briers[[3]])[, 1]), ]
briers2.grridge <- do.call(rbind, results1$briers[[4]])[
  order(do.call(rbind, results1$briers[[4]])[, 1]), ]
briers3.grridge <- do.call(rbind, results1$briers[[5]])[
  order(do.call(rbind, results1$briers[[5]])[, 1]), ]
briers4.grridge <- do.call(rbind, results1$briers[[6]])[
  order(do.call(rbind, results1$briers[[6]])[, 1]), ]

briers1.ridge <- do.call(rbind, results1$briers[[1]])[
  order(do.call(rbind, results1$briers[[1]])[, 1]), ]
briers2.ridge <- do.call(rbind, results1$briers[[2]])[
  order(do.call(rbind, results1$briers[[2]])[, 1]), ]

lbriers1.enet <- lowess(briers1.enet[, 1], briers1.enet[, 2])
lbriers2.enet <- lowess(briers2.enet[, 1], briers2.enet[, 2])
lbriers3.enet <- lowess(briers3.enet[, 1], briers3.enet[, 2])
lbriers4.enet <- lowess(briers4.enet[, 1], briers4.enet[, 2])

lbriers1.greben <- lowess(briers1.greben[, 1], briers1.greben[, 2])
lbriers2.greben <- lowess(briers2.greben[, 1], briers2.greben[, 2])
lbriers3.greben <- lowess(briers3.greben[, 1], briers3.greben[, 2])
lbriers4.greben <- lowess(briers4.greben[, 1], briers4.greben[, 2])

lbriers1.grridge <- lowess(briers1.grridge[, 1], briers1.grridge[, 2])
lbriers2.grridge <- lowess(briers2.grridge[, 1], briers2.grridge[, 2])
lbriers3.grridge <- lowess(briers3.grridge[, 1], briers3.grridge[, 2])
lbriers4.grridge <- lowess(briers4.grridge[, 1], briers4.grridge[, 2])

lbriers1.ridge <- lowess(briers1.ridge[, 1], briers1.ridge[, 2])
lbriers2.ridge <- lowess(briers2.ridge[, 1], briers2.ridge[, 2])

mse1.enet <- do.call(rbind, results1$mse[[7]])[
  order(do.call(rbind, results1$mse[[7]])[, 1]), ]
mse2.enet <- do.call(rbind, results1$mse[[8]])[
  order(do.call(rbind, results1$mse[[8]])[, 1]), ]
mse3.enet <- do.call(rbind, results1$mse[[9]])[
  order(do.call(rbind, results1$mse[[9]])[, 1]), ]
mse4.enet <- do.call(rbind, results1$mse[[10]])[
  order(do.call(rbind, results1$mse[[10]])[, 1]), ]

mse1.greben <- do.call(rbind, results1$mse[[11]])[
  order(do.call(rbind, results1$mse[[11]])[, 1]), ]
mse2.greben <- do.call(rbind, results1$mse[[12]])[
  order(do.call(rbind, results1$mse[[12]])[, 1]), ]
mse3.greben <- do.call(rbind, results1$mse[[13]])[
  order(do.call(rbind, results1$mse[[13]])[, 1]), ]
mse4.greben <- do.call(rbind, results1$mse[[14]])[
  order(do.call(rbind, results1$mse[[14]])[, 1]), ]

mse1.grridge <- do.call(rbind, results1$mse[[3]])[
  order(do.call(rbind, results1$mse[[3]])[, 1]), ]
mse2.grridge <- do.call(rbind, results1$mse[[4]])[
  order(do.call(rbind, results1$mse[[4]])[, 1]), ]
mse3.grridge <- do.call(rbind, results1$mse[[5]])[
  order(do.call(rbind, results1$mse[[5]])[, 1]), ]
mse4.grridge <- do.call(rbind, results1$mse[[6]])[
  order(do.call(rbind, results1$mse[[6]])[, 1]), ]

mse1.ridge <- do.call(rbind, results1$mse[[1]])[
  order(do.call(rbind, results1$mse[[1]])[, 1]), ]
mse2.ridge <- do.call(rbind, results1$mse[[2]])[
  order(do.call(rbind, results1$mse[[2]])[, 1]), ]

lmse1.enet <- lowess(mse1.enet[, 1], mse1.enet[, 2])
lmse2.enet <- lowess(mse2.enet[, 1], mse2.enet[, 2])
lmse3.enet <- lowess(mse3.enet[, 1], mse3.enet[, 2])
lmse4.enet <- lowess(mse4.enet[, 1], mse4.enet[, 2])

lmse1.greben <- lowess(mse1.greben[, 1], mse1.greben[, 2])
lmse2.greben <- lowess(mse2.greben[, 1], mse2.greben[, 2])
lmse3.greben <- lowess(mse3.greben[, 1], mse3.greben[, 2])
lmse4.greben <- lowess(mse4.greben[, 1], mse4.greben[, 2])

lmse1.grridge <- lowess(mse1.grridge[, 1], mse1.grridge[, 2])
lmse2.grridge <- lowess(mse2.grridge[, 1], mse2.grridge[, 2])
lmse3.grridge <- lowess(mse3.grridge[, 1], mse3.grridge[, 2])
lmse4.grridge <- lowess(mse4.grridge[, 1], mse4.grridge[, 2])

lmse1.ridge <- lowess(mse1.ridge[, 1], mse1.ridge[, 2])
lmse2.ridge <- lowess(mse2.ridge[, 1], mse2.ridge[, 2])

### diagnostics: checking lowess fits
# AUC
xlim2.auc <- range(auc1.enet[, 1], auc2.enet[, 1], auc3.enet[, 1], auc4.enet[, 1],
                   auc1.greben[, 1], auc2.greben[, 1], auc3.greben[, 1],
                   auc4.greben[, 1], auc3.grridge[, 1], auc4.grridge[, 1])
ylim2.auc <- range(auc1.enet[, 2], auc2.enet[, 2], auc3.enet[, 2], auc4.enet[, 2],
                   auc1.greben[, 2], auc2.greben[, 2], auc3.greben[, 2],
                   auc4.greben[, 2], auc3.grridge[, 2], auc4.grridge[, 2])

png(paste(path.graph, "grEBEN_test_res1_lowess_auc.png", sep=""),
    units="in", width=12, height=6, res=120)
par(mfrow=c(2, 5))
plot(auc3.grridge[, 1], auc3.grridge[, 2], col=1, ylab="AUC", main="a)",
     xlab="Number of selected variables", xlim=xlim2.auc, ylim=ylim2.auc)
lines(lauc3.grridge, col=2)

plot(auc1.greben[, 1], auc1.greben[, 2], col=1, ylab="AUC", main="b)",
     xlab="Number of selected variables", xlim=xlim2.auc, ylim=ylim2.auc)
lines(lauc1.greben, col=4)

plot(auc2.greben[, 1], auc2.greben[, 2], col=1, ylab="AUC", main="c)",
     xlab="Number of selected variables", xlim=xlim2.auc, ylim=ylim2.auc)
lines(lauc2.greben, col=5)

plot(auc3.greben[, 1], auc4.greben[, 2], col=1, ylab="AUC", main="d)",
     xlab="Number of selected variables", xlim=xlim2.auc, ylim=ylim2.auc)
lines(lauc3.greben, col=6)

plot(auc4.greben[, 1], auc4.greben[, 2], col=1, ylab="AUC", main="e)",
     xlab="Number of selected variables", xlim=xlim2.auc, ylim=ylim2.auc)
lines(lauc4.greben, col=7)

plot(auc4.grridge[, 1], auc4.grridge[, 2], col=1, ylab="AUC", main="f)",
     xlab="Number of selected variables", xlim=xlim2.auc, ylim=ylim2.auc)
lines(lauc4.grridge, col=3)

plot(auc1.enet[, 1], auc1.enet[, 2], col=1, ylab="AUC", main="g)",
     xlab="Number of selected variables", xlim=xlim2.auc, ylim=ylim2.auc)
lines(lauc1.enet, col=4)

plot(auc2.enet[, 1], auc2.enet[, 2], col=1, ylab="AUC", main="h)",
     xlab="Number of selected variables", xlim=xlim2.auc, ylim=ylim2.auc)
lines(lauc2.enet, col=5)

plot(auc3.enet[, 1], auc4.enet[, 2], col=1, ylab="AUC", main="i)",
     xlab="Number of selected variables", xlim=xlim2.auc, ylim=ylim2.auc)
lines(lauc3.enet, col=6)

plot(auc4.enet[, 1], auc4.enet[, 2], col=1, ylab="AUC", main="j)",
     xlab="Number of selected variables", xlim=xlim2.auc, ylim=ylim2.auc)
lines(lauc4.enet, col=7)
dev.off()

# briers
xlim2.briers <- range(briers1.enet[, 1], briers2.enet[, 1], briers3.enet[, 1], briers4.enet[, 1],
                   briers1.greben[, 1], briers2.greben[, 1], briers3.greben[, 1],
                   briers4.greben[, 1], briers3.grridge[, 1], briers4.grridge[, 1])
ylim2.briers <- range(briers1.enet[, 2], briers2.enet[, 2], briers3.enet[, 2], briers4.enet[, 2],
                   briers1.greben[, 2], briers2.greben[, 2], briers3.greben[, 2],
                   briers4.greben[, 2], briers3.grridge[, 2], briers4.grridge[, 2])

png(paste(path.graph, "grEBEN_test_res1_lowess_briers.png", sep=""),
    units="in", width=12, height=6, res=120)
par(mfrow=c(2, 5))
plot(briers3.grridge[, 1], briers3.grridge[, 2], col=1, ylab="Brier skill score", main="a)",
     xlab="Number of selected variables", xlim=xlim2.briers, ylim=ylim2.briers)
lines(lbriers3.grridge, col=2)

plot(briers1.greben[, 1], briers1.greben[, 2], col=1, ylab="Brier skill score", main="b)",
     xlab="Number of selected variables", xlim=xlim2.briers, ylim=ylim2.briers)
lines(lbriers1.greben, col=4)

plot(briers2.greben[, 1], briers2.greben[, 2], col=1, ylab="Brier skill score", main="c)",
     xlab="Number of selected variables", xlim=xlim2.briers, ylim=ylim2.briers)
lines(lbriers2.greben, col=5)

plot(briers3.greben[, 1], briers4.greben[, 2], col=1, ylab="Brier skill score", main="d)",
     xlab="Number of selected variables", xlim=xlim2.briers, ylim=ylim2.briers)
lines(lbriers3.greben, col=6)

plot(briers4.greben[, 1], briers4.greben[, 2], col=1, ylab="Brier skill score", main="e)",
     xlab="Number of selected variables", xlim=xlim2.briers, ylim=ylim2.briers)
lines(lbriers4.greben, col=7)

plot(briers4.grridge[, 1], briers4.grridge[, 2], col=1, ylab="Brier skill score", main="f)",
     xlab="Number of selected variables", xlim=xlim2.briers, ylim=ylim2.briers)
lines(lbriers4.grridge, col=3)

plot(briers1.enet[, 1], briers1.enet[, 2], col=1, ylab="Brier skill score", main="g)",
     xlab="Number of selected variables", xlim=xlim2.briers, ylim=ylim2.briers)
lines(lbriers1.enet, col=4)

plot(briers2.enet[, 1], briers2.enet[, 2], col=1, ylab="Brier skill score", main="h)",
     xlab="Number of selected variables", xlim=xlim2.briers, ylim=ylim2.briers)
lines(lbriers2.enet, col=5)

plot(briers3.enet[, 1], briers4.enet[, 2], col=1, ylab="Brier skill score", main="i)",
     xlab="Number of selected variables", xlim=xlim2.briers, ylim=ylim2.briers)
lines(lbriers3.enet, col=6)

plot(briers4.enet[, 1], briers4.enet[, 2], col=1, ylab="Brier skill score", main="j)",
     xlab="Number of selected variables", xlim=xlim2.briers, ylim=ylim2.briers)
lines(lbriers4.enet, col=7)
dev.off()

# mse
xlim2.mse <- range(mse1.enet[, 1], mse2.enet[, 1], mse3.enet[, 1], mse4.enet[, 1],
                   mse1.greben[, 1], mse2.greben[, 1], mse3.greben[, 1],
                   mse4.greben[, 1], mse3.grridge[, 1], mse4.grridge[, 1])
ylim2.mse <- range(mse1.enet[, 2], mse2.enet[, 2], mse3.enet[, 2], mse4.enet[, 2],
                   mse1.greben[, 2], mse2.greben[, 2], mse3.greben[, 2],
                   mse4.greben[, 2], mse3.grridge[, 2], mse4.grridge[, 2])

png(paste(path.graph, "grEBEN_test_res1_lowess_mse.png", sep=""),
    units="in", width=12, height=6, res=120)
par(mfrow=c(2, 5))
plot(mse3.grridge[, 1], mse3.grridge[, 2], col=1, ylab="MSE", main="a)",
     xlab="Number of selected variables", xlim=xlim2.mse, ylim=ylim2.mse)
lines(lmse3.grridge, col=2)

plot(mse1.greben[, 1], mse1.greben[, 2], col=1, ylab="MSE", main="b)",
     xlab="Number of selected variables", xlim=xlim2.mse, ylim=ylim2.mse)
lines(lmse1.greben, col=4)

plot(mse2.greben[, 1], mse2.greben[, 2], col=1, ylab="MSE", main="c)",
     xlab="Number of selected variables", xlim=xlim2.mse, ylim=ylim2.mse)
lines(lmse2.greben, col=5)

plot(mse3.greben[, 1], mse4.greben[, 2], col=1, ylab="MSE", main="d)",
     xlab="Number of selected variables", xlim=xlim2.mse, ylim=ylim2.mse)
lines(lmse3.greben, col=6)

plot(mse4.greben[, 1], mse4.greben[, 2], col=1, ylab="MSE", main="e)",
     xlab="Number of selected variables", xlim=xlim2.mse, ylim=ylim2.mse)
lines(lmse4.greben, col=7)

plot(mse4.grridge[, 1], mse4.grridge[, 2], col=1, ylab="MSE", main="f)",
     xlab="Number of selected variables", xlim=xlim2.mse, ylim=ylim2.mse)
lines(lmse4.grridge, col=3)

plot(mse1.enet[, 1], mse1.enet[, 2], col=1, ylab="MSE", main="g)",
     xlab="Number of selected variables", xlim=xlim2.mse, ylim=ylim2.mse)
lines(lmse1.enet, col=4)

plot(mse2.enet[, 1], mse2.enet[, 2], col=1, ylab="MSE", main="h)",
     xlab="Number of selected variables", xlim=xlim2.mse, ylim=ylim2.mse)
lines(lmse2.enet, col=5)

plot(mse3.enet[, 1], mse4.enet[, 2], col=1, ylab="MSE", main="i)",
     xlab="Number of selected variables", xlim=xlim2.mse, ylim=ylim2.mse)
lines(lmse3.enet, col=6)

plot(mse4.enet[, 1], mse4.enet[, 2], col=1, ylab="MSE", main="j)",
     xlab="Number of selected variables", xlim=xlim2.mse, ylim=ylim2.mse)
lines(lmse4.enet, col=7)
dev.off()



### Performance measures
# AUC
png(paste(path.graph, "grEBEN_test_res1_performance.png", sep=""),
    units="in", width=12, height=4, res=120)
par(mfrow=c(1, 3))
xlim1.auc <- range(lauc1.enet$x, lauc2.enet$x, lauc3.enet$x, lauc4.enet$x,
                   lauc1.greben$x, lauc2.greben$x, lauc3.greben$x,
                   lauc4.greben$x, lauc3.grridge$x, lauc4.grridge$x)
ylim1.auc <- range(lauc1.enet$y, lauc2.enet$y, lauc3.enet$y, lauc4.enet$y,
                   lauc1.greben$y, lauc2.greben$y, lauc3.greben$y,
                   lauc4.greben$y, lauc1.grridge$y, lauc2.grridge$y,
                   lauc3.grridge$y, lauc4.grridge$y, lauc1.ridge$y,
                   lauc2.ridge$y)

plot(0, 0, col=2, type="n", ylab="AUC", main="a)",
     xlab="Number of selected variables", xlim=xlim1.auc, ylim=ylim1.auc)

lines(lauc3.grridge, col=2, lty=1)
abline(h=lauc1.grridge$y, col=2, lty=2)
abline(h=lauc1.ridge$y, col=2, lty=3)

lines(lauc4.grridge, col=3, lty=1)
abline(h=lauc2.grridge$y, col=3, lty=2)
abline(h=lauc2.ridge$y, col=3, lty=3)

lines(lauc1.greben, col=4, lty=1)
lines(lauc1.enet, col=4, lty=3)

lines(lauc2.greben, col=5, lty=1)
lines(lauc2.enet, col=5, lty=3)

lines(lauc3.greben, col=6, lty=1)
lines(lauc3.enet, col=6, lty=3)

lines(lauc4.greben, col=7, lty=1)
lines(lauc4.enet, col=7, lty=3)

# Brier skilll
xlim1.briers <- range(lbriers1.enet$x, lbriers2.enet$x, lbriers3.enet$x,
                      lbriers4.enet$x, lbriers1.greben$x, lbriers2.greben$x,
                      lbriers3.greben$x, lbriers4.greben$x, lbriers3.grridge$x,
                      lbriers4.grridge$x)
ylim1.briers <- range(lbriers1.enet$y, lbriers2.enet$y, lbriers3.enet$y,
                      lbriers4.enet$y, lbriers1.greben$y, lbriers2.greben$y,
                      lbriers3.greben$y, lbriers4.greben$y, lbriers1.grridge$y,
                      lbriers2.grridge$y, lbriers3.grridge$y, lbriers4.grridge$y,
                      lbriers1.ridge$y, lbriers2.ridge$y)

plot(0, 0, col=2, type="n", ylab="Brier skill score", main="b)",
     xlab="Number of selected variables", xlim=xlim1.briers, ylim=ylim1.briers)

lines(lbriers3.grridge, col=2, lty=1)
abline(h=lbriers1.grridge$y, col=2, lty=2)
abline(h=lbriers1.ridge$y, col=2, lty=3)

lines(lbriers4.grridge, col=3, lty=1)
abline(h=lbriers2.grridge$y, col=3, lty=2)
abline(h=lbriers2.ridge$y, col=3, lty=3)

lines(lbriers1.greben, col=4, lty=1)
lines(lbriers1.enet, col=4, lty=3)

lines(lbriers2.greben, col=5, lty=1)
lines(lbriers2.enet, col=5, lty=3)

lines(lbriers3.greben, col=6, lty=1)
lines(lbriers3.enet, col=6, lty=3)

lines(lbriers4.greben, col=7, lty=1)
lines(lbriers4.enet, col=7, lty=3)

leglabels <- c(expression(paste("ridge, 'true' ", lambda)),
               "ridge", expression(paste("enet, true ", lambda)),
               expression(paste("enet, ", alpha==0.05)),
               expression(paste("enet, ", alpha==0.5)),
               expression(paste("enet, ", alpha==0.95)),
               "group-regularized + selection", "group-regularized",
               "not group-regularized")
legend("bottomright", legend=leglabels, fill=c(2:7, 0, 0, 0),
       lty=c(rep(NA, 6), 1, 2, 3), border=c(rep(1, 6), 0 ,0, 0), merge=TRUE,
       seg.len=1)

# MSE
xlim1.mse <- range(lmse1.enet$x, lmse2.enet$x, lmse3.enet$x,
                      lmse4.enet$x, lmse1.greben$x, lmse2.greben$x,
                      lmse3.greben$x, lmse4.greben$x, lmse3.grridge$x,
                      lmse4.grridge$x)
ylim1.mse <- range(lmse1.enet$y, lmse2.enet$y, lmse3.enet$y,
                      lmse4.enet$y, lmse1.greben$y, lmse2.greben$y,
                      lmse3.greben$y, lmse4.greben$y, lmse1.grridge$y,
                      lmse2.grridge$y, lmse3.grridge$y, lmse4.grridge$y,
                      lmse1.ridge$y, lmse2.ridge$y)

plot(0, 0, col=2, type="n", ylab="MSE", main="c)",
     xlab="Number of selected variables", xlim=xlim1.mse, ylim=ylim1.mse)

lines(lmse3.grridge, col=2, lty=1)
abline(h=lmse1.grridge$y, col=2, lty=2)
abline(h=lmse1.ridge$y, col=2, lty=3)

lines(lmse4.grridge, col=3, lty=1)
abline(h=lmse2.grridge$y, col=3, lty=2)
abline(h=lmse2.ridge$y, col=3, lty=3)

lines(lmse1.greben, col=4, lty=1)
lines(lmse1.enet, col=4, lty=3)

lines(lmse2.greben, col=5, lty=1)
lines(lmse2.enet, col=5, lty=3)

lines(lmse3.greben, col=6, lty=1)
lines(lmse3.enet, col=6, lty=3)

lines(lmse4.greben, col=7, lty=1)
lines(lmse4.enet, col=7, lty=3)
dev.off()

## penalty parameters
lambdag <- exp(seq(-1, 1, length.out=ncol(results1$lambdag[[1]])))
png(paste(path.graph, "grEBEN_test_res1_penalties.png", sep=""),
    units="in", width=8, height=6, res=120)
par(mfrow=c(2, 3))
boxplot(results1$lambdag[[1]], ylim=range(results1$lambdag[[1]], lambdag), 
        xlab="", ylab="", main="a)")
title(ylab=expression(paste(lambda[g], "'")), 
      xlab="Groups in decreasing effect size", line=2.5)
points(1:length(lambdag), lambdag, col=2, pch="x")
boxplot(results1$lambdag[[2]], ylim=range(results1$lambdag[[2]], lambdag), 
        xlab="", ylab="", main="b)")
title(ylab=expression(paste(lambda[g], "'")), 
      xlab="Groups in decreasing effect size", line=2.5)
points(1:length(lambdag), lambdag, col=2, pch="x")
boxplot(results1$lambdag[[3]], ylim=range(results1$lambdag[[3]], lambdag), 
        xlab="", ylab="", main="c)")
title(ylab=expression(paste(lambda[g], "'")), 
      xlab="Groups in decreasing effect size", line=2.5)
points(1:length(lambdag), lambdag, col=2, pch="x")
boxplot(results1$lambdag[[4]], ylim=range(results1$lambdag[[4]], lambdag), 
        xlab="", ylab="", main="d)")
title(ylab=expression(paste(lambda[g], "'")), 
      xlab="Groups in decreasing effect size", line=2.5)
points(1:length(lambdag), lambdag, col=2, pch="x")
boxplot(results1$lambdag[[5]], ylim=range(results1$lambdag[[5]], lambdag), 
        xlab="", ylab="", main="e)")
title(ylab=expression(paste(lambda[g], "'")), 
      xlab="Groups in decreasing effect size", line=2.5)
points(1:length(lambdag), lambdag, col=2, pch="x")
boxplot(results1$lambdag[[6]], ylim=range(results1$lambdag[[6]], lambdag), 
        xlab="", ylab="", main="f)")
title(ylab=expression(paste(lambda[g], "'")), 
      xlab="Groups in decreasing effect size", line=2.5)
points(1:length(lambdag), lambdag, col=2, pch="x")
dev.off()




## simulation 2
# create data
n <- 200
p <- 1000
G <- 5
pblock <- 20
rho <- 0.7
sigma <- matrix(rho, ncol=pblock, nrow=pblock); diag(sigma) <- 1
lambda <- 0.02
alpha <- 0.5
lambdag <- exp(seq(-1, 1, length.out=G))
m <- rep(1, n)
part.greben <- list(groups=rep(1:G, each=p/G))
part.grridge <- list(groups=CreatePartition(as.factor(part.greben$groups)))
ntest <- 1000

methods <- c("ridge+truelambda", "ridge", "GRridge+truelambda", "GRridge",
             "GRridge+truelambda+sel", "GRridge+sel",
             "enet+truelambda", "enet+a=0.05", "enet+a=0.5", "enet+a=0.95",
             "grEBEN+truelambda", "grEBEN+a=0.05", "grEBEN+a=0.5","grEBEN+a=0.95")
nreps <- 50
auc2 <- briers2 <- mse2 <- rep(list(vector(mode="list", length=nreps)),
                               length(methods))
names(auc2) <- names(briers2) <- names(mse2) <- methods
lambdag2 <- rep(list(matrix(NA, ncol=G, nrow=nreps)), 6)
names(lambdag2) <- methods[c(3, 4, 11, 12, 13, 14)]
varbeta2 <- matrix(nrow=nreps, ncol=G)
# the simulations
for(r in 1:nreps) {

  set.seed(200 + r)
  x <- do.call(cbind, replicate(p/pblock, rmvnorm(n, mean=rep(0, pblock), sigma=sigma),
                                simplify=FALSE))
  beta <- as.numeric(sapply(1:G, function(g) {
    renbeta(p/G, 2*n*lambda*alpha*sqrt(lambdag[g]), n*lambda*(1 - alpha)*lambdag[g])}))
  prob <- as.numeric(exp(x %*% beta)/(1 + exp(x %*% beta)))
  y <- rbinom(n, 1, prob)
  optl <- n*lambda*alpha*sum(abs(beta))/sum(beta^2) + 0.5*n*lambda*(1 - alpha)

  xtest <- do.call(cbind, replicate(p/pblock, rmvnorm(ntest, mean=rep(0, pblock),
                                                      sigma=sigma), simplify=FALSE))
  probtest <- as.numeric(exp(xtest %*% beta)/(1 + exp(xtest %*% beta)))
  ytest <- rbinom(ntest, 1, probtest)

  fit1.greben <- grEBEN3(x, y, m, partitions=part.greben, alpha=alpha,
                         lambda=lambda, psel=TRUE)
  fit2.greben <- grEBEN3(x, y, m, partitions=part.greben, alpha=0.05, psel=TRUE)
  fit3.greben <- grEBEN3(x, y, m, partitions=part.greben, alpha=0.5, psel=TRUE)
  fit4.greben <- grEBEN3(x, y, m, partitions=part.greben, alpha=0.95, psel=TRUE)

  # number of selected variables
  psel1.enet <- apply(fit1.greben$beta.nogroups, 2, function(b) {sum(b!=0) - 1})
  psel2.enet <- apply(fit2.greben$beta.nogroups, 2, function(b) {sum(b!=0) - 1})
  psel3.enet <- apply(fit3.greben$beta.nogroups, 2, function(b) {sum(b!=0) - 1})
  psel4.enet <- apply(fit4.greben$beta.nogroups, 2, function(b) {sum(b!=0) - 1})

  psel1.greben <- apply(fit1.greben$beta, 2, function(b) {sum(b!=0) - 1})
  psel2.greben <- apply(fit2.greben$beta, 2, function(b) {sum(b!=0) - 1})
  psel3.greben <- apply(fit3.greben$beta, 2, function(b) {sum(b!=0) - 1})
  psel4.greben <- apply(fit4.greben$beta, 2, function(b) {sum(b!=0) - 1})

  psel.all <- unique(c(psel1.enet, psel2.enet, psel3.enet, psel4.enet,
    psel1.greben, psel2.greben, psel3.greben, psel4.greben))
  psel.in <- floor(quantile(psel.all[psel.all!=0],
    prob=c(seq(0.01, 0.05, 0.01), seq(0.07, 0.25, 0.03), seq(0.3, 1, 0.1))))

  fit1.grridge <- vector(mode="list", length=length(psel.in))
  for(s in 1:length(psel.in)) {
    fit1.grridge[[s]] <- grridge(t(x), y, part.grridge, selectionEN=TRUE,
                                 optl=optl, maxsel=psel.in[s])
  }

  fit2.grridge <- vector(mode="list", length=length(psel.in))
  fit2.grridge[[1]] <- grridge(t(x), y, part.grridge, selectionEN=TRUE,
                               maxsel=psel.in[1])
  for(s in 2:length(psel.in)) {
    fit2.grridge[[s]] <- grridge(t(x), y, part.grridge, selectionEN=TRUE,
                                 optl=fit2.grridge[[1]]$optl,
                                 maxsel=psel.in[s])
  }

  psel1.grridge <- sapply(1:length(fit1.grridge), function(s) {
    return(length(fit1.grridge[[s]]$resEN$whichEN))})
  psel2.grridge <- sapply(1:length(fit2.grridge), function(s) {
    return(length(fit2.grridge[[s]]$resEN$whichEN))})

  # estimates
  est1.ridge <- coef(fit1.grridge[[1]]$predobj$NoGroups, "all")
  est2.ridge <- coef(fit2.grridge[[1]]$predobj$NoGroups, "all")

  est1.grridge <- coef(fit1.grridge[[1]]$predobj$GroupRegul, "all")
  est2.grridge <- coef(fit2.grridge[[1]]$predobj$GroupRegul, "all")

  est3.grridge <- sapply(fit1.grridge, function(s) {
    replace(rep(0, p + 1), c(1, s$resEN$whichEN + 1), coef(s$predobj$EN, "all"))})
  est4.grridge <- sapply(fit2.grridge, function(s) {
    replace(rep(0, p + 1), c(1, s$resEN$whichEN + 1), coef(s$predobj$EN, "all"))})

  est1.enet <- fit1.greben$beta.nogroups
  est2.enet <- fit2.greben$beta.nogroups
  est3.enet <- fit3.greben$beta.nogroups
  est4.enet <- fit4.greben$beta.nogroups

  est1.greben <- fit1.greben$beta
  est2.greben <- fit2.greben$beta
  est3.greben <- fit3.greben$beta
  est4.greben <- fit4.greben$beta

  # predictions on fit data
  pred1.ridge <- predict.grridge(fit1.grridge[[1]], t(xtest))[, 1]
  pred2.ridge <- predict.grridge(fit2.grridge[[1]], t(xtest))[, 1]

  pred1.grridge <- predict.grridge(fit1.grridge[[1]], t(xtest))[, 2]
  pred2.grridge <- predict.grridge(fit2.grridge[[1]], t(xtest))[, 2]
  pred3.grridge <- sapply(fit1.grridge, function(s) {
    predict.grridge(s, t(xtest))[, 3]})
  pred4.grridge <- sapply(fit2.grridge, function(s) {
    predict.grridge(s, t(xtest))[, 3]})

  pred1.enet <- 1/(1 + exp(-xtest %*% est1.enet[-1, ]))
  pred2.enet <- 1/(1 + exp(-xtest %*% est2.enet[-1, ]))
  pred3.enet <- 1/(1 + exp(-xtest %*% est3.enet[-1, ]))
  pred4.enet <- 1/(1 + exp(-xtest %*% est4.enet[-1, ]))

  pred1.greben <- 1/(1 + exp(-xtest %*% est1.greben[-1, ]))
  pred2.greben <- 1/(1 + exp(-xtest %*% est2.greben[-1, ]))
  pred3.greben <- 1/(1 + exp(-xtest %*% est3.greben[-1, ]))
  pred4.greben <- 1/(1 + exp(-xtest %*% est4.greben[-1, ]))

   # AUCs
  auc.true <- pROC::roc(ytest, probtest)$auc

  auc1.ridge <- pROC::roc(ytest, pred1.ridge)$auc
  auc2.ridge <- pROC::roc(ytest, pred2.ridge)$auc

  auc1.grridge <- pROC::roc(ytest, pred1.grridge)$auc
  auc2.grridge <- pROC::roc(ytest, pred2.grridge)$auc
  auc3.grridge <- apply(pred3.grridge, 2, function(r) {pROC::roc(ytest, r)$auc})
  auc4.grridge <- apply(pred4.grridge, 2, function(r) {pROC::roc(ytest, r)$auc})

  auc1.enet <- apply(pred1.enet, 2, function(r) {pROC::roc(ytest, r)$auc})
  auc2.enet <- apply(pred2.enet, 2, function(r) {pROC::roc(ytest, r)$auc})
  auc3.enet <- apply(pred3.enet, 2, function(r) {pROC::roc(ytest, r)$auc})
  auc4.enet <- apply(pred4.enet, 2, function(r) {pROC::roc(ytest, r)$auc})

  auc1.greben <- apply(pred1.greben, 2, function(r) {pROC::roc(ytest, r)$auc})
  auc2.greben <- apply(pred2.greben, 2, function(r) {pROC::roc(ytest, r)$auc})
  auc3.greben <- apply(pred3.greben, 2, function(r) {pROC::roc(ytest, r)$auc})
  auc4.greben <- apply(pred4.greben, 2, function(r) {pROC::roc(ytest, r)$auc})

  # Brier scores
  brier.null <- sum((ytest - mean(ytest))^2)
  briers.true <- 1 - sum((ytest - probtest)^2)/brier.null

  briers1.ridge <- 1 - sum((ytest - pred1.ridge)^2)/brier.null
  briers2.ridge <- 1 - sum((ytest - pred2.ridge)^2)/brier.null

  briers1.grridge <- 1 - sum((ytest - pred1.grridge)^2)/brier.null
  briers2.grridge <- 1 - sum((ytest - pred2.grridge)^2)/brier.null
  briers3.grridge <- apply(pred3.grridge, 2, function(pred) {
    1 - sum((ytest - pred)^2)/brier.null})
  briers4.grridge <- apply(pred4.grridge, 2, function(pred) {
    1 - sum((ytest - pred)^2)/brier.null})

  briers1.enet <- apply(pred1.enet, 2, function(pred) {
    1 - sum((ytest - pred)^2)/brier.null})
  briers2.enet <- apply(pred2.enet, 2, function(pred) {
    1 - sum((ytest - pred)^2)/brier.null})
  briers3.enet <- apply(pred3.enet, 2, function(pred) {
    1 - sum((ytest - pred)^2)/brier.null})
  briers4.enet <- apply(pred4.enet, 2, function(pred) {
    1 - sum((ytest - pred)^2)/brier.null})

  briers1.greben <- apply(pred1.greben, 2, function(pred) {
    1 - sum((ytest - pred)^2)/brier.null})
  briers2.greben <- apply(pred2.greben, 2, function(pred) {
    1 - sum((ytest - pred)^2)/brier.null})
  briers3.greben <- apply(pred3.greben, 2, function(pred) {
    1 - sum((ytest - pred)^2)/brier.null})
  briers4.greben <- apply(pred4.greben, 2, function(pred) {
    1 - sum((ytest - pred)^2)/brier.null})

  # MSE
  mse.true <- 0

  mse1.ridge <- mean((c(0, beta) - est1.ridge)^2)
  mse2.ridge <- mean((c(0, beta) - est2.ridge)^2)

  mse1.grridge <- mean((c(0, beta) - est1.grridge)^2)
  mse2.grridge <- mean((c(0, beta) - est2.grridge)^2)
  mse3.grridge <- apply(est3.grridge, 2, function(b) {mean((c(0, beta) - b)^2)})
  mse4.grridge <- apply(est4.grridge, 2, function(b) {mean((c(0, beta) - b)^2)})

  mse1.enet <- apply(est1.enet, 2, function(b) {mean((c(0, beta) - b)^2)})
  mse2.enet <- apply(est2.enet, 2, function(b) {mean((c(0, beta) - b)^2)})
  mse3.enet <- apply(est3.enet, 2, function(b) {mean((c(0, beta) - b)^2)})
  mse4.enet <- apply(est4.enet, 2, function(b) {mean((c(0, beta) - b)^2)})

  mse1.greben <- apply(est1.greben, 2, function(b) {mean((c(0, beta) - b)^2)})
  mse2.greben <- apply(est2.greben, 2, function(b) {mean((c(0, beta) - b)^2)})
  mse3.greben <- apply(est3.greben, 2, function(b) {mean((c(0, beta) - b)^2)})
  mse4.greben <- apply(est4.greben, 2, function(b) {mean((c(0, beta) - b)^2)})

  auc2[[1]][[r]] <- cbind(psel=p, auc=auc1.ridge)
  auc2[[2]][[r]] <- cbind(psel=p, auc=auc2.ridge)
  auc2[[3]][[r]] <- cbind(psel=p, auc=auc1.grridge)
  auc2[[4]][[r]] <- cbind(psel=p, auc=auc2.grridge)
  auc2[[5]][[r]] <- cbind(psel=psel1.grridge, auc=auc3.grridge)
  auc2[[6]][[r]] <- cbind(psel=psel2.grridge, auc=auc4.grridge)
  auc2[[7]][[r]] <- cbind(psel=psel1.enet, auc=auc1.enet)
  auc2[[8]][[r]] <- cbind(psel=psel2.enet, auc=auc2.enet)
  auc2[[9]][[r]] <- cbind(psel=psel3.enet, auc=auc3.enet)
  auc2[[10]][[r]] <- cbind(psel=psel4.enet, auc=auc4.enet)
  auc2[[11]][[r]] <- cbind(psel=psel1.greben, auc=auc1.greben)
  auc2[[12]][[r]] <- cbind(psel=psel2.greben, auc=auc2.greben)
  auc2[[13]][[r]] <- cbind(psel=psel3.greben, auc=auc3.greben)
  auc2[[14]][[r]] <- cbind(psel=psel4.greben, auc=auc4.greben)

  briers2[[1]][[r]] <- cbind(psel=p, briers=briers1.ridge)
  briers2[[2]][[r]] <- cbind(psel=p, briers=briers2.ridge)
  briers2[[3]][[r]] <- cbind(psel=p, briers=briers1.grridge)
  briers2[[4]][[r]] <- cbind(psel=p, briers=briers2.grridge)
  briers2[[5]][[r]] <- cbind(psel=psel1.grridge, briers=briers3.grridge)
  briers2[[6]][[r]] <- cbind(psel=psel2.grridge, briers=briers4.grridge)
  briers2[[7]][[r]] <- cbind(psel=psel1.enet, briers=briers1.enet)
  briers2[[8]][[r]] <- cbind(psel=psel2.enet, briers=briers2.enet)
  briers2[[9]][[r]] <- cbind(psel=psel3.enet, briers=briers3.enet)
  briers2[[10]][[r]] <- cbind(psel=psel4.enet, briers=briers4.enet)
  briers2[[11]][[r]] <- cbind(psel=psel1.greben, briers=briers1.greben)
  briers2[[12]][[r]] <- cbind(psel=psel2.greben, briers=briers2.greben)
  briers2[[13]][[r]] <- cbind(psel=psel3.greben, briers=briers3.greben)
  briers2[[14]][[r]] <- cbind(psel=psel4.greben, briers=briers4.greben)

  mse2[[1]][[r]] <- cbind(psel=p, mse=mse1.ridge)
  mse2[[2]][[r]] <- cbind(psel=p, mse=mse2.ridge)
  mse2[[3]][[r]] <- cbind(psel=p, mse=mse1.grridge)
  mse2[[4]][[r]] <- cbind(psel=p, mse=mse2.grridge)
  mse2[[5]][[r]] <- cbind(psel=psel1.grridge, mse=mse3.grridge)
  mse2[[6]][[r]] <- cbind(psel=psel2.grridge, mse=mse4.grridge)
  mse2[[7]][[r]] <- cbind(psel=psel1.enet, mse=mse1.enet)
  mse2[[8]][[r]] <- cbind(psel=psel2.enet, mse=mse2.enet)
  mse2[[9]][[r]] <- cbind(psel=psel3.enet, mse=mse3.enet)
  mse2[[10]][[r]] <- cbind(psel=psel4.enet, mse=mse4.enet)
  mse2[[11]][[r]] <- cbind(psel=psel1.greben, mse=mse1.greben)
  mse2[[12]][[r]] <- cbind(psel=psel2.greben, mse=mse2.greben)
  mse2[[13]][[r]] <- cbind(psel=psel3.greben, mse=mse3.greben)
  mse2[[14]][[r]] <- cbind(psel=psel4.greben, mse=mse4.greben)


  lambdag2[[1]][r, ] <- fit1.grridge[[1]]$lambdamults$groups
  lambdag2[[2]][r, ] <- fit2.grridge[[1]]$lambdamults$groups
  lambdag2[[3]][r, ] <- fit1.greben$lambdag$groups[, fit1.greben$nouteriter + 1]
  lambdag2[[4]][r, ] <- fit2.greben$lambdag$groups[, fit2.greben$nouteriter + 1]
  lambdag2[[5]][r, ] <- fit3.greben$lambdag$groups[, fit3.greben$nouteriter + 1]
  lambdag2[[6]][r, ] <- fit4.greben$lambdag$groups[, fit4.greben$nouteriter + 1]

  varbeta2[r, ] <- sapply(1:G, function(g) {var(beta[(p*(g - 1)/G + 1):(p*g/G)])})

  results2 <- list(auc=auc2, briers=briers2, mse=mse2, lambdag=lambdag2,
                   varbeta=varbeta2)
  save(results2, file=paste(path.res, "grEBEN_test_res2.Rdata", sep=""))

}

### plots
load(paste(path.res, "grEBEN_test_res2.Rdata", sep=""))

# preparing data
auc1.enet <- do.call(rbind, results2$auc[[7]])[
  order(do.call(rbind, results2$auc[[7]])[, 1]), ]
auc2.enet <- do.call(rbind, results2$auc[[8]])[
  order(do.call(rbind, results2$auc[[8]])[, 1]), ]
auc3.enet <- do.call(rbind, results2$auc[[9]])[
  order(do.call(rbind, results2$auc[[9]])[, 1]), ]
auc4.enet <- do.call(rbind, results2$auc[[10]])[
  order(do.call(rbind, results2$auc[[10]])[, 1]), ]

auc1.greben <- do.call(rbind, results2$auc[[11]])[
  order(do.call(rbind, results2$auc[[11]])[, 1]), ]
auc2.greben <- do.call(rbind, results2$auc[[12]])[
  order(do.call(rbind, results2$auc[[12]])[, 1]), ]
auc3.greben <- do.call(rbind, results2$auc[[13]])[
  order(do.call(rbind, results2$auc[[13]])[, 1]), ]
auc4.greben <- do.call(rbind, results2$auc[[14]])[
  order(do.call(rbind, results2$auc[[14]])[, 1]), ]

auc1.grridge <- do.call(rbind, results2$auc[[3]])[
  order(do.call(rbind, results2$auc[[3]])[, 1]), ]
auc2.grridge <- do.call(rbind, results2$auc[[4]])[
  order(do.call(rbind, results2$auc[[4]])[, 1]), ]
auc3.grridge <- do.call(rbind, results2$auc[[5]])[
  order(do.call(rbind, results2$auc[[5]])[, 1]), ]
auc4.grridge <- do.call(rbind, results2$auc[[6]])[
  order(do.call(rbind, results2$auc[[6]])[, 1]), ]

auc1.ridge <- do.call(rbind, results2$auc[[1]])[
  order(do.call(rbind, results2$auc[[1]])[, 1]), ]
auc2.ridge <- do.call(rbind, results2$auc[[2]])[
  order(do.call(rbind, results2$auc[[2]])[, 1]), ]

lauc1.enet <- lowess(auc1.enet[, 1], auc1.enet[, 2])
lauc2.enet <- lowess(auc2.enet[, 1], auc2.enet[, 2])
lauc3.enet <- lowess(auc3.enet[, 1], auc3.enet[, 2])
lauc4.enet <- lowess(auc4.enet[, 1], auc4.enet[, 2])

lauc1.greben <- lowess(auc1.greben[, 1], auc1.greben[, 2])
lauc2.greben <- lowess(auc2.greben[, 1], auc2.greben[, 2])
lauc3.greben <- lowess(auc3.greben[, 1], auc3.greben[, 2])
lauc4.greben <- lowess(auc4.greben[, 1], auc4.greben[, 2])

lauc1.grridge <- lowess(auc1.grridge[, 1], auc1.grridge[, 2])
lauc2.grridge <- lowess(auc2.grridge[, 1], auc2.grridge[, 2])
lauc3.grridge <- lowess(auc3.grridge[, 1], auc3.grridge[, 2])
lauc4.grridge <- lowess(auc4.grridge[, 1], auc4.grridge[, 2])

lauc1.ridge <- lowess(auc1.ridge[, 1], auc1.ridge[, 2])
lauc2.ridge <- lowess(auc2.ridge[, 1], auc2.ridge[, 2])

briers1.enet <- do.call(rbind, results2$briers[[7]])[
  order(do.call(rbind, results2$briers[[7]])[, 1]), ]
briers2.enet <- do.call(rbind, results2$briers[[8]])[
  order(do.call(rbind, results2$briers[[8]])[, 1]), ]
briers3.enet <- do.call(rbind, results2$briers[[9]])[
  order(do.call(rbind, results2$briers[[9]])[, 1]), ]
briers4.enet <- do.call(rbind, results2$briers[[10]])[
  order(do.call(rbind, results2$briers[[10]])[, 1]), ]

briers1.greben <- do.call(rbind, results2$briers[[11]])[
  order(do.call(rbind, results2$briers[[11]])[, 1]), ]
briers2.greben <- do.call(rbind, results2$briers[[12]])[
  order(do.call(rbind, results2$briers[[12]])[, 1]), ]
briers3.greben <- do.call(rbind, results2$briers[[13]])[
  order(do.call(rbind, results2$briers[[13]])[, 1]), ]
briers4.greben <- do.call(rbind, results2$briers[[14]])[
  order(do.call(rbind, results2$briers[[14]])[, 1]), ]

briers1.grridge <- do.call(rbind, results2$briers[[3]])[
  order(do.call(rbind, results2$briers[[3]])[, 1]), ]
briers2.grridge <- do.call(rbind, results2$briers[[4]])[
  order(do.call(rbind, results2$briers[[4]])[, 1]), ]
briers3.grridge <- do.call(rbind, results2$briers[[5]])[
  order(do.call(rbind, results2$briers[[5]])[, 1]), ]
briers4.grridge <- do.call(rbind, results2$briers[[6]])[
  order(do.call(rbind, results2$briers[[6]])[, 1]), ]

briers1.ridge <- do.call(rbind, results2$briers[[1]])[
  order(do.call(rbind, results2$briers[[1]])[, 1]), ]
briers2.ridge <- do.call(rbind, results2$briers[[2]])[
  order(do.call(rbind, results2$briers[[2]])[, 1]), ]

lbriers1.enet <- lowess(briers1.enet[, 1], briers1.enet[, 2])
lbriers2.enet <- lowess(briers2.enet[, 1], briers2.enet[, 2])
lbriers3.enet <- lowess(briers3.enet[, 1], briers3.enet[, 2])
lbriers4.enet <- lowess(briers4.enet[, 1], briers4.enet[, 2])

lbriers1.greben <- lowess(briers1.greben[, 1], briers1.greben[, 2])
lbriers2.greben <- lowess(briers2.greben[, 1], briers2.greben[, 2])
lbriers3.greben <- lowess(briers3.greben[, 1], briers3.greben[, 2])
lbriers4.greben <- lowess(briers4.greben[, 1], briers4.greben[, 2])

lbriers1.grridge <- lowess(briers1.grridge[, 1], briers1.grridge[, 2])
lbriers2.grridge <- lowess(briers2.grridge[, 1], briers2.grridge[, 2])
lbriers3.grridge <- lowess(briers3.grridge[, 1], briers3.grridge[, 2])
lbriers4.grridge <- lowess(briers4.grridge[, 1], briers4.grridge[, 2])

lbriers1.ridge <- lowess(briers1.ridge[, 1], briers1.ridge[, 2])
lbriers2.ridge <- lowess(briers2.ridge[, 1], briers2.ridge[, 2])

mse1.enet <- do.call(rbind, results2$mse[[7]])[
  order(do.call(rbind, results2$mse[[7]])[, 1]), ]
mse2.enet <- do.call(rbind, results2$mse[[8]])[
  order(do.call(rbind, results2$mse[[8]])[, 1]), ]
mse3.enet <- do.call(rbind, results2$mse[[9]])[
  order(do.call(rbind, results2$mse[[9]])[, 1]), ]
mse4.enet <- do.call(rbind, results2$mse[[10]])[
  order(do.call(rbind, results2$mse[[10]])[, 1]), ]

mse1.greben <- do.call(rbind, results2$mse[[11]])[
  order(do.call(rbind, results2$mse[[11]])[, 1]), ]
mse2.greben <- do.call(rbind, results2$mse[[12]])[
  order(do.call(rbind, results2$mse[[12]])[, 1]), ]
mse3.greben <- do.call(rbind, results2$mse[[13]])[
  order(do.call(rbind, results2$mse[[13]])[, 1]), ]
mse4.greben <- do.call(rbind, results2$mse[[14]])[
  order(do.call(rbind, results2$mse[[14]])[, 1]), ]

mse1.grridge <- do.call(rbind, results2$mse[[3]])[
  order(do.call(rbind, results2$mse[[3]])[, 1]), ]
mse2.grridge <- do.call(rbind, results2$mse[[4]])[
  order(do.call(rbind, results2$mse[[4]])[, 1]), ]
mse3.grridge <- do.call(rbind, results2$mse[[5]])[
  order(do.call(rbind, results2$mse[[5]])[, 1]), ]
mse4.grridge <- do.call(rbind, results2$mse[[6]])[
  order(do.call(rbind, results2$mse[[6]])[, 1]), ]

mse1.ridge <- do.call(rbind, results2$mse[[1]])[
  order(do.call(rbind, results2$mse[[1]])[, 1]), ]
mse2.ridge <- do.call(rbind, results2$mse[[2]])[
  order(do.call(rbind, results2$mse[[2]])[, 1]), ]

lmse1.enet <- lowess(mse1.enet[, 1], mse1.enet[, 2])
lmse2.enet <- lowess(mse2.enet[, 1], mse2.enet[, 2])
lmse3.enet <- lowess(mse3.enet[, 1], mse3.enet[, 2])
lmse4.enet <- lowess(mse4.enet[, 1], mse4.enet[, 2])

lmse1.greben <- lowess(mse1.greben[, 1], mse1.greben[, 2])
lmse2.greben <- lowess(mse2.greben[, 1], mse2.greben[, 2])
lmse3.greben <- lowess(mse3.greben[, 1], mse3.greben[, 2])
lmse4.greben <- lowess(mse4.greben[, 1], mse4.greben[, 2])

lmse1.grridge <- lowess(mse1.grridge[, 1], mse1.grridge[, 2])
lmse2.grridge <- lowess(mse2.grridge[, 1], mse2.grridge[, 2])
lmse3.grridge <- lowess(mse3.grridge[, 1], mse3.grridge[, 2])
lmse4.grridge <- lowess(mse4.grridge[, 1], mse4.grridge[, 2])

lmse1.ridge <- lowess(mse1.ridge[, 1], mse1.ridge[, 2])
lmse2.ridge <- lowess(mse2.ridge[, 1], mse2.ridge[, 2])

### diagnostics: checking lowess fits
# AUC
xlim2.auc <- range(auc1.enet[, 1], auc2.enet[, 1], auc3.enet[, 1], auc4.enet[, 1],
                   auc1.greben[, 1], auc2.greben[, 1], auc3.greben[, 1],
                   auc4.greben[, 1], auc3.grridge[, 1], auc4.grridge[, 1])
ylim2.auc <- range(auc1.enet[, 2], auc2.enet[, 2], auc3.enet[, 2], auc4.enet[, 2],
                   auc1.greben[, 2], auc2.greben[, 2], auc3.greben[, 2],
                   auc4.greben[, 2], auc3.grridge[, 2], auc4.grridge[, 2])

png(paste(path.graph, "grEBEN_test_res2_lowess_auc.png", sep=""),
    units="in", width=12, height=6, res=120)
par(mfrow=c(2, 5))
plot(auc3.grridge[, 1], auc3.grridge[, 2], col=1, ylab="AUC", main="a)",
     xlab="Number of selected variables", xlim=xlim2.auc, ylim=ylim2.auc)
lines(lauc3.grridge, col=2)

plot(auc1.greben[, 1], auc1.greben[, 2], col=1, ylab="AUC", main="b)",
     xlab="Number of selected variables", xlim=xlim2.auc, ylim=ylim2.auc)
lines(lauc1.greben, col=4)

plot(auc2.greben[, 1], auc2.greben[, 2], col=1, ylab="AUC", main="c)",
     xlab="Number of selected variables", xlim=xlim2.auc, ylim=ylim2.auc)
lines(lauc2.greben, col=5)

plot(auc3.greben[, 1], auc4.greben[, 2], col=1, ylab="AUC", main="d)",
     xlab="Number of selected variables", xlim=xlim2.auc, ylim=ylim2.auc)
lines(lauc3.greben, col=6)

plot(auc4.greben[, 1], auc4.greben[, 2], col=1, ylab="AUC", main="e)",
     xlab="Number of selected variables", xlim=xlim2.auc, ylim=ylim2.auc)
lines(lauc4.greben, col=7)

plot(auc4.grridge[, 1], auc4.grridge[, 2], col=1, ylab="AUC", main="f)",
     xlab="Number of selected variables", xlim=xlim2.auc, ylim=ylim2.auc)
lines(lauc4.grridge, col=3)

plot(auc1.enet[, 1], auc1.enet[, 2], col=1, ylab="AUC", main="g)",
     xlab="Number of selected variables", xlim=xlim2.auc, ylim=ylim2.auc)
lines(lauc1.enet, col=4)

plot(auc2.enet[, 1], auc2.enet[, 2], col=1, ylab="AUC", main="h)",
     xlab="Number of selected variables", xlim=xlim2.auc, ylim=ylim2.auc)
lines(lauc2.enet, col=5)

plot(auc3.enet[, 1], auc4.enet[, 2], col=1, ylab="AUC", main="i)",
     xlab="Number of selected variables", xlim=xlim2.auc, ylim=ylim2.auc)
lines(lauc3.enet, col=6)

plot(auc4.enet[, 1], auc4.enet[, 2], col=1, ylab="AUC", main="j)",
     xlab="Number of selected variables", xlim=xlim2.auc, ylim=ylim2.auc)
lines(lauc4.enet, col=7)
dev.off()

# briers
xlim2.briers <- range(briers1.enet[, 1], briers2.enet[, 1], briers3.enet[, 1], briers4.enet[, 1],
                   briers1.greben[, 1], briers2.greben[, 1], briers3.greben[, 1],
                   briers4.greben[, 1], briers3.grridge[, 1], briers4.grridge[, 1])
ylim2.briers <- range(briers1.enet[, 2], briers2.enet[, 2], briers3.enet[, 2], briers4.enet[, 2],
                   briers1.greben[, 2], briers2.greben[, 2], briers3.greben[, 2],
                   briers4.greben[, 2], briers3.grridge[, 2], briers4.grridge[, 2])

png(paste(path.graph, "grEBEN_test_res2_lowess_briers.png", sep=""),
    units="in", width=12, height=6, res=120)
par(mfrow=c(2, 5))
plot(briers3.grridge[, 1], briers3.grridge[, 2], col=1, ylab="Brier skill score", main="a)",
     xlab="Number of selected variables", xlim=xlim2.briers, ylim=ylim2.briers)
lines(lbriers3.grridge, col=2)

plot(briers1.greben[, 1], briers1.greben[, 2], col=1, ylab="Brier skill score", main="b)",
     xlab="Number of selected variables", xlim=xlim2.briers, ylim=ylim2.briers)
lines(lbriers1.greben, col=4)

plot(briers2.greben[, 1], briers2.greben[, 2], col=1, ylab="Brier skill score", main="c)",
     xlab="Number of selected variables", xlim=xlim2.briers, ylim=ylim2.briers)
lines(lbriers2.greben, col=5)

plot(briers3.greben[, 1], briers4.greben[, 2], col=1, ylab="Brier skill score", main="d)",
     xlab="Number of selected variables", xlim=xlim2.briers, ylim=ylim2.briers)
lines(lbriers3.greben, col=6)

plot(briers4.greben[, 1], briers4.greben[, 2], col=1, ylab="Brier skill score", main="e)",
     xlab="Number of selected variables", xlim=xlim2.briers, ylim=ylim2.briers)
lines(lbriers4.greben, col=7)

plot(briers4.grridge[, 1], briers4.grridge[, 2], col=1, ylab="Brier skill score", main="f)",
     xlab="Number of selected variables", xlim=xlim2.briers, ylim=ylim2.briers)
lines(lbriers4.grridge, col=3)

plot(briers1.enet[, 1], briers1.enet[, 2], col=1, ylab="Brier skill score", main="g)",
     xlab="Number of selected variables", xlim=xlim2.briers, ylim=ylim2.briers)
lines(lbriers1.enet, col=4)

plot(briers2.enet[, 1], briers2.enet[, 2], col=1, ylab="Brier skill score", main="h)",
     xlab="Number of selected variables", xlim=xlim2.briers, ylim=ylim2.briers)
lines(lbriers2.enet, col=5)

plot(briers3.enet[, 1], briers4.enet[, 2], col=1, ylab="Brier skill score", main="i)",
     xlab="Number of selected variables", xlim=xlim2.briers, ylim=ylim2.briers)
lines(lbriers3.enet, col=6)

plot(briers4.enet[, 1], briers4.enet[, 2], col=1, ylab="Brier skill score", main="j)",
     xlab="Number of selected variables", xlim=xlim2.briers, ylim=ylim2.briers)
lines(lbriers4.enet, col=7)
dev.off()

# mse
xlim2.mse <- range(mse1.enet[, 1], mse2.enet[, 1], mse3.enet[, 1], mse4.enet[, 1],
                   mse1.greben[, 1], mse2.greben[, 1], mse3.greben[, 1],
                   mse4.greben[, 1], mse3.grridge[, 1], mse4.grridge[, 1])
ylim2.mse <- range(mse1.enet[, 2], mse2.enet[, 2], mse3.enet[, 2], mse4.enet[, 2],
                   mse1.greben[, 2], mse2.greben[, 2], mse3.greben[, 2],
                   mse4.greben[, 2], mse3.grridge[, 2], mse4.grridge[, 2])

png(paste(path.graph, "grEBEN_test_res2_lowess_mse.png", sep=""),
    units="in", width=12, height=6, res=120)
par(mfrow=c(2, 5))
plot(mse3.grridge[, 1], mse3.grridge[, 2], col=1, ylab="MSE", main="a)",
     xlab="Number of selected variables", xlim=xlim2.mse, ylim=ylim2.mse)
lines(lmse3.grridge, col=2)

plot(mse1.greben[, 1], mse1.greben[, 2], col=1, ylab="MSE", main="b)",
     xlab="Number of selected variables", xlim=xlim2.mse, ylim=ylim2.mse)
lines(lmse1.greben, col=4)

plot(mse2.greben[, 1], mse2.greben[, 2], col=1, ylab="MSE", main="c)",
     xlab="Number of selected variables", xlim=xlim2.mse, ylim=ylim2.mse)
lines(lmse2.greben, col=5)

plot(mse3.greben[, 1], mse4.greben[, 2], col=1, ylab="MSE", main="d)",
     xlab="Number of selected variables", xlim=xlim2.mse, ylim=ylim2.mse)
lines(lmse3.greben, col=6)

plot(mse4.greben[, 1], mse4.greben[, 2], col=1, ylab="MSE", main="e)",
     xlab="Number of selected variables", xlim=xlim2.mse, ylim=ylim2.mse)
lines(lmse4.greben, col=7)

plot(mse4.grridge[, 1], mse4.grridge[, 2], col=1, ylab="MSE", main="f)",
     xlab="Number of selected variables", xlim=xlim2.mse, ylim=ylim2.mse)
lines(lmse4.grridge, col=3)

plot(mse1.enet[, 1], mse1.enet[, 2], col=1, ylab="MSE", main="g)",
     xlab="Number of selected variables", xlim=xlim2.mse, ylim=ylim2.mse)
lines(lmse1.enet, col=4)

plot(mse2.enet[, 1], mse2.enet[, 2], col=1, ylab="MSE", main="h)",
     xlab="Number of selected variables", xlim=xlim2.mse, ylim=ylim2.mse)
lines(lmse2.enet, col=5)

plot(mse3.enet[, 1], mse4.enet[, 2], col=1, ylab="MSE", main="i)",
     xlab="Number of selected variables", xlim=xlim2.mse, ylim=ylim2.mse)
lines(lmse3.enet, col=6)

plot(mse4.enet[, 1], mse4.enet[, 2], col=1, ylab="MSE", main="j)",
     xlab="Number of selected variables", xlim=xlim2.mse, ylim=ylim2.mse)
lines(lmse4.enet, col=7)
dev.off()



### Performance measures
# AUC
png(paste(path.graph, "grEBEN_test_res2_performance.png", sep=""),
    units="in", width=12, height=4, res=120)
par(mfrow=c(1, 3))
xlim1.auc <- range(lauc1.enet$x, lauc2.enet$x, lauc3.enet$x, lauc4.enet$x,
                   lauc1.greben$x, lauc2.greben$x, lauc3.greben$x,
                   lauc4.greben$x, lauc3.grridge$x, lauc4.grridge$x)
ylim1.auc <- range(lauc1.enet$y, lauc2.enet$y, lauc3.enet$y, lauc4.enet$y,
                   lauc1.greben$y, lauc2.greben$y, lauc3.greben$y,
                   lauc4.greben$y, lauc1.grridge$y, lauc2.grridge$y,
                   lauc3.grridge$y, lauc4.grridge$y, lauc1.ridge$y,
                   lauc2.ridge$y)

plot(0, 0, col=2, type="n", ylab="AUC", main="a)",
     xlab="Number of selected variables", xlim=xlim1.auc, ylim=ylim1.auc)

lines(lauc3.grridge, col=2, lty=1)
abline(h=lauc1.grridge$y, col=2, lty=2)
abline(h=lauc1.ridge$y, col=2, lty=3)

lines(lauc4.grridge, col=3, lty=1)
abline(h=lauc2.grridge$y, col=3, lty=2)
abline(h=lauc2.ridge$y, col=3, lty=3)

lines(lauc1.greben, col=4, lty=1)
lines(lauc1.enet, col=4, lty=3)

lines(lauc2.greben, col=5, lty=1)
lines(lauc2.enet, col=5, lty=3)

lines(lauc3.greben, col=6, lty=1)
lines(lauc3.enet, col=6, lty=3)

lines(lauc4.greben, col=7, lty=1)
lines(lauc4.enet, col=7, lty=3)

# Brier skilll
xlim1.briers <- range(lbriers1.enet$x, lbriers2.enet$x, lbriers3.enet$x,
                      lbriers4.enet$x, lbriers1.greben$x, lbriers2.greben$x,
                      lbriers3.greben$x, lbriers4.greben$x, lbriers3.grridge$x,
                      lbriers4.grridge$x)
ylim1.briers <- range(lbriers1.enet$y, lbriers2.enet$y, lbriers3.enet$y,
                      lbriers4.enet$y, lbriers1.greben$y, lbriers2.greben$y,
                      lbriers3.greben$y, lbriers4.greben$y, lbriers1.grridge$y,
                      lbriers2.grridge$y, lbriers3.grridge$y, lbriers4.grridge$y,
                      lbriers1.ridge$y, lbriers2.ridge$y)

plot(0, 0, col=2, type="n", ylab="Brier skill score", main="b)",
     xlab="Number of selected variables", xlim=xlim1.briers, ylim=ylim1.briers)

lines(lbriers3.grridge, col=2, lty=1)
abline(h=lbriers1.grridge$y, col=2, lty=2)
abline(h=lbriers1.ridge$y, col=2, lty=3)

lines(lbriers4.grridge, col=3, lty=1)
abline(h=lbriers2.grridge$y, col=3, lty=2)
abline(h=lbriers2.ridge$y, col=3, lty=3)

lines(lbriers1.greben, col=4, lty=1)
lines(lbriers1.enet, col=4, lty=3)

lines(lbriers2.greben, col=5, lty=1)
lines(lbriers2.enet, col=5, lty=3)

lines(lbriers3.greben, col=6, lty=1)
lines(lbriers3.enet, col=6, lty=3)

lines(lbriers4.greben, col=7, lty=1)
lines(lbriers4.enet, col=7, lty=3)

leglabels <- c(expression(paste("ridge, 'true' ", lambda)),
               "ridge", expression(paste("enet, true ", lambda)),
               expression(paste("enet, ", alpha==0.05)),
               expression(paste("enet, ", alpha==0.5)),
               expression(paste("enet, ", alpha==0.95)),
               "group-regularized + selection", "group-regularized",
               "not group-regularized")
legend("bottomright", legend=leglabels, fill=c(2:7, 0, 0, 0),
       lty=c(rep(NA, 6), 1, 2, 3), border=c(rep(1, 6), 0 ,0, 0), merge=TRUE,
       seg.len=1)

# MSE
xlim1.mse <- range(lmse1.enet$x, lmse2.enet$x, lmse3.enet$x,
                      lmse4.enet$x, lmse1.greben$x, lmse2.greben$x,
                      lmse3.greben$x, lmse4.greben$x, lmse3.grridge$x,
                      lmse4.grridge$x)
ylim1.mse <- range(lmse1.enet$y, lmse2.enet$y, lmse3.enet$y,
                      lmse4.enet$y, lmse1.greben$y, lmse2.greben$y,
                      lmse3.greben$y, lmse4.greben$y, lmse1.grridge$y,
                      lmse2.grridge$y, lmse3.grridge$y, lmse4.grridge$y,
                      lmse1.ridge$y, lmse2.ridge$y)

plot(0, 0, col=2, type="n", ylab="MSE", main="c)",
     xlab="Number of selected variables", xlim=xlim1.mse, ylim=ylim1.mse)

lines(lmse3.grridge, col=2, lty=1)
abline(h=lmse1.grridge$y, col=2, lty=2)
abline(h=lmse1.ridge$y, col=2, lty=3)

lines(lmse4.grridge, col=3, lty=1)
abline(h=lmse2.grridge$y, col=3, lty=2)
abline(h=lmse2.ridge$y, col=3, lty=3)

lines(lmse1.greben, col=4, lty=1)
lines(lmse1.enet, col=4, lty=3)

lines(lmse2.greben, col=5, lty=1)
lines(lmse2.enet, col=5, lty=3)

lines(lmse3.greben, col=6, lty=1)
lines(lmse3.enet, col=6, lty=3)

lines(lmse4.greben, col=7, lty=1)
lines(lmse4.enet, col=7, lty=3)

dev.off()

## penalty parameters
lambdag <- exp(seq(-1, 1, length.out=ncol(results2$lambdag[[1]])))
png(paste(path.graph, "grEBEN_test_res2_penalties.png", sep=""),
    units="in", width=8, height=6, res=120)
par(mfrow=c(2, 3))
boxplot(results2$lambdag[[1]], ylim=range(results2$lambdag[[1]], lambdag), 
        xlab="", ylab="", main="a)")
title(ylab=expression(paste(lambda[g], "'")), 
      xlab="Groups in decreasing effect size", line=2.5)
points(1:length(lambdag), lambdag, col=2, pch="x")
boxplot(results2$lambdag[[2]], ylim=range(results2$lambdag[[2]], lambdag), 
        xlab="", ylab="", main="b)")
title(ylab=expression(paste(lambda[g], "'")), 
      xlab="Groups in decreasing effect size", line=2.5)
points(1:length(lambdag), lambdag, col=2, pch="x")
boxplot(results2$lambdag[[3]], ylim=range(results2$lambdag[[3]], lambdag), 
        xlab="", ylab="", main="c)")
title(ylab=expression(paste(lambda[g], "'")), 
      xlab="Groups in decreasing effect size", line=2.5)
points(1:length(lambdag), lambdag, col=2, pch="x")
boxplot(results2$lambdag[[4]], ylim=range(results2$lambdag[[4]], lambdag), 
        xlab="", ylab="", main="d)")
title(ylab=expression(paste(lambda[g], "'")), 
      xlab="Groups in decreasing effect size", line=2.5)
points(1:length(lambdag), lambdag, col=2, pch="x")
boxplot(results2$lambdag[[5]], ylim=range(results2$lambdag[[5]], lambdag), 
        xlab="", ylab="", main="e)")
title(ylab=expression(paste(lambda[g], "'")), 
      xlab="Groups in decreasing effect size", line=2.5)
points(1:length(lambdag), lambdag, col=2, pch="x")
boxplot(results2$lambdag[[6]], ylim=range(results2$lambdag[[6]], lambdag), 
        xlab="", ylab="", main="f)")
title(ylab=expression(paste(lambda[g], "'")), 
      xlab="Groups in decreasing effect size", line=2.5)
points(1:length(lambdag), lambdag, col=2, pch="x")
dev.off()






# ## simulations
# ## simulation 3
# # create data
# n <- 200
# p <- 1000
# G <- 5
# pblock <- 20
# rho <- 0.7
# sigma <- matrix(rho, ncol=pblock, nrow=pblock); diag(sigma) <- 1
# lambda <- 0.02
# alpha <- 0
# lambdag <- exp(seq(-1, 1, length.out=G))
# m <- rep(1, n)
# part.greben <- list(groups=rep(1:G, each=p/G))
# part.grridge <- list(groups=CreatePartition(as.factor(part.greben$groups)))
# ntest <- 1000
# 
# methods <- c("ridge+truelambda", "ridge", "GRridge+truelambda", "GRridge",
#              "GRridge+truelambda+sel", "GRridge+sel",
#              "enet+truelambda", "enet+a=0.05", "enet+a=0.5", "enet+a=0.95",
#              "grEBEN+truelambda", "grEBEN+a=0.05", "grEBEN+a=0.5","grEBEN+a=0.95")
# nreps <- 50
# auc3 <- briers3 <- mse3 <- rep(list(vector(mode="list", length=nreps)),
#                                length(methods))
# names(auc3) <- names(briers3) <- names(mse3) <- methods
# lambdag3 <- rep(list(matrix(NA, ncol=G, nrow=nreps)), 6)
# names(lambdag3) <- methods[c(3, 4, 11, 12, 13, 14)]
# varbeta3 <- matrix(nrow=nreps, ncol=G)
# # the simulations
# for(r in 1:1) {
# 
#   print(paste("Simulation 3, repetition", r))
#   set.seed(300 + r)
#   x <- do.call(cbind, replicate(p/pblock, rmvnorm(n, mean=rep(0, pblock), sigma=sigma),
#                                 simplify=FALSE))
#   beta <- as.numeric(sapply(1:G, function(g) {
#     renbeta(p/G, 2*n*lambda*alpha*sqrt(lambdag[g]), n*lambda*(1 - alpha)*lambdag[g])}))
#   prob <- as.numeric(exp(x %*% beta)/(1 + exp(x %*% beta)))
#   y <- rbinom(n, 1, prob)
#   optl <- n*lambda*alpha*sum(abs(beta))/sum(beta^2) + 0.5*n*lambda*(1 - alpha)
# 
#   xtest <- do.call(cbind, replicate(p/pblock, rmvnorm(ntest, mean=rep(0, pblock),
#                                                       sigma=sigma), simplify=FALSE))
#   probtest <- as.numeric(exp(xtest %*% beta)/(1 + exp(xtest %*% beta)))
#   ytest <- rbinom(ntest, 1, probtest)
# 
#   fit1.greben <- grEBEN3(x, y, m, partitions=part.greben, alpha=alpha,
#                          lambda=lambda, psel=TRUE)
#   fit2.greben <- grEBEN3(x, y, m, partitions=part.greben, alpha=0.05, psel=TRUE)
#   fit3.greben <- grEBEN3(x, y, m, partitions=part.greben, alpha=0.5, psel=TRUE)
#   fit4.greben <- grEBEN3(x, y, m, partitions=part.greben, alpha=0.95, psel=TRUE)
# 
#   # number of selected variables
#   psel1.enet <- apply(fit1.greben$beta.nogroups, 2, function(b) {sum(b!=0) - 1})
#   psel2.enet <- apply(fit2.greben$beta.nogroups, 2, function(b) {sum(b!=0) - 1})
#   psel3.enet <- apply(fit3.greben$beta.nogroups, 2, function(b) {sum(b!=0) - 1})
#   psel4.enet <- apply(fit4.greben$beta.nogroups, 2, function(b) {sum(b!=0) - 1})
# 
#   psel1.greben <- apply(fit1.greben$beta, 2, function(b) {sum(b!=0) - 1})
#   psel2.greben <- apply(fit2.greben$beta, 2, function(b) {sum(b!=0) - 1})
#   psel3.greben <- apply(fit3.greben$beta, 2, function(b) {sum(b!=0) - 1})
#   psel4.greben <- apply(fit4.greben$beta, 2, function(b) {sum(b!=0) - 1})
# 
#   psel.all <- unique(c(psel1.enet, psel2.enet, psel3.enet, psel4.enet,
#     psel1.greben, psel2.greben, psel3.greben, psel4.greben))
#   psel.in <- floor(quantile(psel.all[psel.all!=0],
#     prob=c(seq(0.01, 0.05, 0.01), seq(0.07, 0.25, 0.03), seq(0.3, 1, 0.1))))
# 
#   fit1.grridge <- vector(mode="list", length=length(psel.in))
#   for(s in 1:length(psel.in)) {
#     fit1.grridge[[s]] <- grridge(t(x), y, part.grridge, selectionEN=TRUE,
#                                  optl=optl, maxsel=psel.in[s])
#   }
# 
#   fit2.grridge <- vector(mode="list", length=length(psel.in))
#   fit2.grridge[[1]] <- grridge(t(x), y, part.grridge, selectionEN=TRUE,
#                                maxsel=psel.in[1])
#   for(s in 2:length(psel.in)) {
#     fit2.grridge[[s]] <- grridge(t(x), y, part.grridge, selectionEN=TRUE,
#                                  optl=fit2.grridge[[1]]$optl,
#                                  maxsel=psel.in[s])
#   }
# 
#   psel1.grridge <- sapply(1:length(fit1.grridge), function(s) {
#     return(length(fit1.grridge[[s]]$resEN$whichEN))})
#   psel2.grridge <- sapply(1:length(fit2.grridge), function(s) {
#     return(length(fit2.grridge[[s]]$resEN$whichEN))})
# 
#   # estimates
#   est1.ridge <- coef(fit1.grridge[[1]]$predobj$NoGroups, "all")
#   est2.ridge <- coef(fit2.grridge[[1]]$predobj$NoGroups, "all")
# 
#   est1.grridge <- coef(fit1.grridge[[1]]$predobj$GroupRegul, "all")
#   est2.grridge <- coef(fit2.grridge[[1]]$predobj$GroupRegul, "all")
# 
#   est3.grridge <- sapply(fit1.grridge, function(s) {
#     replace(rep(0, p + 1), c(1, s$resEN$whichEN + 1), coef(s$predobj$EN, "all"))})
#   est4.grridge <- sapply(fit2.grridge, function(s) {
#     replace(rep(0, p + 1), c(1, s$resEN$whichEN + 1), coef(s$predobj$EN, "all"))})
# 
#   est1.enet <- fit1.greben$beta.nogroups
#   est2.enet <- fit2.greben$beta.nogroups
#   est3.enet <- fit3.greben$beta.nogroups
#   est4.enet <- fit4.greben$beta.nogroups
# 
#   est1.greben <- fit1.greben$beta
#   est2.greben <- fit2.greben$beta
#   est3.greben <- fit3.greben$beta
#   est4.greben <- fit4.greben$beta
# 
#   # predictions on fit data
#   pred1.ridge <- predict.grridge(fit1.grridge[[1]], t(xtest))[, 1]
#   pred2.ridge <- predict.grridge(fit2.grridge[[1]], t(xtest))[, 1]
# 
#   pred1.grridge <- predict.grridge(fit1.grridge[[1]], t(xtest))[, 2]
#   pred2.grridge <- predict.grridge(fit2.grridge[[1]], t(xtest))[, 2]
#   pred3.grridge <- sapply(fit1.grridge, function(s) {
#     predict.grridge(s, t(xtest))[, 3]})
#   pred4.grridge <- sapply(fit2.grridge, function(s) {
#     predict.grridge(s, t(xtest))[, 3]})
# 
#   pred1.enet <- 1/(1 + exp(-xtest %*% est1.enet[-1, ]))
#   pred2.enet <- 1/(1 + exp(-xtest %*% est2.enet[-1, ]))
#   pred3.enet <- 1/(1 + exp(-xtest %*% est3.enet[-1, ]))
#   pred4.enet <- 1/(1 + exp(-xtest %*% est4.enet[-1, ]))
# 
#   pred1.greben <- 1/(1 + exp(-xtest %*% est1.greben[-1, ]))
#   pred2.greben <- 1/(1 + exp(-xtest %*% est2.greben[-1, ]))
#   pred3.greben <- 1/(1 + exp(-xtest %*% est3.greben[-1, ]))
#   pred4.greben <- 1/(1 + exp(-xtest %*% est4.greben[-1, ]))
# 
#    # AUCs
#   auc.true <- pROC::roc(ytest, probtest)$auc
# 
#   auc1.ridge <- pROC::roc(ytest, pred1.ridge)$auc
#   auc2.ridge <- pROC::roc(ytest, pred2.ridge)$auc
# 
#   auc1.grridge <- pROC::roc(ytest, pred1.grridge)$auc
#   auc2.grridge <- pROC::roc(ytest, pred2.grridge)$auc
#   auc3.grridge <- apply(pred3.grridge, 2, function(r) {pROC::roc(ytest, r)$auc})
#   auc4.grridge <- apply(pred4.grridge, 2, function(r) {pROC::roc(ytest, r)$auc})
# 
#   auc1.enet <- apply(pred1.enet, 2, function(r) {pROC::roc(ytest, r)$auc})
#   auc2.enet <- apply(pred2.enet, 2, function(r) {pROC::roc(ytest, r)$auc})
#   auc3.enet <- apply(pred3.enet, 2, function(r) {pROC::roc(ytest, r)$auc})
#   auc4.enet <- apply(pred4.enet, 2, function(r) {pROC::roc(ytest, r)$auc})
# 
#   auc1.greben <- apply(pred1.greben, 2, function(r) {pROC::roc(ytest, r)$auc})
#   auc2.greben <- apply(pred2.greben, 2, function(r) {pROC::roc(ytest, r)$auc})
#   auc3.greben <- apply(pred3.greben, 2, function(r) {pROC::roc(ytest, r)$auc})
#   auc4.greben <- apply(pred4.greben, 2, function(r) {pROC::roc(ytest, r)$auc})
# 
#   # Brier scores
#   brier.null <- sum((ytest - mean(ytest))^2)
#   briers.true <- 1 - sum((ytest - probtest)^2)/brier.null
# 
#   briers1.ridge <- 1 - sum((ytest - pred1.ridge)^2)/brier.null
#   briers2.ridge <- 1 - sum((ytest - pred2.ridge)^2)/brier.null
# 
#   briers1.grridge <- 1 - sum((ytest - pred1.grridge)^2)/brier.null
#   briers2.grridge <- 1 - sum((ytest - pred2.grridge)^2)/brier.null
#   briers3.grridge <- apply(pred3.grridge, 2, function(pred) {
#     1 - sum((ytest - pred)^2)/brier.null})
#   briers4.grridge <- apply(pred4.grridge, 2, function(pred) {
#     1 - sum((ytest - pred)^2)/brier.null})
# 
#   briers1.enet <- apply(pred1.enet, 2, function(pred) {
#     1 - sum((ytest - pred)^2)/brier.null})
#   briers2.enet <- apply(pred2.enet, 2, function(pred) {
#     1 - sum((ytest - pred)^2)/brier.null})
#   briers3.enet <- apply(pred3.enet, 2, function(pred) {
#     1 - sum((ytest - pred)^2)/brier.null})
#   briers4.enet <- apply(pred4.enet, 2, function(pred) {
#     1 - sum((ytest - pred)^2)/brier.null})
# 
#   briers1.greben <- apply(pred1.greben, 2, function(pred) {
#     1 - sum((ytest - pred)^2)/brier.null})
#   briers2.greben <- apply(pred2.greben, 2, function(pred) {
#     1 - sum((ytest - pred)^2)/brier.null})
#   briers3.greben <- apply(pred3.greben, 2, function(pred) {
#     1 - sum((ytest - pred)^2)/brier.null})
#   briers4.greben <- apply(pred4.greben, 2, function(pred) {
#     1 - sum((ytest - pred)^2)/brier.null})
# 
#   # MSE
#   mse.true <- 0
# 
#   mse1.ridge <- mean((c(0, beta) - est1.ridge)^2)
#   mse2.ridge <- mean((c(0, beta) - est2.ridge)^2)
# 
#   mse1.grridge <- mean((c(0, beta) - est1.grridge)^2)
#   mse2.grridge <- mean((c(0, beta) - est2.grridge)^2)
#   mse3.grridge <- apply(est3.grridge, 2, function(b) {mean((c(0, beta) - b)^2)})
#   mse4.grridge <- apply(est4.grridge, 2, function(b) {mean((c(0, beta) - b)^2)})
# 
#   mse1.enet <- apply(est1.enet, 2, function(b) {mean((c(0, beta) - b)^2)})
#   mse2.enet <- apply(est2.enet, 2, function(b) {mean((c(0, beta) - b)^2)})
#   mse3.enet <- apply(est3.enet, 2, function(b) {mean((c(0, beta) - b)^2)})
#   mse4.enet <- apply(est4.enet, 2, function(b) {mean((c(0, beta) - b)^2)})
# 
#   mse1.greben <- apply(est1.greben, 2, function(b) {mean((c(0, beta) - b)^2)})
#   mse2.greben <- apply(est2.greben, 2, function(b) {mean((c(0, beta) - b)^2)})
#   mse3.greben <- apply(est3.greben, 2, function(b) {mean((c(0, beta) - b)^2)})
#   mse4.greben <- apply(est4.greben, 2, function(b) {mean((c(0, beta) - b)^2)})
# 
#   auc3[[1]][[r]] <- cbind(psel=p, auc=auc1.ridge)
#   auc3[[2]][[r]] <- cbind(psel=p, auc=auc2.ridge)
#   auc3[[3]][[r]] <- cbind(psel=p, auc=auc1.grridge)
#   auc3[[4]][[r]] <- cbind(psel=p, auc=auc2.grridge)
#   auc3[[5]][[r]] <- cbind(psel=psel1.grridge, auc=auc3.grridge)
#   auc3[[6]][[r]] <- cbind(psel=psel2.grridge, auc=auc4.grridge)
#   auc3[[7]][[r]] <- cbind(psel=psel1.enet, auc=auc1.enet)
#   auc3[[8]][[r]] <- cbind(psel=psel2.enet, auc=auc2.enet)
#   auc3[[9]][[r]] <- cbind(psel=psel3.enet, auc=auc3.enet)
#   auc3[[10]][[r]] <- cbind(psel=psel4.enet, auc=auc4.enet)
#   auc3[[11]][[r]] <- cbind(psel=psel1.greben, auc=auc1.greben)
#   auc3[[12]][[r]] <- cbind(psel=psel2.greben, auc=auc2.greben)
#   auc3[[13]][[r]] <- cbind(psel=psel3.greben, auc=auc3.greben)
#   auc3[[14]][[r]] <- cbind(psel=psel4.greben, auc=auc4.greben)
# 
#   briers3[[1]][[r]] <- cbind(psel=p, briers=briers1.ridge)
#   briers3[[2]][[r]] <- cbind(psel=p, briers=briers2.ridge)
#   briers3[[3]][[r]] <- cbind(psel=p, briers=briers1.grridge)
#   briers3[[4]][[r]] <- cbind(psel=p, briers=briers2.grridge)
#   briers3[[5]][[r]] <- cbind(psel=psel1.grridge, briers=briers3.grridge)
#   briers3[[6]][[r]] <- cbind(psel=psel2.grridge, briers=briers4.grridge)
#   briers3[[7]][[r]] <- cbind(psel=psel1.enet, briers=briers1.enet)
#   briers3[[8]][[r]] <- cbind(psel=psel2.enet, briers=briers2.enet)
#   briers3[[9]][[r]] <- cbind(psel=psel3.enet, briers=briers3.enet)
#   briers3[[10]][[r]] <- cbind(psel=psel4.enet, briers=briers4.enet)
#   briers3[[11]][[r]] <- cbind(psel=psel1.greben, briers=briers1.greben)
#   briers3[[12]][[r]] <- cbind(psel=psel2.greben, briers=briers2.greben)
#   briers3[[13]][[r]] <- cbind(psel=psel3.greben, briers=briers3.greben)
#   briers3[[14]][[r]] <- cbind(psel=psel4.greben, briers=briers4.greben)
# 
#   mse3[[1]][[r]] <- cbind(psel=p, mse=mse1.ridge)
#   mse3[[2]][[r]] <- cbind(psel=p, mse=mse2.ridge)
#   mse3[[3]][[r]] <- cbind(psel=p, mse=mse1.grridge)
#   mse3[[4]][[r]] <- cbind(psel=p, mse=mse2.grridge)
#   mse3[[5]][[r]] <- cbind(psel=psel1.grridge, mse=mse3.grridge)
#   mse3[[6]][[r]] <- cbind(psel=psel2.grridge, mse=mse4.grridge)
#   mse3[[7]][[r]] <- cbind(psel=psel1.enet, mse=mse1.enet)
#   mse3[[8]][[r]] <- cbind(psel=psel2.enet, mse=mse2.enet)
#   mse3[[9]][[r]] <- cbind(psel=psel3.enet, mse=mse3.enet)
#   mse3[[10]][[r]] <- cbind(psel=psel4.enet, mse=mse4.enet)
#   mse3[[11]][[r]] <- cbind(psel=psel1.greben, mse=mse1.greben)
#   mse3[[12]][[r]] <- cbind(psel=psel2.greben, mse=mse2.greben)
#   mse3[[13]][[r]] <- cbind(psel=psel3.greben, mse=mse3.greben)
#   mse3[[14]][[r]] <- cbind(psel=psel4.greben, mse=mse4.greben)
# 
# 
#   lambdag3[[1]][r, ] <- fit1.grridge[[1]]$lambdamults$groups
#   lambdag3[[2]][r, ] <- fit2.grridge[[1]]$lambdamults$groups
#   lambdag3[[3]][r, ] <- fit1.greben$lambdag$groups[, fit1.greben$nouteriter + 1]
#   lambdag3[[4]][r, ] <- fit2.greben$lambdag$groups[, fit2.greben$nouteriter + 1]
#   lambdag3[[5]][r, ] <- fit3.greben$lambdag$groups[, fit3.greben$nouteriter + 1]
#   lambdag3[[6]][r, ] <- fit4.greben$lambdag$groups[, fit4.greben$nouteriter + 1]
# 
#   varbeta3[r, ] <- sapply(1:G, function(g) {var(beta[(p*(g - 1)/G + 1):(p*g/G)])})
# 
#   results3 <- list(auc=auc3, briers=briers3, mse=mse3, lambdag=lambdag3,
#                    varbeta=varbeta3)
#   save(results3, file=paste(path.res, "grEBEN_test_res3.Rdata", sep=""))
# 
# }

# ### plots
# load(paste(path.res, "grEBEN_test_res3.Rdata", sep=""))
# 
# # preparing data
# auc1.enet <- do.call(rbind, results3$auc[[7]])[
#   order(do.call(rbind, results3$auc[[7]])[, 1]), ]
# auc2.enet <- do.call(rbind, results3$auc[[8]])[
#   order(do.call(rbind, results3$auc[[8]])[, 1]), ]
# auc3.enet <- do.call(rbind, results3$auc[[9]])[
#   order(do.call(rbind, results3$auc[[9]])[, 1]), ]
# auc4.enet <- do.call(rbind, results3$auc[[10]])[
#   order(do.call(rbind, results3$auc[[10]])[, 1]), ]
# 
# auc1.greben <- do.call(rbind, results3$auc[[11]])[
#   order(do.call(rbind, results3$auc[[11]])[, 1]), ]
# auc2.greben <- do.call(rbind, results3$auc[[12]])[
#   order(do.call(rbind, results3$auc[[12]])[, 1]), ]
# auc3.greben <- do.call(rbind, results3$auc[[13]])[
#   order(do.call(rbind, results3$auc[[13]])[, 1]), ]
# auc4.greben <- do.call(rbind, results3$auc[[14]])[
#   order(do.call(rbind, results3$auc[[14]])[, 1]), ]
# 
# auc1.grridge <- do.call(rbind, results3$auc[[3]])[
#   order(do.call(rbind, results3$auc[[3]])[, 1]), ]
# auc2.grridge <- do.call(rbind, results3$auc[[4]])[
#   order(do.call(rbind, results3$auc[[4]])[, 1]), ]
# auc3.grridge <- do.call(rbind, results3$auc[[5]])[
#   order(do.call(rbind, results3$auc[[5]])[, 1]), ]
# auc4.grridge <- do.call(rbind, results3$auc[[6]])[
#   order(do.call(rbind, results3$auc[[6]])[, 1]), ]
# 
# auc1.ridge <- do.call(rbind, results3$auc[[1]])[
#   order(do.call(rbind, results3$auc[[1]])[, 1]), ]
# auc2.ridge <- do.call(rbind, results3$auc[[2]])[
#   order(do.call(rbind, results3$auc[[2]])[, 1]), ]
# 
# lauc1.enet <- lowess(auc1.enet[, 1], auc1.enet[, 2])
# lauc2.enet <- lowess(auc2.enet[, 1], auc2.enet[, 2])
# lauc3.enet <- lowess(auc3.enet[, 1], auc3.enet[, 2])
# lauc4.enet <- lowess(auc4.enet[, 1], auc4.enet[, 2])
# 
# lauc1.greben <- lowess(auc1.greben[, 1], auc1.greben[, 2])
# lauc2.greben <- lowess(auc2.greben[, 1], auc2.greben[, 2])
# lauc3.greben <- lowess(auc3.greben[, 1], auc3.greben[, 2])
# lauc4.greben <- lowess(auc4.greben[, 1], auc4.greben[, 2])
# 
# lauc1.grridge <- lowess(auc1.grridge[, 1], auc1.grridge[, 2])
# lauc2.grridge <- lowess(auc2.grridge[, 1], auc2.grridge[, 2])
# lauc3.grridge <- lowess(auc3.grridge[, 1], auc3.grridge[, 2])
# lauc4.grridge <- lowess(auc4.grridge[, 1], auc4.grridge[, 2])
# 
# lauc1.ridge <- lowess(auc1.ridge[, 1], auc1.ridge[, 2])
# lauc2.ridge <- lowess(auc2.ridge[, 1], auc2.ridge[, 2])
# 
# briers1.enet <- do.call(rbind, results3$briers[[7]])[
#   order(do.call(rbind, results3$briers[[7]])[, 1]), ]
# briers2.enet <- do.call(rbind, results3$briers[[8]])[
#   order(do.call(rbind, results3$briers[[8]])[, 1]), ]
# briers3.enet <- do.call(rbind, results3$briers[[9]])[
#   order(do.call(rbind, results3$briers[[9]])[, 1]), ]
# briers4.enet <- do.call(rbind, results3$briers[[10]])[
#   order(do.call(rbind, results3$briers[[10]])[, 1]), ]
# 
# briers1.greben <- do.call(rbind, results3$briers[[11]])[
#   order(do.call(rbind, results3$briers[[11]])[, 1]), ]
# briers2.greben <- do.call(rbind, results3$briers[[12]])[
#   order(do.call(rbind, results3$briers[[12]])[, 1]), ]
# briers3.greben <- do.call(rbind, results3$briers[[13]])[
#   order(do.call(rbind, results3$briers[[13]])[, 1]), ]
# briers4.greben <- do.call(rbind, results3$briers[[14]])[
#   order(do.call(rbind, results3$briers[[14]])[, 1]), ]
# 
# briers1.grridge <- do.call(rbind, results3$briers[[3]])[
#   order(do.call(rbind, results3$briers[[3]])[, 1]), ]
# briers2.grridge <- do.call(rbind, results3$briers[[4]])[
#   order(do.call(rbind, results3$briers[[4]])[, 1]), ]
# briers3.grridge <- do.call(rbind, results3$briers[[5]])[
#   order(do.call(rbind, results3$briers[[5]])[, 1]), ]
# briers4.grridge <- do.call(rbind, results3$briers[[6]])[
#   order(do.call(rbind, results3$briers[[6]])[, 1]), ]
# 
# briers1.ridge <- do.call(rbind, results3$briers[[1]])[
#   order(do.call(rbind, results3$briers[[1]])[, 1]), ]
# briers2.ridge <- do.call(rbind, results3$briers[[2]])[
#   order(do.call(rbind, results3$briers[[2]])[, 1]), ]
# 
# lbriers1.enet <- lowess(briers1.enet[, 1], briers1.enet[, 2])
# lbriers2.enet <- lowess(briers2.enet[, 1], briers2.enet[, 2])
# lbriers3.enet <- lowess(briers3.enet[, 1], briers3.enet[, 2])
# lbriers4.enet <- lowess(briers4.enet[, 1], briers4.enet[, 2])
# 
# lbriers1.greben <- lowess(briers1.greben[, 1], briers1.greben[, 2])
# lbriers2.greben <- lowess(briers2.greben[, 1], briers2.greben[, 2])
# lbriers3.greben <- lowess(briers3.greben[, 1], briers3.greben[, 2])
# lbriers4.greben <- lowess(briers4.greben[, 1], briers4.greben[, 2])
# 
# lbriers1.grridge <- lowess(briers1.grridge[, 1], briers1.grridge[, 2])
# lbriers2.grridge <- lowess(briers2.grridge[, 1], briers2.grridge[, 2])
# lbriers3.grridge <- lowess(briers3.grridge[, 1], briers3.grridge[, 2])
# lbriers4.grridge <- lowess(briers4.grridge[, 1], briers4.grridge[, 2])
# 
# lbriers1.ridge <- lowess(briers1.ridge[, 1], briers1.ridge[, 2])
# lbriers2.ridge <- lowess(briers2.ridge[, 1], briers2.ridge[, 2])
# 
# mse1.enet <- do.call(rbind, results3$mse[[7]])[
#   order(do.call(rbind, results3$mse[[7]])[, 1]), ]
# mse2.enet <- do.call(rbind, results3$mse[[8]])[
#   order(do.call(rbind, results3$mse[[8]])[, 1]), ]
# mse3.enet <- do.call(rbind, results3$mse[[9]])[
#   order(do.call(rbind, results3$mse[[9]])[, 1]), ]
# mse4.enet <- do.call(rbind, results3$mse[[10]])[
#   order(do.call(rbind, results3$mse[[10]])[, 1]), ]
# 
# mse1.greben <- do.call(rbind, results3$mse[[11]])[
#   order(do.call(rbind, results3$mse[[11]])[, 1]), ]
# mse2.greben <- do.call(rbind, results3$mse[[12]])[
#   order(do.call(rbind, results3$mse[[12]])[, 1]), ]
# mse3.greben <- do.call(rbind, results3$mse[[13]])[
#   order(do.call(rbind, results3$mse[[13]])[, 1]), ]
# mse4.greben <- do.call(rbind, results3$mse[[14]])[
#   order(do.call(rbind, results3$mse[[14]])[, 1]), ]
# 
# mse1.grridge <- do.call(rbind, results3$mse[[3]])[
#   order(do.call(rbind, results3$mse[[3]])[, 1]), ]
# mse2.grridge <- do.call(rbind, results3$mse[[4]])[
#   order(do.call(rbind, results3$mse[[4]])[, 1]), ]
# mse3.grridge <- do.call(rbind, results3$mse[[5]])[
#   order(do.call(rbind, results3$mse[[5]])[, 1]), ]
# mse4.grridge <- do.call(rbind, results3$mse[[6]])[
#   order(do.call(rbind, results3$mse[[6]])[, 1]), ]
# 
# mse1.ridge <- do.call(rbind, results3$mse[[1]])[
#   order(do.call(rbind, results3$mse[[1]])[, 1]), ]
# mse2.ridge <- do.call(rbind, results3$mse[[2]])[
#   order(do.call(rbind, results3$mse[[2]])[, 1]), ]
# 
# lmse1.enet <- lowess(mse1.enet[, 1], mse1.enet[, 2])
# lmse2.enet <- lowess(mse2.enet[, 1], mse2.enet[, 2])
# lmse3.enet <- lowess(mse3.enet[, 1], mse3.enet[, 2])
# lmse4.enet <- lowess(mse4.enet[, 1], mse4.enet[, 2])
# 
# lmse1.greben <- lowess(mse1.greben[, 1], mse1.greben[, 2])
# lmse2.greben <- lowess(mse2.greben[, 1], mse2.greben[, 2])
# lmse3.greben <- lowess(mse3.greben[, 1], mse3.greben[, 2])
# lmse4.greben <- lowess(mse4.greben[, 1], mse4.greben[, 2])
# 
# lmse1.grridge <- lowess(mse1.grridge[, 1], mse1.grridge[, 2])
# lmse2.grridge <- lowess(mse2.grridge[, 1], mse2.grridge[, 2])
# lmse3.grridge <- lowess(mse3.grridge[, 1], mse3.grridge[, 2])
# lmse4.grridge <- lowess(mse4.grridge[, 1], mse4.grridge[, 2])
# 
# lmse1.ridge <- lowess(mse1.ridge[, 1], mse1.ridge[, 2])
# lmse2.ridge <- lowess(mse2.ridge[, 1], mse2.ridge[, 2])
# 
# ### Performance measures
# # AUC
# png(paste(path.graph, "grEBEN_test_res3_performance.png", sep=""),
#     units="in", width=12, height=4, res=120)
# par(mfrow=c(1, 3))
# xlim1.auc <- range(lauc1.enet$x, lauc2.enet$x, lauc3.enet$x, lauc4.enet$x,
#                    lauc1.greben$x, lauc2.greben$x, lauc3.greben$x,
#                    lauc4.greben$x, lauc3.grridge$x, lauc4.grridge$x)
# ylim1.auc <- range(lauc1.enet$y, lauc2.enet$y, lauc3.enet$y, lauc4.enet$y,
#                    lauc1.greben$y, lauc2.greben$y, lauc3.greben$y,
#                    lauc4.greben$y, lauc1.grridge$y, lauc2.grridge$y,
#                    lauc3.grridge$y, lauc4.grridge$y, lauc1.ridge$y,
#                    lauc2.ridge$y)
# 
# plot(0, 0, col=2, type="n", ylab="AUC", main="a)",
#      xlab="Number of selected variables", xlim=xlim1.auc, ylim=ylim1.auc)
# 
# lines(lauc3.grridge, col=2, lty=1)
# abline(h=lauc1.grridge$y, col=2, lty=2)
# abline(h=lauc1.ridge$y, col=2, lty=3)
# 
# lines(lauc4.grridge, col=3, lty=1)
# abline(h=lauc2.grridge$y, col=3, lty=2)
# abline(h=lauc2.ridge$y, col=3, lty=3)
# 
# lines(lauc1.greben, col=4, lty=1)
# lines(lauc1.enet, col=4, lty=3)
# 
# lines(lauc2.greben, col=5, lty=1)
# lines(lauc2.enet, col=5, lty=3)
# 
# lines(lauc3.greben, col=6, lty=1)
# lines(lauc3.enet, col=6, lty=3)
# 
# lines(lauc4.greben, col=7, lty=1)
# lines(lauc4.enet, col=7, lty=3)
# 
# # Brier skilll
# xlim1.briers <- range(lbriers1.enet$x, lbriers2.enet$x, lbriers3.enet$x,
#                       lbriers4.enet$x, lbriers1.greben$x, lbriers2.greben$x,
#                       lbriers3.greben$x, lbriers4.greben$x, lbriers3.grridge$x,
#                       lbriers4.grridge$x)
# ylim1.briers <- range(lbriers1.enet$y, lbriers2.enet$y, lbriers3.enet$y,
#                       lbriers4.enet$y, lbriers1.greben$y, lbriers2.greben$y,
#                       lbriers3.greben$y, lbriers4.greben$y, lbriers1.grridge$y,
#                       lbriers2.grridge$y, lbriers3.grridge$y, 
#                       lbriers4.grridge$y, lbriers1.ridge$y, lbriers2.ridge$y)
# 
# plot(0, 0, col=2, type="n", ylab="Brier skill score", main="b)",
#      xlab="Number of selected variables", xlim=xlim1.briers, ylim=ylim1.briers)
# 
# lines(lbriers3.grridge, col=2, lty=1)
# abline(h=lbriers1.grridge$y, col=2, lty=2)
# abline(h=lbriers1.ridge$y, col=2, lty=3)
# 
# lines(lbriers4.grridge, col=3, lty=1)
# abline(h=lbriers2.grridge$y, col=3, lty=2)
# abline(h=lbriers2.ridge$y, col=3, lty=3)
# 
# lines(lbriers1.greben, col=4, lty=1)
# lines(lbriers1.enet, col=4, lty=3)
# 
# lines(lbriers2.greben, col=5, lty=1)
# lines(lbriers2.enet, col=5, lty=3)
# 
# lines(lbriers3.greben, col=6, lty=1)
# lines(lbriers3.enet, col=6, lty=3)
# 
# lines(lbriers4.greben, col=7, lty=1)
# lines(lbriers4.enet, col=7, lty=3)
# 
# leglabels <- c(expression(paste("ridge, 'true' ", lambda)),
#                "ridge", expression(paste("enet, true ", lambda)),
#                expression(paste("enet, ", alpha==0.05)),
#                expression(paste("enet, ", alpha==0.5)),
#                expression(paste("enet, ", alpha==0.95)),
#                "group-regularized + selection", "group-regularized",
#                "not group-regularized")
# legend("bottomright", legend=leglabels, fill=c(2:7, 0, 0, 0),
#        lty=c(rep(NA, 6), 1, 2, 3), border=c(rep(1, 6), 0 ,0, 0), merge=TRUE,
#        seg.len=1)
# 
# # MSE
# xlim1.mse <- range(lmse1.enet$x, lmse2.enet$x, lmse3.enet$x,
#                       lmse4.enet$x, lmse1.greben$x, lmse2.greben$x,
#                       lmse3.greben$x, lmse4.greben$x, lmse3.grridge$x,
#                       lmse4.grridge$x)
# ylim1.mse <- range(lmse1.enet$y, lmse2.enet$y, lmse3.enet$y,
#                       lmse4.enet$y, lmse1.greben$y, lmse2.greben$y,
#                       lmse3.greben$y, lmse4.greben$y, lmse1.grridge$y,
#                       lmse2.grridge$y, lmse3.grridge$y, lmse4.grridge$y,
#                       lmse1.ridge$y, lmse2.ridge$y)
# 
# plot(0, 0, col=2, type="n", ylab="MSE", main="c)",
#      xlab="Number of selected variables", xlim=xlim1.mse, ylim=ylim1.mse)
# 
# lines(lmse3.grridge, col=2, lty=1)
# abline(h=lmse1.grridge$y, col=2, lty=2)
# abline(h=lmse1.ridge$y, col=2, lty=3)
# 
# lines(lmse4.grridge, col=3, lty=1)
# abline(h=lmse2.grridge$y, col=3, lty=2)
# abline(h=lmse2.ridge$y, col=3, lty=3)
# 
# lines(lmse1.greben, col=4, lty=1)
# lines(lmse1.enet, col=4, lty=3)
# 
# lines(lmse2.greben, col=5, lty=1)
# lines(lmse2.enet, col=5, lty=3)
# 
# lines(lmse3.greben, col=6, lty=1)
# lines(lmse3.enet, col=6, lty=3)
# 
# lines(lmse4.greben, col=7, lty=1)
# lines(lmse4.enet, col=7, lty=3)
# 
# dev.off()
# 
# ### diagnostics: checking lowess fits
# # AUC
# xlim2.auc <- range(auc1.enet[, 1], auc2.enet[, 1], auc3.enet[, 1], auc4.enet[, 1],
#                    auc1.greben[, 1], auc2.greben[, 1], auc3.greben[, 1],
#                    auc4.greben[, 1], auc3.grridge[, 1], auc4.grridge[, 1])
# ylim2.auc <- range(auc1.enet[, 2], auc2.enet[, 2], auc3.enet[, 2], auc4.enet[, 2],
#                    auc1.greben[, 2], auc2.greben[, 2], auc3.greben[, 2],
#                    auc4.greben[, 2], auc3.grridge[, 2], auc4.grridge[, 2])
# 
# png(paste(path.graph, "grEBEN_test_res3_lowess_auc.png", sep=""),
#     units="in", width=12, height=6, res=120)
# par(mfrow=c(2, 5))
# plot(auc3.grridge[, 1], auc3.grridge[, 2], col=1, ylab="AUC", main="a)",
#      xlab="Number of selected variables", xlim=xlim2.auc, ylim=ylim2.auc)
# lines(lauc3.grridge, col=2)
# 
# plot(auc1.greben[, 1], auc1.greben[, 2], col=1, ylab="AUC", main="b)",
#      xlab="Number of selected variables", xlim=xlim2.auc, ylim=ylim2.auc)
# lines(lauc1.greben, col=4)
# 
# plot(auc2.greben[, 1], auc2.greben[, 2], col=1, ylab="AUC", main="c)",
#      xlab="Number of selected variables", xlim=xlim2.auc, ylim=ylim2.auc)
# lines(lauc2.greben, col=5)
# 
# plot(auc3.greben[, 1], auc4.greben[, 2], col=1, ylab="AUC", main="d)",
#      xlab="Number of selected variables", xlim=xlim2.auc, ylim=ylim2.auc)
# lines(lauc3.greben, col=6)
# 
# plot(auc4.greben[, 1], auc4.greben[, 2], col=1, ylab="AUC", main="e)",
#      xlab="Number of selected variables", xlim=xlim2.auc, ylim=ylim2.auc)
# lines(lauc4.greben, col=7)
# 
# plot(auc4.grridge[, 1], auc4.grridge[, 2], col=1, ylab="AUC", main="f)",
#      xlab="Number of selected variables", xlim=xlim2.auc, ylim=ylim2.auc)
# lines(lauc4.grridge, col=3)
# 
# plot(auc1.enet[, 1], auc1.enet[, 2], col=1, ylab="AUC", main="g)",
#      xlab="Number of selected variables", xlim=xlim2.auc, ylim=ylim2.auc)
# lines(lauc1.enet, col=4)
# 
# plot(auc2.enet[, 1], auc2.enet[, 2], col=1, ylab="AUC", main="h)",
#      xlab="Number of selected variables", xlim=xlim2.auc, ylim=ylim2.auc)
# lines(lauc2.enet, col=5)
# 
# plot(auc3.enet[, 1], auc4.enet[, 2], col=1, ylab="AUC", main="i)",
#      xlab="Number of selected variables", xlim=xlim2.auc, ylim=ylim2.auc)
# lines(lauc3.enet, col=6)
# 
# plot(auc4.enet[, 1], auc4.enet[, 2], col=1, ylab="AUC", main="j)",
#      xlab="Number of selected variables", xlim=xlim2.auc, ylim=ylim2.auc)
# lines(lauc4.enet, col=7)
# dev.off()
# 
# # briers
# xlim2.briers <- range(briers1.enet[, 1], briers2.enet[, 1], briers3.enet[, 1], briers4.enet[, 1],
#                    briers1.greben[, 1], briers2.greben[, 1], briers3.greben[, 1],
#                    briers4.greben[, 1], briers3.grridge[, 1], briers4.grridge[, 1])
# ylim2.briers <- range(briers1.enet[, 2], briers2.enet[, 2], briers3.enet[, 2], briers4.enet[, 2],
#                    briers1.greben[, 2], briers2.greben[, 2], briers3.greben[, 2],
#                    briers4.greben[, 2], briers3.grridge[, 2], briers4.grridge[, 2])
# 
# png(paste(path.graph, "grEBEN_test_res3_lowess_briers.png", sep=""),
#     units="in", width=12, height=6, res=120)
# par(mfrow=c(2, 5))
# plot(briers3.grridge[, 1], briers3.grridge[, 2], col=1, ylab="Brier skill score", main="a)",
#      xlab="Number of selected variables", xlim=xlim2.briers, ylim=ylim2.briers)
# lines(lbriers3.grridge, col=2)
# 
# plot(briers1.greben[, 1], briers1.greben[, 2], col=1, ylab="Brier skill score", main="b)",
#      xlab="Number of selected variables", xlim=xlim2.briers, ylim=ylim2.briers)
# lines(lbriers1.greben, col=4)
# 
# plot(briers2.greben[, 1], briers2.greben[, 2], col=1, ylab="Brier skill score", main="c)",
#      xlab="Number of selected variables", xlim=xlim2.briers, ylim=ylim2.briers)
# lines(lbriers2.greben, col=5)
# 
# plot(briers3.greben[, 1], briers4.greben[, 2], col=1, ylab="Brier skill score", main="d)",
#      xlab="Number of selected variables", xlim=xlim2.briers, ylim=ylim2.briers)
# lines(lbriers3.greben, col=6)
# 
# plot(briers4.greben[, 1], briers4.greben[, 2], col=1, ylab="Brier skill score", main="e)",
#      xlab="Number of selected variables", xlim=xlim2.briers, ylim=ylim2.briers)
# lines(lbriers4.greben, col=7)
# 
# plot(briers4.grridge[, 1], briers4.grridge[, 2], col=1, ylab="Brier skill score", main="f)",
#      xlab="Number of selected variables", xlim=xlim2.briers, ylim=ylim2.briers)
# lines(lbriers4.grridge, col=3)
# 
# plot(briers1.enet[, 1], briers1.enet[, 2], col=1, ylab="Brier skill score", main="g)",
#      xlab="Number of selected variables", xlim=xlim2.briers, ylim=ylim2.briers)
# lines(lbriers1.enet, col=4)
# 
# plot(briers2.enet[, 1], briers2.enet[, 2], col=1, ylab="Brier skill score", main="h)",
#      xlab="Number of selected variables", xlim=xlim2.briers, ylim=ylim2.briers)
# lines(lbriers2.enet, col=5)
# 
# plot(briers3.enet[, 1], briers4.enet[, 2], col=1, ylab="Brier skill score", main="i)",
#      xlab="Number of selected variables", xlim=xlim2.briers, ylim=ylim2.briers)
# lines(lbriers3.enet, col=6)
# 
# plot(briers4.enet[, 1], briers4.enet[, 2], col=1, ylab="Brier skill score", main="j)",
#      xlab="Number of selected variables", xlim=xlim2.briers, ylim=ylim2.briers)
# lines(lbriers4.enet, col=7)
# dev.off()
# 
# # mse
# xlim2.mse <- range(mse1.enet[, 1], mse2.enet[, 1], mse3.enet[, 1], mse4.enet[, 1],
#                    mse1.greben[, 1], mse2.greben[, 1], mse3.greben[, 1],
#                    mse4.greben[, 1], mse3.grridge[, 1], mse4.grridge[, 1])
# ylim2.mse <- range(mse1.enet[, 2], mse2.enet[, 2], mse3.enet[, 2], mse4.enet[, 2],
#                    mse1.greben[, 2], mse2.greben[, 2], mse3.greben[, 2],
#                    mse4.greben[, 2], mse3.grridge[, 2], mse4.grridge[, 2])
# 
# png(paste(path.graph, "grEBEN_test_res3_lowess_mse.png", sep=""),
#     units="in", width=12, height=6, res=120)
# par(mfrow=c(2, 5))
# plot(mse3.grridge[, 1], mse3.grridge[, 2], col=1, ylab="MSE", main="a)",
#      xlab="Number of selected variables", xlim=xlim2.mse, ylim=ylim2.mse)
# lines(lmse3.grridge, col=2)
# 
# plot(mse1.greben[, 1], mse1.greben[, 2], col=1, ylab="MSE", main="b)",
#      xlab="Number of selected variables", xlim=xlim2.mse, ylim=ylim2.mse)
# lines(lmse1.greben, col=4)
# 
# plot(mse2.greben[, 1], mse2.greben[, 2], col=1, ylab="MSE", main="c)",
#      xlab="Number of selected variables", xlim=xlim2.mse, ylim=ylim2.mse)
# lines(lmse2.greben, col=5)
# 
# plot(mse3.greben[, 1], mse4.greben[, 2], col=1, ylab="MSE", main="d)",
#      xlab="Number of selected variables", xlim=xlim2.mse, ylim=ylim2.mse)
# lines(lmse3.greben, col=6)
# 
# plot(mse4.greben[, 1], mse4.greben[, 2], col=1, ylab="MSE", main="e)",
#      xlab="Number of selected variables", xlim=xlim2.mse, ylim=ylim2.mse)
# lines(lmse4.greben, col=7)
# 
# plot(mse4.grridge[, 1], mse4.grridge[, 2], col=1, ylab="MSE", main="f)",
#      xlab="Number of selected variables", xlim=xlim2.mse, ylim=ylim2.mse)
# lines(lmse4.grridge, col=3)
# 
# plot(mse1.enet[, 1], mse1.enet[, 2], col=1, ylab="MSE", main="g)",
#      xlab="Number of selected variables", xlim=xlim2.mse, ylim=ylim2.mse)
# lines(lmse1.enet, col=4)
# 
# plot(mse2.enet[, 1], mse2.enet[, 2], col=1, ylab="MSE", main="h)",
#      xlab="Number of selected variables", xlim=xlim2.mse, ylim=ylim2.mse)
# lines(lmse2.enet, col=5)
# 
# plot(mse3.enet[, 1], mse4.enet[, 2], col=1, ylab="MSE", main="i)",
#      xlab="Number of selected variables", xlim=xlim2.mse, ylim=ylim2.mse)
# lines(lmse3.enet, col=6)
# 
# plot(mse4.enet[, 1], mse4.enet[, 2], col=1, ylab="MSE", main="j)",
#      xlab="Number of selected variables", xlim=xlim2.mse, ylim=ylim2.mse)
# lines(lmse4.enet, col=7)
# dev.off()

## penalty parameters
lambdag <- exp(seq(-1, 1, length.out=ncol(results3$lambdag[[1]])))
png(paste(path.graph, "grEBEN_test_res3_penalties.png", sep=""),
    units="in", width=8, height=6, res=120)
par(mfrow=c(2, 3))
boxplot(results3$lambdag[[1]], ylim=range(results3$lambdag[[1]], lambdag), 
        xlab="", ylab="", main="a)")
title(ylab=expression(paste(lambda[g], "'")), 
      xlab="Groups in decreasing effect size", line=2.5)
points(1:length(lambdag), lambdag, col=2, pch="x")
boxplot(results3$lambdag[[2]], ylim=range(results3$lambdag[[2]], lambdag), 
        xlab="", ylab="", main="b)")
title(ylab=expression(paste(lambda[g], "'")), 
      xlab="Groups in decreasing effect size", line=2.5)
points(1:length(lambdag), lambdag, col=2, pch="x")
boxplot(results3$lambdag[[3]], ylim=range(results3$lambdag[[3]], lambdag), 
        xlab="", ylab="", main="c)")
title(ylab=expression(paste(lambda[g], "'")), 
      xlab="Groups in decreasing effect size", line=2.5)
points(1:length(lambdag), lambdag, col=2, pch="x")
boxplot(results3$lambdag[[4]], ylim=range(results3$lambdag[[4]], lambdag), 
        xlab="", ylab="", main="d)")
title(ylab=expression(paste(lambda[g], "'")), 
      xlab="Groups in decreasing effect size", line=2.5)
points(1:length(lambdag), lambdag, col=2, pch="x")
boxplot(results3$lambdag[[5]], ylim=range(results3$lambdag[[5]], lambdag), 
        xlab="", ylab="", main="e)")
title(ylab=expression(paste(lambda[g], "'")), 
      xlab="Groups in decreasing effect size", line=2.5)
points(1:length(lambdag), lambdag, col=2, pch="x")
boxplot(results3$lambdag[[6]], ylim=range(results3$lambdag[[6]], lambdag), 
        xlab="", ylab="", main="f)")
title(ylab=expression(paste(lambda[g], "'")), 
      xlab="Groups in decreasing effect size", line=2.5)
points(1:length(lambdag), lambdag, col=2, pch="x")
dev.off()





# ## simulations
# ## simulation 4
# # create data
# n <- 200
# p <- 1000
# G <- 5
# pblock <- 20
# rho <- 0.7
# sigma <- matrix(rho, ncol=pblock, nrow=pblock); diag(sigma) <- 1
# lambda <- 0.02
# alpha <- 1
# lambdag <- exp(seq(-1, 1, length.out=G))
# m <- rep(1, n)
# part.greben <- list(groups=rep(1:G, each=p/G))
# part.grridge <- list(groups=CreatePartition(as.factor(part.greben$groups)))
# ntest <- 1000
# 
# methods <- c("ridge+truelambda", "ridge", "GRridge+truelambda", "GRridge",
#              "GRridge+truelambda+sel", "GRridge+sel",
#              "enet+truelambda", "enet+a=0.05", "enet+a=0.5", "enet+a=0.95",
#              "grEBEN+truelambda", "grEBEN+a=0.05", "grEBEN+a=0.5","grEBEN+a=0.95")
# nreps <- 50
# auc4 <- briers4 <- mse4 <- rep(list(vector(mode="list", length=nreps)),
#                                length(methods))
# names(auc4) <- names(briers4) <- names(mse4) <- methods
# lambdag4 <- rep(list(matrix(NA, ncol=G, nrow=nreps)), 6)
# names(lambdag4) <- methods[c(3, 4, 11, 12, 13, 14)]
# varbeta4 <- matrix(nrow=nreps, ncol=G)
# # the simulations
# for(r in 1:nreps) {
#   
#   print(paste("Simulation 3, repetition", r))
#   set.seed(400 + r)
#   x <- do.call(cbind, replicate(p/pblock, rmvnorm(n, mean=rep(0, pblock), sigma=sigma),
#                                 simplify=FALSE))
#   beta <- as.numeric(sapply(1:G, function(g) {
#     renbeta(p/G, 2*n*lambda*alpha*sqrt(lambdag[g]), n*lambda*(1 - alpha)*lambdag[g])}))
#   prob <- as.numeric(exp(x %*% beta)/(1 + exp(x %*% beta)))
#   y <- rbinom(n, 1, prob)
#   optl <- n*lambda*alpha*sum(abs(beta))/sum(beta^2) + 0.5*n*lambda*(1 - alpha)
#   
#   xtest <- do.call(cbind, replicate(p/pblock, rmvnorm(ntest, mean=rep(0, pblock),
#                                                       sigma=sigma), simplify=FALSE))
#   probtest <- as.numeric(exp(xtest %*% beta)/(1 + exp(xtest %*% beta)))
#   ytest <- rbinom(ntest, 1, probtest)
#   
#   fit1.greben <- grEBEN3(x, y, m, partitions=part.greben, alpha=alpha - 0.01,
#                          lambda=lambda, psel=TRUE)
#   fit2.greben <- grEBEN3(x, y, m, partitions=part.greben, alpha=0.05, psel=TRUE)
#   fit3.greben <- grEBEN3(x, y, m, partitions=part.greben, alpha=0.5, psel=TRUE)
#   fit4.greben <- grEBEN3(x, y, m, partitions=part.greben, alpha=0.95, psel=TRUE)
#   
#   # number of selected variables
#   psel1.enet <- apply(fit1.greben$beta.nogroups, 2, function(b) {sum(b!=0) - 1})
#   psel2.enet <- apply(fit2.greben$beta.nogroups, 2, function(b) {sum(b!=0) - 1})
#   psel3.enet <- apply(fit3.greben$beta.nogroups, 2, function(b) {sum(b!=0) - 1})
#   psel4.enet <- apply(fit4.greben$beta.nogroups, 2, function(b) {sum(b!=0) - 1})
#   
#   psel1.greben <- apply(fit1.greben$beta, 2, function(b) {sum(b!=0) - 1})
#   psel2.greben <- apply(fit2.greben$beta, 2, function(b) {sum(b!=0) - 1})
#   psel3.greben <- apply(fit3.greben$beta, 2, function(b) {sum(b!=0) - 1})
#   psel4.greben <- apply(fit4.greben$beta, 2, function(b) {sum(b!=0) - 1})
#   
#   psel.all <- unique(c(psel1.enet, psel2.enet, psel3.enet, psel4.enet,
#                        psel1.greben, psel2.greben, psel3.greben, psel4.greben))
#   psel.in <- floor(quantile(psel.all[psel.all!=0],
#                             prob=c(seq(0.01, 0.05, 0.01), seq(0.07, 0.25, 0.03), seq(0.3, 1, 0.1))))
#   
#   fit1.grridge <- vector(mode="list", length=length(psel.in))
#   for(s in 1:length(psel.in)) {
#     fit1.grridge[[s]] <- grridge(t(x), y, part.grridge, selectionEN=TRUE,
#                                  optl=optl, maxsel=psel.in[s])
#   }
#   
#   fit2.grridge <- vector(mode="list", length=length(psel.in))
#   fit2.grridge[[1]] <- grridge(t(x), y, part.grridge, selectionEN=TRUE,
#                                maxsel=psel.in[1])
#   for(s in 2:length(psel.in)) {
#     fit2.grridge[[s]] <- grridge(t(x), y, part.grridge, selectionEN=TRUE,
#                                  optl=fit2.grridge[[1]]$optl,
#                                  maxsel=psel.in[s])
#   }
#   
#   psel1.grridge <- sapply(1:length(fit1.grridge), function(s) {
#     return(length(fit1.grridge[[s]]$resEN$whichEN))})
#   psel2.grridge <- sapply(1:length(fit2.grridge), function(s) {
#     return(length(fit2.grridge[[s]]$resEN$whichEN))})
#   
#   # estimates
#   est1.ridge <- coef(fit1.grridge[[1]]$predobj$NoGroups, "all")
#   est2.ridge <- coef(fit2.grridge[[1]]$predobj$NoGroups, "all")
#   
#   est1.grridge <- coef(fit1.grridge[[1]]$predobj$GroupRegul, "all")
#   est2.grridge <- coef(fit2.grridge[[1]]$predobj$GroupRegul, "all")
#   
#   est3.grridge <- sapply(fit1.grridge, function(s) {
#     replace(rep(0, p + 1), c(1, s$resEN$whichEN + 1), coef(s$predobj$EN, "all"))})
#   est4.grridge <- sapply(fit2.grridge, function(s) {
#     replace(rep(0, p + 1), c(1, s$resEN$whichEN + 1), coef(s$predobj$EN, "all"))})
#   
#   est1.enet <- fit1.greben$beta.nogroups
#   est2.enet <- fit2.greben$beta.nogroups
#   est3.enet <- fit3.greben$beta.nogroups
#   est4.enet <- fit4.greben$beta.nogroups
#   
#   est1.greben <- fit1.greben$beta
#   est2.greben <- fit2.greben$beta
#   est3.greben <- fit3.greben$beta
#   est4.greben <- fit4.greben$beta
#   
#   # predictions on fit data
#   pred1.ridge <- predict.grridge(fit1.grridge[[1]], t(xtest))[, 1]
#   pred2.ridge <- predict.grridge(fit2.grridge[[1]], t(xtest))[, 1]
#   
#   pred1.grridge <- predict.grridge(fit1.grridge[[1]], t(xtest))[, 2]
#   pred2.grridge <- predict.grridge(fit2.grridge[[1]], t(xtest))[, 2]
#   pred3.grridge <- sapply(fit1.grridge, function(s) {
#     predict.grridge(s, t(xtest))[, 3]})
#   pred4.grridge <- sapply(fit2.grridge, function(s) {
#     predict.grridge(s, t(xtest))[, 3]})
#   
#   pred1.enet <- 1/(1 + exp(-xtest %*% est1.enet[-1, ]))
#   pred2.enet <- 1/(1 + exp(-xtest %*% est2.enet[-1, ]))
#   pred3.enet <- 1/(1 + exp(-xtest %*% est3.enet[-1, ]))
#   pred4.enet <- 1/(1 + exp(-xtest %*% est4.enet[-1, ]))
#   
#   pred1.greben <- 1/(1 + exp(-xtest %*% est1.greben[-1, ]))
#   pred2.greben <- 1/(1 + exp(-xtest %*% est2.greben[-1, ]))
#   pred3.greben <- 1/(1 + exp(-xtest %*% est3.greben[-1, ]))
#   pred4.greben <- 1/(1 + exp(-xtest %*% est4.greben[-1, ]))
#   
#   # AUCs
#   auc.true <- pROC::roc(ytest, probtest)$auc
#   
#   auc1.ridge <- pROC::roc(ytest, pred1.ridge)$auc
#   auc2.ridge <- pROC::roc(ytest, pred2.ridge)$auc
#   
#   auc1.grridge <- pROC::roc(ytest, pred1.grridge)$auc
#   auc2.grridge <- pROC::roc(ytest, pred2.grridge)$auc
#   auc3.grridge <- apply(pred3.grridge, 2, function(r) {pROC::roc(ytest, r)$auc})
#   auc4.grridge <- apply(pred4.grridge, 2, function(r) {pROC::roc(ytest, r)$auc})
#   
#   auc1.enet <- apply(pred1.enet, 2, function(r) {pROC::roc(ytest, r)$auc})
#   auc2.enet <- apply(pred2.enet, 2, function(r) {pROC::roc(ytest, r)$auc})
#   auc3.enet <- apply(pred3.enet, 2, function(r) {pROC::roc(ytest, r)$auc})
#   auc4.enet <- apply(pred4.enet, 2, function(r) {pROC::roc(ytest, r)$auc})
#   
#   auc1.greben <- apply(pred1.greben, 2, function(r) {pROC::roc(ytest, r)$auc})
#   auc2.greben <- apply(pred2.greben, 2, function(r) {pROC::roc(ytest, r)$auc})
#   auc3.greben <- apply(pred3.greben, 2, function(r) {pROC::roc(ytest, r)$auc})
#   auc4.greben <- apply(pred4.greben, 2, function(r) {pROC::roc(ytest, r)$auc})
#   
#   # Brier scores
#   brier.null <- sum((ytest - mean(ytest))^2)
#   briers.true <- 1 - sum((ytest - probtest)^2)/brier.null
#   
#   briers1.ridge <- 1 - sum((ytest - pred1.ridge)^2)/brier.null
#   briers2.ridge <- 1 - sum((ytest - pred2.ridge)^2)/brier.null
#   
#   briers1.grridge <- 1 - sum((ytest - pred1.grridge)^2)/brier.null
#   briers2.grridge <- 1 - sum((ytest - pred2.grridge)^2)/brier.null
#   briers3.grridge <- apply(pred3.grridge, 2, function(pred) {
#     1 - sum((ytest - pred)^2)/brier.null})
#   briers4.grridge <- apply(pred4.grridge, 2, function(pred) {
#     1 - sum((ytest - pred)^2)/brier.null})
#   
#   briers1.enet <- apply(pred1.enet, 2, function(pred) {
#     1 - sum((ytest - pred)^2)/brier.null})
#   briers2.enet <- apply(pred2.enet, 2, function(pred) {
#     1 - sum((ytest - pred)^2)/brier.null})
#   briers3.enet <- apply(pred3.enet, 2, function(pred) {
#     1 - sum((ytest - pred)^2)/brier.null})
#   briers4.enet <- apply(pred4.enet, 2, function(pred) {
#     1 - sum((ytest - pred)^2)/brier.null})
#   
#   briers1.greben <- apply(pred1.greben, 2, function(pred) {
#     1 - sum((ytest - pred)^2)/brier.null})
#   briers2.greben <- apply(pred2.greben, 2, function(pred) {
#     1 - sum((ytest - pred)^2)/brier.null})
#   briers3.greben <- apply(pred3.greben, 2, function(pred) {
#     1 - sum((ytest - pred)^2)/brier.null})
#   briers4.greben <- apply(pred4.greben, 2, function(pred) {
#     1 - sum((ytest - pred)^2)/brier.null})
#   
#   # MSE
#   mse.true <- 0
#   
#   mse1.ridge <- mean((c(0, beta) - est1.ridge)^2)
#   mse2.ridge <- mean((c(0, beta) - est2.ridge)^2)
#   
#   mse1.grridge <- mean((c(0, beta) - est1.grridge)^2)
#   mse2.grridge <- mean((c(0, beta) - est2.grridge)^2)
#   mse3.grridge <- apply(est3.grridge, 2, function(b) {mean((c(0, beta) - b)^2)})
#   mse4.grridge <- apply(est4.grridge, 2, function(b) {mean((c(0, beta) - b)^2)})
#   
#   mse1.enet <- apply(est1.enet, 2, function(b) {mean((c(0, beta) - b)^2)})
#   mse2.enet <- apply(est2.enet, 2, function(b) {mean((c(0, beta) - b)^2)})
#   mse3.enet <- apply(est3.enet, 2, function(b) {mean((c(0, beta) - b)^2)})
#   mse4.enet <- apply(est4.enet, 2, function(b) {mean((c(0, beta) - b)^2)})
#   
#   mse1.greben <- apply(est1.greben, 2, function(b) {mean((c(0, beta) - b)^2)})
#   mse2.greben <- apply(est2.greben, 2, function(b) {mean((c(0, beta) - b)^2)})
#   mse3.greben <- apply(est3.greben, 2, function(b) {mean((c(0, beta) - b)^2)})
#   mse4.greben <- apply(est4.greben, 2, function(b) {mean((c(0, beta) - b)^2)})
#   
#   auc4[[1]][[r]] <- cbind(psel=p, auc=auc1.ridge)
#   auc4[[2]][[r]] <- cbind(psel=p, auc=auc2.ridge)
#   auc4[[3]][[r]] <- cbind(psel=p, auc=auc1.grridge)
#   auc4[[4]][[r]] <- cbind(psel=p, auc=auc2.grridge)
#   auc4[[5]][[r]] <- cbind(psel=psel1.grridge, auc=auc3.grridge)
#   auc4[[6]][[r]] <- cbind(psel=psel2.grridge, auc=auc4.grridge)
#   auc4[[7]][[r]] <- cbind(psel=psel1.enet, auc=auc1.enet)
#   auc4[[8]][[r]] <- cbind(psel=psel2.enet, auc=auc2.enet)
#   auc4[[9]][[r]] <- cbind(psel=psel3.enet, auc=auc3.enet)
#   auc4[[10]][[r]] <- cbind(psel=psel4.enet, auc=auc4.enet)
#   auc4[[11]][[r]] <- cbind(psel=psel1.greben, auc=auc1.greben)
#   auc4[[12]][[r]] <- cbind(psel=psel2.greben, auc=auc2.greben)
#   auc4[[13]][[r]] <- cbind(psel=psel3.greben, auc=auc3.greben)
#   auc4[[14]][[r]] <- cbind(psel=psel4.greben, auc=auc4.greben)
#   
#   briers4[[1]][[r]] <- cbind(psel=p, briers=briers1.ridge)
#   briers4[[2]][[r]] <- cbind(psel=p, briers=briers2.ridge)
#   briers4[[3]][[r]] <- cbind(psel=p, briers=briers1.grridge)
#   briers4[[4]][[r]] <- cbind(psel=p, briers=briers2.grridge)
#   briers4[[5]][[r]] <- cbind(psel=psel1.grridge, briers=briers3.grridge)
#   briers4[[6]][[r]] <- cbind(psel=psel2.grridge, briers=briers4.grridge)
#   briers4[[7]][[r]] <- cbind(psel=psel1.enet, briers=briers1.enet)
#   briers4[[8]][[r]] <- cbind(psel=psel2.enet, briers=briers2.enet)
#   briers4[[9]][[r]] <- cbind(psel=psel3.enet, briers=briers3.enet)
#   briers4[[10]][[r]] <- cbind(psel=psel4.enet, briers=briers4.enet)
#   briers4[[11]][[r]] <- cbind(psel=psel1.greben, briers=briers1.greben)
#   briers4[[12]][[r]] <- cbind(psel=psel2.greben, briers=briers2.greben)
#   briers4[[13]][[r]] <- cbind(psel=psel3.greben, briers=briers3.greben)
#   briers4[[14]][[r]] <- cbind(psel=psel4.greben, briers=briers4.greben)
#   
#   mse4[[1]][[r]] <- cbind(psel=p, mse=mse1.ridge)
#   mse4[[2]][[r]] <- cbind(psel=p, mse=mse2.ridge)
#   mse4[[3]][[r]] <- cbind(psel=p, mse=mse1.grridge)
#   mse4[[4]][[r]] <- cbind(psel=p, mse=mse2.grridge)
#   mse4[[5]][[r]] <- cbind(psel=psel1.grridge, mse=mse3.grridge)
#   mse4[[6]][[r]] <- cbind(psel=psel2.grridge, mse=mse4.grridge)
#   mse4[[7]][[r]] <- cbind(psel=psel1.enet, mse=mse1.enet)
#   mse4[[8]][[r]] <- cbind(psel=psel2.enet, mse=mse2.enet)
#   mse4[[9]][[r]] <- cbind(psel=psel3.enet, mse=mse3.enet)
#   mse4[[10]][[r]] <- cbind(psel=psel4.enet, mse=mse4.enet)
#   mse4[[11]][[r]] <- cbind(psel=psel1.greben, mse=mse1.greben)
#   mse4[[12]][[r]] <- cbind(psel=psel2.greben, mse=mse2.greben)
#   mse4[[13]][[r]] <- cbind(psel=psel3.greben, mse=mse3.greben)
#   mse4[[14]][[r]] <- cbind(psel=psel4.greben, mse=mse4.greben)
#   
#   
#   lambdag4[[1]][r, ] <- fit1.grridge[[1]]$lambdamults$groups
#   lambdag4[[2]][r, ] <- fit2.grridge[[1]]$lambdamults$groups
#   lambdag4[[3]][r, ] <- fit1.greben$lambdag$groups[, fit1.greben$nouteriter + 1]
#   lambdag4[[4]][r, ] <- fit2.greben$lambdag$groups[, fit2.greben$nouteriter + 1]
#   lambdag4[[5]][r, ] <- fit3.greben$lambdag$groups[, fit3.greben$nouteriter + 1]
#   lambdag4[[6]][r, ] <- fit4.greben$lambdag$groups[, fit4.greben$nouteriter + 1]
#   
#   varbeta4[r, ] <- sapply(1:G, function(g) {var(beta[(p*(g - 1)/G + 1):(p*g/G)])})
#   
#   results4 <- list(auc=auc4, briers=briers4, mse=mse4, lambdag=lambdag4,
#                    varbeta=varbeta4)
#   save(results4, file=paste(path.res, "grEBEN_test_res4.Rdata", sep=""))
#   
# }

# ### plots
# load(paste(path.res, "grEBEN_test_res4.Rdata", sep=""))
# 
# # preparing data
# auc1.enet <- do.call(rbind, results4$auc[[7]])[
#   order(do.call(rbind, results4$auc[[7]])[, 1]), ]
# auc2.enet <- do.call(rbind, results4$auc[[8]])[
#   order(do.call(rbind, results4$auc[[8]])[, 1]), ]
# auc3.enet <- do.call(rbind, results4$auc[[9]])[
#   order(do.call(rbind, results4$auc[[9]])[, 1]), ]
# auc4.enet <- do.call(rbind, results4$auc[[10]])[
#   order(do.call(rbind, results4$auc[[10]])[, 1]), ]
# 
# auc1.greben <- do.call(rbind, results4$auc[[11]])[
#   order(do.call(rbind, results4$auc[[11]])[, 1]), ]
# auc2.greben <- do.call(rbind, results4$auc[[12]])[
#   order(do.call(rbind, results4$auc[[12]])[, 1]), ]
# auc3.greben <- do.call(rbind, results4$auc[[13]])[
#   order(do.call(rbind, results4$auc[[13]])[, 1]), ]
# auc4.greben <- do.call(rbind, results4$auc[[14]])[
#   order(do.call(rbind, results4$auc[[14]])[, 1]), ]
# 
# auc1.grridge <- do.call(rbind, results4$auc[[3]])[
#   order(do.call(rbind, results4$auc[[3]])[, 1]), ]
# auc2.grridge <- do.call(rbind, results4$auc[[4]])[
#   order(do.call(rbind, results4$auc[[4]])[, 1]), ]
# auc3.grridge <- do.call(rbind, results4$auc[[5]])[
#   order(do.call(rbind, results4$auc[[5]])[, 1]), ]
# auc4.grridge <- do.call(rbind, results4$auc[[6]])[
#   order(do.call(rbind, results4$auc[[6]])[, 1]), ]
# 
# auc1.ridge <- do.call(rbind, results4$auc[[1]])[
#   order(do.call(rbind, results4$auc[[1]])[, 1]), ]
# auc2.ridge <- do.call(rbind, results4$auc[[2]])[
#   order(do.call(rbind, results4$auc[[2]])[, 1]), ]
# 
# lauc1.enet <- lowess(auc1.enet[, 1], auc1.enet[, 2])
# lauc2.enet <- lowess(auc2.enet[, 1], auc2.enet[, 2])
# lauc3.enet <- lowess(auc3.enet[, 1], auc3.enet[, 2])
# lauc4.enet <- lowess(auc4.enet[, 1], auc4.enet[, 2])
# 
# lauc1.greben <- lowess(auc1.greben[, 1], auc1.greben[, 2])
# lauc2.greben <- lowess(auc2.greben[, 1], auc2.greben[, 2])
# lauc3.greben <- lowess(auc3.greben[, 1], auc3.greben[, 2])
# lauc4.greben <- lowess(auc4.greben[, 1], auc4.greben[, 2])
# 
# lauc1.grridge <- lowess(auc1.grridge[, 1], auc1.grridge[, 2])
# lauc2.grridge <- lowess(auc2.grridge[, 1], auc2.grridge[, 2])
# lauc3.grridge <- lowess(auc3.grridge[, 1], auc3.grridge[, 2])
# lauc4.grridge <- lowess(auc4.grridge[, 1], auc4.grridge[, 2])
# 
# lauc1.ridge <- lowess(auc1.ridge[, 1], auc1.ridge[, 2])
# lauc2.ridge <- lowess(auc2.ridge[, 1], auc2.ridge[, 2])
# 
# briers1.enet <- do.call(rbind, results4$briers[[7]])[
#   order(do.call(rbind, results4$briers[[7]])[, 1]), ]
# briers2.enet <- do.call(rbind, results4$briers[[8]])[
#   order(do.call(rbind, results4$briers[[8]])[, 1]), ]
# briers3.enet <- do.call(rbind, results4$briers[[9]])[
#   order(do.call(rbind, results4$briers[[9]])[, 1]), ]
# briers4.enet <- do.call(rbind, results4$briers[[10]])[
#   order(do.call(rbind, results4$briers[[10]])[, 1]), ]
# 
# briers1.greben <- do.call(rbind, results4$briers[[11]])[
#   order(do.call(rbind, results4$briers[[11]])[, 1]), ]
# briers2.greben <- do.call(rbind, results4$briers[[12]])[
#   order(do.call(rbind, results4$briers[[12]])[, 1]), ]
# briers3.greben <- do.call(rbind, results4$briers[[13]])[
#   order(do.call(rbind, results4$briers[[13]])[, 1]), ]
# briers4.greben <- do.call(rbind, results4$briers[[14]])[
#   order(do.call(rbind, results4$briers[[14]])[, 1]), ]
# 
# briers1.grridge <- do.call(rbind, results4$briers[[3]])[
#   order(do.call(rbind, results4$briers[[3]])[, 1]), ]
# briers2.grridge <- do.call(rbind, results4$briers[[4]])[
#   order(do.call(rbind, results4$briers[[4]])[, 1]), ]
# briers3.grridge <- do.call(rbind, results4$briers[[5]])[
#   order(do.call(rbind, results4$briers[[5]])[, 1]), ]
# briers4.grridge <- do.call(rbind, results4$briers[[6]])[
#   order(do.call(rbind, results4$briers[[6]])[, 1]), ]
# 
# briers1.ridge <- do.call(rbind, results4$briers[[1]])[
#   order(do.call(rbind, results4$briers[[1]])[, 1]), ]
# briers2.ridge <- do.call(rbind, results4$briers[[2]])[
#   order(do.call(rbind, results4$briers[[2]])[, 1]), ]
# 
# lbriers1.enet <- lowess(briers1.enet[, 1], briers1.enet[, 2])
# lbriers2.enet <- lowess(briers2.enet[, 1], briers2.enet[, 2])
# lbriers3.enet <- lowess(briers3.enet[, 1], briers3.enet[, 2])
# lbriers4.enet <- lowess(briers4.enet[, 1], briers4.enet[, 2])
# 
# lbriers1.greben <- lowess(briers1.greben[, 1], briers1.greben[, 2])
# lbriers2.greben <- lowess(briers2.greben[, 1], briers2.greben[, 2])
# lbriers3.greben <- lowess(briers3.greben[, 1], briers3.greben[, 2])
# lbriers4.greben <- lowess(briers4.greben[, 1], briers4.greben[, 2])
# 
# lbriers1.grridge <- lowess(briers1.grridge[, 1], briers1.grridge[, 2])
# lbriers2.grridge <- lowess(briers2.grridge[, 1], briers2.grridge[, 2])
# lbriers3.grridge <- lowess(briers3.grridge[, 1], briers3.grridge[, 2])
# lbriers4.grridge <- lowess(briers4.grridge[, 1], briers4.grridge[, 2])
# 
# lbriers1.ridge <- lowess(briers1.ridge[, 1], briers1.ridge[, 2])
# lbriers2.ridge <- lowess(briers2.ridge[, 1], briers2.ridge[, 2])
# 
# mse1.enet <- do.call(rbind, results4$mse[[7]])[
#   order(do.call(rbind, results4$mse[[7]])[, 1]), ]
# mse2.enet <- do.call(rbind, results4$mse[[8]])[
#   order(do.call(rbind, results4$mse[[8]])[, 1]), ]
# mse3.enet <- do.call(rbind, results4$mse[[9]])[
#   order(do.call(rbind, results4$mse[[9]])[, 1]), ]
# mse4.enet <- do.call(rbind, results4$mse[[10]])[
#   order(do.call(rbind, results4$mse[[10]])[, 1]), ]
# 
# mse1.greben <- do.call(rbind, results4$mse[[11]])[
#   order(do.call(rbind, results4$mse[[11]])[, 1]), ]
# mse2.greben <- do.call(rbind, results4$mse[[12]])[
#   order(do.call(rbind, results4$mse[[12]])[, 1]), ]
# mse3.greben <- do.call(rbind, results4$mse[[13]])[
#   order(do.call(rbind, results4$mse[[13]])[, 1]), ]
# mse4.greben <- do.call(rbind, results4$mse[[14]])[
#   order(do.call(rbind, results4$mse[[14]])[, 1]), ]
# 
# mse1.grridge <- do.call(rbind, results4$mse[[3]])[
#   order(do.call(rbind, results4$mse[[3]])[, 1]), ]
# mse2.grridge <- do.call(rbind, results4$mse[[4]])[
#   order(do.call(rbind, results4$mse[[4]])[, 1]), ]
# mse3.grridge <- do.call(rbind, results4$mse[[5]])[
#   order(do.call(rbind, results4$mse[[5]])[, 1]), ]
# mse4.grridge <- do.call(rbind, results4$mse[[6]])[
#   order(do.call(rbind, results4$mse[[6]])[, 1]), ]
# 
# mse1.ridge <- do.call(rbind, results4$mse[[1]])[
#   order(do.call(rbind, results4$mse[[1]])[, 1]), ]
# mse2.ridge <- do.call(rbind, results4$mse[[2]])[
#   order(do.call(rbind, results4$mse[[2]])[, 1]), ]
# 
# lmse1.enet <- lowess(mse1.enet[, 1], mse1.enet[, 2])
# lmse2.enet <- lowess(mse2.enet[, 1], mse2.enet[, 2])
# lmse3.enet <- lowess(mse3.enet[, 1], mse3.enet[, 2])
# lmse4.enet <- lowess(mse4.enet[, 1], mse4.enet[, 2])
# 
# lmse1.greben <- lowess(mse1.greben[, 1], mse1.greben[, 2])
# lmse2.greben <- lowess(mse2.greben[, 1], mse2.greben[, 2])
# lmse3.greben <- lowess(mse3.greben[, 1], mse3.greben[, 2])
# lmse4.greben <- lowess(mse4.greben[, 1], mse4.greben[, 2])
# 
# lmse1.grridge <- lowess(mse1.grridge[, 1], mse1.grridge[, 2])
# lmse2.grridge <- lowess(mse2.grridge[, 1], mse2.grridge[, 2])
# lmse3.grridge <- lowess(mse3.grridge[, 1], mse3.grridge[, 2])
# lmse4.grridge <- lowess(mse4.grridge[, 1], mse4.grridge[, 2])
# 
# lmse1.ridge <- lowess(mse1.ridge[, 1], mse1.ridge[, 2])
# lmse2.ridge <- lowess(mse2.ridge[, 1], mse2.ridge[, 2])
# 
# ### Performance measures
# # AUC
# png(paste(path.graph, "grEBEN_test_res4_performance.png", sep=""),
#     units="in", width=12, height=4, res=120)
# par(mfrow=c(1, 3))
# xlim1.auc <- range(lauc1.enet$x, lauc2.enet$x, lauc3.enet$x, lauc4.enet$x,
#                    lauc1.greben$x, lauc2.greben$x, lauc3.greben$x,
#                    lauc4.greben$x, lauc3.grridge$x, lauc4.grridge$x)
# ylim1.auc <- range(lauc1.enet$y, lauc2.enet$y, lauc3.enet$y, lauc4.enet$y,
#                    lauc1.greben$y, lauc2.greben$y, lauc3.greben$y,
#                    lauc4.greben$y, lauc1.grridge$y, lauc2.grridge$y,
#                    lauc3.grridge$y, lauc4.grridge$y, lauc1.ridge$y,
#                    lauc2.ridge$y)
# 
# plot(0, 0, col=2, type="n", ylab="AUC", main="a)",
#      xlab="Number of selected variables", xlim=xlim1.auc, ylim=ylim1.auc)
# 
# lines(lauc3.grridge, col=2, lty=1)
# abline(h=lauc1.grridge$y, col=2, lty=2)
# abline(h=lauc1.ridge$y, col=2, lty=3)
# 
# lines(lauc4.grridge, col=3, lty=1)
# abline(h=lauc2.grridge$y, col=3, lty=2)
# abline(h=lauc2.ridge$y, col=3, lty=3)
# 
# lines(lauc1.greben, col=4, lty=1)
# lines(lauc1.enet, col=4, lty=3)
# 
# lines(lauc2.greben, col=5, lty=1)
# lines(lauc2.enet, col=5, lty=3)
# 
# lines(lauc3.greben, col=6, lty=1)
# lines(lauc3.enet, col=6, lty=3)
# 
# lines(lauc4.greben, col=7, lty=1)
# lines(lauc4.enet, col=7, lty=3)
# 
# # Brier skilll
# xlim1.briers <- range(lbriers1.enet$x, lbriers2.enet$x, lbriers3.enet$x,
#                       lbriers4.enet$x, lbriers1.greben$x, lbriers2.greben$x,
#                       lbriers3.greben$x, lbriers4.greben$x, lbriers3.grridge$x,
#                       lbriers4.grridge$x)
# ylim1.briers <- range(lbriers1.enet$y, lbriers2.enet$y, lbriers3.enet$y,
#                       lbriers4.enet$y, lbriers1.greben$y, lbriers2.greben$y,
#                       lbriers3.greben$y, lbriers4.greben$y, lbriers1.grridge$y,
#                       lbriers2.grridge$y, lbriers3.grridge$y, 
#                       lbriers4.grridge$y, lbriers1.ridge$y, lbriers2.ridge$y)
# 
# plot(0, 0, col=2, type="n", ylab="Brier skill score", main="b)",
#      xlab="Number of selected variables", xlim=xlim1.briers, ylim=ylim1.briers)
# 
# lines(lbriers3.grridge, col=2, lty=1)
# abline(h=lbriers1.grridge$y, col=2, lty=2)
# abline(h=lbriers1.ridge$y, col=2, lty=3)
# 
# lines(lbriers4.grridge, col=3, lty=1)
# abline(h=lbriers2.grridge$y, col=3, lty=2)
# abline(h=lbriers2.ridge$y, col=3, lty=3)
# 
# lines(lbriers1.greben, col=4, lty=1)
# lines(lbriers1.enet, col=4, lty=3)
# 
# lines(lbriers2.greben, col=5, lty=1)
# lines(lbriers2.enet, col=5, lty=3)
# 
# lines(lbriers3.greben, col=6, lty=1)
# lines(lbriers3.enet, col=6, lty=3)
# 
# lines(lbriers4.greben, col=7, lty=1)
# lines(lbriers4.enet, col=7, lty=3)
# 
# leglabels <- c(expression(paste("ridge, 'true' ", lambda)),
#                "ridge", expression(paste("enet, true ", lambda)),
#                expression(paste("enet, ", alpha==0.05)),
#                expression(paste("enet, ", alpha==0.5)),
#                expression(paste("enet, ", alpha==0.95)),
#                "group-regularized + selection", "group-regularized",
#                "not group-regularized")
# legend("bottomright", legend=leglabels, fill=c(2:7, 0, 0, 0),
#        lty=c(rep(NA, 6), 1, 2, 3), border=c(rep(1, 6), 0 ,0, 0), merge=TRUE,
#        seg.len=1)
# 
# # MSE
# xlim1.mse <- range(lmse1.enet$x, lmse2.enet$x, lmse3.enet$x,
#                    lmse4.enet$x, lmse1.greben$x, lmse2.greben$x,
#                    lmse3.greben$x, lmse4.greben$x, lmse3.grridge$x,
#                    lmse4.grridge$x)
# ylim1.mse <- range(lmse1.enet$y, lmse2.enet$y, lmse3.enet$y,
#                    lmse4.enet$y, lmse1.greben$y, lmse2.greben$y,
#                    lmse3.greben$y, lmse4.greben$y, lmse1.grridge$y,
#                    lmse2.grridge$y, lmse3.grridge$y, lmse4.grridge$y,
#                    lmse1.ridge$y, lmse2.ridge$y)
# 
# plot(0, 0, col=2, type="n", ylab="MSE", main="c)",
#      xlab="Number of selected variables", xlim=xlim1.mse, ylim=ylim1.mse)
# 
# lines(lmse3.grridge, col=2, lty=1)
# abline(h=lmse1.grridge$y, col=2, lty=2)
# abline(h=lmse1.ridge$y, col=2, lty=3)
# 
# lines(lmse4.grridge, col=3, lty=1)
# abline(h=lmse2.grridge$y, col=3, lty=2)
# abline(h=lmse2.ridge$y, col=3, lty=3)
# 
# lines(lmse1.greben, col=4, lty=1)
# lines(lmse1.enet, col=4, lty=3)
# 
# lines(lmse2.greben, col=5, lty=1)
# lines(lmse2.enet, col=5, lty=3)
# 
# lines(lmse3.greben, col=6, lty=1)
# lines(lmse3.enet, col=6, lty=3)
# 
# lines(lmse4.greben, col=7, lty=1)
# lines(lmse4.enet, col=7, lty=3)
# 
# dev.off()
# 
# ### diagnostics: checking lowess fits
# # AUC
# xlim2.auc <- range(auc1.enet[, 1], auc2.enet[, 1], auc3.enet[, 1], auc4.enet[, 1],
#                    auc1.greben[, 1], auc2.greben[, 1], auc3.greben[, 1],
#                    auc4.greben[, 1], auc3.grridge[, 1], auc4.grridge[, 1])
# ylim2.auc <- range(auc1.enet[, 2], auc2.enet[, 2], auc3.enet[, 2], auc4.enet[, 2],
#                    auc1.greben[, 2], auc2.greben[, 2], auc3.greben[, 2],
#                    auc4.greben[, 2], auc3.grridge[, 2], auc4.grridge[, 2])
# 
# png(paste(path.graph, "grEBEN_test_res4_lowess_auc.png", sep=""),
#     units="in", width=12, height=6, res=120)
# par(mfrow=c(2, 5))
# plot(auc3.grridge[, 1], auc3.grridge[, 2], col=1, ylab="AUC", main="a)",
#      xlab="Number of selected variables", xlim=xlim2.auc, ylim=ylim2.auc)
# lines(lauc3.grridge, col=2)
# 
# plot(auc1.greben[, 1], auc1.greben[, 2], col=1, ylab="AUC", main="b)",
#      xlab="Number of selected variables", xlim=xlim2.auc, ylim=ylim2.auc)
# lines(lauc1.greben, col=4)
# 
# plot(auc2.greben[, 1], auc2.greben[, 2], col=1, ylab="AUC", main="c)",
#      xlab="Number of selected variables", xlim=xlim2.auc, ylim=ylim2.auc)
# lines(lauc2.greben, col=5)
# 
# plot(auc3.greben[, 1], auc4.greben[, 2], col=1, ylab="AUC", main="d)",
#      xlab="Number of selected variables", xlim=xlim2.auc, ylim=ylim2.auc)
# lines(lauc3.greben, col=6)
# 
# plot(auc4.greben[, 1], auc4.greben[, 2], col=1, ylab="AUC", main="e)",
#      xlab="Number of selected variables", xlim=xlim2.auc, ylim=ylim2.auc)
# lines(lauc4.greben, col=7)
# 
# plot(auc4.grridge[, 1], auc4.grridge[, 2], col=1, ylab="AUC", main="f)",
#      xlab="Number of selected variables", xlim=xlim2.auc, ylim=ylim2.auc)
# lines(lauc4.grridge, col=3)
# 
# plot(auc1.enet[, 1], auc1.enet[, 2], col=1, ylab="AUC", main="g)",
#      xlab="Number of selected variables", xlim=xlim2.auc, ylim=ylim2.auc)
# lines(lauc1.enet, col=4)
# 
# plot(auc2.enet[, 1], auc2.enet[, 2], col=1, ylab="AUC", main="h)",
#      xlab="Number of selected variables", xlim=xlim2.auc, ylim=ylim2.auc)
# lines(lauc2.enet, col=5)
# 
# plot(auc3.enet[, 1], auc4.enet[, 2], col=1, ylab="AUC", main="i)",
#      xlab="Number of selected variables", xlim=xlim2.auc, ylim=ylim2.auc)
# lines(lauc3.enet, col=6)
# 
# plot(auc4.enet[, 1], auc4.enet[, 2], col=1, ylab="AUC", main="j)",
#      xlab="Number of selected variables", xlim=xlim2.auc, ylim=ylim2.auc)
# lines(lauc4.enet, col=7)
# dev.off()
# 
# # briers
# xlim2.briers <- range(briers1.enet[, 1], briers2.enet[, 1], briers3.enet[, 1], briers4.enet[, 1],
#                       briers1.greben[, 1], briers2.greben[, 1], briers3.greben[, 1],
#                       briers4.greben[, 1], briers3.grridge[, 1], briers4.grridge[, 1])
# ylim2.briers <- range(briers1.enet[, 2], briers2.enet[, 2], briers3.enet[, 2], briers4.enet[, 2],
#                       briers1.greben[, 2], briers2.greben[, 2], briers3.greben[, 2],
#                       briers4.greben[, 2], briers3.grridge[, 2], briers4.grridge[, 2])
# 
# png(paste(path.graph, "grEBEN_test_res4_lowess_briers.png", sep=""),
#     units="in", width=12, height=6, res=120)
# par(mfrow=c(2, 5))
# plot(briers3.grridge[, 1], briers3.grridge[, 2], col=1, ylab="Brier skill score", main="a)",
#      xlab="Number of selected variables", xlim=xlim2.briers, ylim=ylim2.briers)
# lines(lbriers3.grridge, col=2)
# 
# plot(briers1.greben[, 1], briers1.greben[, 2], col=1, ylab="Brier skill score", main="b)",
#      xlab="Number of selected variables", xlim=xlim2.briers, ylim=ylim2.briers)
# lines(lbriers1.greben, col=4)
# 
# plot(briers2.greben[, 1], briers2.greben[, 2], col=1, ylab="Brier skill score", main="c)",
#      xlab="Number of selected variables", xlim=xlim2.briers, ylim=ylim2.briers)
# lines(lbriers2.greben, col=5)
# 
# plot(briers3.greben[, 1], briers4.greben[, 2], col=1, ylab="Brier skill score", main="d)",
#      xlab="Number of selected variables", xlim=xlim2.briers, ylim=ylim2.briers)
# lines(lbriers3.greben, col=6)
# 
# plot(briers4.greben[, 1], briers4.greben[, 2], col=1, ylab="Brier skill score", main="e)",
#      xlab="Number of selected variables", xlim=xlim2.briers, ylim=ylim2.briers)
# lines(lbriers4.greben, col=7)
# 
# plot(briers4.grridge[, 1], briers4.grridge[, 2], col=1, ylab="Brier skill score", main="f)",
#      xlab="Number of selected variables", xlim=xlim2.briers, ylim=ylim2.briers)
# lines(lbriers4.grridge, col=3)
# 
# plot(briers1.enet[, 1], briers1.enet[, 2], col=1, ylab="Brier skill score", main="g)",
#      xlab="Number of selected variables", xlim=xlim2.briers, ylim=ylim2.briers)
# lines(lbriers1.enet, col=4)
# 
# plot(briers2.enet[, 1], briers2.enet[, 2], col=1, ylab="Brier skill score", main="h)",
#      xlab="Number of selected variables", xlim=xlim2.briers, ylim=ylim2.briers)
# lines(lbriers2.enet, col=5)
# 
# plot(briers3.enet[, 1], briers4.enet[, 2], col=1, ylab="Brier skill score", main="i)",
#      xlab="Number of selected variables", xlim=xlim2.briers, ylim=ylim2.briers)
# lines(lbriers3.enet, col=6)
# 
# plot(briers4.enet[, 1], briers4.enet[, 2], col=1, ylab="Brier skill score", main="j)",
#      xlab="Number of selected variables", xlim=xlim2.briers, ylim=ylim2.briers)
# lines(lbriers4.enet, col=7)
# dev.off()
# 
# # mse
# xlim2.mse <- range(mse1.enet[, 1], mse2.enet[, 1], mse3.enet[, 1], mse4.enet[, 1],
#                    mse1.greben[, 1], mse2.greben[, 1], mse3.greben[, 1],
#                    mse4.greben[, 1], mse3.grridge[, 1], mse4.grridge[, 1])
# ylim2.mse <- range(mse1.enet[, 2], mse2.enet[, 2], mse3.enet[, 2], mse4.enet[, 2],
#                    mse1.greben[, 2], mse2.greben[, 2], mse3.greben[, 2],
#                    mse4.greben[, 2], mse3.grridge[, 2], mse4.grridge[, 2])
# 
# png(paste(path.graph, "grEBEN_test_res4_lowess_mse.png", sep=""),
#     units="in", width=12, height=6, res=120)
# par(mfrow=c(2, 5))
# plot(mse3.grridge[, 1], mse3.grridge[, 2], col=1, ylab="MSE", main="a)",
#      xlab="Number of selected variables", xlim=xlim2.mse, ylim=ylim2.mse)
# lines(lmse3.grridge, col=2)
# 
# plot(mse1.greben[, 1], mse1.greben[, 2], col=1, ylab="MSE", main="b)",
#      xlab="Number of selected variables", xlim=xlim2.mse, ylim=ylim2.mse)
# lines(lmse1.greben, col=4)
# 
# plot(mse2.greben[, 1], mse2.greben[, 2], col=1, ylab="MSE", main="c)",
#      xlab="Number of selected variables", xlim=xlim2.mse, ylim=ylim2.mse)
# lines(lmse2.greben, col=5)
# 
# plot(mse3.greben[, 1], mse4.greben[, 2], col=1, ylab="MSE", main="d)",
#      xlab="Number of selected variables", xlim=xlim2.mse, ylim=ylim2.mse)
# lines(lmse3.greben, col=6)
# 
# plot(mse4.greben[, 1], mse4.greben[, 2], col=1, ylab="MSE", main="e)",
#      xlab="Number of selected variables", xlim=xlim2.mse, ylim=ylim2.mse)
# lines(lmse4.greben, col=7)
# 
# plot(mse4.grridge[, 1], mse4.grridge[, 2], col=1, ylab="MSE", main="f)",
#      xlab="Number of selected variables", xlim=xlim2.mse, ylim=ylim2.mse)
# lines(lmse4.grridge, col=3)
# 
# plot(mse1.enet[, 1], mse1.enet[, 2], col=1, ylab="MSE", main="g)",
#      xlab="Number of selected variables", xlim=xlim2.mse, ylim=ylim2.mse)
# lines(lmse1.enet, col=4)
# 
# plot(mse2.enet[, 1], mse2.enet[, 2], col=1, ylab="MSE", main="h)",
#      xlab="Number of selected variables", xlim=xlim2.mse, ylim=ylim2.mse)
# lines(lmse2.enet, col=5)
# 
# plot(mse3.enet[, 1], mse4.enet[, 2], col=1, ylab="MSE", main="i)",
#      xlab="Number of selected variables", xlim=xlim2.mse, ylim=ylim2.mse)
# lines(lmse3.enet, col=6)
# 
# plot(mse4.enet[, 1], mse4.enet[, 2], col=1, ylab="MSE", main="j)",
#      xlab="Number of selected variables", xlim=xlim2.mse, ylim=ylim2.mse)
# lines(lmse4.enet, col=7)
# dev.off()

## penalty parameters
lambdag <- exp(seq(-1, 1, length.out=ncol(results4$lambdag[[1]])))
png(paste(path.graph, "grEBEN_test_res4_penalties.png", sep=""),
    units="in", width=8, height=6, res=120)
par(mfrow=c(2, 3))
boxplot(results4$lambdag[[1]], ylim=range(results4$lambdag[[1]], lambdag), 
        xlab="", ylab="", main="a)")
title(ylab=expression(paste(lambda[g], "'")), 
      xlab="Groups in decreasing effect size", line=2.5)
points(1:length(lambdag), lambdag, col=2, pch="x")
boxplot(results4$lambdag[[2]], ylim=range(results4$lambdag[[2]], lambdag), 
        xlab="", ylab="", main="b)")
title(ylab=expression(paste(lambda[g], "'")), 
      xlab="Groups in decreasing effect size", line=2.5)
points(1:length(lambdag), lambdag, col=2, pch="x")
boxplot(results4$lambdag[[3]], ylim=range(results4$lambdag[[3]], lambdag), 
        xlab="", ylab="", main="c)")
title(ylab=expression(paste(lambda[g], "'")), 
      xlab="Groups in decreasing effect size", line=2.5)
points(1:length(lambdag), lambdag, col=2, pch="x")
boxplot(results4$lambdag[[4]], ylim=range(results4$lambdag[[4]], lambdag), 
        xlab="", ylab="", main="d)")
title(ylab=expression(paste(lambda[g], "'")), 
      xlab="Groups in decreasing effect size", line=2.5)
points(1:length(lambdag), lambdag, col=2, pch="x")
boxplot(results4$lambdag[[5]], ylim=range(results4$lambdag[[5]], lambdag), 
        xlab="", ylab="", main="e)")
title(ylab=expression(paste(lambda[g], "'")), 
      xlab="Groups in decreasing effect size", line=2.5)
points(1:length(lambdag), lambdag, col=2, pch="x")
boxplot(results4$lambdag[[6]], ylim=range(results4$lambdag[[6]], lambdag), 
        xlab="", ylab="", main="f)")
title(ylab=expression(paste(lambda[g], "'")), 
      xlab="Groups in decreasing effect size", line=2.5)
points(1:length(lambdag), lambdag, col=2, pch="x")
dev.off()





# ## simulation 5
# # create data
# n <- 200
# p <- 1000
# G <- 5
# pblock <- 20
# rho <- 0.7
# sigma <- matrix(rho, ncol=pblock, nrow=pblock); diag(sigma) <- 1
# lambda <- 0.02
# alpha <- 0.95
# lambdag <- exp(seq(-1, 1, length.out=G))
# m <- rep(1, n)
# part.greben <- list(groups=rep(1:G, each=p/G))
# part.grridge <- list(groups=CreatePartition(as.factor(part.greben$groups)))
# ntest <- 1000
# 
# methods <- c("ridge+truelambda", "ridge", "GRridge+truelambda", "GRridge",
#              "GRridge+truelambda+sel", "GRridge+sel",
#              "enet+truelambda", "enet+a=0.05", "enet+a=0.5", "enet+a=0.95",
#              "grEBEN+truelambda", "grEBEN+a=0.05", "grEBEN+a=0.5","grEBEN+a=0.95")
# nreps <- 50
# auc5 <- briers5 <- mse5 <- rep(list(vector(mode="list", length=nreps)),
#                                length(methods))
# names(auc5) <- names(briers5) <- names(mse5) <- methods
# lambdag5 <- rep(list(matrix(NA, ncol=G, nrow=nreps)), 6)
# names(lambdag5) <- methods[c(3, 4, 11, 12, 13, 14)]
# varbeta5 <- matrix(nrow=nreps, ncol=G)
# # the simulations
# for(r in 1:nreps) {
#   
#   print(paste("Simulation 3, repetition", r))
#   set.seed(600 + r)
#   x <- do.call(cbind, replicate(p/pblock, rmvnorm(n, mean=rep(0, pblock), sigma=sigma),
#                                 simplify=FALSE))
#   beta <- as.numeric(sapply(1:G, function(g) {
#     renbeta(p/G, 2*n*lambda*alpha*sqrt(lambdag[g]), n*lambda*(1 - alpha)*lambdag[g])}))
#   prob <- as.numeric(exp(x %*% beta)/(1 + exp(x %*% beta)))
#   y <- rbinom(n, 1, prob)
#   optl <- n*lambda*alpha*sum(abs(beta))/sum(beta^2) + 0.5*n*lambda*(1 - alpha)
#   
#   xtest <- do.call(cbind, replicate(p/pblock, rmvnorm(ntest, mean=rep(0, pblock),
#                                                       sigma=sigma), simplify=FALSE))
#   probtest <- as.numeric(exp(xtest %*% beta)/(1 + exp(xtest %*% beta)))
#   ytest <- rbinom(ntest, 1, probtest)
#   
#   fit1.greben <- grEBEN3(x, y, m, partitions=part.greben, alpha=alpha - 0.01,
#                          lambda=lambda, psel=TRUE)
#   fit2.greben <- grEBEN3(x, y, m, partitions=part.greben, alpha=0.05, psel=TRUE)
#   fit3.greben <- grEBEN3(x, y, m, partitions=part.greben, alpha=0.5, psel=TRUE)
#   fit4.greben <- grEBEN3(x, y, m, partitions=part.greben, alpha=0.95, psel=TRUE)
#   
#   # number of selected variables
#   psel1.enet <- apply(fit1.greben$beta.nogroups, 2, function(b) {sum(b!=0) - 1})
#   psel2.enet <- apply(fit2.greben$beta.nogroups, 2, function(b) {sum(b!=0) - 1})
#   psel3.enet <- apply(fit3.greben$beta.nogroups, 2, function(b) {sum(b!=0) - 1})
#   psel4.enet <- apply(fit4.greben$beta.nogroups, 2, function(b) {sum(b!=0) - 1})
#   
#   psel1.greben <- apply(fit1.greben$beta, 2, function(b) {sum(b!=0) - 1})
#   psel2.greben <- apply(fit2.greben$beta, 2, function(b) {sum(b!=0) - 1})
#   psel3.greben <- apply(fit3.greben$beta, 2, function(b) {sum(b!=0) - 1})
#   psel4.greben <- apply(fit4.greben$beta, 2, function(b) {sum(b!=0) - 1})
#   
#   psel.all <- unique(c(psel1.enet, psel2.enet, psel3.enet, psel4.enet,
#                        psel1.greben, psel2.greben, psel3.greben, psel4.greben))
#   psel.in <- floor(quantile(psel.all[psel.all!=0],
#                             prob=c(seq(0.01, 0.05, 0.01), seq(0.07, 0.25, 0.03), seq(0.3, 1, 0.1))))
#   
#   fit1.grridge <- vector(mode="list", length=length(psel.in))
#   for(s in 1:length(psel.in)) {
#     fit1.grridge[[s]] <- grridge(t(x), y, part.grridge, selectionEN=TRUE,
#                                  optl=optl, maxsel=psel.in[s])
#   }
#   
#   fit2.grridge <- vector(mode="list", length=length(psel.in))
#   fit2.grridge[[1]] <- grridge(t(x), y, part.grridge, selectionEN=TRUE,
#                                maxsel=psel.in[1])
#   for(s in 2:length(psel.in)) {
#     fit2.grridge[[s]] <- grridge(t(x), y, part.grridge, selectionEN=TRUE,
#                                  optl=fit2.grridge[[1]]$optl,
#                                  maxsel=psel.in[s])
#   }
#   
#   psel1.grridge <- sapply(1:length(fit1.grridge), function(s) {
#     return(length(fit1.grridge[[s]]$resEN$whichEN))})
#   psel2.grridge <- sapply(1:length(fit2.grridge), function(s) {
#     return(length(fit2.grridge[[s]]$resEN$whichEN))})
#   
#   # estimates
#   est1.ridge <- coef(fit1.grridge[[1]]$predobj$NoGroups, "all")
#   est2.ridge <- coef(fit2.grridge[[1]]$predobj$NoGroups, "all")
#   
#   est1.grridge <- coef(fit1.grridge[[1]]$predobj$GroupRegul, "all")
#   est2.grridge <- coef(fit2.grridge[[1]]$predobj$GroupRegul, "all")
#   
#   est3.grridge <- sapply(fit1.grridge, function(s) {
#     replace(rep(0, p + 1), c(1, s$resEN$whichEN + 1), coef(s$predobj$EN, "all"))})
#   est4.grridge <- sapply(fit2.grridge, function(s) {
#     replace(rep(0, p + 1), c(1, s$resEN$whichEN + 1), coef(s$predobj$EN, "all"))})
#   
#   est1.enet <- fit1.greben$beta.nogroups
#   est2.enet <- fit2.greben$beta.nogroups
#   est3.enet <- fit3.greben$beta.nogroups
#   est4.enet <- fit4.greben$beta.nogroups
#   
#   est1.greben <- fit1.greben$beta
#   est2.greben <- fit2.greben$beta
#   est3.greben <- fit3.greben$beta
#   est4.greben <- fit4.greben$beta
#   
#   # predictions on fit data
#   pred1.ridge <- predict.grridge(fit1.grridge[[1]], t(xtest))[, 1]
#   pred2.ridge <- predict.grridge(fit2.grridge[[1]], t(xtest))[, 1]
#   
#   pred1.grridge <- predict.grridge(fit1.grridge[[1]], t(xtest))[, 2]
#   pred2.grridge <- predict.grridge(fit2.grridge[[1]], t(xtest))[, 2]
#   pred3.grridge <- sapply(fit1.grridge, function(s) {
#     predict.grridge(s, t(xtest))[, 3]})
#   pred4.grridge <- sapply(fit2.grridge, function(s) {
#     predict.grridge(s, t(xtest))[, 3]})
#   
#   pred1.enet <- 1/(1 + exp(-xtest %*% est1.enet[-1, ]))
#   pred2.enet <- 1/(1 + exp(-xtest %*% est2.enet[-1, ]))
#   pred3.enet <- 1/(1 + exp(-xtest %*% est3.enet[-1, ]))
#   pred4.enet <- 1/(1 + exp(-xtest %*% est4.enet[-1, ]))
#   
#   pred1.greben <- 1/(1 + exp(-xtest %*% est1.greben[-1, ]))
#   pred2.greben <- 1/(1 + exp(-xtest %*% est2.greben[-1, ]))
#   pred3.greben <- 1/(1 + exp(-xtest %*% est3.greben[-1, ]))
#   pred4.greben <- 1/(1 + exp(-xtest %*% est4.greben[-1, ]))
#   
#   # AUCs
#   auc.true <- pROC::roc(ytest, probtest)$auc
#   
#   auc1.ridge <- pROC::roc(ytest, pred1.ridge)$auc
#   auc2.ridge <- pROC::roc(ytest, pred2.ridge)$auc
#   
#   auc1.grridge <- pROC::roc(ytest, pred1.grridge)$auc
#   auc2.grridge <- pROC::roc(ytest, pred2.grridge)$auc
#   auc3.grridge <- apply(pred3.grridge, 2, function(r) {pROC::roc(ytest, r)$auc})
#   auc4.grridge <- apply(pred4.grridge, 2, function(r) {pROC::roc(ytest, r)$auc})
#   
#   auc1.enet <- apply(pred1.enet, 2, function(r) {pROC::roc(ytest, r)$auc})
#   auc2.enet <- apply(pred2.enet, 2, function(r) {pROC::roc(ytest, r)$auc})
#   auc3.enet <- apply(pred3.enet, 2, function(r) {pROC::roc(ytest, r)$auc})
#   auc4.enet <- apply(pred4.enet, 2, function(r) {pROC::roc(ytest, r)$auc})
#   
#   auc1.greben <- apply(pred1.greben, 2, function(r) {pROC::roc(ytest, r)$auc})
#   auc2.greben <- apply(pred2.greben, 2, function(r) {pROC::roc(ytest, r)$auc})
#   auc3.greben <- apply(pred3.greben, 2, function(r) {pROC::roc(ytest, r)$auc})
#   auc4.greben <- apply(pred4.greben, 2, function(r) {pROC::roc(ytest, r)$auc})
#   
#   # Brier scores
#   brier.null <- sum((ytest - mean(ytest))^2)
#   briers.true <- 1 - sum((ytest - probtest)^2)/brier.null
#   
#   briers1.ridge <- 1 - sum((ytest - pred1.ridge)^2)/brier.null
#   briers2.ridge <- 1 - sum((ytest - pred2.ridge)^2)/brier.null
#   
#   briers1.grridge <- 1 - sum((ytest - pred1.grridge)^2)/brier.null
#   briers2.grridge <- 1 - sum((ytest - pred2.grridge)^2)/brier.null
#   briers3.grridge <- apply(pred3.grridge, 2, function(pred) {
#     1 - sum((ytest - pred)^2)/brier.null})
#   briers4.grridge <- apply(pred4.grridge, 2, function(pred) {
#     1 - sum((ytest - pred)^2)/brier.null})
#   
#   briers1.enet <- apply(pred1.enet, 2, function(pred) {
#     1 - sum((ytest - pred)^2)/brier.null})
#   briers2.enet <- apply(pred2.enet, 2, function(pred) {
#     1 - sum((ytest - pred)^2)/brier.null})
#   briers3.enet <- apply(pred3.enet, 2, function(pred) {
#     1 - sum((ytest - pred)^2)/brier.null})
#   briers4.enet <- apply(pred4.enet, 2, function(pred) {
#     1 - sum((ytest - pred)^2)/brier.null})
#   
#   briers1.greben <- apply(pred1.greben, 2, function(pred) {
#     1 - sum((ytest - pred)^2)/brier.null})
#   briers2.greben <- apply(pred2.greben, 2, function(pred) {
#     1 - sum((ytest - pred)^2)/brier.null})
#   briers3.greben <- apply(pred3.greben, 2, function(pred) {
#     1 - sum((ytest - pred)^2)/brier.null})
#   briers4.greben <- apply(pred4.greben, 2, function(pred) {
#     1 - sum((ytest - pred)^2)/brier.null})
#   
#   # MSE
#   mse.true <- 0
#   
#   mse1.ridge <- mean((c(0, beta) - est1.ridge)^2)
#   mse2.ridge <- mean((c(0, beta) - est2.ridge)^2)
#   
#   mse1.grridge <- mean((c(0, beta) - est1.grridge)^2)
#   mse2.grridge <- mean((c(0, beta) - est2.grridge)^2)
#   mse3.grridge <- apply(est3.grridge, 2, function(b) {mean((c(0, beta) - b)^2)})
#   mse4.grridge <- apply(est4.grridge, 2, function(b) {mean((c(0, beta) - b)^2)})
#   
#   mse1.enet <- apply(est1.enet, 2, function(b) {mean((c(0, beta) - b)^2)})
#   mse2.enet <- apply(est2.enet, 2, function(b) {mean((c(0, beta) - b)^2)})
#   mse3.enet <- apply(est3.enet, 2, function(b) {mean((c(0, beta) - b)^2)})
#   mse4.enet <- apply(est4.enet, 2, function(b) {mean((c(0, beta) - b)^2)})
#   
#   mse1.greben <- apply(est1.greben, 2, function(b) {mean((c(0, beta) - b)^2)})
#   mse2.greben <- apply(est2.greben, 2, function(b) {mean((c(0, beta) - b)^2)})
#   mse3.greben <- apply(est3.greben, 2, function(b) {mean((c(0, beta) - b)^2)})
#   mse4.greben <- apply(est4.greben, 2, function(b) {mean((c(0, beta) - b)^2)})
#   
#   auc5[[1]][[r]] <- cbind(psel=p, auc=auc1.ridge)
#   auc5[[2]][[r]] <- cbind(psel=p, auc=auc2.ridge)
#   auc5[[3]][[r]] <- cbind(psel=p, auc=auc1.grridge)
#   auc5[[4]][[r]] <- cbind(psel=p, auc=auc2.grridge)
#   auc5[[5]][[r]] <- cbind(psel=psel1.grridge, auc=auc3.grridge)
#   auc5[[6]][[r]] <- cbind(psel=psel2.grridge, auc=auc4.grridge)
#   auc5[[7]][[r]] <- cbind(psel=psel1.enet, auc=auc1.enet)
#   auc5[[8]][[r]] <- cbind(psel=psel2.enet, auc=auc2.enet)
#   auc5[[9]][[r]] <- cbind(psel=psel3.enet, auc=auc3.enet)
#   auc5[[10]][[r]] <- cbind(psel=psel4.enet, auc=auc4.enet)
#   auc5[[11]][[r]] <- cbind(psel=psel1.greben, auc=auc1.greben)
#   auc5[[12]][[r]] <- cbind(psel=psel2.greben, auc=auc2.greben)
#   auc5[[13]][[r]] <- cbind(psel=psel3.greben, auc=auc3.greben)
#   auc5[[14]][[r]] <- cbind(psel=psel4.greben, auc=auc4.greben)
#   
#   briers5[[1]][[r]] <- cbind(psel=p, briers=briers1.ridge)
#   briers5[[2]][[r]] <- cbind(psel=p, briers=briers2.ridge)
#   briers5[[3]][[r]] <- cbind(psel=p, briers=briers1.grridge)
#   briers5[[4]][[r]] <- cbind(psel=p, briers=briers2.grridge)
#   briers5[[5]][[r]] <- cbind(psel=psel1.grridge, briers=briers3.grridge)
#   briers5[[6]][[r]] <- cbind(psel=psel2.grridge, briers=briers4.grridge)
#   briers5[[7]][[r]] <- cbind(psel=psel1.enet, briers=briers1.enet)
#   briers5[[8]][[r]] <- cbind(psel=psel2.enet, briers=briers2.enet)
#   briers5[[9]][[r]] <- cbind(psel=psel3.enet, briers=briers3.enet)
#   briers5[[10]][[r]] <- cbind(psel=psel4.enet, briers=briers4.enet)
#   briers5[[11]][[r]] <- cbind(psel=psel1.greben, briers=briers1.greben)
#   briers5[[12]][[r]] <- cbind(psel=psel2.greben, briers=briers2.greben)
#   briers5[[13]][[r]] <- cbind(psel=psel3.greben, briers=briers3.greben)
#   briers5[[14]][[r]] <- cbind(psel=psel4.greben, briers=briers4.greben)
#   
#   mse5[[1]][[r]] <- cbind(psel=p, mse=mse1.ridge)
#   mse5[[2]][[r]] <- cbind(psel=p, mse=mse2.ridge)
#   mse5[[3]][[r]] <- cbind(psel=p, mse=mse1.grridge)
#   mse5[[4]][[r]] <- cbind(psel=p, mse=mse2.grridge)
#   mse5[[5]][[r]] <- cbind(psel=psel1.grridge, mse=mse3.grridge)
#   mse5[[6]][[r]] <- cbind(psel=psel2.grridge, mse=mse4.grridge)
#   mse5[[7]][[r]] <- cbind(psel=psel1.enet, mse=mse1.enet)
#   mse5[[8]][[r]] <- cbind(psel=psel2.enet, mse=mse2.enet)
#   mse5[[9]][[r]] <- cbind(psel=psel3.enet, mse=mse3.enet)
#   mse5[[10]][[r]] <- cbind(psel=psel4.enet, mse=mse4.enet)
#   mse5[[11]][[r]] <- cbind(psel=psel1.greben, mse=mse1.greben)
#   mse5[[12]][[r]] <- cbind(psel=psel2.greben, mse=mse2.greben)
#   mse5[[13]][[r]] <- cbind(psel=psel3.greben, mse=mse3.greben)
#   mse5[[14]][[r]] <- cbind(psel=psel4.greben, mse=mse4.greben)
#   
#   
#   lambdag5[[1]][r, ] <- fit1.grridge[[1]]$lambdamults$groups
#   lambdag5[[2]][r, ] <- fit2.grridge[[1]]$lambdamults$groups
#   lambdag5[[3]][r, ] <- fit1.greben$lambdag$groups[, fit1.greben$nouteriter + 1]
#   lambdag5[[4]][r, ] <- fit2.greben$lambdag$groups[, fit2.greben$nouteriter + 1]
#   lambdag5[[5]][r, ] <- fit3.greben$lambdag$groups[, fit3.greben$nouteriter + 1]
#   lambdag5[[6]][r, ] <- fit4.greben$lambdag$groups[, fit4.greben$nouteriter + 1]
#   
#   varbeta5[r, ] <- sapply(1:G, function(g) {var(beta[(p*(g - 1)/G + 1):(p*g/G)])})
#   
#   results5 <- list(auc=auc5, briers=briers5, mse=mse5, lambdag=lambdag5,
#                    varbeta=varbeta5)
#   save(results5, file=paste(path.res, "grEBEN_test_res5.Rdata", sep=""))
#   
# }

### plots
load(paste(path.res, "grEBEN_test_res5.Rdata", sep=""))

# preparing data
auc1.enet <- do.call(rbind, results5$auc[[7]])[
  order(do.call(rbind, results5$auc[[7]])[, 1]), ]
auc2.enet <- do.call(rbind, results5$auc[[8]])[
  order(do.call(rbind, results5$auc[[8]])[, 1]), ]
auc3.enet <- do.call(rbind, results5$auc[[9]])[
  order(do.call(rbind, results5$auc[[9]])[, 1]), ]
auc4.enet <- do.call(rbind, results5$auc[[10]])[
  order(do.call(rbind, results5$auc[[10]])[, 1]), ]

auc1.greben <- do.call(rbind, results5$auc[[11]])[
  order(do.call(rbind, results5$auc[[11]])[, 1]), ]
auc2.greben <- do.call(rbind, results5$auc[[12]])[
  order(do.call(rbind, results5$auc[[12]])[, 1]), ]
auc3.greben <- do.call(rbind, results5$auc[[13]])[
  order(do.call(rbind, results5$auc[[13]])[, 1]), ]
auc4.greben <- do.call(rbind, results5$auc[[14]])[
  order(do.call(rbind, results5$auc[[14]])[, 1]), ]

auc1.grridge <- do.call(rbind, results5$auc[[3]])[
  order(do.call(rbind, results5$auc[[3]])[, 1]), ]
auc2.grridge <- do.call(rbind, results5$auc[[4]])[
  order(do.call(rbind, results5$auc[[4]])[, 1]), ]
auc3.grridge <- do.call(rbind, results5$auc[[5]])[
  order(do.call(rbind, results5$auc[[5]])[, 1]), ]
auc4.grridge <- do.call(rbind, results5$auc[[6]])[
  order(do.call(rbind, results5$auc[[6]])[, 1]), ]

auc1.ridge <- do.call(rbind, results5$auc[[1]])[
  order(do.call(rbind, results5$auc[[1]])[, 1]), ]
auc2.ridge <- do.call(rbind, results5$auc[[2]])[
  order(do.call(rbind, results5$auc[[2]])[, 1]), ]

lauc1.enet <- lowess(auc1.enet[, 1], auc1.enet[, 2])
lauc2.enet <- lowess(auc2.enet[, 1], auc2.enet[, 2])
lauc3.enet <- lowess(auc3.enet[, 1], auc3.enet[, 2])
lauc4.enet <- lowess(auc4.enet[, 1], auc4.enet[, 2])

lauc1.greben <- lowess(auc1.greben[, 1], auc1.greben[, 2])
lauc2.greben <- lowess(auc2.greben[, 1], auc2.greben[, 2])
lauc3.greben <- lowess(auc3.greben[, 1], auc3.greben[, 2])
lauc4.greben <- lowess(auc4.greben[, 1], auc4.greben[, 2])

lauc1.grridge <- lowess(auc1.grridge[, 1], auc1.grridge[, 2])
lauc2.grridge <- lowess(auc2.grridge[, 1], auc2.grridge[, 2])
lauc3.grridge <- lowess(auc3.grridge[, 1], auc3.grridge[, 2])
lauc4.grridge <- lowess(auc4.grridge[, 1], auc4.grridge[, 2])

lauc1.ridge <- lowess(auc1.ridge[, 1], auc1.ridge[, 2])
lauc2.ridge <- lowess(auc2.ridge[, 1], auc2.ridge[, 2])

briers1.enet <- do.call(rbind, results5$briers[[7]])[
  order(do.call(rbind, results5$briers[[7]])[, 1]), ]
briers2.enet <- do.call(rbind, results5$briers[[8]])[
  order(do.call(rbind, results5$briers[[8]])[, 1]), ]
briers3.enet <- do.call(rbind, results5$briers[[9]])[
  order(do.call(rbind, results5$briers[[9]])[, 1]), ]
briers4.enet <- do.call(rbind, results5$briers[[10]])[
  order(do.call(rbind, results5$briers[[10]])[, 1]), ]

briers1.greben <- do.call(rbind, results5$briers[[11]])[
  order(do.call(rbind, results5$briers[[11]])[, 1]), ]
briers2.greben <- do.call(rbind, results5$briers[[12]])[
  order(do.call(rbind, results5$briers[[12]])[, 1]), ]
briers3.greben <- do.call(rbind, results5$briers[[13]])[
  order(do.call(rbind, results5$briers[[13]])[, 1]), ]
briers4.greben <- do.call(rbind, results5$briers[[14]])[
  order(do.call(rbind, results5$briers[[14]])[, 1]), ]

briers1.grridge <- do.call(rbind, results5$briers[[3]])[
  order(do.call(rbind, results5$briers[[3]])[, 1]), ]
briers2.grridge <- do.call(rbind, results5$briers[[4]])[
  order(do.call(rbind, results5$briers[[4]])[, 1]), ]
briers3.grridge <- do.call(rbind, results5$briers[[5]])[
  order(do.call(rbind, results5$briers[[5]])[, 1]), ]
briers4.grridge <- do.call(rbind, results5$briers[[6]])[
  order(do.call(rbind, results5$briers[[6]])[, 1]), ]

briers1.ridge <- do.call(rbind, results5$briers[[1]])[
  order(do.call(rbind, results5$briers[[1]])[, 1]), ]
briers2.ridge <- do.call(rbind, results5$briers[[2]])[
  order(do.call(rbind, results5$briers[[2]])[, 1]), ]

lbriers1.enet <- lowess(briers1.enet[, 1], briers1.enet[, 2])
lbriers2.enet <- lowess(briers2.enet[, 1], briers2.enet[, 2])
lbriers3.enet <- lowess(briers3.enet[, 1], briers3.enet[, 2])
lbriers4.enet <- lowess(briers4.enet[, 1], briers4.enet[, 2])

lbriers1.greben <- lowess(briers1.greben[, 1], briers1.greben[, 2])
lbriers2.greben <- lowess(briers2.greben[, 1], briers2.greben[, 2])
lbriers3.greben <- lowess(briers3.greben[, 1], briers3.greben[, 2])
lbriers4.greben <- lowess(briers4.greben[, 1], briers4.greben[, 2])

lbriers1.grridge <- lowess(briers1.grridge[, 1], briers1.grridge[, 2])
lbriers2.grridge <- lowess(briers2.grridge[, 1], briers2.grridge[, 2])
lbriers3.grridge <- lowess(briers3.grridge[, 1], briers3.grridge[, 2])
lbriers4.grridge <- lowess(briers4.grridge[, 1], briers4.grridge[, 2])

lbriers1.ridge <- lowess(briers1.ridge[, 1], briers1.ridge[, 2])
lbriers2.ridge <- lowess(briers2.ridge[, 1], briers2.ridge[, 2])

mse1.enet <- do.call(rbind, results5$mse[[7]])[
  order(do.call(rbind, results5$mse[[7]])[, 1]), ]
mse2.enet <- do.call(rbind, results5$mse[[8]])[
  order(do.call(rbind, results5$mse[[8]])[, 1]), ]
mse3.enet <- do.call(rbind, results5$mse[[9]])[
  order(do.call(rbind, results5$mse[[9]])[, 1]), ]
mse4.enet <- do.call(rbind, results5$mse[[10]])[
  order(do.call(rbind, results5$mse[[10]])[, 1]), ]

mse1.greben <- do.call(rbind, results5$mse[[11]])[
  order(do.call(rbind, results5$mse[[11]])[, 1]), ]
mse2.greben <- do.call(rbind, results5$mse[[12]])[
  order(do.call(rbind, results5$mse[[12]])[, 1]), ]
mse3.greben <- do.call(rbind, results5$mse[[13]])[
  order(do.call(rbind, results5$mse[[13]])[, 1]), ]
mse4.greben <- do.call(rbind, results5$mse[[14]])[
  order(do.call(rbind, results5$mse[[14]])[, 1]), ]

mse1.grridge <- do.call(rbind, results5$mse[[3]])[
  order(do.call(rbind, results5$mse[[3]])[, 1]), ]
mse2.grridge <- do.call(rbind, results5$mse[[4]])[
  order(do.call(rbind, results5$mse[[4]])[, 1]), ]
mse3.grridge <- do.call(rbind, results5$mse[[5]])[
  order(do.call(rbind, results5$mse[[5]])[, 1]), ]
mse4.grridge <- do.call(rbind, results5$mse[[6]])[
  order(do.call(rbind, results5$mse[[6]])[, 1]), ]

mse1.ridge <- do.call(rbind, results5$mse[[1]])[
  order(do.call(rbind, results5$mse[[1]])[, 1]), ]
mse2.ridge <- do.call(rbind, results5$mse[[2]])[
  order(do.call(rbind, results5$mse[[2]])[, 1]), ]

lmse1.enet <- lowess(mse1.enet[, 1], mse1.enet[, 2])
lmse2.enet <- lowess(mse2.enet[, 1], mse2.enet[, 2])
lmse3.enet <- lowess(mse3.enet[, 1], mse3.enet[, 2])
lmse4.enet <- lowess(mse4.enet[, 1], mse4.enet[, 2])

lmse1.greben <- lowess(mse1.greben[, 1], mse1.greben[, 2])
lmse2.greben <- lowess(mse2.greben[, 1], mse2.greben[, 2])
lmse3.greben <- lowess(mse3.greben[, 1], mse3.greben[, 2])
lmse4.greben <- lowess(mse4.greben[, 1], mse4.greben[, 2])

lmse1.grridge <- lowess(mse1.grridge[, 1], mse1.grridge[, 2])
lmse2.grridge <- lowess(mse2.grridge[, 1], mse2.grridge[, 2])
lmse3.grridge <- lowess(mse3.grridge[, 1], mse3.grridge[, 2])
lmse4.grridge <- lowess(mse4.grridge[, 1], mse4.grridge[, 2])

lmse1.ridge <- lowess(mse1.ridge[, 1], mse1.ridge[, 2])
lmse2.ridge <- lowess(mse2.ridge[, 1], mse2.ridge[, 2])

### diagnostics: checking lowess fits
# AUC
xlim2.auc <- range(auc1.enet[, 1], auc2.enet[, 1], auc3.enet[, 1], auc4.enet[, 1],
                   auc1.greben[, 1], auc2.greben[, 1], auc3.greben[, 1],
                   auc4.greben[, 1], auc3.grridge[, 1], auc4.grridge[, 1])
ylim2.auc <- range(auc1.enet[, 2], auc2.enet[, 2], auc3.enet[, 2], auc4.enet[, 2],
                   auc1.greben[, 2], auc2.greben[, 2], auc3.greben[, 2],
                   auc4.greben[, 2], auc3.grridge[, 2], auc4.grridge[, 2])

png(paste(path.graph, "grEBEN_test_res5_lowess_auc.png", sep=""),
    units="in", width=12, height=6, res=120)
par(mfrow=c(2, 5))
plot(auc3.grridge[, 1], auc3.grridge[, 2], col=1, ylab="AUC", main="a)",
     xlab="Number of selected variables", xlim=xlim2.auc, ylim=ylim2.auc)
lines(lauc3.grridge, col=2)

plot(auc1.greben[, 1], auc1.greben[, 2], col=1, ylab="AUC", main="b)",
     xlab="Number of selected variables", xlim=xlim2.auc, ylim=ylim2.auc)
lines(lauc1.greben, col=4)

plot(auc2.greben[, 1], auc2.greben[, 2], col=1, ylab="AUC", main="c)",
     xlab="Number of selected variables", xlim=xlim2.auc, ylim=ylim2.auc)
lines(lauc2.greben, col=5)

plot(auc3.greben[, 1], auc4.greben[, 2], col=1, ylab="AUC", main="d)",
     xlab="Number of selected variables", xlim=xlim2.auc, ylim=ylim2.auc)
lines(lauc3.greben, col=6)

plot(auc4.greben[, 1], auc4.greben[, 2], col=1, ylab="AUC", main="e)",
     xlab="Number of selected variables", xlim=xlim2.auc, ylim=ylim2.auc)
lines(lauc4.greben, col=7)

plot(auc4.grridge[, 1], auc4.grridge[, 2], col=1, ylab="AUC", main="f)",
     xlab="Number of selected variables", xlim=xlim2.auc, ylim=ylim2.auc)
lines(lauc4.grridge, col=3)

plot(auc1.enet[, 1], auc1.enet[, 2], col=1, ylab="AUC", main="g)",
     xlab="Number of selected variables", xlim=xlim2.auc, ylim=ylim2.auc)
lines(lauc1.enet, col=4)

plot(auc2.enet[, 1], auc2.enet[, 2], col=1, ylab="AUC", main="h)",
     xlab="Number of selected variables", xlim=xlim2.auc, ylim=ylim2.auc)
lines(lauc2.enet, col=5)

plot(auc3.enet[, 1], auc4.enet[, 2], col=1, ylab="AUC", main="i)",
     xlab="Number of selected variables", xlim=xlim2.auc, ylim=ylim2.auc)
lines(lauc3.enet, col=6)

plot(auc4.enet[, 1], auc4.enet[, 2], col=1, ylab="AUC", main="j)",
     xlab="Number of selected variables", xlim=xlim2.auc, ylim=ylim2.auc)
lines(lauc4.enet, col=7)
dev.off()

# briers
xlim2.briers <- range(briers1.enet[, 1], briers2.enet[, 1], briers3.enet[, 1], briers4.enet[, 1],
                      briers1.greben[, 1], briers2.greben[, 1], briers3.greben[, 1],
                      briers4.greben[, 1], briers3.grridge[, 1], briers4.grridge[, 1])
ylim2.briers <- range(briers1.enet[, 2], briers2.enet[, 2], briers3.enet[, 2], briers4.enet[, 2],
                      briers1.greben[, 2], briers2.greben[, 2], briers3.greben[, 2],
                      briers4.greben[, 2], briers3.grridge[, 2], briers4.grridge[, 2])

png(paste(path.graph, "grEBEN_test_res5_lowess_briers.png", sep=""),
    units="in", width=12, height=6, res=120)
par(mfrow=c(2, 5))
plot(briers3.grridge[, 1], briers3.grridge[, 2], col=1, ylab="Brier skill score", main="a)",
     xlab="Number of selected variables", xlim=xlim2.briers, ylim=ylim2.briers)
lines(lbriers3.grridge, col=2)

plot(briers1.greben[, 1], briers1.greben[, 2], col=1, ylab="Brier skill score", main="b)",
     xlab="Number of selected variables", xlim=xlim2.briers, ylim=ylim2.briers)
lines(lbriers1.greben, col=4)

plot(briers2.greben[, 1], briers2.greben[, 2], col=1, ylab="Brier skill score", main="c)",
     xlab="Number of selected variables", xlim=xlim2.briers, ylim=ylim2.briers)
lines(lbriers2.greben, col=5)

plot(briers3.greben[, 1], briers4.greben[, 2], col=1, ylab="Brier skill score", main="d)",
     xlab="Number of selected variables", xlim=xlim2.briers, ylim=ylim2.briers)
lines(lbriers3.greben, col=6)

plot(briers4.greben[, 1], briers4.greben[, 2], col=1, ylab="Brier skill score", main="e)",
     xlab="Number of selected variables", xlim=xlim2.briers, ylim=ylim2.briers)
lines(lbriers4.greben, col=7)

plot(briers4.grridge[, 1], briers4.grridge[, 2], col=1, ylab="Brier skill score", main="f)",
     xlab="Number of selected variables", xlim=xlim2.briers, ylim=ylim2.briers)
lines(lbriers4.grridge, col=3)

plot(briers1.enet[, 1], briers1.enet[, 2], col=1, ylab="Brier skill score", main="g)",
     xlab="Number of selected variables", xlim=xlim2.briers, ylim=ylim2.briers)
lines(lbriers1.enet, col=4)

plot(briers2.enet[, 1], briers2.enet[, 2], col=1, ylab="Brier skill score", main="h)",
     xlab="Number of selected variables", xlim=xlim2.briers, ylim=ylim2.briers)
lines(lbriers2.enet, col=5)

plot(briers3.enet[, 1], briers4.enet[, 2], col=1, ylab="Brier skill score", main="i)",
     xlab="Number of selected variables", xlim=xlim2.briers, ylim=ylim2.briers)
lines(lbriers3.enet, col=6)

plot(briers4.enet[, 1], briers4.enet[, 2], col=1, ylab="Brier skill score", main="j)",
     xlab="Number of selected variables", xlim=xlim2.briers, ylim=ylim2.briers)
lines(lbriers4.enet, col=7)
dev.off()

# mse
xlim2.mse <- range(mse1.enet[, 1], mse2.enet[, 1], mse3.enet[, 1], mse4.enet[, 1],
                   mse1.greben[, 1], mse2.greben[, 1], mse3.greben[, 1],
                   mse4.greben[, 1], mse3.grridge[, 1], mse4.grridge[, 1])
ylim2.mse <- range(mse1.enet[, 2], mse2.enet[, 2], mse3.enet[, 2], mse4.enet[, 2],
                   mse1.greben[, 2], mse2.greben[, 2], mse3.greben[, 2],
                   mse4.greben[, 2], mse3.grridge[, 2], mse4.grridge[, 2])

png(paste(path.graph, "grEBEN_test_res5_lowess_mse.png", sep=""),
    units="in", width=12, height=6, res=120)
par(mfrow=c(2, 5))
plot(mse3.grridge[, 1], mse3.grridge[, 2], col=1, ylab="MSE", main="a)",
     xlab="Number of selected variables", xlim=xlim2.mse, ylim=ylim2.mse)
lines(lmse3.grridge, col=2)

plot(mse1.greben[, 1], mse1.greben[, 2], col=1, ylab="MSE", main="b)",
     xlab="Number of selected variables", xlim=xlim2.mse, ylim=ylim2.mse)
lines(lmse1.greben, col=4)

plot(mse2.greben[, 1], mse2.greben[, 2], col=1, ylab="MSE", main="c)",
     xlab="Number of selected variables", xlim=xlim2.mse, ylim=ylim2.mse)
lines(lmse2.greben, col=5)

plot(mse3.greben[, 1], mse4.greben[, 2], col=1, ylab="MSE", main="d)",
     xlab="Number of selected variables", xlim=xlim2.mse, ylim=ylim2.mse)
lines(lmse3.greben, col=6)

plot(mse4.greben[, 1], mse4.greben[, 2], col=1, ylab="MSE", main="e)",
     xlab="Number of selected variables", xlim=xlim2.mse, ylim=ylim2.mse)
lines(lmse4.greben, col=7)

plot(mse4.grridge[, 1], mse4.grridge[, 2], col=1, ylab="MSE", main="f)",
     xlab="Number of selected variables", xlim=xlim2.mse, ylim=ylim2.mse)
lines(lmse4.grridge, col=3)

plot(mse1.enet[, 1], mse1.enet[, 2], col=1, ylab="MSE", main="g)",
     xlab="Number of selected variables", xlim=xlim2.mse, ylim=ylim2.mse)
lines(lmse1.enet, col=4)

plot(mse2.enet[, 1], mse2.enet[, 2], col=1, ylab="MSE", main="h)",
     xlab="Number of selected variables", xlim=xlim2.mse, ylim=ylim2.mse)
lines(lmse2.enet, col=5)

plot(mse3.enet[, 1], mse4.enet[, 2], col=1, ylab="MSE", main="i)",
     xlab="Number of selected variables", xlim=xlim2.mse, ylim=ylim2.mse)
lines(lmse3.enet, col=6)

plot(mse4.enet[, 1], mse4.enet[, 2], col=1, ylab="MSE", main="j)",
     xlab="Number of selected variables", xlim=xlim2.mse, ylim=ylim2.mse)
lines(lmse4.enet, col=7)
dev.off()



### Performance measures
# AUC
png(paste(path.graph, "grEBEN_test_res5_performance.png", sep=""),
    units="in", width=12, height=4, res=120)
par(mfrow=c(1, 3))
xlim1.auc <- range(lauc1.enet$x, lauc2.enet$x, lauc3.enet$x, lauc4.enet$x,
                   lauc1.greben$x, lauc2.greben$x, lauc3.greben$x,
                   lauc4.greben$x, lauc3.grridge$x, lauc4.grridge$x)
ylim1.auc <- range(lauc1.enet$y, lauc2.enet$y, lauc3.enet$y, lauc4.enet$y,
                   lauc1.greben$y, lauc2.greben$y, lauc3.greben$y,
                   lauc4.greben$y, lauc1.grridge$y, lauc2.grridge$y,
                   lauc3.grridge$y, lauc4.grridge$y, lauc1.ridge$y,
                   lauc2.ridge$y)

plot(0, 0, col=2, type="n", ylab="AUC", main="a)",
     xlab="Number of selected variables", xlim=xlim1.auc, ylim=ylim1.auc)

lines(lauc3.grridge, col=2, lty=1)
abline(h=lauc1.grridge$y, col=2, lty=2)
abline(h=lauc1.ridge$y, col=2, lty=3)

lines(lauc4.grridge, col=3, lty=1)
abline(h=lauc2.grridge$y, col=3, lty=2)
abline(h=lauc2.ridge$y, col=3, lty=3)

lines(lauc1.greben, col=4, lty=1)
lines(lauc1.enet, col=4, lty=3)

lines(lauc2.greben, col=5, lty=1)
lines(lauc2.enet, col=5, lty=3)

lines(lauc3.greben, col=6, lty=1)
lines(lauc3.enet, col=6, lty=3)

lines(lauc4.greben, col=7, lty=1)
lines(lauc4.enet, col=7, lty=3)

# Brier skilll
xlim1.briers <- range(lbriers1.enet$x, lbriers2.enet$x, lbriers3.enet$x,
                      lbriers4.enet$x, lbriers1.greben$x, lbriers2.greben$x,
                      lbriers3.greben$x, lbriers4.greben$x, lbriers3.grridge$x,
                      lbriers4.grridge$x)
ylim1.briers <- range(lbriers1.enet$y, lbriers2.enet$y, lbriers3.enet$y,
                      lbriers4.enet$y, lbriers1.greben$y, lbriers2.greben$y,
                      lbriers3.greben$y, lbriers4.greben$y, lbriers1.grridge$y,
                      lbriers2.grridge$y, lbriers3.grridge$y, lbriers4.grridge$y,
                      lbriers1.ridge$y, lbriers2.ridge$y)

plot(0, 0, col=2, type="n", ylab="Brier skill score", main="b)",
     xlab="Number of selected variables", xlim=xlim1.briers, ylim=ylim1.briers)

lines(lbriers3.grridge, col=2, lty=1)
abline(h=lbriers1.grridge$y, col=2, lty=2)
abline(h=lbriers1.ridge$y, col=2, lty=3)

lines(lbriers4.grridge, col=3, lty=1)
abline(h=lbriers2.grridge$y, col=3, lty=2)
abline(h=lbriers2.ridge$y, col=3, lty=3)

lines(lbriers1.greben, col=4, lty=1)
lines(lbriers1.enet, col=4, lty=3)

lines(lbriers2.greben, col=5, lty=1)
lines(lbriers2.enet, col=5, lty=3)

lines(lbriers3.greben, col=6, lty=1)
lines(lbriers3.enet, col=6, lty=3)

lines(lbriers4.greben, col=7, lty=1)
lines(lbriers4.enet, col=7, lty=3)

leglabels <- c(expression(paste("ridge, 'true' ", lambda)),
               "ridge", expression(paste("enet, true ", lambda)),
               expression(paste("enet, ", alpha==0.05)),
               expression(paste("enet, ", alpha==0.5)),
               expression(paste("enet, ", alpha==0.95)),
               "group-regularized + selection", "group-regularized",
               "not group-regularized")
legend("bottomright", legend=leglabels, fill=c(2:7, 0, 0, 0),
       lty=c(rep(NA, 6), 1, 2, 3), border=c(rep(1, 6), 0 ,0, 0), merge=TRUE,
       seg.len=1)

# MSE
xlim1.mse <- range(lmse1.enet$x, lmse2.enet$x, lmse3.enet$x,
                   lmse4.enet$x, lmse1.greben$x, lmse2.greben$x,
                   lmse3.greben$x, lmse4.greben$x, lmse3.grridge$x,
                   lmse4.grridge$x)
ylim1.mse <- range(lmse1.enet$y, lmse2.enet$y, lmse3.enet$y,
                   lmse4.enet$y, lmse1.greben$y, lmse2.greben$y,
                   lmse3.greben$y, lmse4.greben$y, lmse1.grridge$y,
                   lmse2.grridge$y, lmse3.grridge$y, lmse4.grridge$y,
                   lmse1.ridge$y, lmse2.ridge$y)

plot(0, 0, col=2, type="n", ylab="MSE", main="c)",
     xlab="Number of selected variables", xlim=xlim1.mse, ylim=ylim1.mse)

lines(lmse3.grridge, col=2, lty=1)
abline(h=lmse1.grridge$y, col=2, lty=2)
abline(h=lmse1.ridge$y, col=2, lty=3)

lines(lmse4.grridge, col=3, lty=1)
abline(h=lmse2.grridge$y, col=3, lty=2)
abline(h=lmse2.ridge$y, col=3, lty=3)

lines(lmse1.greben, col=4, lty=1)
lines(lmse1.enet, col=4, lty=3)

lines(lmse2.greben, col=5, lty=1)
lines(lmse2.enet, col=5, lty=3)

lines(lmse3.greben, col=6, lty=1)
lines(lmse3.enet, col=6, lty=3)

lines(lmse4.greben, col=7, lty=1)
lines(lmse4.enet, col=7, lty=3)

dev.off()

### penalty parameters
lambdag <- exp(seq(-1, 1, length.out=ncol(results5$lambdag[[1]])))
png(paste(path.graph, "grEBEN_test_res5_penalties.png", sep=""),
    units="in", width=8, height=6, res=120)
par(mfrow=c(2, 3))
boxplot(results5$lambdag[[1]], ylim=range(results5$lambdag[[1]], lambdag), 
        xlab="", ylab="", main="a)")
title(ylab=expression(paste(lambda[g], "'")), 
      xlab="Groups in decreasing effect size", line=2.5)
points(1:length(lambdag), lambdag, col=2, pch="x")
boxplot(results5$lambdag[[2]], ylim=range(results5$lambdag[[2]], lambdag), 
        xlab="", ylab="", main="b)")
title(ylab=expression(paste(lambda[g], "'")), 
      xlab="Groups in decreasing effect size", line=2.5)
points(1:length(lambdag), lambdag, col=2, pch="x")
boxplot(results5$lambdag[[3]], ylim=range(results5$lambdag[[3]], lambdag), 
        xlab="", ylab="", main="c)")
title(ylab=expression(paste(lambda[g], "'")), 
      xlab="Groups in decreasing effect size", line=2.5)
points(1:length(lambdag), lambdag, col=2, pch="x")
boxplot(results5$lambdag[[4]], ylim=range(results5$lambdag[[4]], lambdag), 
        xlab="", ylab="", main="d)")
title(ylab=expression(paste(lambda[g], "'")), 
      xlab="Groups in decreasing effect size", line=2.5)
points(1:length(lambdag), lambdag, col=2, pch="x")
boxplot(results5$lambdag[[5]], ylim=range(results5$lambdag[[5]], lambdag), 
        xlab="", ylab="", main="e)")
title(ylab=expression(paste(lambda[g], "'")), 
      xlab="Groups in decreasing effect size", line=2.5)
points(1:length(lambdag), lambdag, col=2, pch="x")
boxplot(results5$lambdag[[6]], ylim=range(results5$lambdag[[6]], lambdag), 
        xlab="", ylab="", main="f)")
title(ylab=expression(paste(lambda[g], "'")), 
      xlab="Groups in decreasing effect size", line=2.5)
points(1:length(lambdag), lambdag, col=2, pch="x")
dev.off()










