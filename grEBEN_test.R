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
library(pROC)

### source grEBEN
source(paste(path.code, "grVBEM.R", sep=""))
source(paste(path.code, "mygrridge.R", sep=""))

### simulations
## simulation 1
# create data
n <- 200
p <- 1000
G <- 5
pblock <- 20
rho <- 0.9
sigma <- matrix(rho, ncol=pblock, nrow=pblock); diag(sigma) <- 1
x <- do.call(cbind, replicate(p/pblock, rmvnorm(n, mean=rep(0, pblock), sigma=sigma), 
                              simplify=FALSE))
lambda <- 0.02
alpha <- 0.5
lambdag <- exp(seq(-2, 2, length.out=G))
beta <- as.numeric(sapply(1:G, function(g) {
  renbeta(p/G, 2*n*lambda*alpha*sqrt(lambdag[g]), n*lambda*(1 - alpha)*lambdag[g])}))
prob <- as.numeric(exp(x %*% beta)/(1 + exp(x %*% beta)))
y <- rbinom(n, 1, prob)
m <- rep(1, n)
partitions <- list(groups=rep(1:G, each=p/G))

# fit models
test1.greben <- grEBEN(x, y, m, partitions=partitions, alpha=alpha, lambda=lambda)
test2.greben <- grEBEN(x, y, m, partitions=partitions, alpha=0.05)
test3.greben <- grEBEN(x, y, m, partitions=partitions, alpha=0.5)
test4.greben <- grEBEN(x, y, m, partitions=partitions, alpha=0.95)
maxsel <- c(sum(test1.greben$beta[-1]!=0), sum(test2.greben$beta[-1]!=0), 
            sum(test3.greben$beta[-1]!=0), sum(test4.greben$beta[-1]!=0))
test1.grridge <- grridge(t(x), y, partitions=list(groups=CreatePartition(
  as.factor(partitions$groups))), selectionEN=TRUE, maxsel=maxsel[1],
  optl=n*lambda*(1 - alpha)/2 + n*lambda*alpha*sum(abs(beta))/sum(beta^2))
test2.grridge <- grridge(t(x), y, partitions=list(groups=CreatePartition(
  as.factor(partitions$groups))), selectionEN=TRUE, maxsel=maxsel[2],
  optl=n*lambda*(1 - alpha)/2 + n*lambda*alpha*sum(abs(beta))/sum(beta^2))
test3.grridge <- grridge(t(x), y, partitions=list(groups=CreatePartition(
  as.factor(partitions$groups))), selectionEN=TRUE, maxsel=maxsel[3],
  optl=n*lambda*(1 - alpha)/2 + n*lambda*alpha*sum(abs(beta))/sum(beta^2))
test4.grridge <- grridge(t(x), y, partitions=list(groups=CreatePartition(
  as.factor(partitions$groups))), selectionEN=TRUE, maxsel=maxsel[4],
  optl=n*lambda*(1 - alpha)/2 + n*lambda*alpha*sum(abs(beta))/sum(beta^2))
test5.grridge <- grridge(t(x), y, partitions=list(groups=CreatePartition(
  as.factor(partitions$groups))), selectionEN=TRUE, maxsel=maxsel[1])
test6.grridge <- grridge(t(x), y, partitions=list(groups=CreatePartition(
  as.factor(partitions$groups))), optl=test5.grridge$optl, selectionEN=TRUE, 
  maxsel=maxsel[2])
test7.grridge <- grridge(t(x), y, partitions=list(groups=CreatePartition(
  as.factor(partitions$groups))), optl=test5.grridge$optl, selectionEN=TRUE, 
  maxsel=maxsel[3])
test8.grridge <- grridge(t(x), y, partitions=list(groups=CreatePartition(
  as.factor(partitions$groups))), optl=test5.grridge$optl, selectionEN=TRUE, 
  maxsel=maxsel[4])

# create test data
ntest <- 1000
xtest <- do.call(cbind, replicate(p/pblock, rmvnorm(ntest, mean=rep(0, pblock),
                                                    sigma=sigma), simplify=FALSE))
probtest <- as.numeric(exp(xtest %*% beta)/(1 + exp(xtest %*% beta)))
ytest <- rbinom(ntest, 1, probtest)

# AUCs
auc.true <- pROC::roc(ytest, probtest)$auc

auc1.ridge <- pROC::roc(ytest, predict.grridge(test1.grridge, t(xtest))[, 1])$auc
auc2.ridge <- pROC::roc(ytest, predict.grridge(test5.grridge, t(xtest))[, 1])$auc

auc1.grridge <- pROC::roc(ytest, predict.grridge(test1.grridge, t(xtest))[, 2])$auc
auc2.grridge <- pROC::roc(ytest, predict.grridge(test5.grridge, t(xtest))[, 2])$auc

auc3.grridge <- pROC::roc(ytest, predict.grridge(test1.grridge, t(xtest))[, 3])$auc
auc4.grridge <- pROC::roc(ytest, predict.grridge(test2.grridge, t(xtest))[, 3])$auc
auc5.grridge <- pROC::roc(ytest, predict.grridge(test3.grridge, t(xtest))[, 3])$auc
auc6.grridge <- pROC::roc(ytest, predict.grridge(test4.grridge, t(xtest))[, 3])$auc
auc7.grridge <- pROC::roc(ytest, predict.grridge(test5.grridge, t(xtest))[, 3])$auc
auc8.grridge <- pROC::roc(ytest, predict.grridge(test6.grridge, t(xtest))[, 3])$auc
auc9.grridge <- pROC::roc(ytest, predict.grridge(test7.grridge, t(xtest))[, 3])$auc
auc10.grridge <- pROC::roc(ytest, predict.grridge(test8.grridge, t(xtest))[, 3])$auc

auc1.enet <- pROC::roc(ytest, predict.grEBEN(test1.greben, xtest, type="nogroups"))$auc
auc2.enet <- pROC::roc(ytest, predict.grEBEN(test2.greben, xtest, type="nogroups"))$auc
auc3.enet <- pROC::roc(ytest, predict.grEBEN(test3.greben, xtest, type="nogroups"))$auc
auc4.enet <- pROC::roc(ytest, predict.grEBEN(test4.greben, xtest, type="nogroups"))$auc

auc1.greben <- pROC::roc(ytest, predict.grEBEN(test1.greben, xtest, type="penalized"))$auc
auc2.greben <- pROC::roc(ytest, predict.grEBEN(test2.greben, xtest, type="penalized"))$auc
auc3.greben <- pROC::roc(ytest, predict.grEBEN(test3.greben, xtest, type="penalized"))$auc
auc4.greben <- pROC::roc(ytest, predict.grEBEN(test4.greben, xtest, type="penalized"))$auc

auc1.vbgreben <- pROC::roc(ytest, predict.grEBEN(test1.greben, xtest, type="VB"))$auc
auc2.vbgreben <- pROC::roc(ytest, predict.grEBEN(test2.greben, xtest, type="VB"))$auc
auc3.vbgreben <- pROC::roc(ytest, predict.grEBEN(test3.greben, xtest, type="VB"))$auc
auc4.vbgreben <- pROC::roc(ytest, predict.grEBEN(test4.greben, xtest, type="VB"))$auc

# Brier scores
briers.true <- 1 - sum((ytest - probtest)^2)/sum((ytest - mean(ytest))^2)

briers1.ridge <- 1 - sum((ytest - predict.grridge(test1.grridge, t(xtest))[, 1])^2)/sum((ytest - mean(ytest))^2)
briers2.ridge <- 1 - sum((ytest - predict.grridge(test5.grridge, t(xtest))[, 1])^2)/sum((ytest - mean(ytest))^2)

briers1.grridge <- 1 - sum((ytest - predict.grridge(test1.grridge, t(xtest))[, 2])^2)/sum((ytest - mean(ytest))^2)
briers2.grridge <- 1 - sum((ytest - predict.grridge(test5.grridge, t(xtest))[, 2])^2)/sum((ytest - mean(ytest))^2)

briers3.grridge <- 1 - sum((ytest - predict.grridge(test1.grridge, t(xtest))[, 3])^2)/sum((ytest - mean(ytest))^2)
briers4.grridge <- 1 - sum((ytest - predict.grridge(test2.grridge, t(xtest))[, 3])^2)/sum((ytest - mean(ytest))^2)
briers5.grridge <- 1 - sum((ytest - predict.grridge(test3.grridge, t(xtest))[, 3])^2)/sum((ytest - mean(ytest))^2)
briers6.grridge <- 1 - sum((ytest - predict.grridge(test4.grridge, t(xtest))[, 3])^2)/sum((ytest - mean(ytest))^2)
briers7.grridge <- 1 - sum((ytest - predict.grridge(test5.grridge, t(xtest))[, 3])^2)/sum((ytest - mean(ytest))^2)
briers8.grridge <- 1 - sum((ytest - predict.grridge(test6.grridge, t(xtest))[, 3])^2)/sum((ytest - mean(ytest))^2)
briers9.grridge <- 1 - sum((ytest - predict.grridge(test7.grridge, t(xtest))[, 3])^2)/sum((ytest - mean(ytest))^2)
briers10.grridge <- 1 - sum((ytest - predict.grridge(test8.grridge, t(xtest))[, 3])^2)/sum((ytest - mean(ytest))^2)

briers1.enet <- 1 - sum((ytest - predict.grEBEN(test1.greben, xtest, type="nogroups"))^2)/sum((ytest - mean(ytest))^2)
briers2.enet <- 1 - sum((ytest - predict.grEBEN(test2.greben, xtest, type="nogroups"))^2)/sum((ytest - mean(ytest))^2)
briers3.enet <- 1 - sum((ytest - predict.grEBEN(test3.greben, xtest, type="nogroups"))^2)/sum((ytest - mean(ytest))^2)
briers4.enet <- 1 - sum((ytest - predict.grEBEN(test4.greben, xtest, type="nogroups"))^2)/sum((ytest - mean(ytest))^2)

briers1.greben <- 1 - sum((ytest - predict.grEBEN(test1.greben, xtest, type="penalized"))^2)/sum((ytest - mean(ytest))^2)
briers2.greben <- 1 - sum((ytest - predict.grEBEN(test2.greben, xtest, type="penalized"))^2)/sum((ytest - mean(ytest))^2)
briers3.greben <- 1 - sum((ytest - predict.grEBEN(test3.greben, xtest, type="penalized"))^2)/sum((ytest - mean(ytest))^2)
briers4.greben <- 1 - sum((ytest - predict.grEBEN(test4.greben, xtest, type="penalized"))^2)/sum((ytest - mean(ytest))^2)

briers1.vbgreben <- 1 - sum((ytest - predict.grEBEN(test1.greben, xtest, type="VB"))^2)/sum((ytest - mean(ytest))^2)
briers2.vbgreben <- 1 - sum((ytest - predict.grEBEN(test2.greben, xtest, type="VB"))^2)/sum((ytest - mean(ytest))^2)
briers3.vbgreben <- 1 - sum((ytest - predict.grEBEN(test3.greben, xtest, type="VB"))^2)/sum((ytest - mean(ytest))^2)
briers4.vbgreben <- 1 - sum((ytest - predict.grEBEN(test4.greben, xtest, type="VB"))^2)/sum((ytest - mean(ytest))^2)

# MSE
mse.true <- 0

mse1.ridge <- mean((c(0, beta) - coef(test1.grridge$predobj$NoGroups, which="all"))^2)
mse2.ridge <- mean((c(0, beta) - coef(test5.grridge$predobj$NoGroups, which="all"))^2)

mse1.grridge <- mean((c(0, beta) - coef(test1.grridge$predobj$GroupRegul, which="all"))^2)
mse2.grridge <- mean((c(0, beta) - coef(test5.grridge$predobj$GroupRegul, which="all"))^2)

mse3.grridge <- mean((c(0, beta) - replace(rep(0, p + 1), c(1, test1.grridge$resEN$whichEN + 1), 
                                           coef(test1.grridge$predobj$EN, which="all")))^2)
mse4.grridge <- mean((c(0, beta) - replace(rep(0, p + 1), c(1, test2.grridge$resEN$whichEN + 1), 
                                           coef(test2.grridge$predobj$EN, which="all")))^2)
mse5.grridge <- mean((c(0, beta) - replace(rep(0, p + 1), c(1, test3.grridge$resEN$whichEN + 1), 
                                           coef(test3.grridge$predobj$EN, which="all")))^2)
mse6.grridge <- mean((c(0, beta) - replace(rep(0, p + 1), c(1, test4.grridge$resEN$whichEN + 1), 
                                           coef(test4.grridge$predobj$EN, which="all")))^2)
mse7.grridge <- mean((c(0, beta) - replace(rep(0, p + 1), c(1, test5.grridge$resEN$whichEN + 1), 
                                           coef(test5.grridge$predobj$EN, which="all")))^2)
mse8.grridge <- mean((c(0, beta) - replace(rep(0, p + 1), c(1, test6.grridge$resEN$whichEN + 1), 
                                           coef(test6.grridge$predobj$EN, which="all")))^2)
mse9.grridge <- mean((c(0, beta) - replace(rep(0, p + 1), c(1, test7.grridge$resEN$whichEN + 1), 
                                           coef(test7.grridge$predobj$EN, which="all")))^2)
mse10.grridge <- mean((c(0, beta) - replace(rep(0, p + 1), c(1, test8.grridge$resEN$whichEN + 1), 
                                            coef(test8.grridge$predobj$EN, which="all")))^2)

mse1.enet <- mean((c(0, beta) - test1.greben$beta.nogroups)^2)
mse2.enet <- mean((c(0, beta) - test2.greben$beta.nogroups)^2)
mse3.enet <- mean((c(0, beta) - test3.greben$beta.nogroups)^2)
mse4.enet <- mean((c(0, beta) - test4.greben$beta.nogroups)^2)

mse1.greben <- mean((c(0, beta) - test1.greben$beta)^2)
mse2.greben <- mean((c(0, beta) - test2.greben$beta)^2)
mse3.greben <- mean((c(0, beta) - test3.greben$beta)^2)
mse4.greben <- mean((c(0, beta) - test4.greben$beta)^2)

mse1.vbgreben <- mean((c(0, beta) - test1.greben$mu)^2)
mse2.vbgreben <- mean((c(0, beta) - test2.greben$mu)^2)
mse3.vbgreben <- mean((c(0, beta) - test3.greben$mu)^2)
mse4.vbgreben <- mean((c(0, beta) - test4.greben$mu)^2)

# combining everything
methods <- c("true", "ridge+true", "ridge", "grridge+true", "grridge", "grridge+true+sel1", 
             "grridge+true+sel2", "grridge+true+sel3", "grridge+true+sel4",
             "grridge+sel1", "grridge+sel2", "grridge+sel3", "grridge+sel4", "enet+true", 
             "enet+a=0.05", "enet+a=0.5", "enet+a=0.95", 
             "greben+true", "greben+a=0.05", "greben+a=0.5", "greben+a=0.95", 
             "vbgreben+true", "vbgreben+a=0.05", "vbgreven+a=0.5", "vbgreven+a=0.95")
results1 <- list(psel=setNames(c(rep(p, 5), 
                                 length(test1.grridge$resEN$whichEN), 
                                 length(test2.grridge$resEN$whichEN), 
                                 length(test3.grridge$resEN$whichEN), 
                                 length(test4.grridge$resEN$whichEN), 
                                 length(test5.grridge$resEN$whichEN), 
                                 length(test6.grridge$resEN$whichEN), 
                                 length(test7.grridge$resEN$whichEN), 
                                 length(test8.grridge$resEN$whichEN), 
                                 sum(test1.greben$beta.nogroups[-1]!=0), 
                                 sum(test2.greben$beta.nogroups[-1]!=0), sum(test3.greben$beta.nogroups[-1]!=0),
                                 sum(test4.greben$beta.nogroups[-1]!=0), sum(test1.greben$beta[-1]!=0), 
                                 sum(test2.greben$beta[-1]!=0), sum(test3.greben$beta[-1]!=0), 
                                 sum(test4.greben$beta[-1]!=0), rep(p, 4)), methods),
                 auc=setNames(c(auc.true, auc1.ridge, auc2.ridge, auc1.grridge, auc2.grridge, 
                                auc3.grridge, auc4.grridge, auc5.grridge, auc6.grridge, 
                                auc7.grridge, auc8.grridge, auc9.grridge, auc10.grridge, 
                                auc1.enet, auc2.enet, auc3.enet, auc4.enet, auc1.greben, 
                                auc2.greben, auc3.greben, auc4.greben, auc1.vbgreben, 
                                auc2.vbgreben, auc3.vbgreben, auc4.vbgreben), methods),
                 briers=setNames(c(briers.true, briers1.ridge, briers2.ridge, briers1.grridge, 
                                   briers2.grridge, briers3.grridge, briers4.grridge, 
                                   briers5.grridge, briers6.grridge, briers7.grridge, 
                                   briers8.grridge, briers9.grridge, briers10.grridge, 
                                   briers1.enet, briers2.enet, briers3.enet, briers4.enet, briers1.greben, 
                                   briers2.greben, briers3.greben, briers4.greben, briers1.vbgreben, 
                                   briers2.vbgreben, briers3.vbgreben, briers4.vbgreben), methods),
                 mse=setNames(c(mse.true, mse1.ridge, mse2.ridge, mse1.grridge, 
                                mse2.grridge, mse3.grridge, mse4.grridge, 
                                mse5.grridge, mse6.grridge, mse7.grridge, 
                                mse8.grridge, mse9.grridge, mse10.grridge, 
                                mse1.enet, mse2.enet, mse3.enet, mse4.enet, mse1.greben, 
                                mse2.greben, mse3.greben, mse4.greben, mse1.vbgreben, 
                                mse2.vbgreben, mse3.vbgreben, mse4.vbgreben), methods))

do.call("cbind", results1)

varbeta <- sapply(1:G, function(g) {var(beta[(100*(g - 1) + 1):(100*g)])})

### graphs
# graph simulation 1
leglabels <- c(expression(paste("grridge, true ", lambda)), "grridge", 
               expression(paste("true ", lambda[g])), 
               expression(paste("greben, true ", lambda)), 
               expression(paste("greben, true ", alpha)),
               expression(paste("greben, ", alpha==0.95)),
               expression(paste("greben, ", alpha==0.5)),
               "ridge and enet")
xlabels <- paste("Group ", c(1:G), sep="")

bar1 <- barplot(rbind(test1.grridge$lambdamults$groups, test2.grridge$lambdamults$groups, lambdag,
                      test1.greben$lambdag$groups[, test1.greben$nouteriter + 1],
                      test2.greben$lambdag$groups[, test2.greben$nouteriter + 1],
                      test3.greben$lambdag$groups[, test3.greben$nouteriter + 1],
                      test4.greben$lambdag$groups[, test4.greben$nouteriter + 1]), 
                beside=TRUE, ylab=expression(paste(lambda[g], "'")), axes=FALSE, axisnames=FALSE, 
                legend.text=leglabels,
                args.legend=list(x="topleft", fill=c(gray.colors(7), 0), lty=c(rep(NA, 7), 2),
                                 border=c(rep(1, 7), 0), merge=TRUE, seg.len=1))
axis(2)
text(colMeans(bar1), par("usr")[3] - 0.25, srt=45, adj=1, labels=xlabels, xpd=TRUE)
abline(h=1, lty=2)




n <- 100
p <- 10
x <- matrix(rnorm(n*p), nrow=n, ncol=p)
beta <- c(1:10)
y <- as.numeric(x %*% beta + rnorm(n))

R <- cor(x)
r <- t(cor(y, x))
S <- cov(x)
D <- diag(apply(x, 2, sd))
beta <- solve(S) %*% t(cov(y, x))








