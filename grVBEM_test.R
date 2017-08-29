##############################  preamble  #############################
# testing grMCEM                                                      #
# version: 01                                                         #
# author: Magnus M?nch                                                #
# created: 15-03-2017                                                 #
# last edited: 15-03-2017                                             #
#######################################################################

###############################  notes  ###############################
#                                                                     #
#######################################################################

### paths
path.code <- as.character(ifelse(Sys.info()[1]=="Darwin","/Users/magnusmunch/Documents/PhD/EBEN/code/" ,"~/EBEN/code/"))
path.graph <- "/Users/magnusmunch/Documents/PhD/EBEN/graphs/"
path.data <- as.character(ifelse(Sys.info()[1]=="Darwin","/Users/magnusmunch/Documents/PhD/EBEN/data/" ,"~/EBEN/data/"))
path.res <- as.character(ifelse(Sys.info()[1]=="Darwin","/Users/magnusmunch/Documents/PhD/EBEN/results/" ,"~/EBEN/results/"))

## libraries
library(GRridge)
library(mvtnorm)
library(grpreg)
library(SGL)
library(psych)
library(glmnet)

# source grENVB functions
source(paste(path.code, "grVBEM.R", sep=""))

### simulation 4 (one partition)
set.seed(123)
n <- 200
p <- 150
G <- 3
groups <- rep(1:G, each=p/G)
x <- matrix(rnorm(n*p), ncol=p, nrow=n)
lambda1 <- 1
lambda2 <- 1
lambdag <- c(0.2, 3, 1/(0.2*3))
m <- rep(1, n)
b0 <- rnorm(p + 1)
sigma0 <- diag(rchisq(p + 1, 1))
beta <- c(renbeta(p/G, lambda1*sqrt(lambdag[1]), lambda2*lambdag[1]), 
          renbeta(p/G, lambda1*sqrt(lambdag[2]), lambda2*lambdag[2]),
          renbeta(p/G, lambda1*sqrt(lambdag[3]), lambda2*lambdag[3]))
y <- rbinom(n, m, exp(x %*% beta)/(1 + exp(x %*% beta)))

test4.grVBEM1 <- grVBEM2(x, y, m, list(groups=groups), lambda1=NULL, lambda2=NULL, intercept=TRUE, 
                         monotone=FALSE, posterior=FALSE, ELBO=TRUE, eps=0.00001, maxiter=100, trace=TRUE, 
                         alphastart=NULL)
brew test4.grVBEM2 <- grVBEM(x, y, m, groups, lambda1=NULL, lambda2=NULL, intercept=TRUE, 
                        eps=0.00001, maxiter=100, trace=TRUE, QNacc=TRUE)
test4.grridge <- grridge(t(x), y, list(CreatePartition(as.factor(groups))), unpenal=~1)
test4.grpreg <- cv.grpreg(x, y, group=groups, penalty="grLasso", family="binomial")
test4.SGL <- cvSGL(list(x=x, y=y), type="logit")

ntest <- 1000
xtest <- matrix(rnorm(ntest*p), ncol=p, nrow=ntest)
ytest <- rbinom(ntest, m, exp(xtest %*% beta)/(1 + exp(xtest %*% beta)))

pred4.grVBEM1 <- as.numeric(exp(cbind(1, xtest) %*% test4.grVBEM1$mu)/
                             (1 + exp(cbind(1, xtest) %*% test4.grVBEM1 $mu)))
pred4.grVBEM2 <- as.numeric(exp(cbind(1, xtest) %*% test4.grVBEM2$mu)/
                              (1 + exp(cbind(1, xtest) %*% test4.grVBEM2$mu)))
pred4.grridge <- predict.grridge(test4.grridge, t(xtest))[, 2]
pred4.grpreg <- predict(test4.grpreg, xtest, type="response")
pred4.SGL <- as.numeric(exp(cbind(1, xtest) %*% c(test4.SGL$fit$intercepts[which.min(test4.SGL$lldiff)], test4.SGL$fit$beta[, which.min(test4.SGL$lldiff)]))/
                          (1 + exp(cbind(1, xtest) %*% c(test4.SGL$fit$intercepts[which.min(test4.SGL$lldiff)], test4.SGL$fit$beta[, which.min(test4.SGL$lldiff)]))))
pred4.truth <- as.numeric(exp(xtest %*% beta)/(1 + exp(xtest %*% beta)))

auc4.grVBEM1 <- pROC::roc(ytest, pred4.grVBEM1)$auc
auc4.grVBEM2 <- pROC::roc(ytest, pred4.grVBEM2)$auc
auc4.grridge <- pROC::roc(ytest, pred4.grridge)$auc
auc4.grpreg <- pROC::roc(ytest, pred4.grpreg)$auc
auc4.SGL <- pROC::roc(ytest, pred4.SGL)$auc
auc4.truth <- pROC::roc(ytest, pred4.truth)$auc

barplot(rbind(test4.grridge$lambdamults[[1]], lambdag, 
              test4.grVBEM1$lambdag[[1]][, test4.grVBEM1$nouteriter + 1],
              test4.grVBEM2$lambdag[[1]][, test4.grVBEM2$nouteriter + 1]), beside=TRUE,
        names.arg=c(expression(lambda[1]), expression(lambda[2]), expression(lambda[3])),
        legend.text=paste(c("GRridge", "truth", "VBEM", "VBEM + QN", "grlasso", "SGL"), 
                          paste(", AUC=", round(c(auc4.grridge, auc4.truth, auc4.grVBEM1, 
                                                  auc4.grVBEM2, auc4.grpreg, auc4.SGL), 2)), 
                          sep=""), args.legend=list(x="topleft"))


### simulation 5 (one noise partition)
set.seed(123)
n <- 200
p <- 150
G <- 3
groups <- list(rep(1:G, each=p/G), rep(1:5, each=p/5))
x <- matrix(rnorm(n*p), ncol=p, nrow=n)
lambda1 <- 1
lambda2 <- 1
lambdag <- c(0.2, 3, 1/(0.2*3))
m <- rep(1, n)
b0 <- rnorm(p + 1)
sigma0 <- diag(rchisq(p + 1, 1))
beta <- c(renbeta(p/G, lambda1*sqrt(lambdag[1]), lambda2*lambdag[1]), 
          renbeta(p/G, lambda1*sqrt(lambdag[2]), lambda2*lambdag[2]),
          renbeta(p/G, lambda1*sqrt(lambdag[3]), lambda2*lambdag[3]))
y <- rbinom(n, m, exp(x %*% beta)/(1 + exp(x %*% beta)))


test5.grVBEM <- grVBEM(x, y, m, groups, lambda1=NULL, lambda2=NULL, intercept=TRUE, 
                       eps=0.001, maxiter=100)
test5.grridge <- grridge(t(x), y, list(CreatePartition(as.factor(groups[[1]])), 
                                       CreatePartition(as.factor(groups[[2]]))), unpenal=~1)
test5.grpreg <- cv.grpreg(x, y, group=groups, penalty="grLasso", family="binomial")
test5.SGL <- cvSGL(list(x=x, y=y), type="logit")

ntest <- 1000
xtest <- matrix(rnorm(ntest*p), ncol=p, nrow=ntest)
ytest <- rbinom(ntest, m, exp(xtest %*% beta)/(1 + exp(xtest %*% beta)))

pred5.grVBEM <- as.numeric(exp(cbind(1, xtest) %*% test5.grVBEM$mu)/
                             (1 + exp(cbind(1, xtest) %*% test5.grVBEM$mu)))
pred5.grridge <- predict.grridge(test5.grridge, t(xtest))[, 2]
pred5.grpreg <- predict(test5.grpreg, xtest, type="response")
pred5.SGL <- as.numeric(exp(cbind(1, xtest) %*% c(test5.SGL$fit$intercepts[which.min(test5.SGL$lldiff)], test5.SGL$fit$beta[, which.min(test5.SGL$lldiff)]))/
                          (1 + exp(cbind(1, xtest) %*% c(test5.SGL$fit$intercepts[which.min(test5.SGL$lldiff)], test5.SGL$fit$beta[, which.min(test5.SGL$lldiff)]))))
pred5.truth <- as.numeric(exp(xtest %*% beta)/(1 + exp(xtest %*% beta)))

auc5.grVBEM <- pROC::roc(ytest, pred5.grVBEM)$auc
auc5.grridge <- pROC::roc(ytest, pred5.grridge)$auc
auc5.grpreg <- pROC::roc(ytest, pred5.grpreg)$auc
auc5.SGL <- pROC::roc(ytest, pred5.SGL)$auc
auc5.truth <- pROC::roc(ytest, pred5.truth)$auc

barplot(rbind(test5.grridge$lambdamults[[1]], lambdag, 
              test5.grVBEM$lambdag[[1]][, test5.grVBEM$nouteriter + 1]), beside=TRUE,
        names.arg=c(expression(lambda[1]), expression(lambda[2]), expression(lambda[3])),
        legend.text=paste(c("GRridge", "truth", "VBEM", "grlasso", "SGL"), 
                          paste(", AUC=", round(c(auc5.grridge, auc5.truth, auc5.grVBEM, 
                                                  auc5.grpreg, auc5.SGL), 2)), sep=""),
        args.legend=list(x="topleft"))

barplot(rbind(test5.grridge$lambdamults[[2]], rep(1, 5), 
              test5.grVBEM$lambdag[[2]][, test5.grVBEM$nouteriter + 1]), beside=TRUE,
        names.arg=c(expression(lambda[1]), expression(lambda[2]), expression(lambda[3]), expression(lambda[4]), 
                    expression(lambda[5])),
        legend.text=paste(c("GRridge", "truth", "VBEM", "grlasso", "SGL"), 
                          paste(", AUC=", round(c(auc5.grridge, auc5.truth, auc5.grVBEM, 
                                                  auc5.grpreg, auc5.SGL), 2)), sep=""),
        args.legend=list(x="topright"))

plot(beta, test5.grVBEM$mu[-1])
plot(beta, test5.grridge$betas)
plot(beta, coef(test5.grpreg)[-1])
plot(beta, test5.SGL$fit$beta[, which.min(test5.SGL$lldiff)])

plot(test5.grVBEM$lambdag[[1]][1, ], ylim=range(test5.grVBEM$lambdag[[1]]),
     type="l")
lines(test5.grVBEM$lambdag[[1]][2, ], col=2)
lines(test5.grVBEM$lambdag[[1]][3, ], col=3)

plot(test5.grVBEM$lambdag[[2]][1, ], ylim=range(test5.grVBEM$lambdag[[2]]),
     type="l")
lines(test5.grVBEM$lambdag[[2]][2, ], col=2)
lines(test5.grVBEM$lambdag[[2]][3, ], col=3)
lines(test5.grVBEM$lambdag[[2]][4, ], col=4)
lines(test5.grVBEM$lambdag[[2]][5, ], col=5)

plot(test5.grVBEM$lowermll, type="l")

### simulation 6 (group wise l1 and l2 penalization)
set.seed(123)
n <- 200
p <- 150
G <- 3
groups <- rep(1:G, each=p/G)
x <- matrix(rnorm(n*p), ncol=p, nrow=n)
m <- rep(1, n)
b0 <- rnorm(p + 1)
sigma0 <- diag(rchisq(p + 1, 1))
beta <- c(rep(0.5, p/G), rep(0.8, p/G), rep(1.5, p/G))
y <- rbinom(n, m, exp(x %*% beta)/(1 + exp(x %*% beta)))

test6.grVBEM <- grVBEM(x, y, m, groups, lambda1=NULL, lambda2=NULL, intercept=TRUE, 
                       eps=0.001, maxiter=100)
test6.grridge <- grridge(t(x), y, list(group=CreatePartition(as.factor(groups))), unpenal=~1)
test6.grpreg <- cv.grpreg(x, y, group=groups, penalty="grLasso", family="binomial")

ntest <- 1000
xtest <- matrix(rnorm(ntest*p), ncol=p, nrow=ntest)
ytest <- rbinom(ntest, m, exp(xtest %*% beta)/(1 + exp(xtest %*% beta)))

pred6.grVBEM <- as.numeric(exp(cbind(1, xtest) %*% test6.grVBEM$mu)/
                             (1 + exp(cbind(1, xtest) %*% test6.grVBEM$mu)))
pred6.grridge <- predict.grridge(test6.grridge, t(xtest))[, 2]
pred6.grpreg <- predict(test6.grpreg, xtest, type="response")
pred6.truth <- as.numeric(exp(xtest %*% beta)/(1 + xtest %*% beta))

auc6.grVBEM <- pROC::roc(ytest, pred6.grVBEM)$auc
auc6.grridge <- pROC::roc(ytest, pred6.grridge)$auc
auc6.grpreg <- pROC::roc(ytest, pred6.grpreg)$auc
auc6.truth <- pROC::roc(ytest, pred6.truth)$auc

barplot(rbind(test6.grridge$lambdamults$group,
              test6.grVBEM$lambdag[, test6.grVBEM$nouteriter + 1]), 
        beside=TRUE, names.arg=c(expression(lambda[1]), expression(lambda[2]), expression(lambda[3])),
        legend.text=paste(c("GRridge", "VBEM", "truth", "grlasso"), 
                          paste(", AUC=", round(c(auc6.grridge, auc6.grVBEM, auc6.truth,
                                                  auc6.grpreg), 2)), sep=""))
plot(beta, test6.grVBEM$mu[-1])
plot(beta, test6.grridge$betas)
plot(beta, coef(test6.grpreg)[-1])

### simulation 7 (group wise l1 and l2 penalization)
set.seed(123)
n <- 200
p <- 150
G <- 3
groups <- rep(1:G, each=p/G)
x <- matrix(rnorm(n*p), ncol=p, nrow=n)
m <- rep(1, n)
b0 <- rnorm(p + 1)
sigma0 <- diag(rchisq(p + 1, 1))
beta <- c(rep(0, p/G), rep(0.5, p/G), rep(1.5, p/G))
y <- rbinom(n, m, exp(x %*% beta)/(1 + exp(x %*% beta)))

test7.grVBEM <- grVBEM(x, y, m, groups, lambda1=NULL, lambda2=NULL, sigma0, intercept=TRUE, 
                       eps=0.001, maxiter=100)
test7.grridge <- grridge(t(x), y, list(group=CreatePartition(as.factor(groups))), unpenal=~1)

ntest <- 1000
xtest <- matrix(rnorm(ntest*p), ncol=p, nrow=ntest)
ytest <- rbinom(ntest, m, exp(xtest %*% beta)/(1 + exp(xtest %*% beta)))

pred7.grVBEM <- as.numeric(exp(cbind(1, xtest) %*% test7.grVBEM$mu)/
                             (1 + exp(cbind(1, xtest) %*% test7.grVBEM$mu)))
pred7.grridge <- predict.grridge(test7.grridge, t(xtest))[, 2]
pred7.truth <- as.numeric(exp(xtest %*% beta)/(1 + exp(xtest %*% beta)))

auc7.grVBEM <- pROC::roc(ytest, pred7.grVBEM)$auc
auc7.grridge <- pROC::roc(ytest, pred7.grridge)$auc
auc7.truth <- pROC::roc(ytest, pred7.truth)$auc

barplot(rbind(test7.grridge$lambdamults$group, 
              test7.grVBEM$lambdag[, test7.grVBEM$nouteriter + 1]), 
        beside=TRUE, names.arg=c(expression(lambda[1]), expression(lambda[2]), expression(lambda[3])),
        legend.text=paste(c("GRridge", "VBEM", "truth"), 
                          paste(", AUC=", round(c(auc7.grridge, auc7.grVBEM, auc7.truth), 2)), sep=""))
plot(beta, test7.grVBEM$mu[-1])
plot(beta, test7.grridge$betas)

### simulation 9 (group wise l1 and l2 penalization)
set.seed(123)
set.seed(2345)
n <- 200
p <- 150
G <- 3
groups <- rep(1:G, each=p/G)
rho <- 0.5
sigma <- matrix(rho, ncol=p, nrow=p)
diag(sigma) <- 1
x <- rmvnorm(n, mean=rep(0, p), sigma=sigma)
m <- rep(1, n)
b0 <- rnorm(p + 1)
sigma0 <- diag(rchisq(p + 1, 1))
beta <- c(rep(0, p/G), rep(0.5, p/G), rep(1.5, p/G))
y <- rbinom(n, m, exp(x %*% beta)/(1 + exp(x %*% beta)))

test9.grVBEM <- grVBEM(x, y, m, groups, lambda1=NULL, lambda2=NULL, sigma0, intercept=TRUE, 
                       eps=0.001, maxiter=100)
test9.grridge <- grridge(t(x), y, list(group=CreatePartition(as.factor(groups))), unpenal=~1)

ntest <- 1000
xtest <- matrix(rnorm(ntest*p), ncol=p, nrow=ntest)
ytest <- rbinom(ntest, m, exp(xtest %*% beta)/(1 + exp(xtest %*% beta)))

pred9.grVBEM <- as.numeric(exp(cbind(1, xtest) %*% test9.grVBEM$mu)/
                             (1 + exp(cbind(1, xtest) %*% test9.grVBEM$mu)))
pred9.grridge <- predict.grridge(test9.grridge, t(xtest))[, 2]
pred9.truth <- as.numeric(exp(xtest %*% beta)/(1 + exp(xtest %*% beta)))

auc9.grVBEM <- pROC::roc(ytest, pred9.grVBEM)$auc
auc9.grridge <- pROC::roc(ytest, pred9.grridge)$auc
auc9.truth <- pROC::roc(ytest, pred9.truth)$auc

barplot(rbind(test9.grridge$lambdamults$group, 
              test9.grVBEM$lambdag[, test9.grVBEM$nouteriter + 1]), 
        beside=TRUE, names.arg=c(expression(lambda[1]), expression(lambda[2]), expression(lambda[3])),
        legend.text=paste(c("GRridge", "VBEM", "truth"), 
                          paste(", AUC=", round(c(auc9.grridge, auc9.grVBEM, auc9.truth), 2)), sep=""))
plot(beta, test9.grVBEM$mu[-1])
plot(beta, test9.grridge$betas)

### simulation 10 (group wise l1 and l2 penalization)
set.seed(123)
n <- 100
p <- 150
G <- 3
groups <- rep(1:G, each=p/G)
rho <- 0.3
sigma <- matrix(rho, ncol=p, nrow=p)
diag(sigma) <- 1
x <- rmvnorm(n, mean=rep(0, p), sigma=sigma)
m <- rep(1, n)
b0 <- rnorm(p + 1)
sigma0 <- diag(rchisq(p + 1, 1))
beta <- c(rep(0, p/G), rep(0.5, p/G), rep(1.5, p/G))
y <- rbinom(n, m, exp(x %*% beta)/(1 + exp(x %*% beta)))

test10.grVBEM <- grVBEM(x, y, m, groups, lambda1=NULL, lambda2=NULL, sigma0, intercept=TRUE, 
                        eps=0.001, maxiter=100)
test10.grridge <- grridge(t(x), y, list(group=CreatePartition(as.factor(groups))), unpenal=~1)

ntest <- 1000
xtest <- matrix(rnorm(ntest*p), ncol=p, nrow=ntest)
ytest <- rbinom(ntest, m, exp(xtest %*% beta)/(1 + exp(xtest %*% beta)))

pred10.grVBEM <- as.numeric(exp(cbind(1, xtest) %*% test10.grVBEM$mu)/
                              (1 + exp(cbind(1, xtest) %*% test10.grVBEM$mu)))
pred10.grridge <- predict.grridge(test10.grridge, t(xtest))[, 2]
pred10.truth <- as.numeric(exp(xtest %*% beta)/(1 + exp(xtest %*% beta)))

auc10.grVBEM <- pROC::roc(ytest, pred10.grVBEM)$auc
auc10.grridge <- pROC::roc(ytest, pred10.grridge)$auc
auc10.truth <- pROC::roc(ytest, pred10.truth)$auc

barplot(rbind(test10.grridge$lambdamults$group, 
              test10.grVBEM$lambdag[, test10.grVBEM$nouteriter + 1]), 
        beside=TRUE, names.arg=c(expression(lambda[1]), expression(lambda[2]), expression(lambda[3])),
        legend.text=paste(c("GRridge", "VBEM", "truth"), 
                          paste(", AUC=", round(c(auc10.grridge, auc10.grVBEM, auc10.truth), 2)), sep=""))
plot(beta, test10.grVBEM$mu[-1])
plot(beta, test10.grridge$betas)

### simulation 11 (group wise l1 and l2 penalization)
set.seed(123)
n <- 100
p <- 150
G <- 3
groups <- rep(1:G, each=p/G)
rho <- 0.3
sigma <- matrix(rho, ncol=p, nrow=p)
diag(sigma) <- 1
x <- rmvnorm(n, mean=rep(0, p), sigma=sigma)
m <- rep(1, n)
b0 <- rnorm(p + 1)
sigma0 <- diag(rchisq(p + 1, 1))
beta <- rep(0, p)
y <- rbinom(n, m, exp(x %*% beta)/(1 + exp(x %*% beta)))

test11.grVBEM <- grVBEM(x, y, m, groups, lambda1=NULL, lambda2=NULL, sigma0, intercept=TRUE, 
                        eps=0.001, maxiter=100)
test11.grridge <- grridge(t(x), y, list(group=CreatePartition(as.factor(groups))), unpenal=~1)

ntest <- 1000
xtest <- matrix(rnorm(ntest*p), ncol=p, nrow=ntest)
ytest <- rbinom(ntest, m, exp(xtest %*% beta)/(1 + exp(xtest %*% beta)))

pred11.grVBEM <- as.numeric(exp(cbind(1, xtest) %*% test11.grVBEM$mu)/
                              (1 + exp(cbind(1, xtest) %*% test11.grVBEM$mu)))
pred11.grridge <- predict.grridge(test11.grridge, t(xtest))[, 2]
pred11.truth <- as.numeric(exp(xtest %*% beta)/(1 + exp(xtest %*% beta)))

auc11.grVBEM <- pROC::roc(ytest, pred11.grVBEM)$auc
auc11.grridge <- pROC::roc(ytest, pred11.grridge)$auc
auc11.truth <- pROC::roc(ytest, pred11.truth)$auc

barplot(rbind(test11.grridge$lambdamults$group, 
              test11.grVBEM$lambdag[, test11.grVBEM$nouteriter + 1]), 
        beside=TRUE, names.arg=c(expression(lambda[1]), expression(lambda[2]), expression(lambda[3])),
        legend.text=paste(c("GRridge", "VBEM", "truth"), 
                          paste(", AUC=", round(c(auc11.grridge, auc11.grVBEM, auc11.truth), 2)), sep=""))
plot(beta, test11.grVBEM$mu[-1])
plot(beta, test11.grridge$betas)

### simulation 12 (group wise l1 and l2 penalization)
set.seed(123)
n <- 100
p <- 150
G <- 3
groups <- rep(1:G, each=p/G)
rho <- 0.3
sigma <- matrix(rho, ncol=p, nrow=p)
diag(sigma) <- 1
x <- rmvnorm(n, mean=rep(0, p), sigma=sigma)
m <- rep(1, n)
b0 <- rnorm(p + 1)
sigma0 <- diag(rchisq(p + 1, 1))
beta <- rep(1, p)
y <- rbinom(n, m, exp(x %*% beta)/(1 + exp(x %*% beta)))

test12.grVBEM <- grVBEM(x, y, m, groups, lambda1=NULL, lambda2=NULL, sigma0, intercept=TRUE, 
                        eps=0.001, maxiter=100)
test12.grVBEM2 <- penalized(y, x, unpenalized=~1, lambda1=0.5*test12.grVBEM$lambda1*sqrt(rep(test12.grVBEM$lambdag[, test12.grVBEM$nouteriter + 1], times=rle(groups)$lengths)),
                            lambda2=0.5*test12.grVBEM$lambda2*rep(test12.grVBEM$lambdag[, test12.grVBEM$nouteriter + 1], times=rle(groups)$lengths),
                            model="logistic")
test12.grridge <- grridge(t(x), y, list(group=CreatePartition(as.factor(groups))), unpenal=~1)

ntest <- 1000
xtest <- matrix(rnorm(ntest*p), ncol=p, nrow=ntest)
ytest <- rbinom(ntest, m, exp(xtest %*% beta)/(1 + exp(xtest %*% beta)))

pred12.grVBEM <- as.numeric(exp(cbind(1, xtest) %*% test12.grVBEM$mu)/
                              (1 + exp(cbind(1, xtest) %*% test12.grVBEM$mu)))
pred12.grVBEM2 <- predict(test12.grVBEM2, xtest)
pred12.grridge <- predict.grridge(test12.grridge, t(xtest))[, 2]
pred12.truth <- as.numeric(exp(xtest %*% beta)/(1 + exp(xtest %*% beta)))

auc12.grVBEM <- pROC::roc(ytest, pred12.grVBEM)$auc
auc12.grVBEM2 <- pROC::roc(ytest, pred12.grVBEM2)$auc
auc12.grridge <- pROC::roc(ytest, pred12.grridge)$auc
auc12.truth <- pROC::roc(ytest, pred12.truth)$auc

barplot(rbind(test12.grridge$lambdamults$group, 
              test12.grVBEM$lambdag[, test12.grVBEM$nouteriter + 1]), 
        beside=TRUE, names.arg=c(expression(lambda[1]), expression(lambda[2]), expression(lambda[3])),
        legend.text=paste(c("GRridge", "VBEM", "VBEM + EN", "truth"), 
                          paste(", AUC=", round(c(auc12.grridge, auc12.grVBEM, auc12.grVBEM2, auc12.truth), 2)), sep=""))
plot(beta, test12.grVBEM$mu[-1])
plot(beta, test12.grridge$betas)

# simulation 13
set.seed(456)
n <- 100
p <- 900
G <- 3
groups <- rep(1:G, each=p/G)
rho <- 0.3
sigma <- matrix(rho, ncol=p, nrow=p)
diag(sigma) <- 1
x <- rmvnorm(n, mean=rep(0, p), sigma=sigma)
m <- rep(1, n)
b0 <- rnorm(p + 1)
sigma0 <- diag(rchisq(p + 1, 1))
beta <- rep(0, p)
y <- rbinom(n, m, exp(x %*% beta)/(1 + exp(x %*% beta)))

test13.grVBEM1 <- grVBEM(x, y, m, groups, lambda1=NULL, lambda2=NULL, intercept=TRUE, 
                         eps=0.001, maxiter=100, trace=TRUE, QNacc=FALSE)
test13.grVBEM2 <- grVBEM(x, y, m, groups, lambda1=NULL, lambda2=NULL, intercept=TRUE, 
                         eps=0.001, maxiter=100, trace=TRUE, QNacc=TRUE)
test13.grridge <- grridge(t(x), y, unpenal=~1, innfold=10,
                          list(group=CreatePartition(as.factor(groups))))

ntest <- 1000
xtest <- matrix(rnorm(ntest*p), ncol=p, nrow=ntest)
ytest <- rbinom(ntest, m, exp(xtest %*% beta)/(1 + exp(xtest %*% beta)))

pred13.grVBEM <- as.numeric(exp(cbind(1, xtest) %*% test13.grVBEM$mu)/
                              (1 + exp(cbind(1, xtest) %*% test13.grVBEM$mu)))
pred13.grridge <- predict.grridge(test13.grridge, t(xtest))[, 2]
pred13.truth <- as.numeric(exp(xtest %*% beta)/(1 + exp(xtest %*% beta)))

auc13.grVBEM <- pROC::roc(ytest, pred13.grVBEM)$auc
auc13.grridge <- pROC::roc(ytest, pred13.grridge)$auc
auc13.truth <- pROC::roc(ytest, pred13.truth)$auc

barplot(rbind(test13.grridge$lambdamults$group, 
              test13.grVBEM$lambdag[, test13.grVBEM$nouteriter + 1]), 
        beside=TRUE, names.arg=c(expression(lambda[1]), expression(lambda[2]), expression(lambda[3])),
        legend.text=paste(c("GRridge", "VBEM", "truth"), 
                          paste(", AUC=", round(c(auc13.grridge, auc13.grVBEM, auc13.truth), 2)), sep=""))
plot(beta, test13.grVBEM$mu[-1])
plot(beta, test13.grridge$betas)

# simulation 14
set.seed(456)
n <- 100
p <- 900
G <- 3
groups <- rep(1:G, each=p/G)
rho <- 0.3
sigma <- matrix(rho, ncol=p, nrow=p)
diag(sigma) <- 1
x <- rmvnorm(n, mean=rep(0, p), sigma=sigma)
m <- rep(1, n)
b0 <- rnorm(p + 1)
sigma0 <- diag(rchisq(p + 1, 1))
beta <- rep(0.005, p)
y <- rbinom(n, m, exp(x %*% beta)/(1 + exp(x %*% beta)))

test14.grVBEM <- grVBEM(x, y, m, groups, lambda1=NULL, lambda2=NULL, sigma0, intercept=TRUE, 
                        eps=0.001, maxiter=100)
test14.grridge <- grridge(t(x), y, unpenal=~1, #innfold=10,
                          list(group=CreatePartition(as.factor(groups))))

ntest <- 1000
xtest <- matrix(rnorm(ntest*p), ncol=p, nrow=ntest)
ytest <- rbinom(ntest, m, exp(xtest %*% beta)/(1 + exp(xtest %*% beta)))

pred14.grVBEM <- as.numeric(exp(cbind(1, xtest) %*% test14.grVBEM$mu)/
                              (1 + exp(cbind(1, xtest) %*% test14.grVBEM$mu)))
pred14.grridge <- predict.grridge(test14.grridge, t(xtest))[, 2]
pred14.truth <- as.numeric(exp(xtest %*% beta)/(1 + exp(xtest %*% beta)))

auc14.grVBEM <- pROC::roc(ytest, pred14.grVBEM)$auc
auc14.grridge <- pROC::roc(ytest, pred14.grridge)$auc
auc14.truth <- pROC::roc(ytest, pred14.truth)$auc

barplot(rbind(test14.grridge$lambdamults$group, 
              test14.grVBEM$lambdag[, test14.grVBEM$nouteriter + 1]), 
        beside=TRUE, names.arg=c(expression(lambda[1]), expression(lambda[2]), expression(lambda[3])),
        legend.text=paste(c("GRridge", "VBEM", "truth"), 
                          paste(", AUC=", round(c(auc14.grridge, auc14.grVBEM, auc14.truth), 2)), sep=""))
plot(beta, test14.grVBEM$mu[-1])
plot(beta, test14.grridge$betas)

# simulation 15
set.seed(456)
n <- 100
p <- 900
G <- 3
groups <- rep(1:G, each=p/G)
rho <- 0.3
sigma <- matrix(rho, ncol=p, nrow=p)
diag(sigma) <- 1
x <- rmvnorm(n, mean=rep(0, p), sigma=sigma)
m <- rep(1, n)
b0 <- rnorm(p + 1)
sigma0 <- diag(rchisq(p + 1, 1))
beta <- c(rep(0, p/G), rep(0.004, p/G), rep(0.008, p/G))
y <- rbinom(n, m, exp(x %*% beta)/(1 + exp(x %*% beta)))

test15.grVBEM <- grVBEM(x, y, m, groups, lambda1=NULL, lambda2=NULL, sigma0, intercept=TRUE, 
                        eps=0.001, maxiter=100)
test15.grridge <- grridge(t(x), y, unpenal=~1, #innfold=10,
                          list(group=CreatePartition(as.factor(groups))))

ntest <- 1000
xtest <- matrix(rnorm(ntest*p), ncol=p, nrow=ntest)
ytest <- rbinom(ntest, m, exp(xtest %*% beta)/(1 + exp(xtest %*% beta)))

pred15.grVBEM <- as.numeric(exp(cbind(1, xtest) %*% test15.grVBEM$mu)/
                              (1 + exp(cbind(1, xtest) %*% test15.grVBEM$mu)))
pred15.grridge <- predict.grridge(test15.grridge, t(xtest))[, 2]
pred15.truth <- as.numeric(exp(xtest %*% beta)/(1 + exp(xtest %*% beta)))

auc15.grVBEM <- pROC::roc(ytest, pred15.grVBEM)$auc
auc15.grridge <- pROC::roc(ytest, pred15.grridge)$auc
auc15.truth <- pROC::roc(ytest, pred15.truth)$auc

barplot(rbind(test15.grridge$lambdamults$group, 
              test15.grVBEM$lambdag[, test15.grVBEM$nouteriter + 1]), 
        beside=TRUE, names.arg=c(expression(lambda[1]), expression(lambda[2]), expression(lambda[3])),
        legend.text=paste(c("GRridge", "VBEM", "truth"), 
                          paste(", AUC=", round(c(auc15.grridge, auc15.grVBEM, auc15.truth), 2)), sep=""))
plot(beta, test15.grVBEM$mu[-1])
plot(beta, test15.grridge$betas)

### data 1 (group wise l1 and l2 penalization)
load(paste(path.data, "mirsData.RData", sep=""))
parCons <- CreatePartition(mirsData$conservation) # using conservation status as grouping
set.seed(123)
groups <- rep(1:length(parCons), unlist(lapply(parCons, length)))
x <- apply(t(as.matrix(mirsData$transformedData))[, unlist(parCons)], 2, 
           function(x) {(x - mean(x))/sd(x)})       
y <- as.numeric(mirsData$response) - 1
n <- nrow(x)
p <- ncol(x)
m <- rep(1, n)

test16.grVBEM <- grVBEM(x, y, m, groups, lambda1=NULL, lambda2=NULL, intercept=TRUE, eps=0.001, 
                        maxiter=500)
test16.grridge <- grridge(t(x), y, list(group=CreatePartition(as.factor(groups))), unpenal=~1)

barplot(rbind(test16.grridge$lambdamults$group, 
              test16.grVBEM$lambdag[, test16.grVBEM$nouteriter + 1]), 
        beside=TRUE, names.arg=c(expression(lambda[1]), expression(lambda[2]), expression(lambda[3])),
        legend.text=c("GRridge", "VBEM"))

pen.fit <- penalized(y, x, unpenalized=~1, model="logistic",
                     lambda1=rep(sqrt(test16.grridge$lambdamults$group)*test16.grVBEM$lambda1, rle(groups)$lengths),
                     lambda2=rep(test16.grridge$lambdamults$group*test16.grVBEM$lambda2, rle(groups)$lengths))


### data 2 (group wise l1 and l2 penalization) permuting the groups (no info)
load(paste(path.data, "mirsData.RData", sep=""))
parCons <- CreatePartition(mirsData$conservation) # using conservation status as grouping
set.seed(123)
n <- ncol(mirsData$transformedData)
p <- nrow(mirsData$transformedData)
groups <- rep(1:length(parCons), unlist(lapply(parCons, length)))
x <- apply(t(as.matrix(mirsData$transformedData))[, unlist(parCons)], 2, 
           function(x) {(x - mean(x))/sd(x)})[, sample(1:p)]       
y <- as.numeric(mirsData$response) - 1
m <- rep(1, n)
fit.optL2 <- optL2(y, x, unpenalized=~1, lambda1=0, model="logistic")
b0 <- coef(fit.optL2$fullfit)
pred.b0 <- as.numeric(exp(cbind(1, x) %*% b0)/(1 + exp(cbind(1, x) %*% b0)))
W <- diag(sqrt(pred.b0*(1 - pred.b0)))
Xw <- W %*% cbind(1, x)
invmat <- solve(t(Xw) %*% Xw + diag(c(0, rep(2*fit.optL2$lambda, p))))
sigma0 <- invmat %*% t(Xw) %*% Xw %*% invmat

test17.grVBEM <- grVBEM(x, y, m, groups, lambda1=NULL, lambda2=NULL, sigma0, 
                        intercept=TRUE, eps=0.001, maxiter=500)
test17.grridge <- grridge(t(x), y, list(group=CreatePartition(as.factor(groups))), unpenal=~1)

barplot(rbind(test17.grridge$lambdamults$group, 
              test17.grVBEM$lambdag[, test17.grVBEM$nouteriter + 1]), 
        beside=TRUE, names.arg=c(expression(lambda[1]), expression(lambda[2]), expression(lambda[3])),
        legend.text=c("GRridge", "VBEM"))

### simulation 18 (more groups)
set.seed(789)
n <- 200
p <- 400
G <- 4
groups <- rep(1:G, each=p/G)
rho <- 0.3
sigma <- matrix(rho, ncol=p, nrow=p)
diag(sigma) <- 1
x <- rmvnorm(n, mean=rep(0, p), sigma=sigma)
m <- rep(1, n)
b0 <- rnorm(p + 1)
sigma0 <- diag(rchisq(p + 1, 1))
beta <- c(rep(0, p/G), rep(0, p/G), rep(0.06, p/G), rep(0.06, p/G))
y <- rbinom(n, m, exp(x %*% beta)/(1 + exp(x %*% beta)))

test18.grVBEM <- grVBEM(x, y, m, groups, lambda1=NULL, lambda2=NULL, sigma0, 
                        intercept=TRUE, eps=0.001, maxiter=500)
test18.grridge <- grridge(t(x), y, list(group=CreatePartition(as.factor(groups))), unpenal=~1)

ntest <- 1000
xtest <- matrix(rnorm(ntest*p), ncol=p, nrow=ntest)
ytest <- rbinom(ntest, m, exp(xtest %*% beta)/(1 + exp(xtest %*% beta)))

pred18.grVBEM <- as.numeric(exp(cbind(1, xtest) %*% test18.grVBEM$mu)/
                              (1 + exp(cbind(1, xtest) %*% test18.grVBEM$mu)))
pred18.grridge <- predict.grridge(test18.grridge, t(xtest))[, 2]
pred18.truth <- as.numeric(exp(xtest %*% beta)/(1 + exp(xtest %*% beta)))

auc18.grVBEM <- pROC::roc(ytest, pred18.grVBEM)$auc
auc18.grridge <- pROC::roc(ytest, pred18.grridge)$auc
auc18.truth <- pROC::roc(ytest, pred18.truth)$auc

barplot(rbind(test18.grridge$lambdamults$group, 
              test18.grVBEM$lambdag[, test18.grVBEM$nouteriter + 1]), 
        beside=TRUE, names.arg=c(expression(lambda[1]), expression(lambda[2]), 
                                 expression(lambda[3]), expression(lambda[4])),
        legend.text=paste(c("GRridge", "VBEM", "truth"), 
                          paste(", AUC=", round(c(auc18.grridge, auc18.grVBEM, auc18.truth), 2)), sep=""))
plot(beta, test18.grVBEM$mu[-1])
plot(beta, test18.grridge$betas)

### data 3 (group wise l1 and l2 penalization) permuting the groups (no info)
load(paste(path.data, "mirsData.RData", sep=""))
parAbund <- CreatePartition(rowSums(mirsData$countData), mingr=25, ngroup=10, decreasing=TRUE) # using abundance as grouping
set.seed(123)
n <- ncol(mirsData$transformedData)
p <- nrow(mirsData$transformedData)
groups <- rep(1:length(parAbund), unlist(lapply(parAbund, length)))
G <- length(unique(groups))
x <- apply(t(as.matrix(mirsData$transformedData))[, unlist(parAbund)], 2, 
           function(x) {(x - mean(x))/sd(x)})[, sample(1:p)]       
y <- as.numeric(mirsData$response) - 1
m <- rep(1, n)
fit.optL2 <- optL2(y, x, unpenalized=~1, lambda1=0, model="logistic")
b0 <- coef(fit.optL2$fullfit)
pred.b0 <- as.numeric(exp(cbind(1, x) %*% b0)/(1 + exp(cbind(1, x) %*% b0)))
W <- diag(sqrt(pred.b0*(1 - pred.b0)))
Xw <- W %*% cbind(1, x)
invmat <- solve(t(Xw) %*% Xw + diag(c(0, rep(2*fit.optL2$lambda, p))))
sigma0 <- invmat %*% t(Xw) %*% Xw %*% invmat

test19.grVBEM <- grVBEM(x, y, m, groups, lambda1=NULL, lambda2=NULL, sigma0, 
                        intercept=TRUE, eps=0.001, maxiter=500)
test19.grridge <- grridge(t(x), y, list(group=CreatePartition(as.factor(groups))), unpenal=~1)

names <- vector(mode="expression", G)
for(g in 1:G) {
  names[g] <- substitute(expression(lambda[i]), list(i=g))[2]
}
barplot(rbind(test19.grridge$lambdamults$group, 
              test19.grVBEM$lambdag[, test19.grVBEM$nouteriter + 1]), 
        beside=TRUE, names.arg=names,
        legend.text=c("GRridge", "VBEM"))

# data 4 (permuted data (no info))
load(paste(path.data, "mirsData.RData", sep=""))
parCons <- CreatePartition(mirsData$conservation) # using conservation status as grouping
set.seed(123)
n <- ncol(mirsData$transformedData)
p <- nrow(mirsData$transformedData)
groups <- rep(1:length(parCons), unlist(lapply(parCons, length)))
x <- apply(t(as.matrix(mirsData$transformedData))[, unlist(parCons)], 2, 
           function(x) {(x - mean(x))/sd(x)})[, sample(1:p)]       
y <- as.numeric(mirsData$response) - 1
m <- rep(1, n)
fit.optL2 <- optL2(y, x, unpenalized=~1, lambda1=0, model="logistic")
b0 <- coef(fit.optL2$fullfit)
pred.b0 <- as.numeric(exp(cbind(1, x) %*% b0)/(1 + exp(cbind(1, x) %*% b0)))
W <- diag(sqrt(pred.b0*(1 - pred.b0)))
Xw <- W %*% cbind(1, x)
invmat <- solve(t(Xw) %*% Xw + diag(c(0, rep(2*fit.optL2$lambda, p))))
sigma0 <- invmat %*% t(Xw) %*% Xw %*% invmat
nfolds <- 10
rest <- n %% nfolds
foldsize <- c(rep(n %/% nfolds + as.numeric(rest!=0), times=rest),
              rep(n %/% nfolds, times=nfolds - rest))
foldid <- sample(rep(1:nfolds, times=foldsize))

pred20.grVBEM <- numeric(n)
for(k in 1:nfolds) {
  print(paste("Fold", k, sep=" "))
  xtrain <- x[foldid!=k, ]
  xtest <- x[foldid==k, ]
  ytrain <- y[foldid!=k]
  mtrain <- rep(1, times=sum(foldid!=k))
  
  fit.grVBEM <- grVBEM(xtrain, ytrain, m=mtrain, groups=groups, lambda1=NULL, 
                       lambda2=NULL, sigma0=sigma0, intercept=TRUE, 
                       eps=0.001, maxiter=500)
  best <- fit.grVBEM$mu
  pred20.grVBEM[foldid==k] <- as.numeric(exp(cbind(1, xtest) %*% best)/
                                           (1 + exp(cbind(1, xtest) %*% best)))
  
}

# test20.grVBEM <- cv.grVBEM(x, y, m, groups, lambda1=NULL, lambda2=NULL, sigma0, 
#                            intercept=TRUE, eps=0.001, maxiter=500, nfolds=nfolds,
#                            foldid=foldid)
test20a.grridge <- grridge(t(x), y, list(group=CreatePartition(as.factor(groups))), unpenal=~1)
test20b.grridge <- grridgeCV(test20a.grridge, t(x), y, outerfold=foldid, fixedfolds=TRUE)

auc20.grVBEM <- pROC::roc(y, pred20.grVBEM)$auc
auc20.grridge <- pROC::roc(y, test20b.grridge[, 3])$auc

# data 5 (checking predictive performance by CV)
load(paste(path.data, "mirsData.RData", sep=""))
parCons <- CreatePartition(mirsData$conservation) # using conservation status as grouping
set.seed(123)
groups <- rep(1:length(parCons), unlist(lapply(parCons, length)))
x <- apply(t(as.matrix(mirsData$transformedData))[, unlist(parCons)], 2, 
           function(x) {(x - mean(x))/sd(x)})       
y <- as.numeric(mirsData$response) - 1
n <- nrow(x)
p <- ncol(x)
m <- rep(1, n)
fit.optL2 <- optL2(y, x, unpenalized=~1, lambda1=0, model="logistic")
b0 <- coef(fit.optL2$fullfit)
pred.b0 <- as.numeric(exp(cbind(1, x) %*% b0)/(1 + exp(cbind(1, x) %*% b0)))
W <- diag(sqrt(pred.b0*(1 - pred.b0)))
Xw <- W %*% cbind(1, x)
invmat <- solve(t(Xw) %*% Xw + diag(c(0, rep(2*fit.optL2$lambda, p))))
sigma0 <- invmat %*% t(Xw) %*% Xw %*% invmat
nfolds <- 10
rest <- n %% nfolds
foldsize <- c(rep(n %/% nfolds + as.numeric(rest!=0), times=rest),
              rep(n %/% nfolds, times=nfolds - rest))
foldid <- sample(rep(1:nfolds, times=foldsize))

pred21.grVBEM <- numeric(n)
for(k in 1:nfolds) {
  print(paste("Fold", k, sep=" "))
  xtrain <- x[foldid!=k, ]
  xtest <- x[foldid==k, ]
  ytrain <- y[foldid!=k]
  mtrain <- rep(1, times=sum(foldid!=k))
  
  fit.grVBEM <- grVBEM(xtrain, ytrain, m=mtrain, groups=groups, lambda1=NULL, 
                       lambda2=NULL, sigma0=sigma0, intercept=TRUE, 
                       eps=0.001, maxiter=500)
  best <- fit.grVBEM$mu
  pred21.grVBEM[foldid==k] <- as.numeric(exp(cbind(1, xtest) %*% best)/
                                           (1 + exp(cbind(1, xtest) %*% best)))
  
}

test21a.grridge <- grridge(t(x), y, list(group=CreatePartition(as.factor(groups))), unpenal=~1)
test21b.grridge <- grridgeCV(test21a.grridge, t(x), y, outerfold=foldid, fixedfolds=TRUE)

pred21.en <- numeric(n)
for(k in 1:nfolds) {
  print(paste("Fold", k, sep=" "))
  xtrain <- x[foldid!=k, ]
  xtest <- x[foldid==k, ]
  ytrain <- y[foldid!=k]
  
  fit.cvpen <- cv.pen(xtrain, ytrain, intercept=TRUE)
  fit.en <- 
    best <- fit.grVBEM$mu
  pred21.grVBEM[foldid==k] <- as.numeric(exp(cbind(1, xtest) %*% best)/
                                           (1 + exp(cbind(1, xtest) %*% best)))
  
}

auc21.grVBEM <- pROC::roc(y, pred21.grVBEM)$auc
auc21.grridge <- pROC::roc(y, test21b.grridge[, 3])$auc
cbind(grVBEM=round(pred21.grVBEM, 2), grridge=round(test21b.grridge[, 3], 2))

### simulation 19 (variable selection)
set.seed(123)
n <- 100
p <- 150
G <- 3
groups <- rep(1:G, each=p/G)
rho <- 0.3
sigma <- matrix(rho, ncol=p, nrow=p)
diag(sigma) <- 1
x <- rmvnorm(n, mean=rep(0, p), sigma=sigma)
m <- rep(1, n)
q <- 0.3 # proportion set to zero per group
beta <- c(rep(0, q*p/G), rep(0.03, (1 - q)*p/G), rep(0, q*p/G), rep(0.05, (1 - q)*p/G), rep(0, q*p/G), 
          rep(0.08, (1 - q)*p/G))
y <- rbinom(n, m, exp(x %*% beta)/(1 + exp(x %*% beta)))

test22.grVBEM <- grVBEM(x, y, m, groups, lambda1=NULL, lambda2=NULL, intercept=TRUE, eps=0.001, maxiter=500)
test22.grVBEMsel <- penalized(y, x, unpenalized=~1, model="logistic",
                              lambda1=rep(0.5*test22.grVBEM$lambda1*sqrt(test22.grVBEM$lambdag[, test22.grVBEM$nouteriter + 1]), each=p/G),
                              lambda2=rep(0.5*test22.grVBEM$lambda2*test22.grVBEM$lambdag[, test22.grVBEM$nouteriter + 1], each=p/G))
test22.grridge <- grridge(t(x), y, list(group=CreatePartition(as.factor(groups))), unpenal=~1)
test22.grridgesel <- grridge(t(x), y, list(group=CreatePartition(as.factor(groups))), unpenal=~1, selectionEN=TRUE,
                             maxsel=sum(test22.grVBEMsel@penalized!=0))

ntest <- 1000
xtest <- matrix(rnorm(ntest*p), ncol=p, nrow=ntest)
ytest <- rbinom(ntest, m, exp(xtest %*% beta)/(1 + exp(xtest %*% beta)))

pred22.grVBEM <- as.numeric(exp(cbind(1, xtest) %*% test22.grVBEM$mu)/
                              (1 + exp(cbind(1, xtest) %*% test22.grVBEM$mu)))
pred22.grVBEMsel <- predict(test22.grVBEMsel, xtest)
pred22.grridge <- predict.grridge(test22.grridge, t(xtest))[, 2]
pred22.grridgesel <- predict.grridge(test22.grridgesel, t(xtest))[, 2]
pred22.truth <- as.numeric(exp(xtest %*% beta)/(1 + exp(xtest %*% beta)))

kappa22.grVBEM <- cohen.kappa(cbind(as.numeric(beta!=0), as.numeric(test22.grVBEMsel@penalized!=0)))$kappa
kappa22.grridge <- cohen.kappa(cbind(as.numeric(beta!=0), 
                                     sapply(1:p, function(j) {ifelse(j %in% test22.grridgesel$resEN$whichEN, 1, 0)})))$kappa

auc22.grVBEM <- pROC::roc(ytest, pred22.grVBEM)$auc
auc22.grVBEMsel <- pROC::roc(ytest, pred22.grVBEMsel)$auc
auc22.grridge <- pROC::roc(ytest, pred22.grridge)$auc
auc22.grridgesel <- pROC::roc(ytest, pred22.grridgesel)$auc
auc22.truth <- pROC::roc(ytest, pred22.truth)$auc

barplot(rbind(test22.grridge$lambdamults$group, 
              test22.grVBEM$lambdag[, test22.grVBEM$nouteriter + 1]), 
        beside=TRUE, names.arg=c(expression(lambda[1]), expression(lambda[2]), 
                                 expression(lambda[3])),
        legend.text=paste(c("GRridge", "GRridge + sel", "VBEM", "VBEM + sel", "truth"), 
                          paste(", AUC=", round(c(auc22.grridge, auc22.grridgesel, auc22.grVBEM, auc22.grVBEMsel, 
                                                  auc22.truth), 2)), sep=""))
plot(beta, test22.grVBEM$mu[-1])
plot(beta, test22.grridge$betas)



# data 6 (Verlaat data)
set.seed(2017)
data(dataVerlaat)
partitions <- list(annotation=sort(as.numeric(CpGann)))
x <- apply(datcenVerlaat, 1, function(x) {(x - mean(x))/sd(x)})[, order(as.numeric(CpGann))]
y <- respVerlaat
n <- nrow(x)
p <- ncol(x)
m <- rep(1, n)

test23.grVBEM <- grVBEM2(x, y, m, partitions, lambda1=0.8936099, lambda2=0.2624093, 
                         intercept=TRUE, posterior=FALSE, eps=0.001, maxiter=500, trace=TRUE)
test23.grridge <- grridge(t(x), y, list(annotation=CreatePartition(as.factor(partitions[[1]]))), 
                          unpenal=~1, optl=22033)

barplot(rbind(test23.grridge$lambdamults$annotation, 
              test23.grVBEM$lambdag$annotation[, test23.grVBEM$nouteriter + 1]), 
        beside=TRUE, names.arg=levels(CpGann), legend.text=c("GRridge", "VBEM"),
        args.legend=list(x="topleft", fill=c(gray.colors(2), 0, 0), lty=c(NA, NA, 2, 2),
                         border=c(rep(1, 2), 0, 0), merge=TRUE, seg.len=1))
abline(h=1, lty=2)


# simulation 15
set.seed(234)
n <- 100
p <- 900
G <- 3
groups <- rep(1:G, each=p/G)
rho <- 0.3
sigma <- matrix(rho, ncol=p, nrow=p)
diag(sigma) <- 1
x <- rmvnorm(n, mean=rep(0, p), sigma=sigma)
m <- rep(1, n)
b0 <- rnorm(p + 1)
sigma0 <- diag(rchisq(p + 1, 1))
beta <- c(rep(0, p/G), rep(0.004, p/G), rep(0.008, p/G))
y <- rbinom(n, m, exp(x %*% beta)/(1 + exp(x %*% beta)))

test24.grVBEM <- grVBEM2(x, y, m, groups, lambda1=4.586998, lambda2=0.1726289, 
                         intercept=TRUE, posterior=TRUE, eps=0.001, maxiter=100, trace=TRUE)
test24.grridge <- grridge(t(x), y, unpenal=~1, #innfold=10,
                          list(group=CreatePartition(as.factor(groups))))

ntest <- 1000
xtest <- matrix(rnorm(ntest*p), ncol=p, nrow=ntest)
ytest <- rbinom(ntest, m, exp(xtest %*% beta)/(1 + exp(xtest %*% beta)))

pred24.grVBEM <- as.numeric(exp(cbind(1, xtest) %*% test24.grVBEM$mu)/
                              (1 + exp(cbind(1, xtest) %*% test24.grVBEM$mu)))
pred24.grridge <- predict.grridge(test24.grridge, t(xtest))[, 2]
pred24.truth <- as.numeric(exp(xtest %*% beta)/(1 + exp(xtest %*% beta)))

auc24.grVBEM <- pROC::roc(ytest, pred24.grVBEM)$auc
auc24.grridge <- pROC::roc(ytest, pred24.grridge)$auc
auc24.truth <- pROC::roc(ytest, pred24.truth)$auc

barplot(rbind(test24.grridge$lambdamults$group, 
              test24.grVBEM$lambdag$partition1[, test24.grVBEM$nouteriter + 1]), 
        beside=TRUE, names.arg=c(expression(lambda[1]), expression(lambda[2]), 
                                 expression(lambda[3])),
        legend.text=paste(c("GRridge", "VBEM", "truth"), 
                          paste(", AUC=", round(c(auc24.grridge, auc24.grVBEM, auc24.truth), 2)), 
                          sep=""),
        args.legend=list(x="bottomright"))
abline(h=1, lty=2)
plot(beta, test24.grVBEM$mu[-1])
plot(beta, test24.grridge$betas)
plot(as.numeric(test24.grVBEM$mu)[-1], test24.grridge$betas)

# simulation 16
set.seed(234)
n <- 100
p <- 900
G <- 9
partitions <- list(mono=rep(1:G, each=p/G))
rho <- 0.3
sigma <- matrix(rho, ncol=p, nrow=p)
diag(sigma) <- 1
x <- rmvnorm(n, mean=rep(0, p), sigma=sigma)
m <- rep(1, n)
b0 <- rnorm(p + 1)
sigma0 <- diag(rchisq(p + 1, 1))
beta <- c(rep(0.1, p/G), rep(0, p - p/G))
y <- rbinom(n, m, exp(x %*% beta)/(1 + exp(x %*% beta)))

test25.grVBEM <- grVBEM2(x, y, m, partitions, lambda1=5.220753, lambda2=0.02636744, intercept=TRUE, 
                         monotone=TRUE, posterior=FALSE, eps=0.00001, maxiter=500, trace=TRUE,
                         alphastart=c(0.3, exp(seq(0.1, 0.3, length.out=G-1))))


plot(test25.grVBEM$lambdag$mono[1, ], type="l", ylim=range(test25.grVBEM$lambdag$mono))
lines(test25.grVBEM$lambdag$mono[2, ], col=2)
lines(test25.grVBEM$lambdag$mono[3, ], col=3)
lines(test25.grVBEM$lambdag$mono[4, ], col=4)
lines(test25.grVBEM$lambdag$mono[5, ], col=5)
lines(test25.grVBEM$lambdag$mono[6, ], col=6)
lines(test25.grVBEM$lambdag$mono[7, ], col=7)
lines(test25.grVBEM$lambdag$mono[8, ], col=8)
lines(test25.grVBEM$lambdag$mono[9, ], col=9)





