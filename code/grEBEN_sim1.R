### paths
path.data <- ifelse(as.character(Sys.info()[1])!="Darwin", "~/EBEN/data/",
                    "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/data/")
path.code <- ifelse(as.character(Sys.info()[1])!="Darwin", "~/EBEN/code/",
                    "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/code/")
path.res <- ifelse(as.character(Sys.info()[1])!="Darwin", "~/EBEN/results/",
                   "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/results/")
path.graph <- "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/graphs/"

### loading libraries
library(GRridge)
library(glmnet)
library(pROC)
library(irr)
library(gglasso)
library(sp)
library(SGL)

### set color scheme for graphs
colors <- bpy.colors(8)[-2]
plot(c(1:length(colors)), rep(1, length(colors)), col=colors)
# colors <- bpy.colors(16)[-c(1, 2, 4, 6, 8, 10, 12, 14, 16)]

### load functions
source(paste(path.code, "grVBEM.R", sep=""))

### the settings
n <- 100
p <- 1000
alpha <- 0.5
lambda <- 100
G <- 4
# lambdag <- exp(c(-1.0625, -0.0625, 0.4375, 0.6875))
# lambdag <- exp(log(seq(0.3, 1.8, length.out=G)) -
#                  mean(log(seq(0.3, 1.8, length.out=G))))
lambdag <- exp(seq(-2, 2, length.out=4))
part <- list(groups=rep(c(1:G), each=p/G))
g <- 4
q <- 0.5
csel <- c(seq(1, 10, 2), seq(15, 50, 5), seq(60, 140, 10), seq(160, 200, 20))
nreps <- 50

## the simulation
psel <- kappa <- list(grridge=matrix(NA, nrow=nreps, ncol=length(csel)),
                      enet1=matrix(NA, nrow=nreps, ncol=100),
                      enet2=matrix(NA, nrow=nreps, ncol=100),
                      enet3=matrix(NA, nrow=nreps, ncol=100),
                      greben1=matrix(NA, nrow=nreps, ncol=100),
                      greben2=matrix(NA, nrow=nreps, ncol=100),
                      greben3=matrix(NA, nrow=nreps, ncol=100),
                      grplasso=matrix(NA, nrow=nreps, ncol=100),
                      sglasso1=matrix(NA, nrow=nreps, ncol=100),
                      sglasso2=matrix(NA, nrow=nreps, ncol=100),
                      sglasso3=matrix(NA, nrow=nreps, ncol=100))
auc <- mse <- briers <- list(ridge=numeric(nreps),
                             grridge=matrix(NA, nrow=nreps, ncol=length(csel)),
                             enet1=matrix(NA, nrow=nreps, ncol=100),
                             enet2=matrix(NA, nrow=nreps, ncol=100),
                             enet3=matrix(NA, nrow=nreps, ncol=100),
                             greben1=matrix(NA, nrow=nreps, ncol=100),
                             greben2=matrix(NA, nrow=nreps, ncol=100),
                             greben3=matrix(NA, nrow=nreps, ncol=100),
                             grplasso=matrix(NA, nrow=nreps, ncol=100),
                             sglasso1=matrix(NA, nrow=nreps, ncol=100),
                             sglasso2=matrix(NA, nrow=nreps, ncol=100),
                             sglasso3=matrix(NA, nrow=nreps, ncol=100))
lambdagest <- list(grridge=matrix(NA, nrow=nreps, ncol=4),
                   greben1=matrix(NA, nrow=nreps, ncol=4),
                   greben2=matrix(NA, nrow=nreps, ncol=4),
                   greben3=matrix(NA, nrow=nreps, ncol=4))
for(r in 1:nreps) {
  set.seed(2018 + r - 1)
  print(paste("rep", r))
  # beta <- as.numeric(sapply(1:G, function(g) {
  #   c(rep(0, q*p/G), renbeta(ceiling((1 - q)*p/G), lambda*lambdag[g]*alpha,
  #                                  0.5*(1 - alpha)*lambda*lambdag[g]))}))
  beta <- as.numeric(sapply(1:G, function(g) {
    b <- renbeta(p/G, lambda*lambdag[g]*alpha, 0.5*(1 - alpha)*lambda*lambdag[g]);
    b[abs(b)<=quantile(abs(b), q)] <- 0
    return(b)}))
  pblock <- 25
  rho <- 0.7
  sigma <- matrix(rho, ncol=pblock, nrow=pblock); diag(sigma) <- 1
  x <- do.call(cbind, replicate(p/pblock, rmvnorm(n, mean=rep(0, pblock),
                                                  sigma=sigma), simplify=FALSE))
  prob <- as.numeric(exp(x %*% beta)/(1 + exp(x %*% beta)))
  y <- rbinom(n, 1, prob)

  ntest <- 1000
  xtest <- do.call(cbind, replicate(p/pblock, rmvnorm(ntest, mean=rep(
    0, pblock), sigma=sigma), simplify=FALSE))
  probtest <- as.numeric(exp(xtest %*% beta)/(1 + exp(xtest %*% beta)))
  ytest <- rbinom(ntest, 1, probtest)
  
  # fit1.enet <- glmnet(x, y, alpha=0.05, standardize=FALSE)
  # fit2.enet <- glmnet(x, y, alpha=0.5, standardize=FALSE)
  # fit3.enet <- glmnet(x, y, alpha=0.95, standardize=FALSE)
  # fit.ridge <- cv.glmnet(x, y, alpha=0, standardize=FALSE)
  # fit1.greben <- grEBEN3(x, y, rep(1, length(y)), partitions=part, alpha=0.05,
  #                        trace=FALSE)
  # fit2.greben <- grEBEN3(x, y, rep(1, length(y)), partitions=part, alpha=0.5,
  #                        trace=FALSE)
  # fit3.greben <- grEBEN3(x, y, rep(1, length(y)), partitions=part, alpha=0.95,
  #                        trace=FALSE)
  # fit1.greben2 <- glmnet(x, y, alpha=0.05, standardize=FALSE,
  #                        penalty.factor=rep(fit1.greben$lambdag$groups[
  #                          , fit1.greben$nouteriter + 1], each=p/G))
  # fit2.greben2 <- glmnet(x, y, alpha=0.5, standardize=FALSE,
  #                        penalty.factor=rep(fit2.greben$lambdag$groups[
  #                          , fit2.greben$nouteriter + 1], each=p/G))
  # fit3.greben2 <- glmnet(x, y, alpha=0.95, standardize=FALSE,
  #                        penalty.factor=rep(fit3.greben$lambdag$groups[
  #                          , fit3.greben$nouteriter + 1], each=p/G))
  # 
  # fit.grridge <- vector("list", length(psel))
  # fit.grridge[[1]] <- grridge(t(x), y, partitions=list(groups=CreatePartition(
  #   as.factor(part$groups))), selection=TRUE, maxsel=csel[1], trace=FALSE)
  # for(s in 2:length(csel)) {
  #   fit.grridge[[s]] <- grridge(t(x), y, partitions=list(groups=CreatePartition(
  #     as.factor(part$groups))), selection=TRUE, maxsel=csel[s],
  #     optl=fit.grridge[[1]]$optl, trace=FALSE)
  # }
  # 
  # fit.grplasso <- gglasso(x, (y - 0.5)*2, group=part$groups, loss="logit")
  
  fit1.sglasso <- SGL(list(x=x, y=y), part$groups, type="logit", nlam=100,
                      alpha=0.05)
  fit2.sglasso <- SGL(list(x=x, y=y), part$groups, type="logit", nlam=100, 
                      alpha=0.5)
  fit3.sglasso <- SGL(list(x=x, y=y), part$groups, type="logit", nlam=100, 
                      alpha=0.95)
  
  # auc$ridge[r] <- pROC::roc(ytest, as.numeric(predict(fit.ridge, xtest,
  #                                                     "lambda.min")))$auc
  # auc$grridge[r, ] <- sapply(fit.grridge, function(s) {
  #   pROC::roc(ytest, predict.grridge(s, t(xtest))[, 3])$auc})
  # auc$enet1[r, ] <- apply(predict(fit1.enet, xtest, type="response"), 2,
  #                         function(s) {pROC::roc(ytest, s)$auc})
  # auc$enet2[r, ] <- apply(predict(fit2.enet, xtest, type="response"), 2,
  #                         function(s) {pROC::roc(ytest, s)$auc})
  # auc$enet3[r, ] <- apply(predict(fit3.enet, xtest, type="response"), 2,
  #                         function(s) {pROC::roc(ytest, s)$auc})
  # auc$greben1[r, ] <- apply(predict(fit1.greben2, xtest, type="response"), 2,
  #                         function(s) {pROC::roc(ytest, s)$auc})
  # auc$greben2[r, ] <- apply(predict(fit2.greben2, xtest, type="response"), 2,
  #                           function(s) {pROC::roc(ytest, s)$auc})
  # auc$greben3[r, ] <- apply(predict(fit3.greben2, xtest, type="response"), 2,
  #                           function(s) {pROC::roc(ytest, s)$auc})
  # pred.grplasso <- 1/(1 + exp(-predict(fit.grplasso, xtest, type="link")))
  # auc$grplasso[r, ] <- apply(pred.grplasso, 2,
  #                            function(s) {pROC::roc(ytest, s)$auc})
  pred1.sglasso <- 1/(1 + exp(-xtest %*% fit1.sglasso$beta))
  auc$sglasso1[r, ] <- apply(pred1.sglasso, 2, 
                             function(s) {pROC::roc(ytest, s)$auc})
  pred2.sglasso <- 1/(1 + exp(-xtest %*% fit2.sglasso$beta))
  auc$sglasso2[r, ] <- apply(pred2.sglasso, 2, 
                             function(s) {pROC::roc(ytest, s)$auc})
  pred3.sglasso <- 1/(1 + exp(-xtest %*% fit3.sglasso$beta))
  auc$sglasso3[r, ] <- apply(pred3.sglasso, 2, 
                             function(s) {pROC::roc(ytest, s)$auc})
  
  const <- sum((ytest - mean(ytest))^2)
  # briers$ridge[r] <- 1 - sum((ytest - as.numeric(
  #   predict(fit.ridge, xtest, "lambda.min")))^2)/const
  # briers$grridge[r, ] <- sapply(fit.grridge, function(s) {
  #   1 - sum((ytest - predict.grridge(s, t(xtest))[, 3])^2)/const})
  # briers$enet1[r, ] <- apply(predict(fit1.enet, xtest, type="response"), 2,
  #                            function(pred) {1 - sum((ytest - pred)^2)/const})
  # briers$enet2[r, ] <- apply(predict(fit2.enet, xtest, type="response"), 2,
  #                            function(pred) {1 - sum((ytest - pred)^2)/const})
  # briers$enet3[r, ] <- apply(predict(fit3.enet, xtest, type="response"), 2,
  #                            function(pred) {1 - sum((ytest - pred)^2)/const})
  # briers$greben1[r, ] <- apply(predict(fit1.greben2, xtest, type="response"), 2,
  #                            function(pred) {1 - sum((ytest - pred)^2)/const})
  # briers$greben2[r, ] <- apply(predict(fit2.greben2, xtest, type="response"), 2,
  #                            function(pred) {1 - sum((ytest - pred)^2)/const})
  # briers$greben3[r, ] <- apply(predict(fit3.greben2, xtest, type="response"), 2,
  #                            function(pred) {1 - sum((ytest - pred)^2)/const})
  # briers$grplasso[r, ] <- apply(pred.grplasso, 2, function(pred) {
  #   1 - sum((ytest - pred)^2)/const})
  briers$sglasso1[r, ] <- apply(pred1.sglasso, 2, function(pred) {
    1 - sum((ytest - pred)^2)/const})
  briers$sglasso2[r, ] <- apply(pred2.sglasso, 2, function(pred) {
    1 - sum((ytest - pred)^2)/const})
  briers$sglasso3[r, ] <- apply(pred3.sglasso, 2, function(pred) {
    1 - sum((ytest - pred)^2)/const})

  # mse$ridge[r] <- mean((coef(fit.ridge)[-1] - beta)^2)
  # mse$grridge[r, ] <- sapply(fit.grridge, function(s) {
  #   mean((replace(rep(0, p), s$resEN$whichEN, s$resEN$betasEN) - beta)^2)})
  # mse$enet1[r, ] <- apply(fit1.enet$beta, 2, function(b) {mean((b - beta)^2)})
  # mse$enet2[r, ] <- apply(fit2.enet$beta, 2, function(b) {mean((b - beta)^2)})
  # mse$enet3[r, ] <- apply(fit3.enet$beta, 2, function(b) {mean((b - beta)^2)})
  # mse$greben1[r, ] <- apply(fit1.greben2$beta, 2, function(b) {
  #   mean((b - beta)^2)})
  # mse$greben2[r, ] <- apply(fit2.greben2$beta, 2, function(b) {
  #   mean((b - beta)^2)})
  # mse$greben3[r, ] <- apply(fit3.greben2$beta, 2, function(b) {
  #   mean((b - beta)^2)})
  # mse$grplasso[r, ] <- apply(fit.grplasso$beta, 2, function(b) {
  #   mean((b - beta)^2)})
  mse$sglasso1[r, ] <- apply(fit1.sglasso$beta, 2, function(b) {
    mean((b - beta)^2)})
  mse$sglasso2[r, ] <- apply(fit2.sglasso$beta, 2, function(b) {
    mean((b - beta)^2)})
  mse$sglasso3[r, ] <- apply(fit3.sglasso$beta, 2, function(b) {
    mean((b - beta)^2)})
  
  # psel$grridge[r, ] <- sapply(fit.grridge, function(s) {
  #   length(s$resEN$whichEN)})
  # psel$enet1[r, ] <- fit1.enet$df
  # psel$enet2[r, ] <- fit2.enet$df
  # psel$enet3[r, ] <- fit3.enet$df
  # psel$greben1[r, ] <- fit1.greben2$df
  # psel$greben2[r, ] <- fit2.greben2$df
  # psel$greben3[r, ] <- fit3.greben2$df
  # psel$grplasso[r, ] <- fit.grplasso$df
  psel$sglasso1[r, ] <- apply(fit1.sglasso$beta, 2, function(b) {sum(b!=0)})
  psel$sglasso2[r, ] <- apply(fit2.sglasso$beta, 2, function(b) {sum(b!=0)})
  psel$sglasso3[r, ] <- apply(fit3.sglasso$beta, 2, function(b) {sum(b!=0)})
    
  # kappa$grridge[r, ] <- sapply(fit.grridge, function(s) {
  #   kappa2(cbind(beta!=0, replace(rep(FALSE, p), s$resEN$whichEN,
  #                                 TRUE)))$value})
  # kappa$enet1[r, ] <- apply(fit1.enet$beta, 2, function(b) {
  #   kappa2(cbind(beta!=0, b!=0))$value})
  # kappa$enet2[r, ] <- apply(fit2.enet$beta, 2, function(b) {
  #   kappa2(cbind(beta!=0, b!=0))$value})
  # kappa$enet3[r, ] <- apply(fit3.enet$beta, 2, function(b) {
  #   kappa2(cbind(beta!=0, b!=0))$value})
  # kappa$greben1[r, ] <- apply(fit1.greben2$beta, 2, function(b) {
  #   kappa2(cbind(beta!=0, b!=0))$value})
  # kappa$greben2[r, ] <- apply(fit2.greben2$beta, 2, function(b) {
  #   kappa2(cbind(beta!=0, b!=0))$value})
  # kappa$greben3[r, ] <- apply(fit3.greben2$beta, 2, function(b) {
  #   kappa2(cbind(beta!=0, b!=0))$value})
  # kappa$grplasso[r, ] <- apply(fit.grplasso$beta, 2, function(b) {
  #   kappa2(cbind(beta!=0, b!=0))$value})
  kappa$sglasso1[r, ] <- apply(fit1.sglasso$beta, 2, function(b) {
    kappa2(cbind(beta!=0, b!=0))$value})
  kappa$sglasso2[r, ] <- apply(fit2.sglasso$beta, 2, function(b) {
    kappa2(cbind(beta!=0, b!=0))$value})
  kappa$sglasso3[r, ] <- apply(fit3.sglasso$beta, 2, function(b) {
    kappa2(cbind(beta!=0, b!=0))$value})

  # lambdagest$grridge[r, ] <- fit.grridge[[1]]$lambdamults$groups
  # lambdagest$greben1[r, ] <- fit1.greben$lambdag$groups[
  #   , fit1.greben$nouteriter + 1]
  # lambdagest$greben2[r, ] <- fit2.greben$lambdag$groups[
  #   , fit2.greben$nouteriter + 1]
  # lambdagest$greben3[r, ] <- fit3.greben$lambdag$groups[
  #   , fit3.greben$nouteriter + 1]

}

# save(results1, file=paste(path.res, "grEBEN_sim1_res1.Rdata", sep=""))

### plots
load(paste(path.res, "grEBEN_sim1_res1.Rdata", sep=""))
# results1$psel <- c(results1$psel, sglasso1=list(psel$sglasso1), 
#                    sglasso2=list(psel$sglasso2),
#                    sglasso3=list(psel$sglasso3))
# results1$kappa <- c(results1$kappa, sglasso1=list(kappa$sglasso1), 
#                    sglasso2=list(kappa$sglasso2),
#                    sglasso3=list(kappa$sglasso3))
# results1$briers <- c(results1$briers, sglasso1=list(briers$sglasso1), 
#                    sglasso2=list(briers$sglasso2),
#                    sglasso3=list(briers$sglasso3))
# results1$auc <- c(results1$auc, sglasso1=list(auc$sglasso1), 
#                    sglasso2=list(auc$sglasso2),
#                    sglasso3=list(auc$sglasso3))
# results1$mse <- c(results1$mse, sglasso1=list(mse$sglasso1), 
#                   sglasso2=list(mse$sglasso2),
#                   sglasso3=list(mse$sglasso3))

namauc <- names(results1$auc)
results1$auc <- c(list(results1$auc$ridge),
                  sapply(2:length(results1$auc), function(s) {
                    results1$auc[[s]][results1$psel[[s - 1]]!=0]}))
names(results1$auc) <- namauc

nambriers <- names(results1$briers)
results1$briers <- c(list(results1$briers$ridge),
                  sapply(2:length(results1$briers), function(s) {
                    results1$briers[[s]][results1$psel[[s - 1]]!=0]}))
names(results1$briers) <- nambriers

namkappa <- names(results1$kappa)
results1$kappa <- sapply(1:length(results1$kappa), function(s) {
                    results1$kappa[[s]][results1$psel[[s]]!=0]})
names(results1$kappa) <- namkappa

nammse <- names(results1$mse)
results1$mse <- c(list(results1$mse$ridge),
                  sapply(2:length(results1$mse), function(s) {
                    results1$mse[[s]][results1$psel[[s - 1]]!=0]}))
names(results1$mse) <- nammse

nampsel <- names(results1$psel)
results1$psel <- sapply(results1$psel, function(s) {s[s!=0]})
names(results1$psel) <- nampsel

# combined
leglabels <- c("ridge", expression(paste("enet, ", alpha==0.05)),
               expression(paste("enet, ", alpha==0.5)),
               expression(paste("enet, ", alpha==0.95)),
               expression(paste("sparse group lasso, ", alpha==0.05)),
               expression(paste("sparse group lasso, ", alpha==0.5)),
               expression(paste("sparse group lasso, ", alpha==0.95)),
               "group-regularized", "not group-regularized")

png(paste(path.graph, "grEBEN_sim1_res1_performance3.png", sep=""),
    units="in", width=14, height=10, res=120)
par(mfrow=c(2, 2), mar=c(5.1, 5.1, 4.1, 2.1))
pselm.grridge <- unique(sort(results1$psel$grridge))
pselm.greben1 <- unique(sort(results1$psel$greben1))
pselm.greben2 <- unique(sort(results1$psel$greben2))
pselm.greben3 <- unique(sort(results1$psel$greben3))
pselm.enet1 <- unique(sort(results1$psel$enet1))
pselm.enet2 <- unique(sort(results1$psel$enet2))
pselm.enet3 <- unique(sort(results1$psel$enet3))
pselm.grplasso <- unique(sort(results1$psel$grplasso))
pselm.sglasso1 <- unique(sort(results1$psel$sglasso1))
pselm.sglasso2 <- unique(sort(results1$psel$sglasso2))
pselm.sglasso3 <- unique(sort(results1$psel$sglasso3))
plot(unique(sort(results1$psel$grridge)),
     unique(predict(loess(results1$kappa$grridge[order(results1$psel$grridge)] ~
                            sort(results1$psel$grridge)))), type="l",
     ylim=range(sapply(1:length(results1$kappa), function(s) {
       predict(loess(as.numeric(results1$kappa[[s]]) ~
                       as.numeric(results1$psel[[s]])))})),
     xlim=range(pselm.grridge, pselm.greben1, pselm.greben2, pselm.greben3,
                pselm.enet1, pselm.enet2, pselm.enet3, pselm.sglasso1, 
                pselm.sglasso2, pselm.sglasso3),
     col=colors[1], main="a)", xlab="Number of selected features",
     ylab="Cohen's kappa", lwd=1.5, cex.axis=1.5, cex.lab=2, cex.main=2)
lines(unique(sort(results1$psel$greben1)),
      unique(predict(loess(results1$kappa$greben1[order(results1$psel$greben1)] ~
                             sort(results1$psel$greben1)))), col=colors[2], 
      lwd=1.5)
lines(unique(sort(results1$psel$greben2)),
      unique(predict(loess(results1$kappa$greben2[order(results1$psel$greben2)] ~
                             sort(results1$psel$greben2)))), col=colors[3], 
      lwd=1.5)
lines(unique(sort(results1$psel$greben3)),
      unique(predict(loess(results1$kappa$greben3[order(results1$psel$greben3)] ~
                             sort(results1$psel$greben3)))), col=colors[4], 
      lwd=1.5)
lines(unique(sort(results1$psel$enet1)),
      unique(predict(loess(results1$kappa$enet1[order(results1$psel$enet1)] ~
                             sort(results1$psel$enet1)))), col=colors[2], lty=2, 
      lwd=1.5)
lines(unique(sort(results1$psel$enet2)),
      unique(predict(loess(results1$kappa$enet2[order(results1$psel$enet2)] ~
                             sort(results1$psel$enet2)))), col=colors[3], lty=2, 
      lwd=1.5)
lines(unique(sort(results1$psel$enet3)),
      unique(predict(loess(results1$kappa$enet3[order(results1$psel$enet3)] ~
                             sort(results1$psel$enet3)))), col=colors[4], lty=2, 
      lwd=1.5)
# points(unique(sort(results1$psel$grplasso)), sapply(
#   sort(unique(results1$psel$grplasso)), function(psel) {
#     mean(results1$kappa$grplasso[results1$psel$grplasso==psel])}), 
#   col=colors[5])
lines(unique(sort(results1$psel$sglasso1)),
      unique(predict(loess(results1$kappa$sglasso1[order(results1$psel$sglasso1)] ~
                             sort(results1$psel$sglasso1)))), col=colors[5], 
      lty=1, lwd=1.5)
lines(unique(sort(results1$psel$sglasso2)),
      unique(predict(loess(results1$kappa$sglasso2[order(results1$psel$sglasso2)] ~
                             sort(results1$psel$sglasso2)))), col=colors[6], 
      lty=1, lwd=1.5)
lines(unique(sort(results1$psel$sglasso3)),
      unique(predict(loess(results1$kappa$sglasso3[order(results1$psel$sglasso3)] ~
                             sort(results1$psel$sglasso3)))), col=colors[7], 
      lty=1, lwd=1.5)

pselm.ridge <- c(0, max(unlist(results1$psel)))
pselm.grridge <- unique(sort(results1$psel$grridge))
pselm.greben1 <- unique(sort(results1$psel$greben1))
pselm.greben2 <- unique(sort(results1$psel$greben2))
pselm.greben3 <- unique(sort(results1$psel$greben3))
pselm.enet1 <- unique(sort(results1$psel$enet1))
pselm.enet2 <- unique(sort(results1$psel$enet2))
pselm.enet3 <- unique(sort(results1$psel$enet3))
pselm.grplasso <- unique(sort(results1$psel$grplasso))
pselm.sglasso1 <- unique(sort(results1$psel$sglasso1))
pselm.sglasso2 <- unique(sort(results1$psel$sglasso2))
pselm.sglasso3 <- unique(sort(results1$psel$sglasso3))
plot(unique(sort(results1$psel$grridge)),
     unique(predict(loess(results1$mse$grridge[order(results1$psel$grridge)] ~
                            sort(results1$psel$grridge)))), type="l",
     ylim=range(sapply(2:(length(results1$mse) - 3), function(s) {
       predict(loess(as.numeric(results1$mse[[s]]) ~
                       as.numeric(results1$psel[[s - 1]])))})),
     xlim=range(pselm.ridge, pselm.grridge, pselm.greben1, pselm.greben2,
                pselm.greben3, pselm.enet1, pselm.enet2, pselm.enet3),
     col=colors[1], main="b)", xlab="Number of selected features", ylab="MSE", 
     lwd=1.5, cex.axis=1.5, cex.lab=2, cex.main=2)
lines(range(results1$psel), rep(median(results1$mse$ridge), 2), col=colors[1], 
      lty=2, lwd=1.5)
lines(unique(sort(results1$psel$greben1)),
      unique(predict(loess(results1$mse$greben1[order(results1$psel$greben1)] ~
                             sort(results1$psel$greben1)))), col=colors[2], 
      lwd=1.5)
lines(unique(sort(results1$psel$greben2)),
      unique(predict(loess(results1$mse$greben2[order(results1$psel$greben2)] ~
                             sort(results1$psel$greben2)))), col=colors[3], 
      lwd=1.5)
lines(unique(sort(results1$psel$greben3)),
      unique(predict(loess(results1$mse$greben3[order(results1$psel$greben3)] ~
                             sort(results1$psel$greben3)))), col=colors[4], 
      lwd=1.5)
lines(unique(sort(results1$psel$enet1)),
      unique(predict(loess(results1$mse$enet1[order(results1$psel$enet1)] ~
                             sort(results1$psel$enet1)))), col=colors[2], lty=2, 
      lwd=1.5)
lines(unique(sort(results1$psel$enet2)),
      unique(predict(loess(results1$mse$enet2[order(results1$psel$enet2)] ~
                             sort(results1$psel$enet2)))), col=colors[3], lty=2, 
      lwd=1.5)
lines(unique(sort(results1$psel$enet3)),
      unique(predict(loess(results1$mse$enet3[order(results1$psel$enet3)] ~
                             sort(results1$psel$enet3)))), col=colors[4], lty=2, 
      lwd=1.5)
# points(unique(sort(results1$psel$grplasso)), sapply(
#   sort(unique(results1$psel$grplasso)), function(psel) {
#     mean(results1$mse$grplasso[results1$psel$grplasso==psel])}), col=colors[5])
lines(unique(sort(results1$psel$sglasso1)),
      unique(predict(loess(results1$mse$sglasso1[order(results1$psel$sglasso1)] ~
                             sort(results1$psel$sglasso1)))), col=colors[5], 
      lty=1, lwd=1.5)
lines(unique(sort(results1$psel$sglasso2)),
      unique(predict(loess(results1$mse$sglasso2[order(results1$psel$sglasso2)] ~
                             sort(results1$psel$sglasso2)))), col=colors[6], 
      lty=1, lwd=1.5)
lines(unique(sort(results1$psel$sglasso3)),
      unique(predict(loess(results1$mse$sglasso3[order(results1$psel$sglasso3)] ~
                             sort(results1$psel$sglasso3)))), col=colors[7], 
      lty=1, lwd=1.5)
legend("bottomright", legend=leglabels, fill=c(colors, 0, 0),
       lty=c(rep(NA, 7), 1, 2), lwd=c(rep(NA, 7), 1.5, 1.5), 
       border=c(rep(1, 7), 0, 0), merge=TRUE, seg.len=1, cex=1.3,
       col=c(rep(NA, 7), 1, 1))
plot(unique(sort(results1$psel$grridge)),
     unique(predict(loess(results1$auc$grridge[order(results1$psel$grridge)] ~
                            sort(results1$psel$grridge)))), type="l",
     ylim=range(sapply(2:length(results1$auc), function(s) {
       predict(loess(as.numeric(results1$auc[[s]]) ~
                       as.numeric(results1$psel[[s - 1]])))})),
     xlim=range(pselm.ridge, pselm.grridge, pselm.greben1, pselm.greben2,
                pselm.greben3, pselm.enet1, pselm.enet2, pselm.enet3, 
                pselm.sglasso1, pselm.sglasso2, pselm.sglasso3),
     col=colors[1], main="c)", xlab="Number of selected features", ylab="AUC", 
     lwd=1.5, cex.axis=1.5, cex.lab=2, cex.main=2)
lines(range(results1$psel), rep(mean(results1$auc$ridge), 2), col=colors[1], 
      lty=2, lwd=1.5)
lines(unique(sort(results1$psel$greben1)),
      unique(predict(loess(results1$auc$greben1[order(results1$psel$greben1)] ~
                             sort(results1$psel$greben1)))), col=colors[2], 
      lwd=1.5)
lines(unique(sort(results1$psel$greben2)),
      unique(predict(loess(results1$auc$greben2[order(results1$psel$greben2)] ~
                             sort(results1$psel$greben2)))), col=colors[3], 
      lwd=1.5)
lines(unique(sort(results1$psel$greben3)),
      unique(predict(loess(results1$auc$greben3[order(results1$psel$greben3)] ~
                             sort(results1$psel$greben3)))), col=colors[4], 
      lwd=1.5)
lines(unique(sort(results1$psel$enet1)),
      unique(predict(loess(results1$auc$enet1[order(results1$psel$enet1)] ~
                             sort(results1$psel$enet1)))), col=colors[2], lty=2, 
      lwd=1.5)
lines(unique(sort(results1$psel$enet2)),
      unique(predict(loess(results1$auc$enet2[order(results1$psel$enet2)] ~
                             sort(results1$psel$enet2)))), col=colors[3], lty=2, 
      lwd=1.5)
lines(unique(sort(results1$psel$enet3)),
      unique(predict(loess(results1$auc$enet3[order(results1$psel$enet3)] ~
                             sort(results1$psel$enet3)))), col=colors[4], lty=2, 
      lwd=1.5)
# points(unique(sort(results1$psel$grplasso)), sapply(
#   sort(unique(results1$psel$grplasso)), function(psel) {
#     mean(results1$auc$grplasso[results1$psel$grplasso==psel])}), col=colors[5])
lines(unique(sort(results1$psel$sglasso1)),
      unique(predict(loess(results1$auc$sglasso1[order(results1$psel$sglasso1)] ~
                             sort(results1$psel$sglasso1)))), col=colors[5], 
      lty=1, lwd=1.5)
lines(unique(sort(results1$psel$sglasso2)),
      unique(predict(loess(results1$auc$sglasso2[order(results1$psel$sglasso2)] ~
                             sort(results1$psel$sglasso2)))), col=colors[6], 
      lty=1, lwd=1.5)
lines(unique(sort(results1$psel$sglasso3)),
      unique(predict(loess(results1$auc$sglasso3[order(results1$psel$sglasso3)] ~
                             sort(results1$psel$sglasso3)))), col=colors[7], 
      lty=1, lwd=1.5)
plot(unique(sort(results1$psel$grridge)),
     unique(predict(loess(results1$briers$grridge[order(results1$psel$grridge)] ~
                            sort(results1$psel$grridge)))), type="l",
     ylim=range(sapply(2:length(results1$briers), function(s) {
       predict(loess(as.numeric(results1$briers[[s]]) ~
                       as.numeric(results1$psel[[s - 1]])))})),
     xlim=range(pselm.ridge, pselm.grridge, pselm.greben1, pselm.greben2,
                pselm.greben3, pselm.enet1, pselm.enet2, pselm.enet3, 
                pselm.sglasso1, pselm.sglasso2, pselm.sglasso3),
     col=colors[1], main="d)", xlab="Number of selected features", 
     ylab="Brier skill score", lwd=1.5, cex.axis=1.5, cex.lab=2, cex.main=2)
lines(range(results1$psel), rep(mean(results1$briers$ridge), 2), col=colors[1], 
      lty=2, lwd=1.5)
lines(unique(sort(results1$psel$greben1)),
      unique(predict(loess(results1$briers$greben1[order(results1$psel$greben1)] ~
                             sort(results1$psel$greben1)))), col=colors[2], 
      lwd=1.5)
lines(unique(sort(results1$psel$greben2)),
      unique(predict(loess(results1$briers$greben2[order(results1$psel$greben2)] ~
                             sort(results1$psel$greben2)))), col=colors[3], 
      lwd=1.5)
lines(unique(sort(results1$psel$greben3)),
      unique(predict(loess(results1$briers$greben3[order(results1$psel$greben3)] ~
                             sort(results1$psel$greben3)))), col=colors[4], 
      lwd=1.5)
lines(unique(sort(results1$psel$enet1)),
      unique(predict(loess(results1$briers$enet1[order(results1$psel$enet1)] ~
                             sort(results1$psel$enet1)))), col=colors[2], lty=2, 
      lwd=1.5)
lines(unique(sort(results1$psel$enet2)),
      unique(predict(loess(results1$briers$enet2[order(results1$psel$enet2)] ~
                             sort(results1$psel$enet2)))), col=colors[3], lty=2, 
      lwd=1.5)
lines(unique(sort(results1$psel$enet3)),
      unique(predict(loess(results1$briers$enet3[order(results1$psel$enet3)] ~
                             sort(results1$psel$enet3)))), col=colors[4], lty=2, 
      lwd=1.5)
# points(unique(sort(results1$psel$grplasso)), sapply(
#   sort(unique(results1$psel$grplasso)), function(psel) {
#     mean(results1$briers$grplasso[results1$psel$grplasso==psel])}), 
#   col=colors[5])
lines(unique(sort(results1$psel$sglasso1)),
      unique(predict(loess(results1$briers$sglasso1[order(results1$psel$sglasso1)] ~
                             sort(results1$psel$sglasso1)))), col=colors[5], 
      lty=1, lwd=1.5)
lines(unique(sort(results1$psel$sglasso2)),
      unique(predict(loess(results1$briers$sglasso2[order(results1$psel$sglasso2)] ~
                             sort(results1$psel$sglasso2)))), col=colors[6], 
      lty=1, lwd=1.5)
lines(unique(sort(results1$psel$sglasso3)),
      unique(predict(loess(results1$briers$sglasso3[order(results1$psel$sglasso3)] ~
                             sort(results1$psel$sglasso3)))), col=colors[7], 
      lty=1, lwd=1.5)
dev.off()

# # auc
# aucm.ridge <- rep(mean(results1$auc$ridge), 2)
# aucm.grridge <- tapply(results1$auc$grridge[order(results1$psel$grridge)],
#                        sort(results1$psel$grridge), function(a) {mean(a)})
# aucm.greben1 <- tapply(results1$auc$greben1[order(results1$psel$greben1)],
#                        sort(results1$psel$greben1), function(a) {mean(a)})
# aucm.greben2 <- tapply(results1$auc$greben2[order(results1$psel$greben2)],
#                        sort(results1$psel$greben2), function(a) {mean(a)})
# aucm.greben3 <- tapply(results1$auc$greben3[order(results1$psel$greben3)],
#                        sort(results1$psel$greben3), function(a) {mean(a)})
# aucm.enet1 <- tapply(results1$auc$enet1[order(results1$psel$enet1)],
#                      sort(results1$psel$enet1), function(a) {mean(a)})
# aucm.enet2 <- tapply(results1$auc$enet2[order(results1$psel$enet2)],
#                      sort(results1$psel$enet2), function(a) {mean(a)})
# aucm.enet3 <- tapply(results1$auc$enet3[order(results1$psel$enet3)],
#                      sort(results1$psel$enet3), function(a) {mean(a)})
# 
# pselm.ridge <- c(0, max(unlist(results1$psel)))
# pselm.grridge <- unique(sort(results1$psel$grridge))
# pselm.greben1 <- unique(sort(results1$psel$greben1))
# pselm.greben2 <- unique(sort(results1$psel$greben2))
# pselm.greben3 <- unique(sort(results1$psel$greben3))
# pselm.enet1 <- unique(sort(results1$psel$enet1))
# pselm.enet2 <- unique(sort(results1$psel$enet2))
# pselm.enet3 <- unique(sort(results1$psel$enet3))
# 
# png(paste(path.graph, "grEBEN_sim1_res1_auc.png", sep=""),
#     units="in", width=12, height=6, res=120)
# par(mfrow=c(1, 2))
# plot(pselm.grridge, aucm.grridge, type="l",
#      ylim=range(aucm.ridge, aucm.grridge, aucm.greben1, aucm.greben2,
#                 aucm.greben3, aucm.enet1, aucm.enet2, aucm.enet3),
#      xlim=range(pselm.ridge, pselm.grridge, pselm.greben1, pselm.greben2,
#                 pselm.greben3, pselm.enet1, pselm.enet2, pselm.enet3),
#      col=2, main="a)", xlab="Number of selected variables", ylab="AUC")
# lines(pselm.ridge, aucm.ridge, col=2, lty=2)
# lines(pselm.greben1, aucm.greben1, col=3)
# lines(pselm.greben2, aucm.greben2, col=4)
# lines(pselm.greben3, aucm.greben3, col=5)
# lines(pselm.enet1, aucm.enet1, col=3, lty=2)
# lines(pselm.enet2, aucm.enet2, col=4, lty=2)
# lines(pselm.enet3, aucm.enet3, col=5, lty=2)
# 
# plot(unique(sort(results1$psel$grridge)),
#      unique(predict(loess(results1$auc$grridge[order(results1$psel$grridge)] ~
#                             sort(results1$psel$grridge)))), type="l",
#      ylim=range(sapply(2:length(results1$auc), function(s) {
#        predict(loess(as.numeric(results1$auc[[s]]) ~
#                        as.numeric(results1$psel[[s - 1]])))})),
#      xlim=range(pselm.ridge, pselm.grridge, pselm.greben1, pselm.greben2,
#                 pselm.greben3, pselm.enet1, pselm.enet2, pselm.enet3),
#      col=2, main="b)", xlab="Number of selected variables", ylab="AUC")
# lines(range(results1$psel), rep(mean(results1$auc$ridge), 2), col=2, lty=2)
# lines(unique(sort(results1$psel$greben1)),
#       unique(predict(loess(results1$auc$greben1[order(results1$psel$greben1)] ~
#                              sort(results1$psel$greben1)))), col=3)
# lines(unique(sort(results1$psel$greben2)),
#       unique(predict(loess(results1$auc$greben2[order(results1$psel$greben2)] ~
#                              sort(results1$psel$greben2)))), col=4)
# lines(unique(sort(results1$psel$greben3)),
#       unique(predict(loess(results1$auc$greben3[order(results1$psel$greben3)] ~
#                              sort(results1$psel$greben3)))), col=5)
# lines(unique(sort(results1$psel$enet1)),
#       unique(predict(loess(results1$auc$enet1[order(results1$psel$enet1)] ~
#                              sort(results1$psel$enet1)))), col=3, lty=2)
# lines(unique(sort(results1$psel$enet2)),
#       unique(predict(loess(results1$auc$enet2[order(results1$psel$enet2)] ~
#                              sort(results1$psel$enet2)))), col=4, lty=2)
# lines(unique(sort(results1$psel$enet3)),
#       unique(predict(loess(results1$auc$enet3[order(results1$psel$enet3)] ~
#                              sort(results1$psel$enet3)))), col=5, lty=2)
# legend("bottomright", lty=c(rep(NA, 4), 1, 2), border=c(rep(1, 4), 0, 0),
#        merge=TRUE, seg.len=1, fill=c(2:5, 0, 0),
#        legend=c("ridge", expression(paste("enet, ", alpha==0.05)),
#                 expression(paste("enet, ", alpha==0.5)),
#                 expression(paste("enet, ", alpha==0.95)),
#                 "group-regularized", "not group-regularized"))
# dev.off()
# 
# # Brier skill score
# briersm.ridge <- rep(mean(results1$briers$ridge), 2)
# briersm.grridge <- tapply(results1$briers$grridge[order(results1$psel$grridge)],
#                        sort(results1$psel$grridge), function(a) {mean(a)})
# briersm.greben1 <- tapply(results1$briers$greben1[order(results1$psel$greben1)],
#                        sort(results1$psel$greben1), function(a) {mean(a)})
# briersm.greben2 <- tapply(results1$briers$greben2[order(results1$psel$greben2)],
#                        sort(results1$psel$greben2), function(a) {mean(a)})
# briersm.greben3 <- tapply(results1$briers$greben3[order(results1$psel$greben3)],
#                        sort(results1$psel$greben3), function(a) {mean(a)})
# briersm.enet1 <- tapply(results1$briers$enet1[order(results1$psel$enet1)],
#                      sort(results1$psel$enet1), function(a) {mean(a)})
# briersm.enet2 <- tapply(results1$briers$enet2[order(results1$psel$enet2)],
#                      sort(results1$psel$enet2), function(a) {mean(a)})
# briersm.enet3 <- tapply(results1$briers$enet3[order(results1$psel$enet3)],
#                      sort(results1$psel$enet3), function(a) {mean(a)})
# 
# pselm.ridge <- c(0, max(unlist(results1$psel)))
# pselm.grridge <- unique(sort(results1$psel$grridge))
# pselm.greben1 <- unique(sort(results1$psel$greben1))
# pselm.greben2 <- unique(sort(results1$psel$greben2))
# pselm.greben3 <- unique(sort(results1$psel$greben3))
# pselm.enet1 <- unique(sort(results1$psel$enet1))
# pselm.enet2 <- unique(sort(results1$psel$enet2))
# pselm.enet3 <- unique(sort(results1$psel$enet3))
# 
# png(paste(path.graph, "grEBEN_sim1_res1_briers.png", sep=""),
#     units="in", width=12, height=6, res=120)
# par(mfrow=c(1, 2))
# plot(pselm.grridge, briersm.grridge, type="l",
#      ylim=range(briersm.ridge, briersm.grridge, briersm.greben1, briersm.greben2,
#                 briersm.greben3, briersm.enet1, briersm.enet2, briersm.enet3),
#      xlim=range(pselm.ridge, pselm.grridge, pselm.greben1, pselm.greben2,
#                 pselm.greben3, pselm.enet1, pselm.enet2, pselm.enet3),
#      col=2, main="a)", xlab="Number of selected variables", 
#      ylab="Brier skill score")
# lines(pselm.ridge, briersm.ridge, col=2, lty=2)
# lines(pselm.greben1, briersm.greben1, col=3)
# lines(pselm.greben2, briersm.greben2, col=4)
# lines(pselm.greben3, briersm.greben3, col=5)
# lines(pselm.enet1, briersm.enet1, col=3, lty=2)
# lines(pselm.enet2, briersm.enet2, col=4, lty=2)
# lines(pselm.enet3, briersm.enet3, col=5, lty=2)
# 
# plot(unique(sort(results1$psel$grridge)),
#      unique(predict(loess(results1$briers$grridge[order(results1$psel$grridge)] ~
#                             sort(results1$psel$grridge)))), type="l",
#      ylim=range(sapply(2:length(results1$briers), function(s) {
#        predict(loess(as.numeric(results1$briers[[s]]) ~
#                        as.numeric(results1$psel[[s - 1]])))})),
#      xlim=range(pselm.ridge, pselm.grridge, pselm.greben1, pselm.greben2,
#                 pselm.greben3, pselm.enet1, pselm.enet2, pselm.enet3),
#      col=2, main="b)", xlab="Number of selected variables", 
#      ylab="Brier skill score")
# lines(range(results1$psel), rep(mean(results1$briers$ridge), 2), col=2, lty=2)
# lines(unique(sort(results1$psel$greben1)),
#       unique(predict(loess(results1$briers$greben1[order(results1$psel$greben1)] ~
#                              sort(results1$psel$greben1)))), col=3)
# lines(unique(sort(results1$psel$greben2)),
#       unique(predict(loess(results1$briers$greben2[order(results1$psel$greben2)] ~
#                              sort(results1$psel$greben2)))), col=4)
# lines(unique(sort(results1$psel$greben3)),
#       unique(predict(loess(results1$briers$greben3[order(results1$psel$greben3)] ~
#                              sort(results1$psel$greben3)))), col=5)
# lines(unique(sort(results1$psel$enet1)),
#       unique(predict(loess(results1$briers$enet1[order(results1$psel$enet1)] ~
#                              sort(results1$psel$enet1)))), col=3, lty=2)
# lines(unique(sort(results1$psel$enet2)),
#       unique(predict(loess(results1$briers$enet2[order(results1$psel$enet2)] ~
#                              sort(results1$psel$enet2)))), col=4, lty=2)
# lines(unique(sort(results1$psel$enet3)),
#       unique(predict(loess(results1$briers$enet3[order(results1$psel$enet3)] ~
#                              sort(results1$psel$enet3)))), col=5, lty=2)
# legend("bottomright", lty=c(rep(NA, 4), 1, 2), border=c(rep(1, 4), 0, 0),
#        merge=TRUE, seg.len=1, fill=c(2:5, 0, 0),
#        legend=c("ridge", expression(paste("enet, ", alpha==0.05)),
#                 expression(paste("enet, ", alpha==0.5)),
#                 expression(paste("enet, ", alpha==0.95)),
#                 "group-regularized", "not group-regularized"))
# par(mfrow=c(1, 1))
# dev.off()
# 
# 
# # mse
# msem.ridge <- rep(mean(results1$mse$ridge), 2)
# msem.grridge <- tapply(results1$mse$grridge[order(results1$psel$grridge)],
#                        sort(results1$psel$grridge), function(a) {mean(a)})
# msem.greben1 <- tapply(results1$mse$greben1[order(results1$psel$greben1)],
#                        sort(results1$psel$greben1), function(a) {mean(a)})
# msem.greben2 <- tapply(results1$mse$greben2[order(results1$psel$greben2)],
#                        sort(results1$psel$greben2), function(a) {mean(a)})
# msem.greben3 <- tapply(results1$mse$greben3[order(results1$psel$greben3)],
#                        sort(results1$psel$greben3), function(a) {mean(a)})
# msem.enet1 <- tapply(results1$mse$enet1[order(results1$psel$enet1)],
#                      sort(results1$psel$enet1), function(a) {mean(a)})
# msem.enet2 <- tapply(results1$mse$enet2[order(results1$psel$enet2)],
#                      sort(results1$psel$enet2), function(a) {mean(a)})
# msem.enet3 <- tapply(results1$mse$enet3[order(results1$psel$enet3)],
#                      sort(results1$psel$enet3), function(a) {mean(a)})
# 
# pselm.ridge <- c(0, max(unlist(results1$psel)))
# pselm.grridge <- unique(sort(results1$psel$grridge))
# pselm.greben1 <- unique(sort(results1$psel$greben1))
# pselm.greben2 <- unique(sort(results1$psel$greben2))
# pselm.greben3 <- unique(sort(results1$psel$greben3))
# pselm.enet1 <- unique(sort(results1$psel$enet1))
# pselm.enet2 <- unique(sort(results1$psel$enet2))
# pselm.enet3 <- unique(sort(results1$psel$enet3))
# 
# png(paste(path.graph, "grEBEN_sim1_res1_mse.png", sep=""),
#     units="in", width=12, height=6, res=120)
# par(mfrow=c(1, 2))
# plot(pselm.grridge, msem.grridge, type="l",
#      ylim=range(msem.ridge, msem.grridge, msem.greben1, msem.greben2,
#                 msem.greben3, msem.enet1, msem.enet2, msem.enet3),
#      xlim=range(pselm.ridge, pselm.grridge, pselm.greben1, pselm.greben2,
#                 pselm.greben3, pselm.enet1, pselm.enet2, pselm.enet3),
#      col=2, main="a)", xlab="Number of selected variables", ylab="MSE")
# lines(pselm.ridge, msem.ridge, col=2, lty=2)
# lines(pselm.greben1, msem.greben1, col=3)
# lines(pselm.greben2, msem.greben2, col=4)
# lines(pselm.greben3, msem.greben3, col=5)
# lines(pselm.enet1, msem.enet1, col=3, lty=2)
# lines(pselm.enet2, msem.enet2, col=4, lty=2)
# lines(pselm.enet3, msem.enet3, col=5, lty=2)
# 
# plot(unique(sort(results1$psel$grridge)),
#      unique(predict(loess(results1$mse$grridge[order(results1$psel$grridge)] ~
#                             sort(results1$psel$grridge)))), type="l",
#      ylim=range(sapply(2:length(results1$mse), function(s) {
#        predict(loess(as.numeric(results1$mse[[s]]) ~
#                        as.numeric(results1$psel[[s - 1]])))})),
#      xlim=range(pselm.ridge, pselm.grridge, pselm.greben1, pselm.greben2,
#                 pselm.greben3, pselm.enet1, pselm.enet2, pselm.enet3),
#      col=2, main="b)", xlab="Number of selected variables", ylab="MSE")
# lines(range(results1$psel), rep(median(results1$mse$ridge), 2), col=2, lty=2)
# lines(unique(sort(results1$psel$greben1)),
#       unique(predict(loess(results1$mse$greben1[order(results1$psel$greben1)] ~
#                              sort(results1$psel$greben1)))), col=3)
# lines(unique(sort(results1$psel$greben2)),
#       unique(predict(loess(results1$mse$greben2[order(results1$psel$greben2)] ~
#                              sort(results1$psel$greben2)))), col=4)
# lines(unique(sort(results1$psel$greben3)),
#       unique(predict(loess(results1$mse$greben3[order(results1$psel$greben3)] ~
#                              sort(results1$psel$greben3)))), col=5)
# lines(unique(sort(results1$psel$enet1)),
#       unique(predict(loess(results1$mse$enet1[order(results1$psel$enet1)] ~
#                              sort(results1$psel$enet1)))), col=3, lty=2)
# lines(unique(sort(results1$psel$enet2)),
#       unique(predict(loess(results1$mse$enet2[order(results1$psel$enet2)] ~
#                              sort(results1$psel$enet2)))), col=4, lty=2)
# lines(unique(sort(results1$psel$enet3)),
#       unique(predict(loess(results1$mse$enet3[order(results1$psel$enet3)] ~
#                              sort(results1$psel$enet3)))), col=5, lty=2)
# legend("bottomleft", lty=c(rep(NA, 4), 1, 2), border=c(rep(1, 4), 0, 0),
#        merge=TRUE, seg.len=1, fill=c(2:5, 0, 0),
#        legend=c("ridge", expression(paste("enet, ", alpha==0.05)),
#                 expression(paste("enet, ", alpha==0.5)),
#                 expression(paste("enet, ", alpha==0.95)),
#                 "group-regularized", "not group-regularized"))
# par(mfrow=c(1, 1))
# dev.off()
# 
# 
# # kappa
# kappam.grridge <- tapply(results1$kappa$grridge[order(results1$psel$grridge)],
#                        sort(results1$psel$grridge), function(a) {mean(a)})
# kappam.greben1 <- tapply(results1$kappa$greben1[order(results1$psel$greben1)],
#                        sort(results1$psel$greben1), function(a) {mean(a)})
# kappam.greben2 <- tapply(results1$kappa$greben2[order(results1$psel$greben2)],
#                        sort(results1$psel$greben2), function(a) {mean(a)})
# kappam.greben3 <- tapply(results1$kappa$greben3[order(results1$psel$greben3)],
#                        sort(results1$psel$greben3), function(a) {mean(a)})
# kappam.enet1 <- tapply(results1$kappa$enet1[order(results1$psel$enet1)],
#                      sort(results1$psel$enet1), function(a) {mean(a)})
# kappam.enet2 <- tapply(results1$kappa$enet2[order(results1$psel$enet2)],
#                      sort(results1$psel$enet2), function(a) {mean(a)})
# kappam.enet3 <- tapply(results1$kappa$enet3[order(results1$psel$enet3)],
#                      sort(results1$psel$enet3), function(a) {mean(a)})
# 
# pselm.grridge <- unique(sort(results1$psel$grridge))
# pselm.greben1 <- unique(sort(results1$psel$greben1))
# pselm.greben2 <- unique(sort(results1$psel$greben2))
# pselm.greben3 <- unique(sort(results1$psel$greben3))
# pselm.enet1 <- unique(sort(results1$psel$enet1))
# pselm.enet2 <- unique(sort(results1$psel$enet2))
# pselm.enet3 <- unique(sort(results1$psel$enet3))
# 
# png(paste(path.graph, "grEBEN_sim1_res1_kappa.png", sep=""),
#     units="in", width=12, height=6, res=120)
# par(mfrow=c(1, 2))
# plot(pselm.grridge, kappam.grridge, type="l",
#      ylim=range(kappam.grridge, kappam.greben1, kappam.greben2, kappam.greben3,
#                 kappam.enet1, kappam.enet2, kappam.enet3),
#      xlim=range(pselm.grridge, pselm.greben1, pselm.greben2, pselm.greben3,
#                 pselm.enet1, pselm.enet2, pselm.enet3),
#      col=2, main="a)", xlab="Number of selected variables",
#      ylab="Cohen's kappa")
# lines(pselm.greben1, kappam.greben1, col=3)
# lines(pselm.greben2, kappam.greben2, col=4)
# lines(pselm.greben3, kappam.greben3, col=5)
# lines(pselm.enet1, kappam.enet1, col=3, lty=2)
# lines(pselm.enet2, kappam.enet2, col=4, lty=2)
# lines(pselm.enet3, kappam.enet3, col=5, lty=2)
# 
# plot(unique(sort(results1$psel$grridge)),
#      unique(predict(loess(results1$kappa$grridge[order(results1$psel$grridge)] ~
#                             sort(results1$psel$grridge)))), type="l",
#      ylim=range(sapply(1:length(results1$kappa), function(s) {
#        predict(loess(as.numeric(results1$kappa[[s]]) ~
#                        as.numeric(results1$psel[[s]])))})),
#      xlim=range(pselm.grridge, pselm.greben1, pselm.greben2, pselm.greben3,
#                 pselm.enet1, pselm.enet2, pselm.enet3),
#      col=2, main="b)", xlab="Number of selected variables",
#      ylab="Cohen's kappa")
# lines(unique(sort(results1$psel$greben1)),
#       unique(predict(loess(results1$kappa$greben1[order(results1$psel$greben1)] ~
#                              sort(results1$psel$greben1)))), col=3)
# lines(unique(sort(results1$psel$greben2)),
#       unique(predict(loess(results1$kappa$greben2[order(results1$psel$greben2)] ~
#                              sort(results1$psel$greben2)))), col=4)
# lines(unique(sort(results1$psel$greben3)),
#       unique(predict(loess(results1$kappa$greben3[order(results1$psel$greben3)] ~
#                              sort(results1$psel$greben3)))), col=5)
# lines(unique(sort(results1$psel$enet1)),
#       unique(predict(loess(results1$kappa$enet1[order(results1$psel$enet1)] ~
#                              sort(results1$psel$enet1)))), col=3, lty=2)
# lines(unique(sort(results1$psel$enet2)),
#       unique(predict(loess(results1$kappa$enet2[order(results1$psel$enet2)] ~
#                              sort(results1$psel$enet2)))), col=4, lty=2)
# lines(unique(sort(results1$psel$enet3)),
#       unique(predict(loess(results1$kappa$enet3[order(results1$psel$enet3)] ~
#                              sort(results1$psel$enet3)))), col=5, lty=2)
# legend("bottomright", lty=c(rep(NA, 4), 1, 2), border=c(rep(1, 4), 0, 0),
#        merge=TRUE, seg.len=1, fill=c(2:5, 0, 0),
#        legend=c("ridge", expression(paste("enet, ", alpha==0.05)),
#                 expression(paste("enet, ", alpha==0.5)),
#                 expression(paste("enet, ", alpha==0.95)),
#                 "group-regularized", "not group-regularized"))
# par(mfrow=c(1, 1))
# dev.off()


# penalties
load(paste(path.res, "grEBEN_sim1_res1.Rdata", sep=""))

boxdata <- data.frame(mults=c(as.numeric(results1$lambdag[[1]]), 
                              as.numeric(results1$lambdag[[2]]), 
                              as.numeric(results1$lambdag[[3]]), 
                              as.numeric(results1$lambdag[[4]])),
                      group=rep(rep(c(1:4), each=50), 4),
                      method=rep(c(1:4), each=200))

leglabels <- c("GRridge", expression(paste("gren, ", alpha==0.05)),
               expression(paste("gren, ", alpha==0.5)),
               expression(paste("gren, ", alpha==0.95)))

png(paste(path.graph, "grEBEN_sim1_res1_penalties.png", sep=""),
    units="in", width=12, height=6, res=120)
par(mar=c(5.1, 5.6, 4.1, 2.1))
boxplot(mults ~ interaction(method, group), data=boxdata, 
        at=c(c(1:4), c(6:9), c(11:14), c(16:19)), xaxt="n", col=rep(colors, 4), 
        outline=FALSE, cex.lab=2, cex.names=1.5, 
        ylab=expression({lambda^{"'"}}[g]), 
        cex.axis=1.5, boxlwd=0.5)#,  border=rep(colors, 3))
axis(1, c(2.5, 7.5, 12.5, 17.5), c("Group 1", "Group 2", "Group 3", "Group 4"), 
     tick=FALSE, cex.axis=1.5)
abline(h=1, lty=2, lwd=1.5)
legend("topleft", legend=leglabels, fill=colors, border=rep(1, 4),
       seg.len=1, cex=1.3)
dev.off()

