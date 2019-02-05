path.res <- ifelse(as.character(Sys.info()[1])!="Darwin", "~/EBEN/results/",
                   "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/results/")

library(mvtnorm)
library(gren)
library(pROC)
library(GRridge)
library(grpreg)
library(SGL)
library(irr)

##################### scenario 4
n <- 100
p <- 1000
ntest <- 1000

alpha <- 0.5
lambda <- 100
G <- 4
lambdag <- exp(seq(-2, 2, length.out=4))
q <- 0.5
pblock <- 25
rho <- 0.7
Sigma <- matrix(rho, ncol=pblock, nrow=pblock); diag(Sigma) <- 1

part1 <- rep(c(1:G), each=p/G)
csel <- c(seq(1, 10, 2), seq(15, 50, 5), seq(60, 140, 10), seq(160, 200, 20))
nreps <- 100
psel4 <- kappa4 <- list(grridge=matrix(NA, nrow=nreps, ncol=length(csel)),
                        enet1=matrix(NA, nrow=nreps, ncol=100),
                        enet2=matrix(NA, nrow=nreps, ncol=100),
                        enet3=matrix(NA, nrow=nreps, ncol=100),
                        gren1=matrix(NA, nrow=nreps, ncol=100),
                        gren2=matrix(NA, nrow=nreps, ncol=100),
                        gren3=matrix(NA, nrow=nreps, ncol=100),
                        sglasso1=matrix(NA, nrow=nreps, ncol=100),
                        sglasso2=matrix(NA, nrow=nreps, ncol=100),
                        sglasso3=matrix(NA, nrow=nreps, ncol=100),
                        cmcp1=matrix(NA, nrow=nreps, ncol=100),
                        cmcp2=matrix(NA, nrow=nreps, ncol=100),
                        cmcp3=matrix(NA, nrow=nreps, ncol=100),
                        gelasso1=matrix(NA, nrow=nreps, ncol=100),
                        gelasso2=matrix(NA, nrow=nreps, ncol=100),
                        gelasso3=matrix(NA, nrow=nreps, ncol=100))
auc4 <- mse4 <- briers4 <- list(ridge=numeric(nreps),
                                grridge=matrix(NA, nrow=nreps, ncol=length(csel)),
                                enet1=matrix(NA, nrow=nreps, ncol=100),
                                enet2=matrix(NA, nrow=nreps, ncol=100),
                                enet3=matrix(NA, nrow=nreps, ncol=100),
                                gren1=matrix(NA, nrow=nreps, ncol=100),
                                gren2=matrix(NA, nrow=nreps, ncol=100),
                                gren3=matrix(NA, nrow=nreps, ncol=100),
                                sglasso1=matrix(NA, nrow=nreps, ncol=100),
                                sglasso2=matrix(NA, nrow=nreps, ncol=100),
                                sglasso3=matrix(NA, nrow=nreps, ncol=100),
                                cmcp1=matrix(NA, nrow=nreps, ncol=100),
                                cmcp2=matrix(NA, nrow=nreps, ncol=100),
                                cmcp3=matrix(NA, nrow=nreps, ncol=100),
                                gelasso1=matrix(NA, nrow=nreps, ncol=100),
                                gelasso2=matrix(NA, nrow=nreps, ncol=100),
                                gelasso3=matrix(NA, nrow=nreps, ncol=100))
lambdagest4 <- list(grridge=matrix(NA, nrow=nreps, ncol=G),
                    gren1=matrix(NA, nrow=nreps, ncol=G),
                    gren2=matrix(NA, nrow=nreps, ncol=G),
                    gren3=matrix(NA, nrow=nreps, ncol=G))
for(r in 1:nreps) {
  set.seed(2018 + r)
  print(paste("rep", r))
  beta <- as.numeric(sapply(1:G, function(g) {
    b <- renet(p/G, lambda*lambdag[g]*alpha, 0.5*(1 - alpha)*lambda*lambdag[g]);
    b[abs(b)<=quantile(abs(b), q)] <- 0
    return(b)}))
  
  x <- do.call(cbind, replicate(p/pblock, rmvnorm(n, mean=rep(0, pblock),
                                                  sigma=Sigma), simplify=FALSE))
  y <- rbinom(n, 1, as.numeric(exp(x %*% beta)/(1 + exp(x %*% beta))))
  
  
  xtest <- do.call(cbind, replicate(p/pblock, rmvnorm(
    ntest, mean=rep(0, pblock), sigma=Sigma), simplify=FALSE))
  ytest <- rbinom(ntest, 1, as.numeric(exp(xtest %*% beta)/
                                         (1 + exp(xtest %*% beta))))
  
  fit4.ridge <- cv.glmnet(x, y, alpha=0, standardize=FALSE)
  
  fit4.grridge <- vector("list", length(csel))
  invisible(capture.output(
    fit4.grridge[[1]] <- grridge(t(x), y, partitions=list(
      groups=CreatePartition(as.factor(part1))), selection=TRUE, maxsel=csel[1],
      trace=FALSE, standardizeX=FALSE)))
  for(s in 2:length(csel)) {
    invisible(capture.output(
      fit4.grridge[[s]] <- grridge(t(x), y, partitions=list(
        groups=CreatePartition(as.factor(part1))), selection=TRUE,
        maxsel=csel[s], optl=fit4.grridge[[1]]$optl, trace=FALSE,
        standardizeX=FALSE)))
  }
  
  fit4.gren1 <- gren(x, y, partitions=list(groups=part1), alpha=0.05,
                     trace=FALSE)
  fit4.gren2 <- gren(x, y, partitions=list(groups=part1), alpha=0.5,
                     trace=FALSE)
  fit4.gren3 <- gren(x, y, partitions=list(groups=part1), alpha=0.95,
                     trace=FALSE)
  
  fit4.sglasso1 <- SGL(list(x=x, y=y), part1, type="logit", alpha=0.05,
                       standardize=FALSE, nlam=100)
  fit4.sglasso2 <- SGL(list(x=x, y=y), part1, type="logit", alpha=0.5,
                       standardize=FALSE, nlam=100)
  fit4.sglasso3 <- SGL(list(x=x, y=y), part1, type="logit", alpha=0.95,
                       standardize=FALSE, nlam=100)
  
  fit4.cmcp1 <- grpreg(x, y, part1, penalty="cMCP", alpha=0.05)
  fit4.cmcp2 <- grpreg(x, y, part1, penalty="cMCP", alpha=0.5)
  fit4.cmcp3 <- grpreg(x, y, part1, penalty="cMCP", alpha=0.95)
  
  fit4.gelasso1 <- grpreg(x, y, part1, penalty="gel", alpha=0.05)
  fit4.gelasso2 <- grpreg(x, y, part1, penalty="gel", alpha=0.5)
  fit4.gelasso3 <- grpreg(x, y, part1, penalty="gel", alpha=0.95)
  
  pred4.ridge <- as.numeric(predict(fit4.ridge, xtest, "lambda.min"))
  
  pred4.grridge <- sapply(fit4.grridge, function(s) {
    predict.grridge(s, t(xtest))[, 3]})
  
  pred4.gren1 <- predict(fit4.gren1, xtest, type="groupreg",
                         s=fit4.gren1$freq.model$groupreg$lambda)
  pred4.gren2 <- predict(fit4.gren2, xtest, type="groupreg",
                         s=fit4.gren2$freq.model$groupreg$lambda)
  pred4.gren3 <- predict(fit4.gren3, xtest, type="groupreg",
                         s=fit4.gren3$freq.model$groupreg$lambda)
  
  pred4.enet1 <- predict(fit4.gren1, xtest, type="regular",
                         s=fit4.gren1$freq.model$regular$lambda)
  pred4.enet2 <- predict(fit4.gren2, xtest, type="regular",
                         s=fit4.gren2$freq.model$regular$lambda)
  pred4.enet3 <- predict(fit4.gren3, xtest, type="regular",
                         s=fit4.gren3$freq.model$regular$lambda)
  
  pred4.sglasso1 <- predictSGL(fit4.sglasso1, xtest)
  pred4.sglasso2 <- predictSGL(fit4.sglasso2, xtest)
  pred4.sglasso3 <- predictSGL(fit4.sglasso3, xtest)
  
  pred4.cmcp1 <- predict(fit4.cmcp1, xtest)
  pred4.cmcp2 <- predict(fit4.cmcp2, xtest)
  pred4.cmcp3 <- predict(fit4.cmcp3, xtest)
  
  pred4.gelasso1 <- predict(fit4.gelasso1, xtest)
  pred4.gelasso2 <- predict(fit4.gelasso2, xtest)
  pred4.gelasso3 <- predict(fit4.gelasso3, xtest)
  
  auc4$ridge[r] <- pROC::auc(ytest, pred4.ridge)
  
  auc4$grridge[r, ] <- apply(pred4.grridge, 2, function(pred) {
    pROC::auc(ytest, pred)})
  
  auc4$gren1[r, ] <- apply(pred4.gren1, 2, function(pred) {
    pROC::auc(ytest, pred)})
  auc4$gren2[r, ] <- apply(pred4.gren2, 2, function(pred) {
    pROC::auc(ytest, pred)})
  auc4$gren3[r, ] <- apply(pred4.gren3, 2, function(pred) {
    pROC::auc(ytest, pred)})
  
  auc4$enet1[r, ] <- apply(pred4.enet1, 2, function(pred) {
    pROC::auc(ytest, pred)})
  auc4$enet2[r, ] <- apply(pred4.enet2, 2, function(pred) {
    pROC::auc(ytest, pred)})
  auc4$enet3[r, ] <- apply(pred4.enet3, 2, function(pred) {
    pROC::auc(ytest, pred)})
  
  auc4$sglasso1[r, ] <- apply(pred4.sglasso1, 2, function(pred) {
    pROC::auc(ytest, pred)})
  auc4$sglasso2[r, ] <- apply(pred4.sglasso2, 2, function(pred) {
    pROC::auc(ytest, pred)})
  auc4$sglasso3[r, ] <- apply(pred4.sglasso3, 2, function(pred) {
    pROC::auc(ytest, pred)})
  
  auc4$cmcp1[r, ] <- apply(pred4.cmcp1, 2, function(pred) {
    pROC::auc(ytest, pred)})
  auc4$cmcp2[r, ] <- apply(pred4.cmcp2, 2, function(pred) {
    pROC::auc(ytest, pred)})
  auc4$cmcp3[r, ] <- apply(pred4.cmcp3, 2, function(pred) {
    pROC::auc(ytest, pred)})
  
  auc4$gelasso1[r, ] <- apply(pred4.gelasso1, 2, function(pred) {
    pROC::auc(ytest, pred)})
  auc4$gelasso2[r, ] <- apply(pred4.gelasso2, 2, function(pred) {
    pROC::auc(ytest, pred)})
  auc4$gelasso3[r, ] <- apply(pred4.gelasso3, 2, function(pred) {
    pROC::auc(ytest, pred)})
  
  const <- sum((ytest - mean(ytest))^2)
  briers4$ridge[r] <- 1 - sum((ytest - pred4.ridge)^2)/const
  
  briers4$grridge[r, ] <- apply(pred4.grridge, 2, function(pred) {
    1 - sum((ytest - pred)^2)/const})
  
  briers4$gren1[r, ] <- apply(pred4.gren1, 2, function(pred) {
    1 - sum((ytest - pred)^2)/const})
  briers4$gren2[r, ] <- apply(pred4.gren2, 2, function(pred) {
    1 - sum((ytest - pred)^2)/const})
  briers4$gren3[r, ] <- apply(pred4.gren3, 2, function(pred) {
    1 - sum((ytest - pred)^2)/const})
  
  briers4$enet1[r, ] <- apply(pred4.enet1, 2, function(pred) {
    1 - sum((ytest - pred)^2)/const})
  briers4$enet2[r, ] <- apply(pred4.enet2, 2, function(pred) {
    1 - sum((ytest - pred)^2)/const})
  briers4$enet3[r, ] <- apply(pred4.enet3, 2, function(pred) {
    1 - sum((ytest - pred)^2)/const})
  
  briers4$sglasso1[r, ] <- apply(pred4.sglasso1, 2, function(pred) {
    1 - sum((ytest - pred)^2)/const})
  briers4$sglasso2[r, ] <- apply(pred4.sglasso2, 2, function(pred) {
    1 - sum((ytest - pred)^2)/const})
  briers4$sglasso3[r, ] <- apply(pred4.sglasso3, 2, function(pred) {
    1 - sum((ytest - pred)^2)/const})
  
  briers4$cmcp1[r, ] <- apply(pred4.cmcp1, 2, function(pred) {
    1 - sum((ytest - pred)^2)/const})
  briers4$cmcp2[r, ] <- apply(pred4.cmcp2, 2, function(pred) {
    1 - sum((ytest - pred)^2)/const})
  briers4$cmcp3[r, ] <- apply(pred4.cmcp3, 2, function(pred) {
    1 - sum((ytest - pred)^2)/const})
  
  briers4$gelasso1[r, ] <- apply(pred4.gelasso1, 2, function(pred) {
    1 - sum((ytest - pred)^2)/const})
  briers4$gelasso2[r, ] <- apply(pred4.gelasso2, 2, function(pred) {
    1 - sum((ytest - pred)^2)/const})
  briers4$gelasso3[r, ] <- apply(pred4.gelasso3, 2, function(pred) {
    1 - sum((ytest - pred)^2)/const})
  
  mse4$ridge[r] <- mean((coef(fit4.ridge)[-1] - beta)^2)
  
  mse4$grridge[r, ] <- sapply(fit4.grridge, function(s) {
    mean((replace(rep(0, p), s$resEN$whichEN, s$resEN$betasEN) - beta)^2)})
  
  mse4$gren1[r, ] <- apply(fit4.gren1$freq.model$groupreg$beta, 2, function(b) {
    mean((b - beta)^2)})
  mse4$gren2[r, ] <- apply(fit4.gren2$freq.model$groupreg$beta, 2, function(b) {
    mean((b - beta)^2)})
  mse4$gren3[r, ] <- apply(fit4.gren3$freq.model$groupreg$beta, 2, function(b) {
    mean((b - beta)^2)})
  
  mse4$enet1[r, ] <- apply(fit4.gren1$freq.model$regular$beta, 2, function(b) {
    mean((b - beta)^2)})
  mse4$enet2[r, ] <- apply(fit4.gren2$freq.model$regular$beta, 2, function(b) {
    mean((b - beta)^2)})
  mse4$enet3[r, ] <- apply(fit4.gren3$freq.model$regular$beta, 2, function(b) {
    mean((b - beta)^2)})
  
  mse4$sglasso1[r, ] <- apply(fit4.sglasso1$beta, 2, function(b) {
    mean((b - beta)^2)})
  mse4$sglasso2[r, ] <- apply(fit4.sglasso2$beta, 2, function(b) {
    mean((b - beta)^2)})
  mse4$sglasso3[r, ] <- apply(fit4.sglasso3$beta, 2, function(b) {
    mean((b - beta)^2)})
  
  mse4$cmcp1[r, ] <- apply(fit4.cmcp1$beta[-1, ], 2, function(b) {
    mean((b - beta)^2)})
  mse4$cmcp2[r, ] <- apply(fit4.cmcp2$beta[-1, ], 2, function(b) {
    mean((b - beta)^2)})
  mse4$cmcp3[r, ] <- apply(fit4.cmcp3$beta[-1, ], 2, function(b) {
    mean((b - beta)^2)})
  
  mse4$gelasso1[r, ] <- apply(fit4.gelasso1$beta[-1, ], 2, function(b) {
    mean((b - beta)^2)})
  mse4$gelasso2[r, ] <- apply(fit4.gelasso2$beta[-1, ], 2, function(b) {
    mean((b - beta)^2)})
  mse4$gelasso3[r, ] <- apply(fit4.gelasso3$beta[-1, ], 2, function(b) {
    mean((b - beta)^2)})
  
  kappa4$grridge[r, ] <- sapply(fit4.grridge, function(s) {
    kappa2(cbind(beta!=0, replace(rep(FALSE, p), s$resEN$whichEN,
                                  TRUE)))$value})
  
  kappa4$gren1[r, ] <- apply(fit4.gren1$freq.model$groupreg$beta, 2,
                             function(b) {kappa2(cbind(beta!=0, b!=0))$value})
  kappa4$gren2[r, ] <- apply(fit4.gren2$freq.model$groupreg$beta, 2,
                             function(b) {kappa2(cbind(beta!=0, b!=0))$value})
  kappa4$gren3[r, ] <- apply(fit4.gren3$freq.model$groupreg$beta, 2,
                             function(b) {kappa2(cbind(beta!=0, b!=0))$value})
  
  kappa4$enet1[r, ] <- apply(fit4.gren1$freq.model$regular$beta, 2,
                             function(b) {kappa2(cbind(beta!=0, b!=0))$value})
  kappa4$enet2[r, ] <- apply(fit4.gren2$freq.model$regular$beta, 2,
                             function(b) {kappa2(cbind(beta!=0, b!=0))$value})
  kappa4$enet3[r, ] <- apply(fit4.gren3$freq.model$regular$beta, 2,
                             function(b) {kappa2(cbind(beta!=0, b!=0))$value})
  
  kappa4$sglasso1[r, ] <- apply(fit4.sglasso1$beta, 2, function(b) {
    kappa2(cbind(beta!=0, b!=0))$value})
  kappa4$sglasso2[r, ] <- apply(fit4.sglasso2$beta, 2, function(b) {
    kappa2(cbind(beta!=0, b!=0))$value})
  kappa4$sglasso3[r, ] <- apply(fit4.sglasso3$beta, 2, function(b) {
    kappa2(cbind(beta!=0, b!=0))$value})
  
  kappa4$cmcp1[r, ] <- apply(fit4.cmcp1$beta[-1, ], 2, function(b) {
    kappa2(cbind(beta!=0, b!=0))$value})
  kappa4$cmcp2[r, ] <- apply(fit4.cmcp2$beta[-1, ], 2, function(b) {
    kappa2(cbind(beta!=0, b!=0))$value})
  kappa4$cmcp3[r, ] <- apply(fit4.cmcp3$beta[-1, ], 2, function(b) {
    kappa2(cbind(beta!=0, b!=0))$value})
  
  kappa4$gelasso1[r, ] <- apply(fit4.gelasso1$beta[-1, ], 2, function(b) {
    kappa2(cbind(beta!=0, b!=0))$value})
  kappa4$gelasso2[r, ] <- apply(fit4.gelasso2$beta[-1, ], 2, function(b) {
    kappa2(cbind(beta!=0, b!=0))$value})
  kappa4$gelasso3[r, ] <- apply(fit4.gelasso3$beta[-1, ], 2, function(b) {
    kappa2(cbind(beta!=0, b!=0))$value})
  
  psel4$grridge[r, ] <- sapply(fit4.grridge, function(s) {
    length(s$resEN$whichEN)})
  
  psel4$gren1[r, ] <- fit4.gren1$freq.model$groupreg$df
  psel4$gren2[r, ] <- fit4.gren2$freq.model$groupreg$df
  psel4$gren3[r, ] <- fit4.gren3$freq.model$groupreg$df
  
  psel4$enet1[r, ] <- fit4.gren1$freq.model$regular$df
  psel4$enet2[r, ] <- fit4.gren2$freq.model$regular$df
  psel4$enet3[r, ] <- fit4.gren3$freq.model$regular$df
  
  psel4$sglasso1[r, ] <- apply(fit4.sglasso1$beta, 2, function(b) {sum(b!=0)})
  psel4$sglasso2[r, ] <- apply(fit4.sglasso2$beta, 2, function(b) {sum(b!=0)})
  psel4$sglasso3[r, ] <- apply(fit4.sglasso3$beta, 2, function(b) {sum(b!=0)})
  
  psel4$cmcp1[r, ] <- apply(fit4.cmcp1$beta, 2, function(b) {sum(b!=0)})
  psel4$cmcp2[r, ] <- apply(fit4.cmcp2$beta, 2, function(b) {sum(b!=0)})
  psel4$cmcp3[r, ] <- apply(fit4.cmcp3$beta, 2, function(b) {sum(b!=0)})
  
  psel4$gelasso1[r, ] <- apply(fit4.gelasso1$beta, 2, function(b) {sum(b!=0)})
  psel4$gelasso2[r, ] <- apply(fit4.gelasso2$beta, 2, function(b) {sum(b!=0)})
  psel4$gelasso3[r, ] <- apply(fit4.gelasso3$beta, 2, function(b) {sum(b!=0)})
  
  lambdagest4$grridge[r, ] <- fit4.grridge[[1]]$lambdamults$groups
  lambdagest4$gren1[r, ] <- fit4.gren1$lambdag$groups
  lambdagest4$gren2[r, ] <- fit4.gren2$lambdag$groups
  lambdagest4$gren3[r, ] <- fit4.gren3$lambdag$groups
  
  results4 <- list(auc=auc4, briers=briers4, mse=mse4, kappa=kappa4, psel=psel4,
                   lambdag=lambdagest4)
  save(results4, file=paste(path.res, "gren_sim4_res1.Rdata", sep=""))
  
}




path.res <- ifelse(as.character(Sys.info()[1])!="Darwin", "~/EBEN/results/",
                   "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/results/")
path.graph <- "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/graphs/"
load(paste(path.res, "gren_sim4_res1.Rdata", sep=""))
library(sp)

pred.loess <- function(method, measure, data) {
  ind <- data[["psel"]][[method]]
  dep <- data[[measure]][[method]]
  out <- tryCatch(list(x=sort(ind), y=predict(loess(dep[order(
    ind)] ~ sort(ind)))), warning=function(w) {
      list(x=unique(sort(ind)), y=sapply(unique(sort(ind)), function(psel) {
        mean(dep[ind==psel], na.rm=TRUE)}))})
  return(list(x=unique(out$x), y=sapply(unique(out$x), function(psel) {
    mean(out$y[out$x==psel])})))
}

### all performances
colors <- bpy.colors(8)[-c(1, 8)]
png(paste(path.graph, "gren_sim4_res1_performance1.png", sep=""),
    units="in", width=14, height=10, res=120)
par(mfrow=c(2, 2), mar=c(5.1, 5.1, 4.1, 2.1))

# kappa
ylim <- range(sapply(names(results4$kappa)[-1], pred.loess, "kappa", results4)[
  seq(2, 2*(length(results4$kappa) - 1), 2)])
plot(pred.loess("gren1", "kappa", results4), ylim=ylim, xlim=range(results4$psel), 
     main="a)", xlab="Number of selected features", ylab="Cohen's kappa", 
     cex.axis=1.5, cex.lab=2, cex.main=2, lwd=1.5, col=colors[1], type="l")
lines(pred.loess("gren2", "kappa", results4), lwd=1.5, col=colors[1], lty=2)
lines(pred.loess("gren3", "kappa", results4), lwd=1.5, col=colors[1], lty=3)

lines(pred.loess("enet1", "kappa", results4), lwd=1.5, col=colors[2], lty=1)
lines(pred.loess("enet2", "kappa", results4), lwd=1.5, col=colors[2], lty=2)
lines(pred.loess("enet3", "kappa", results4), lwd=1.5, col=colors[2], lty=3)

lines(pred.loess("sglasso1", "kappa", results4), lwd=1.5, col=colors[3], lty=1)
lines(pred.loess("sglasso2", "kappa", results4), lwd=1.5, col=colors[3], lty=2)
lines(pred.loess("sglasso3", "kappa", results4), lwd=1.5, col=colors[3], lty=3)

lines(pred.loess("cmcp1", "kappa", results4), lwd=1.5, col=colors[4], lty=1)
lines(pred.loess("cmcp2", "kappa", results4), lwd=1.5, col=colors[4], lty=2)
lines(pred.loess("cmcp3", "kappa", results4), lwd=1.5, col=colors[4], lty=3)

lines(pred.loess("gelasso1", "kappa", results4), lwd=1.5, col=colors[5], lty=1)
lines(pred.loess("gelasso2", "kappa", results4), lwd=1.5, col=colors[5], lty=2)
lines(pred.loess("gelasso3", "kappa", results4), lwd=1.5, col=colors[5], lty=3)

lines(pred.loess("grridge", "kappa", results4), lwd=1.5, col=colors[6], lty=1)

# mse
ylim <- range(sapply(names(results4$mse)[-1], pred.loess, "mse", results4)[
  seq(2, 2*(length(results4$mse) - 1), 2)], mean(results4$mse$ridge))
plot(pred.loess("gren1", "mse", results4), ylim=ylim, xlim=range(results4$psel), 
     main="b)", xlab="Number of selected features", ylab="MSE", cex.axis=1.5, 
     cex.lab=2, cex.main=2, lwd=1.5, col=colors[1], type="l")
lines(pred.loess("gren2", "mse", results4), lwd=1.5, col=colors[1], lty=2)
lines(pred.loess("gren3", "mse", results4), lwd=1.5, col=colors[1], lty=3)

lines(pred.loess("enet1", "mse", results4), lwd=1.5, col=colors[2], lty=1)
lines(pred.loess("enet2", "mse", results4), lwd=1.5, col=colors[2], lty=2)
lines(pred.loess("enet3", "mse", results4), lwd=1.5, col=colors[2], lty=3)

lines(pred.loess("sglasso1", "mse", results4), lwd=1.5, col=colors[3], lty=1)
lines(pred.loess("sglasso2", "mse", results4), lwd=1.5, col=colors[3], lty=2)
lines(pred.loess("sglasso3", "mse", results4), lwd=1.5, col=colors[3], lty=3)

lines(pred.loess("cmcp1", "mse", results4), lwd=1.5, col=colors[4], lty=1)
lines(pred.loess("cmcp2", "mse", results4), lwd=1.5, col=colors[4], lty=2)
lines(pred.loess("cmcp3", "mse", results4), lwd=1.5, col=colors[4], lty=3)

lines(pred.loess("gelasso1", "mse", results4), lwd=1.5, col=colors[5], lty=1)
lines(pred.loess("gelasso2", "mse", results4), lwd=1.5, col=colors[5], lty=2)
lines(pred.loess("gelasso3", "mse", results4), lwd=1.5, col=colors[5], lty=3)

lines(pred.loess("grridge", "mse", results4), lwd=1.5, col=colors[6], lty=1)
lines(range(results4$psel), rep(mean(results4$mse$ridge), 2), lwd=1.5, 
      col=colors[6], lty=2)

# auc
ylim <- range(sapply(names(results4$auc)[-1], pred.loess, "auc", results4)[
  seq(2, 2*(length(results4$auc) - 1), 2)], mean(results4$auc$ridge))
plot(pred.loess("gren1", "auc", results4), ylim=ylim, xlim=range(results4$psel), 
     main="c)", xlab="Number of selected features", ylab="AUC", cex.axis=1.5, 
     cex.lab=2, cex.main=2, lwd=1.5, col=colors[1], type="l")
lines(pred.loess("gren2", "auc", results4), lwd=1.5, col=colors[1], lty=2)
lines(pred.loess("gren3", "auc", results4), lwd=1.5, col=colors[1], lty=3)

lines(pred.loess("enet1", "auc", results4), lwd=1.5, col=colors[2], lty=1)
lines(pred.loess("enet2", "auc", results4), lwd=1.5, col=colors[2], lty=2)
lines(pred.loess("enet3", "auc", results4), lwd=1.5, col=colors[2], lty=3)

lines(pred.loess("sglasso1", "auc", results4), lwd=1.5, col=colors[3], lty=1)
lines(pred.loess("sglasso2", "auc", results4), lwd=1.5, col=colors[3], lty=2)
lines(pred.loess("sglasso3", "auc", results4), lwd=1.5, col=colors[3], lty=3)

lines(pred.loess("cmcp1", "auc", results4), lwd=1.5, col=colors[4], lty=1)
lines(pred.loess("cmcp2", "auc", results4), lwd=1.5, col=colors[4], lty=2)
lines(pred.loess("cmcp3", "auc", results4), lwd=1.5, col=colors[4], lty=3)

lines(pred.loess("gelasso1", "auc", results4), lwd=1.5, col=colors[5], lty=1)
lines(pred.loess("gelasso2", "auc", results4), lwd=1.5, col=colors[5], lty=2)
lines(pred.loess("gelasso3", "auc", results4), lwd=1.5, col=colors[5], lty=3)

lines(pred.loess("grridge", "auc", results4), lwd=1.5, col=colors[6], lty=1)
lines(range(results4$psel), rep(mean(results4$auc$ridge), 2), lwd=1.5, 
      col=colors[6], lty=2)

# briers
ylim <- range(sapply(names(results4$briers)[-1], pred.loess, "briers", results4)[
  seq(2, 2*(length(results4$briers) - 1), 2)], mean(results4$briers$ridge))
plot(pred.loess("gren1", "briers", results4), ylim=ylim, xlim=range(results4$psel), 
     main="d)", xlab="Number of selected features", ylab="Brier skill score", 
     cex.axis=1.5, cex.lab=2, cex.main=2, lwd=1.5, col=colors[1], type="l")
lines(pred.loess("gren2", "briers", results4), lwd=1.5, col=colors[1], lty=2)
lines(pred.loess("gren3", "briers", results4), lwd=1.5, col=colors[1], lty=3)

lines(pred.loess("enet1", "briers", results4), lwd=1.5, col=colors[2], lty=1)
lines(pred.loess("enet2", "briers", results4), lwd=1.5, col=colors[2], lty=2)
lines(pred.loess("enet3", "briers", results4), lwd=1.5, col=colors[2], lty=3)

lines(pred.loess("sglasso1", "briers", results4), lwd=1.5, col=colors[3], lty=1)
lines(pred.loess("sglasso2", "briers", results4), lwd=1.5, col=colors[3], lty=2)
lines(pred.loess("sglasso3", "briers", results4), lwd=1.5, col=colors[3], lty=3)

lines(pred.loess("cmcp1", "briers", results4), lwd=1.5, col=colors[4], lty=1)
lines(pred.loess("cmcp2", "briers", results4), lwd=1.5, col=colors[4], lty=2)
lines(pred.loess("cmcp3", "briers", results4), lwd=1.5, col=colors[4], lty=3)

lines(pred.loess("gelasso1", "briers", results4), lwd=1.5, col=colors[5], lty=1)
lines(pred.loess("gelasso2", "briers", results4), lwd=1.5, col=colors[5], lty=2)
lines(pred.loess("gelasso3", "briers", results4), lwd=1.5, col=colors[5], lty=3)

lines(pred.loess("grridge", "briers", results4), lwd=1.5, col=colors[6], lty=1)
lines(range(results4$psel), rep(mean(results4$briers$ridge), 2), lwd=1.5, 
      col=colors[6], lty=2)

# legend
leglabels <- c("gren", "enet", "sglasso", "cMCP", "gelasso", "ridge",
               expression(alpha==0.05), expression(alpha==0.5),
               expression(alpha==0.95))
legend("bottomright", merge=TRUE, seg.len=1, cex=1.3, fill=c(colors, 0, 0, 0), 
       lty=c(rep(NA, length(colors)), 1, 2, 3), 
       lwd=c(rep(NA, length(colors)), 1.5, 1.5),
       border=c(rep(1, length(colors)), 0, 0, 0), legend=leglabels)
dev.off()


### performance for alpha=0.05
colors <- bpy.colors(7)[-c(1, 7)]
png(paste(path.graph, "gren_sim4_res1_performance1_set1.png", sep=""),
    units="in", width=14, height=10, res=120)
par(mfrow=c(2, 2), mar=c(5.1, 5.1, 4.1, 2.1))

# kappa
ylim <- range(sapply(names(results4$kappa)[c(1, 2, 5, 8, 11, 14)], pred.loess, 
                     "kappa", results4)[seq(2, 12, 2)])
# xlim <- range(results4$psel)
xlim <- c(0, 250)
plot(pred.loess("gren1", "kappa", results4), ylim=ylim, xlim=xlim, 
     main="a)", xlab="Number of selected features", ylab="Cohen's kappa", 
     cex.axis=1.5, cex.lab=2, cex.main=2, lwd=1.5, col=colors[1], type="l")
lines(pred.loess("enet1", "kappa", results4), lwd=1.5, col=colors[1], lty=2)
lines(pred.loess("sglasso1", "kappa", results4), lwd=1.5, col=colors[2], lty=1)
lines(pred.loess("cmcp1", "kappa", results4), lwd=1.5, col=colors[3], lty=1)
lines(pred.loess("gelasso1", "kappa", results4), lwd=1.5, col=colors[4], lty=1)
lines(pred.loess("grridge", "kappa", results4), lwd=1.5, col=colors[5], lty=1)

# mse
ylim <- range(sapply(names(results4$mse)[c(2, 3, 6, 9, 12, 15)], pred.loess, 
                     "mse", results4)[seq(2, 12, 2)], mean(results4$mse$ridge))
xlim <- c(0, 250)
plot(pred.loess("gren1", "mse", results4), ylim=ylim, xlim=xlim, 
     main="b)", xlab="Number of selected features", ylab="MSE", cex.axis=1.5, 
     cex.lab=2, cex.main=2, lwd=1.5, col=colors[1], type="l")
lines(pred.loess("enet1", "mse", results4), lwd=1.5, col=colors[1], lty=2)
lines(pred.loess("sglasso1", "mse", results4), lwd=1.5, col=colors[2], lty=1)
lines(pred.loess("cmcp1", "mse", results4), lwd=1.5, col=colors[3], lty=1)
lines(pred.loess("gelasso1", "mse", results4), lwd=1.5, col=colors[4], lty=1)
lines(pred.loess("grridge", "mse", results4), lwd=1.5, col=colors[5], lty=1)
lines(xlim, rep(mean(results4$mse$ridge), 2), lwd=1.5, 
      col=colors[5], lty=2)

# auc
ylim <- range(sapply(names(results4$auc)[c(2, 3, 6, 9, 12, 15)], pred.loess, 
                     "auc", results4)[seq(2, 12, 2)], mean(results4$auc$ridge))
xlim <- c(0, 250)
plot(pred.loess("gren1", "auc", results4), ylim=ylim, xlim=xlim, 
     main="c)", xlab="Number of selected features", ylab="AUC", cex.axis=1.5, 
     cex.lab=2, cex.main=2, lwd=1.5, col=colors[1], type="l")
lines(pred.loess("enet1", "auc", results4), lwd=1.5, col=colors[1], lty=2)
lines(pred.loess("sglasso1", "auc", results4), lwd=1.5, col=colors[2], lty=1)
lines(pred.loess("cmcp1", "auc", results4), lwd=1.5, col=colors[3], lty=1)
lines(pred.loess("gelasso1", "auc", results4), lwd=1.5, col=colors[4], lty=1)
lines(pred.loess("grridge", "auc", results4), lwd=1.5, col=colors[5], lty=1)
lines(xlim, rep(mean(results4$auc$ridge), 2), lwd=1.5, 
      col=colors[5], lty=2)

# legend
leglabels <- c("enet", "sglasso", "cMCP", "gelasso", "ridge",
               "group-regularized", "regular")
legend("bottomright", merge=TRUE, seg.len=1, cex=1.3, fill=c(colors, 0, 0), 
       lty=c(rep(NA, length(colors)), 1, 2), 
       lwd=c(rep(NA, length(colors)), 1.5, 1.5),
       border=c(rep(1, length(colors)), 0, 0), legend=leglabels)

# briers
ylim <- range(sapply(names(results4$briers)[c(2, 3, 6, 9, 12, 15)], pred.loess, 
                     "briers", results4)[seq(2, 12, 2)], 
              mean(results4$briers$ridge))
xlim <- c(0, 250)
plot(pred.loess("gren1", "briers", results4), ylim=ylim, xlim=xlim, 
     main="d)", xlab="Number of selected features", ylab="Brier skill score", 
     cex.axis=1.5, cex.lab=2, cex.main=2, lwd=1.5, col=colors[1], type="l")
lines(pred.loess("enet1", "briers", results4), lwd=1.5, col=colors[1], lty=2)
lines(pred.loess("sglasso1", "briers", results4), lwd=1.5, col=colors[2], lty=1)
lines(pred.loess("cmcp1", "briers", results4), lwd=1.5, col=colors[3], lty=1)
lines(pred.loess("gelasso1", "briers", results4), lwd=1.5, col=colors[4], lty=1)
lines(pred.loess("grridge", "briers", results4), lwd=1.5, col=colors[5], lty=1)
lines(xlim, rep(mean(results4$briers$ridge), 2), lwd=1.5, 
      col=colors[5], lty=2)
dev.off()



### performance for alpha=0.5
colors <- bpy.colors(7)[-c(1, 7)]
png(paste(path.graph, "gren_sim4_res1_performance1_set2.png", sep=""),
    units="in", width=14, height=10, res=120)
par(mfrow=c(2, 2), mar=c(5.1, 5.1, 4.1, 2.1))

# kappa
ylim <- range(sapply(names(results4$kappa)[c(1, 3, 6, 9, 12, 15)], pred.loess, 
                     "kappa", results4)[seq(2, 12, 2)])
# xlim <- range(results4$psel)
xlim <- c(0, 250)
plot(pred.loess("gren2", "kappa", results4), ylim=ylim, xlim=xlim, 
     main="a)", xlab="Number of selected features", ylab="Cohen's kappa", 
     cex.axis=1.5, cex.lab=2, cex.main=2, lwd=1.5, col=colors[1], type="l")
lines(pred.loess("enet2", "kappa", results4), lwd=1.5, col=colors[1], lty=2)
lines(pred.loess("sglasso2", "kappa", results4), lwd=1.5, col=colors[2], lty=1)
lines(pred.loess("cmcp2", "kappa", results4), lwd=1.5, col=colors[3], lty=1)
lines(pred.loess("gelasso2", "kappa", results4), lwd=1.5, col=colors[4], lty=1)
lines(pred.loess("grridge", "kappa", results4), lwd=1.5, col=colors[5], lty=1)

# mse
ylim <- range(sapply(names(results4$mse)[c(2, 4, 7, 10, 13, 16)], pred.loess, 
                     "mse", results4)[seq(2, 12, 2)], mean(results4$mse$ridge))
xlim <- c(0, 250)
plot(pred.loess("gren2", "mse", results4), ylim=ylim, xlim=xlim, 
     main="b)", xlab="Number of selected features", ylab="MSE", cex.axis=1.5, 
     cex.lab=2, cex.main=2, lwd=1.5, col=colors[1], type="l")
lines(pred.loess("enet2", "mse", results4), lwd=1.5, col=colors[1], lty=2)
lines(pred.loess("sglasso2", "mse", results4), lwd=1.5, col=colors[2], lty=1)
lines(pred.loess("cmcp2", "mse", results4), lwd=1.5, col=colors[3], lty=1)
lines(pred.loess("gelasso2", "mse", results4), lwd=1.5, col=colors[4], lty=1)
lines(pred.loess("grridge", "mse", results4), lwd=1.5, col=colors[5], lty=1)
lines(xlim, rep(mean(results4$mse$ridge), 2), lwd=1.5, 
      col=colors[5], lty=2)

# auc
ylim <- range(sapply(names(results4$auc)[c(2, 4, 7, 10, 13, 16)], pred.loess, 
                     "auc", results4)[seq(2, 12, 2)], mean(results4$auc$ridge))
xlim <- c(0, 250)
plot(pred.loess("gren2", "auc", results4), ylim=ylim, xlim=xlim, 
     main="c)", xlab="Number of selected features", ylab="AUC", cex.axis=1.5, 
     cex.lab=2, cex.main=2, lwd=1.5, col=colors[1], type="l")
lines(pred.loess("enet2", "auc", results4), lwd=1.5, col=colors[1], lty=2)
lines(pred.loess("sglasso2", "auc", results4), lwd=1.5, col=colors[2], lty=1)
lines(pred.loess("cmcp2", "auc", results4), lwd=1.5, col=colors[3], lty=1)
lines(pred.loess("gelasso2", "auc", results4), lwd=1.5, col=colors[4], lty=1)
lines(pred.loess("grridge", "auc", results4), lwd=1.5, col=colors[5], lty=1)
lines(xlim, rep(mean(results4$auc$ridge), 2), lwd=1.5, 
      col=colors[5], lty=2)

# legend
leglabels <- c("enet", "sglasso", "cMCP", "gelasso", "ridge",
               "group-regularized", "regular")
legend("bottomright", merge=TRUE, seg.len=1, cex=1.3, fill=c(colors, 0, 0), 
       lty=c(rep(NA, length(colors)), 1, 2), 
       lwd=c(rep(NA, length(colors)), 1.5, 1.5),
       border=c(rep(1, length(colors)), 0, 0), legend=leglabels)

# briers
ylim <- range(sapply(names(results4$briers)[c(2, 4, 7, 10, 13, 16)], pred.loess, 
                     "briers", results4)[seq(2, 12, 2)], 
              mean(results4$briers$ridge))
xlim <- c(0, 250)
plot(pred.loess("gren2", "briers", results4), ylim=ylim, xlim=xlim, 
     main="d)", xlab="Number of selected features", ylab="Brier skill score", 
     cex.axis=1.5, cex.lab=2, cex.main=2, lwd=1.5, col=colors[1], type="l")
lines(pred.loess("enet2", "briers", results4), lwd=1.5, col=colors[1], lty=2)
lines(pred.loess("sglasso2", "briers", results4), lwd=1.5, col=colors[2], lty=1)
lines(pred.loess("cmcp2", "briers", results4), lwd=1.5, col=colors[3], lty=1)
lines(pred.loess("gelasso2", "briers", results4), lwd=1.5, col=colors[4], lty=1)
lines(pred.loess("grridge", "briers", results4), lwd=1.5, col=colors[5], lty=1)
lines(xlim, rep(mean(results4$briers$ridge), 2), lwd=1.5, 
      col=colors[5], lty=2)
dev.off()



### performance for alpha=0.95
colors <- bpy.colors(7)[-c(1, 7)]
png(paste(path.graph, "gren_sim4_res1_performance1_set3.png", sep=""),
    units="in", width=14, height=10, res=120)
par(mfrow=c(2, 2), mar=c(5.1, 5.1, 4.1, 2.1))

# kappa
ylim <- range(sapply(names(results4$kappa)[c(1, 4, 7, 10, 13, 16)], pred.loess, 
                     "kappa", results4)[seq(2, 12, 2)])
# xlim <- range(results4$psel)
xlim <- c(0, 250)
plot(pred.loess("gren3", "kappa", results4), ylim=ylim, xlim=xlim, 
     main="a)", xlab="Number of selected features", ylab="Cohen's kappa", 
     cex.axis=1.5, cex.lab=2, cex.main=2, lwd=1.5, col=colors[1], type="l")
lines(pred.loess("enet3", "kappa", results4), lwd=1.5, col=colors[1], lty=2)
lines(pred.loess("sglasso3", "kappa", results4), lwd=1.5, col=colors[2], lty=1)
lines(pred.loess("cmcp3", "kappa", results4), lwd=1.5, col=colors[3], lty=1)
lines(pred.loess("gelasso3", "kappa", results4), lwd=1.5, col=colors[4], lty=1)
lines(pred.loess("grridge", "kappa", results4), lwd=1.5, col=colors[5], lty=1)

# mse
ylim <- range(sapply(names(results4$mse)[c(2, 5, 8, 11, 14, 17)], pred.loess, 
                     "mse", results4)[seq(2, 12, 2)], mean(results4$mse$ridge))
xlim <- c(0, 250)
plot(pred.loess("gren3", "mse", results4), ylim=ylim, xlim=xlim, 
     main="b)", xlab="Number of selected features", ylab="MSE", cex.axis=1.5, 
     cex.lab=2, cex.main=2, lwd=1.5, col=colors[1], type="l")
lines(pred.loess("enet3", "mse", results4), lwd=1.5, col=colors[1], lty=2)
lines(pred.loess("sglasso3", "mse", results4), lwd=1.5, col=colors[2], lty=1)
lines(pred.loess("cmcp3", "mse", results4), lwd=1.5, col=colors[3], lty=1)
lines(pred.loess("gelasso3", "mse", results4), lwd=1.5, col=colors[4], lty=1)
lines(pred.loess("grridge", "mse", results4), lwd=1.5, col=colors[5], lty=1)
lines(xlim, rep(mean(results4$mse$ridge), 2), lwd=1.5, 
      col=colors[5], lty=2)

# auc
ylim <- range(sapply(names(results4$auc)[c(2, 5, 8, 11, 14, 17)], pred.loess, 
                     "auc", results4)[seq(2, 12, 2)], mean(results4$auc$ridge))
xlim <- c(0, 250)
plot(pred.loess("gren3", "auc", results4), ylim=ylim, xlim=xlim, 
     main="c)", xlab="Number of selected features", ylab="AUC", cex.axis=1.5, 
     cex.lab=2, cex.main=2, lwd=1.5, col=colors[1], type="l")
lines(pred.loess("enet3", "auc", results4), lwd=1.5, col=colors[1], lty=2)
lines(pred.loess("sglasso3", "auc", results4), lwd=1.5, col=colors[2], lty=1)
lines(pred.loess("cmcp3", "auc", results4), lwd=1.5, col=colors[3], lty=1)
lines(pred.loess("gelasso3", "auc", results4), lwd=1.5, col=colors[4], lty=1)
lines(pred.loess("grridge", "auc", results4), lwd=1.5, col=colors[5], lty=1)
lines(xlim, rep(mean(results4$auc$ridge), 2), lwd=1.5, 
      col=colors[5], lty=2)

# legend
leglabels <- c("enet", "sglasso", "cMCP", "gelasso", "ridge",
               "group-regularized", "regular")
legend("bottomright", merge=TRUE, seg.len=1, cex=1.3, fill=c(colors, 0, 0), 
       lty=c(rep(NA, length(colors)), 1, 2), 
       lwd=c(rep(NA, length(colors)), 1.5, 1.5),
       border=c(rep(1, length(colors)), 0, 0), legend=leglabels)

# briers
ylim <- range(sapply(names(results4$briers)[c(2, 5, 8, 11, 14, 17)], pred.loess, 
                     "briers", results4)[seq(2, 12, 2)], 
              mean(results4$briers$ridge))
xlim <- c(0, 250)
plot(pred.loess("gren3", "briers", results4), ylim=ylim, xlim=xlim, 
     main="d)", xlab="Number of selected features", ylab="Brier skill score", 
     cex.axis=1.5, cex.lab=2, cex.main=2, lwd=1.5, col=colors[1], type="l")
lines(pred.loess("enet3", "briers", results4), lwd=1.5, col=colors[1], lty=2)
lines(pred.loess("sglasso3", "briers", results4), lwd=1.5, col=colors[2], lty=1)
lines(pred.loess("cmcp3", "briers", results4), lwd=1.5, col=colors[3], lty=1)
lines(pred.loess("gelasso3", "briers", results4), lwd=1.5, col=colors[4], lty=1)
lines(pred.loess("grridge", "briers", results4), lwd=1.5, col=colors[5], lty=1)
lines(xlim, rep(mean(results4$briers$ridge), 2), lwd=1.5, 
      col=colors[5], lty=2)
dev.off()