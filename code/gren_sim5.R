path.res <- ifelse(as.character(Sys.info()[1])!="Darwin", "~/EBEN/results/",
                   "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/results/")

library(mvtnorm)
library(gren)
library(pROC)
library(GRridge)
library(grpreg)
library(SGL)
library(irr)

##################### scenario 5
n <- 100
p <- 1000
ntest <- 1000

G <- 10
rho <- 0.5
Sigma <- diag(p)
for(i in 1:p) {
  for(j in 1:p) {
    Sigma[i, j] <- rho^abs(i - j)
  }
}

beta.group <- c(rep(0, 8), 0.2, 0.5)
q.zero <- 0.85

part1 <- rep(c(1:G), each=p/G)
csel <- c(seq(1, 8, 1), seq(10, 20, 2), seq(25, 60, 5), seq(70, 140, 10))
nreps <- 100
psel5 <- kappa5 <- list(grridge=vector("list", nreps),
                        enet1=vector("list", nreps),
                        enet2=vector("list", nreps),
                        enet3=vector("list", nreps),
                        gren1=vector("list", nreps),
                        gren2=vector("list", nreps),
                        gren3=vector("list", nreps),
                        sglasso1=vector("list", nreps),
                        sglasso2=vector("list", nreps),
                        sglasso3=vector("list", nreps),
                        cmcp1=vector("list", nreps),
                        cmcp2=vector("list", nreps),
                        cmcp3=vector("list", nreps),
                        gelasso1=vector("list", nreps),
                        gelasso2=vector("list", nreps),
                        gelasso3=vector("list", nreps))
auc5 <- mse5 <- briers5 <- list(ridge=numeric(nreps),
                                grridge=vector("list", nreps),
                                enet1=vector("list", nreps),
                                enet2=vector("list", nreps),
                                enet3=vector("list", nreps),
                                gren1=vector("list", nreps),
                                gren2=vector("list", nreps),
                                gren3=vector("list", nreps),
                                sglasso1=vector("list", nreps),
                                sglasso2=vector("list", nreps),
                                sglasso3=vector("list", nreps),
                                cmcp1=vector("list", nreps),
                                cmcp2=vector("list", nreps),
                                cmcp3=vector("list", nreps),
                                gelasso1=vector("list", nreps),
                                gelasso2=vector("list", nreps),
                                gelasso3=vector("list", nreps))
lambdagest5 <- list(grridge=matrix(NA, nrow=nreps, ncol=G),
                    gren1=matrix(NA, nrow=nreps, ncol=G),
                    gren2=matrix(NA, nrow=nreps, ncol=G),
                    gren3=matrix(NA, nrow=nreps, ncol=G))
for(r in 1:nreps) {
  set.seed(2018 + r)
  print(paste("rep", r))
  beta <- numeric(p)
  for(g in 1:G) {
    beta[((g - 1)*p/G + 1):(g*p/G)] <- rep(beta.group[g], p/G)
    if(beta.group[g]!=0) {
      beta[((g - 1)*p/G + 1):((g - 1)*p/G + q.zero*p/G)] <- 0
    }
  }
  
  x <- rmvnorm(n, mean=rep(0, p), sigma=Sigma)
  y <- rbinom(n, 1, as.numeric(exp(x %*% beta)/(1 + exp(x %*% beta))))
  
  xtest <- rmvnorm(ntest, mean=rep(0, p), sigma=Sigma)
  ytest <- rbinom(ntest, 1, as.numeric(exp(xtest %*% beta)/
                                         (1 + exp(xtest %*% beta))))
  
  fit5.ridge <- cv.glmnet(x, y, alpha=0, standardize=FALSE)
  
  fit5.grridge <- vector("list", length(csel))
  invisible(capture.output(
    fit5.grridge[[1]] <- grridge(t(x), y, partitions=list(
      groups=CreatePartition(as.factor(part1))), selection=TRUE, maxsel=csel[1],
      trace=FALSE, standardizeX=FALSE)))
  for(s in 2:length(csel)) {
    invisible(capture.output(
      fit5.grridge[[s]] <- grridge(t(x), y, partitions=list(
        groups=CreatePartition(as.factor(part1))), selection=TRUE,
        maxsel=csel[s], optl=fit5.grridge[[1]]$optl, trace=FALSE,
        standardizeX=FALSE)))
  }
  
  fit5.gren1 <- gren(x, y, partitions=list(groups=part1), alpha=0.05,
                     trace=FALSE)
  fit5.gren2 <- gren(x, y, partitions=list(groups=part1), alpha=0.5,
                     trace=FALSE)
  fit5.gren3 <- gren(x, y, partitions=list(groups=part1), alpha=0.95,
                     trace=FALSE)
  
  fit5.sglasso1 <- SGL(list(x=x, y=y), part1, type="logit", alpha=0.05,
                       standardize=FALSE, nlam=100)
  fit5.sglasso2 <- SGL(list(x=x, y=y), part1, type="logit", alpha=0.5,
                       standardize=FALSE, nlam=100)
  fit5.sglasso3 <- SGL(list(x=x, y=y), part1, type="logit", alpha=0.95,
                       standardize=FALSE, nlam=100)
  
  fit5.cmcp1 <- grpreg(x, y, part1, penalty="cMCP", alpha=0.05)
  fit5.cmcp2 <- grpreg(x, y, part1, penalty="cMCP", alpha=0.5)
  fit5.cmcp3 <- grpreg(x, y, part1, penalty="cMCP", alpha=0.95)
  
  fit5.gelasso1 <- grpreg(x, y, part1, penalty="gel", alpha=0.05)
  fit5.gelasso2 <- grpreg(x, y, part1, penalty="gel", alpha=0.5)
  fit5.gelasso3 <- grpreg(x, y, part1, penalty="gel", alpha=0.95)
  
  pred5.ridge <- as.numeric(predict(fit5.ridge, xtest, "lambda.min"))
  
  pred5.grridge <- sapply(fit5.grridge, function(s) {
    predict.grridge(s, t(xtest))[, 3]})
  
  pred5.gren1 <- predict(fit5.gren1, xtest, type="groupreg",
                         s=fit5.gren1$freq.model$groupreg$lambda)
  pred5.gren2 <- predict(fit5.gren2, xtest, type="groupreg",
                         s=fit5.gren2$freq.model$groupreg$lambda)
  pred5.gren3 <- predict(fit5.gren3, xtest, type="groupreg",
                         s=fit5.gren3$freq.model$groupreg$lambda)
  
  pred5.enet1 <- predict(fit5.gren1, xtest, type="regular",
                         s=fit5.gren1$freq.model$regular$lambda)
  pred5.enet2 <- predict(fit5.gren2, xtest, type="regular",
                         s=fit5.gren2$freq.model$regular$lambda)
  pred5.enet3 <- predict(fit5.gren3, xtest, type="regular",
                         s=fit5.gren3$freq.model$regular$lambda)
  
  pred5.sglasso1 <- predictSGL(fit5.sglasso1, xtest)
  pred5.sglasso2 <- predictSGL(fit5.sglasso2, xtest)
  pred5.sglasso3 <- predictSGL(fit5.sglasso3, xtest)
  
  pred5.cmcp1 <- predict(fit5.cmcp1, xtest)
  pred5.cmcp2 <- predict(fit5.cmcp2, xtest)
  pred5.cmcp3 <- predict(fit5.cmcp3, xtest)
  
  pred5.gelasso1 <- predict(fit5.gelasso1, xtest)
  pred5.gelasso2 <- predict(fit5.gelasso2, xtest)
  pred5.gelasso3 <- predict(fit5.gelasso3, xtest)
  
  auc5$ridge[r] <- pROC::auc(ytest, pred5.ridge)
  
  auc5$grridge[[r]] <- apply(pred5.grridge, 2, function(pred) {
    pROC::auc(ytest, pred)})
  
  auc5$gren1[[r]] <- apply(pred5.gren1, 2, function(pred) {
    pROC::auc(ytest, pred)})
  auc5$gren2[[r]] <- apply(pred5.gren2, 2, function(pred) {
    pROC::auc(ytest, pred)})
  auc5$gren3[[r]] <- apply(pred5.gren3, 2, function(pred) {
    pROC::auc(ytest, pred)})
  
  auc5$enet1[[r]] <- apply(pred5.enet1, 2, function(pred) {
    pROC::auc(ytest, pred)})
  auc5$enet2[[r]] <- apply(pred5.enet2, 2, function(pred) {
    pROC::auc(ytest, pred)})
  auc5$enet3[[r]] <- apply(pred5.enet3, 2, function(pred) {
    pROC::auc(ytest, pred)})
  
  auc5$sglasso1[[r]] <- apply(pred5.sglasso1, 2, function(pred) {
    pROC::auc(ytest, pred)})
  auc5$sglasso2[[r]] <- apply(pred5.sglasso2, 2, function(pred) {
    pROC::auc(ytest, pred)})
  auc5$sglasso3[[r]] <- apply(pred5.sglasso3, 2, function(pred) {
    pROC::auc(ytest, pred)})
  
  auc5$cmcp1[[r]] <- apply(pred5.cmcp1, 2, function(pred) {
    pROC::auc(ytest, pred)})
  auc5$cmcp2[[r]] <- apply(pred5.cmcp2, 2, function(pred) {
    pROC::auc(ytest, pred)})
  auc5$cmcp3[[r]] <- apply(pred5.cmcp3, 2, function(pred) {
    pROC::auc(ytest, pred)})
  
  auc5$gelasso1[[r]] <- apply(pred5.gelasso1, 2, function(pred) {
    pROC::auc(ytest, pred)})
  auc5$gelasso2[[r]] <- apply(pred5.gelasso2, 2, function(pred) {
    pROC::auc(ytest, pred)})
  auc5$gelasso3[[r]] <- apply(pred5.gelasso3, 2, function(pred) {
    pROC::auc(ytest, pred)})
  
  const <- sum((ytest - mean(ytest))^2)
  briers5$ridge[r] <- 1 - sum((ytest - pred5.ridge)^2)/const
  
  briers5$grridge[[r]] <- apply(pred5.grridge, 2, function(pred) {
    1 - sum((ytest - pred)^2)/const})
  
  briers5$gren1[[r]] <- apply(pred5.gren1, 2, function(pred) {
    1 - sum((ytest - pred)^2)/const})
  briers5$gren2[[r]] <- apply(pred5.gren2, 2, function(pred) {
    1 - sum((ytest - pred)^2)/const})
  briers5$gren3[[r]] <- apply(pred5.gren3, 2, function(pred) {
    1 - sum((ytest - pred)^2)/const})
  
  briers5$enet1[[r]] <- apply(pred5.enet1, 2, function(pred) {
    1 - sum((ytest - pred)^2)/const})
  briers5$enet2[[r]] <- apply(pred5.enet2, 2, function(pred) {
    1 - sum((ytest - pred)^2)/const})
  briers5$enet3[[r]] <- apply(pred5.enet3, 2, function(pred) {
    1 - sum((ytest - pred)^2)/const})
  
  briers5$sglasso1[[r]] <- apply(pred5.sglasso1, 2, function(pred) {
    1 - sum((ytest - pred)^2)/const})
  briers5$sglasso2[[r]] <- apply(pred5.sglasso2, 2, function(pred) {
    1 - sum((ytest - pred)^2)/const})
  briers5$sglasso3[[r]] <- apply(pred5.sglasso3, 2, function(pred) {
    1 - sum((ytest - pred)^2)/const})
  
  briers5$cmcp1[[r]] <- apply(pred5.cmcp1, 2, function(pred) {
    1 - sum((ytest - pred)^2)/const})
  briers5$cmcp2[[r]] <- apply(pred5.cmcp2, 2, function(pred) {
    1 - sum((ytest - pred)^2)/const})
  briers5$cmcp3[[r]] <- apply(pred5.cmcp3, 2, function(pred) {
    1 - sum((ytest - pred)^2)/const})
  
  briers5$gelasso1[[r]] <- apply(pred5.gelasso1, 2, function(pred) {
    1 - sum((ytest - pred)^2)/const})
  briers5$gelasso2[[r]] <- apply(pred5.gelasso2, 2, function(pred) {
    1 - sum((ytest - pred)^2)/const})
  briers5$gelasso3[[r]] <- apply(pred5.gelasso3, 2, function(pred) {
    1 - sum((ytest - pred)^2)/const})
  
  mse5$ridge[r] <- mean((coef(fit5.ridge)[-1] - beta)^2)
  
  mse5$grridge[[r]] <- sapply(fit5.grridge, function(s) {
    mean((replace(rep(0, p), s$resEN$whichEN, s$resEN$betasEN) - beta)^2)})
  
  mse5$gren1[[r]] <- apply(fit5.gren1$freq.model$groupreg$beta, 2, function(b) {
    mean((b - beta)^2)})
  mse5$gren2[[r]] <- apply(fit5.gren2$freq.model$groupreg$beta, 2, function(b) {
    mean((b - beta)^2)})
  mse5$gren3[[r]] <- apply(fit5.gren3$freq.model$groupreg$beta, 2, function(b) {
    mean((b - beta)^2)})
  
  mse5$enet1[[r]] <- apply(fit5.gren1$freq.model$regular$beta, 2, function(b) {
    mean((b - beta)^2)})
  mse5$enet2[[r]] <- apply(fit5.gren2$freq.model$regular$beta, 2, function(b) {
    mean((b - beta)^2)})
  mse5$enet3[[r]] <- apply(fit5.gren3$freq.model$regular$beta, 2, function(b) {
    mean((b - beta)^2)})
  
  mse5$sglasso1[[r]] <- apply(fit5.sglasso1$beta, 2, function(b) {
    mean((b - beta)^2)})
  mse5$sglasso2[[r]] <- apply(fit5.sglasso2$beta, 2, function(b) {
    mean((b - beta)^2)})
  mse5$sglasso3[[r]] <- apply(fit5.sglasso3$beta, 2, function(b) {
    mean((b - beta)^2)})
  
  mse5$cmcp1[[r]] <- apply(fit5.cmcp1$beta[-1, ], 2, function(b) {
    mean((b - beta)^2)})
  mse5$cmcp2[[r]] <- apply(fit5.cmcp2$beta[-1, ], 2, function(b) {
    mean((b - beta)^2)})
  mse5$cmcp3[[r]] <- apply(fit5.cmcp3$beta[-1, ], 2, function(b) {
    mean((b - beta)^2)})
  
  mse5$gelasso1[[r]] <- apply(fit5.gelasso1$beta[-1, ], 2, function(b) {
    mean((b - beta)^2)})
  mse5$gelasso2[[r]] <- apply(fit5.gelasso2$beta[-1, ], 2, function(b) {
    mean((b - beta)^2)})
  mse5$gelasso3[[r]] <- apply(fit5.gelasso3$beta[-1, ], 2, function(b) {
    mean((b - beta)^2)})
  
  kappa5$grridge[[r]] <- sapply(fit5.grridge, function(s) {
    kappa2(cbind(beta!=0, replace(rep(FALSE, p), s$resEN$whichEN,
                                  TRUE)))$value})
  
  kappa5$gren1[[r]] <- apply(fit5.gren1$freq.model$groupreg$beta, 2,
                             function(b) {kappa2(cbind(beta!=0, b!=0))$value})
  kappa5$gren2[[r]] <- apply(fit5.gren2$freq.model$groupreg$beta, 2,
                             function(b) {kappa2(cbind(beta!=0, b!=0))$value})
  kappa5$gren3[[r]] <- apply(fit5.gren3$freq.model$groupreg$beta, 2,
                             function(b) {kappa2(cbind(beta!=0, b!=0))$value})
  
  kappa5$enet1[[r]] <- apply(fit5.gren1$freq.model$regular$beta, 2,
                             function(b) {kappa2(cbind(beta!=0, b!=0))$value})
  kappa5$enet2[[r]] <- apply(fit5.gren2$freq.model$regular$beta, 2,
                             function(b) {kappa2(cbind(beta!=0, b!=0))$value})
  kappa5$enet3[[r]] <- apply(fit5.gren3$freq.model$regular$beta, 2,
                             function(b) {kappa2(cbind(beta!=0, b!=0))$value})
  
  kappa5$sglasso1[[r]] <- apply(fit5.sglasso1$beta, 2, function(b) {
    kappa2(cbind(beta!=0, b!=0))$value})
  kappa5$sglasso2[[r]] <- apply(fit5.sglasso2$beta, 2, function(b) {
    kappa2(cbind(beta!=0, b!=0))$value})
  kappa5$sglasso3[[r]] <- apply(fit5.sglasso3$beta, 2, function(b) {
    kappa2(cbind(beta!=0, b!=0))$value})
  
  kappa5$cmcp1[[r]] <- apply(fit5.cmcp1$beta[-1, ], 2, function(b) {
    kappa2(cbind(beta!=0, b!=0))$value})
  kappa5$cmcp2[[r]] <- apply(fit5.cmcp2$beta[-1, ], 2, function(b) {
    kappa2(cbind(beta!=0, b!=0))$value})
  kappa5$cmcp3[[r]] <- apply(fit5.cmcp3$beta[-1, ], 2, function(b) {
    kappa2(cbind(beta!=0, b!=0))$value})
  
  kappa5$gelasso1[[r]] <- apply(fit5.gelasso1$beta[-1, ], 2, function(b) {
    kappa2(cbind(beta!=0, b!=0))$value})
  kappa5$gelasso2[[r]] <- apply(fit5.gelasso2$beta[-1, ], 2, function(b) {
    kappa2(cbind(beta!=0, b!=0))$value})
  kappa5$gelasso3[[r]] <- apply(fit5.gelasso3$beta[-1, ], 2, function(b) {
    kappa2(cbind(beta!=0, b!=0))$value})
  
  psel5$grridge[[r]] <- sapply(fit5.grridge, function(s) {
    length(s$resEN$whichEN)})
  
  psel5$gren1[[r]] <- fit5.gren1$freq.model$groupreg$df
  psel5$gren2[[r]] <- fit5.gren2$freq.model$groupreg$df
  psel5$gren3[[r]] <- fit5.gren3$freq.model$groupreg$df
  
  psel5$enet1[[r]] <- fit5.gren1$freq.model$regular$df
  psel5$enet2[[r]] <- fit5.gren2$freq.model$regular$df
  psel5$enet3[[r]] <- fit5.gren3$freq.model$regular$df
  
  psel5$sglasso1[[r]] <- apply(fit5.sglasso1$beta, 2, function(b) {sum(b!=0)})
  psel5$sglasso2[[r]] <- apply(fit5.sglasso2$beta, 2, function(b) {sum(b!=0)})
  psel5$sglasso3[[r]] <- apply(fit5.sglasso3$beta, 2, function(b) {sum(b!=0)})
  
  psel5$cmcp1[[r]] <- apply(fit5.cmcp1$beta, 2, function(b) {sum(b!=0)})
  psel5$cmcp2[[r]] <- apply(fit5.cmcp2$beta, 2, function(b) {sum(b!=0)})
  psel5$cmcp3[[r]] <- apply(fit5.cmcp3$beta, 2, function(b) {sum(b!=0)})
  
  psel5$gelasso1[[r]] <- apply(fit5.gelasso1$beta, 2, function(b) {sum(b!=0)})
  psel5$gelasso2[[r]] <- apply(fit5.gelasso2$beta, 2, function(b) {sum(b!=0)})
  psel5$gelasso3[[r]] <- apply(fit5.gelasso3$beta, 2, function(b) {sum(b!=0)})
  
  lambdagest5$grridge[r, ] <- fit5.grridge[[1]]$lambdamults$groups
  lambdagest5$gren1[r, ] <- fit5.gren1$lambdag$groups
  lambdagest5$gren2[r, ] <- fit5.gren2$lambdag$groups
  lambdagest5$gren3[r, ] <- fit5.gren3$lambdag$groups
  
  results5 <- list(auc=auc5, briers=briers5, mse=mse5, kappa=kappa5, psel=psel5,
                   lambdag=lambdagest5)
  save(results5, file=paste(path.res, "gren_sim5_res1.Rdata", sep=""))
  
}


path.res <- ifelse(as.character(Sys.info()[1])!="Darwin", "~/EBEN/results/",
                   "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/results/")

path.graph <- "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/graphs/"
load(paste(path.res, "gren_sim5_res1.Rdata", sep=""))
library(sp)

pred.loess <- function(method, measure, data) {
  ind <- unlist(data[["psel"]][[method]])
  dep <- unlist(data[[measure]][[method]])
  out <- tryCatch(list(x=sort(ind), y=predict(loess(dep[order(
    ind)] ~ sort(ind)))), warning=function(w) {
      list(x=unique(sort(ind)), y=sapply(unique(sort(ind)), function(psel) {
        mean(dep[ind==psel], na.rm=TRUE)}))})
  return(list(x=unique(out$x), y=sapply(unique(out$x), function(psel) {
    mean(out$y[out$x==psel])})))
}



### performance for alpha=0.05
colors <- bpy.colors(7)[-c(1, 7)]
png(paste(path.graph, "gren_sim5_res1_performance1_set1.png", sep=""),
    units="in", width=14, height=10, res=120)
par(mfrow=c(2, 2), mar=c(5.1, 5.1, 4.1, 2.1))

# kappa
ylim <- range(sapply(names(results5$kappa)[c(1, 2, 5, 8, 11, 14)], pred.loess,
                     "kappa", results5)[seq(2, 12, 2)])
# xlim <- range(results5$psel)
xlim <- c(0, 500)
plot(pred.loess("gren1", "kappa", results5), ylim=ylim, xlim=xlim,
     main="a)", xlab="Number of selected features", ylab="Cohen's kappa",
     cex.axis=1.5, cex.lab=2, cex.main=2, lwd=1.5, col=colors[1], type="l")
lines(pred.loess("enet1", "kappa", results5), lwd=1.5, col=colors[1], lty=2)
lines(pred.loess("sglasso1", "kappa", results5), lwd=1.5, col=colors[2], lty=1)
lines(pred.loess("cmcp1", "kappa", results5), lwd=1.5, col=colors[3], lty=1)
lines(pred.loess("gelasso1", "kappa", results5), lwd=1.5, col=colors[4], lty=1)
lines(pred.loess("grridge", "kappa", results5), lwd=1.5, col=colors[5], lty=1)

# mse
ylim <- range(sapply(names(results5$mse)[c(2, 3, 6, 9, 12, 15)], pred.loess,
                     "mse", results5)[seq(2, 12, 2)], mean(results5$mse$ridge))
xlim <- c(0, 500)
plot(pred.loess("gren1", "mse", results5), ylim=ylim, xlim=xlim,
     main="b)", xlab="Number of selected features", ylab="MSE", cex.axis=1.5,
     cex.lab=2, cex.main=2, lwd=1.5, col=colors[1], type="l")
lines(pred.loess("enet1", "mse", results5), lwd=1.5, col=colors[1], lty=2)
lines(pred.loess("sglasso1", "mse", results5), lwd=1.5, col=colors[2], lty=1)
lines(pred.loess("cmcp1", "mse", results5), lwd=1.5, col=colors[3], lty=1)
lines(pred.loess("gelasso1", "mse", results5), lwd=1.5, col=colors[4], lty=1)
lines(pred.loess("grridge", "mse", results5), lwd=1.5, col=colors[5], lty=1)
lines(xlim, rep(mean(results5$mse$ridge), 2), lwd=1.5,
      col=colors[5], lty=2)

# auc
ylim <- range(sapply(names(results5$auc)[c(2, 3, 6, 9, 12, 15)], pred.loess,
                     "auc", results5)[seq(2, 12, 2)], mean(results5$auc$ridge))
xlim <- c(0, 500)
plot(pred.loess("gren1", "auc", results5), ylim=ylim, xlim=xlim,
     main="c)", xlab="Number of selected features", ylab="AUC", cex.axis=1.5,
     cex.lab=2, cex.main=2, lwd=1.5, col=colors[1], type="l")
lines(pred.loess("enet1", "auc", results5), lwd=1.5, col=colors[1], lty=2)
lines(pred.loess("sglasso1", "auc", results5), lwd=1.5, col=colors[2], lty=1)
lines(pred.loess("cmcp1", "auc", results5), lwd=1.5, col=colors[3], lty=1)
lines(pred.loess("gelasso1", "auc", results5), lwd=1.5, col=colors[4], lty=1)
lines(pred.loess("grridge", "auc", results5), lwd=1.5, col=colors[5], lty=1)
lines(xlim, rep(mean(results5$auc$ridge), 2), lwd=1.5,
      col=colors[5], lty=2)

# legend
leglabels <- c("enet", "sglasso", "cMCP", "gelasso", "ridge",
               "group-regularized", "regular")
legend("bottomright", merge=TRUE, seg.len=1, cex=1.3, fill=c(colors, 0, 0),
       lty=c(rep(NA, length(colors)), 1, 2),
       lwd=c(rep(NA, length(colors)), 1.5, 1.5),
       border=c(rep(1, length(colors)), 0, 0), legend=leglabels)

# briers
ylim <- range(sapply(names(results5$briers)[c(2, 3, 6, 9, 12, 15)], pred.loess,
                     "briers", results5)[seq(2, 12, 2)],
              mean(results5$briers$ridge))
xlim <- c(0, 500)
plot(pred.loess("gren1", "briers", results5), ylim=ylim, xlim=xlim,
     main="d)", xlab="Number of selected features", ylab="Brier skill score",
     cex.axis=1.5, cex.lab=2, cex.main=2, lwd=1.5, col=colors[1], type="l")
lines(pred.loess("enet1", "briers", results5), lwd=1.5, col=colors[1], lty=2)
lines(pred.loess("sglasso1", "briers", results5), lwd=1.5, col=colors[2], lty=1)
lines(pred.loess("cmcp1", "briers", results5), lwd=1.5, col=colors[3], lty=1)
lines(pred.loess("gelasso1", "briers", results5), lwd=1.5, col=colors[4], lty=1)
lines(pred.loess("grridge", "briers", results5), lwd=1.5, col=colors[5], lty=1)
lines(xlim, rep(mean(results5$briers$ridge), 2), lwd=1.5,
      col=colors[5], lty=2)
dev.off()



### performance for alpha=0.5
colors <- bpy.colors(7)[-c(1, 7)]
png(paste(path.graph, "gren_sim5_res1_performance1_set2.png", sep=""),
    units="in", width=14, height=10, res=120)
par(mfrow=c(2, 2), mar=c(5.1, 5.1, 4.1, 2.1))

# kappa
ylim <- range(sapply(names(results5$kappa)[c(1, 3, 6, 9, 12, 15)], pred.loess,
                     "kappa", results5)[seq(2, 12, 2)])
# xlim <- range(results5$psel)
xlim <- c(0, 500)
plot(pred.loess("gren2", "kappa", results5), ylim=ylim, xlim=xlim,
     main="a)", xlab="Number of selected features", ylab="Cohen's kappa",
     cex.axis=1.5, cex.lab=2, cex.main=2, lwd=1.5, col=colors[1], type="l")
lines(pred.loess("enet2", "kappa", results5), lwd=1.5, col=colors[1], lty=2)
lines(pred.loess("sglasso2", "kappa", results5), lwd=1.5, col=colors[2], lty=1)
lines(pred.loess("cmcp2", "kappa", results5), lwd=1.5, col=colors[3], lty=1)
lines(pred.loess("gelasso2", "kappa", results5), lwd=1.5, col=colors[4], lty=1)
lines(pred.loess("grridge", "kappa", results5), lwd=1.5, col=colors[5], lty=1)

# mse
ylim <- range(sapply(names(results5$mse)[c(2, 4, 7, 10, 13, 16)], pred.loess,
                     "mse", results5)[seq(2, 12, 2)], mean(results5$mse$ridge))
xlim <- c(0, 500)
plot(pred.loess("gren2", "mse", results5), ylim=ylim, xlim=xlim,
     main="b)", xlab="Number of selected features", ylab="MSE", cex.axis=1.5,
     cex.lab=2, cex.main=2, lwd=1.5, col=colors[1], type="l")
lines(pred.loess("enet2", "mse", results5), lwd=1.5, col=colors[1], lty=2)
lines(pred.loess("sglasso2", "mse", results5), lwd=1.5, col=colors[2], lty=1)
lines(pred.loess("cmcp2", "mse", results5), lwd=1.5, col=colors[3], lty=1)
lines(pred.loess("gelasso2", "mse", results5), lwd=1.5, col=colors[4], lty=1)
lines(pred.loess("grridge", "mse", results5), lwd=1.5, col=colors[5], lty=1)
lines(xlim, rep(mean(results5$mse$ridge), 2), lwd=1.5,
      col=colors[5], lty=2)

# auc
ylim <- range(sapply(names(results5$auc)[c(2, 4, 7, 10, 13, 16)], pred.loess,
                     "auc", results5)[seq(2, 12, 2)], mean(results5$auc$ridge))
xlim <- c(0, 500)
plot(pred.loess("gren2", "auc", results5), ylim=ylim, xlim=xlim,
     main="c)", xlab="Number of selected features", ylab="AUC", cex.axis=1.5,
     cex.lab=2, cex.main=2, lwd=1.5, col=colors[1], type="l")
lines(pred.loess("enet2", "auc", results5), lwd=1.5, col=colors[1], lty=2)
lines(pred.loess("sglasso2", "auc", results5), lwd=1.5, col=colors[2], lty=1)
lines(pred.loess("cmcp2", "auc", results5), lwd=1.5, col=colors[3], lty=1)
lines(pred.loess("gelasso2", "auc", results5), lwd=1.5, col=colors[4], lty=1)
lines(pred.loess("grridge", "auc", results5), lwd=1.5, col=colors[5], lty=1)
lines(xlim, rep(mean(results5$auc$ridge), 2), lwd=1.5,
      col=colors[5], lty=2)

# legend
leglabels <- c("enet", "sglasso", "cMCP", "gelasso", "ridge",
               "group-regularized", "regular")
legend("bottomright", merge=TRUE, seg.len=1, cex=1.3, fill=c(colors, 0, 0),
       lty=c(rep(NA, length(colors)), 1, 2),
       lwd=c(rep(NA, length(colors)), 1.5, 1.5),
       border=c(rep(1, length(colors)), 0, 0), legend=leglabels)

# briers
ylim <- range(sapply(names(results5$briers)[c(2, 4, 7, 10, 13, 16)], pred.loess,
                     "briers", results5)[seq(2, 12, 2)],
              mean(results5$briers$ridge))
xlim <- c(0, 500)
plot(pred.loess("gren2", "briers", results5), ylim=ylim, xlim=xlim,
     main="d)", xlab="Number of selected features", ylab="Brier skill score",
     cex.axis=1.5, cex.lab=2, cex.main=2, lwd=1.5, col=colors[1], type="l")
lines(pred.loess("enet2", "briers", results5), lwd=1.5, col=colors[1], lty=2)
lines(pred.loess("sglasso2", "briers", results5), lwd=1.5, col=colors[2], lty=1)
lines(pred.loess("cmcp2", "briers", results5), lwd=1.5, col=colors[3], lty=1)
lines(pred.loess("gelasso2", "briers", results5), lwd=1.5, col=colors[4], lty=1)
lines(pred.loess("grridge", "briers", results5), lwd=1.5, col=colors[5], lty=1)
lines(xlim, rep(mean(results5$briers$ridge), 2), lwd=1.5,
      col=colors[5], lty=2)
dev.off()



### performance for alpha=0.95
colors <- bpy.colors(7)[-c(1, 7)]
png(paste(path.graph, "gren_sim5_res1_performance1_set3.png", sep=""),
    units="in", width=14, height=10, res=120)
par(mfrow=c(2, 2), mar=c(5.1, 5.1, 4.1, 2.1))

# kappa
ylim <- range(sapply(names(results5$kappa)[c(1, 4, 7, 10, 13, 16)], pred.loess,
                     "kappa", results5)[seq(2, 12, 2)])
# xlim <- range(results5$psel)
xlim <- c(0, 500)
plot(pred.loess("gren3", "kappa", results5), ylim=ylim, xlim=xlim,
     main="a)", xlab="Number of selected features", ylab="Cohen's kappa",
     cex.axis=1.5, cex.lab=2, cex.main=2, lwd=1.5, col=colors[1], type="l")
lines(pred.loess("enet3", "kappa", results5), lwd=1.5, col=colors[1], lty=2)
lines(pred.loess("sglasso3", "kappa", results5), lwd=1.5, col=colors[2], lty=1)
lines(pred.loess("cmcp3", "kappa", results5), lwd=1.5, col=colors[3], lty=1)
lines(pred.loess("gelasso3", "kappa", results5), lwd=1.5, col=colors[4], lty=1)
lines(pred.loess("grridge", "kappa", results5), lwd=1.5, col=colors[5], lty=1)

# mse
ylim <- range(sapply(names(results5$mse)[c(2, 5, 8, 11, 14, 17)], pred.loess,
                     "mse", results5)[seq(2, 12, 2)], mean(results5$mse$ridge))
xlim <- c(0, 500)
plot(pred.loess("gren3", "mse", results5), ylim=ylim, xlim=xlim,
     main="b)", xlab="Number of selected features", ylab="MSE", cex.axis=1.5,
     cex.lab=2, cex.main=2, lwd=1.5, col=colors[1], type="l")
lines(pred.loess("enet3", "mse", results5), lwd=1.5, col=colors[1], lty=2)
lines(pred.loess("sglasso3", "mse", results5), lwd=1.5, col=colors[2], lty=1)
lines(pred.loess("cmcp3", "mse", results5), lwd=1.5, col=colors[3], lty=1)
lines(pred.loess("gelasso3", "mse", results5), lwd=1.5, col=colors[4], lty=1)
lines(pred.loess("grridge", "mse", results5), lwd=1.5, col=colors[5], lty=1)
lines(xlim, rep(mean(results5$mse$ridge), 2), lwd=1.5,
      col=colors[5], lty=2)

# auc
ylim <- range(sapply(names(results5$auc)[c(2, 5, 8, 11, 14, 17)], pred.loess,
                     "auc", results5)[seq(2, 12, 2)], mean(results5$auc$ridge))
xlim <- c(0, 500)
plot(pred.loess("gren3", "auc", results5), ylim=ylim, xlim=xlim,
     main="c)", xlab="Number of selected features", ylab="AUC", cex.axis=1.5,
     cex.lab=2, cex.main=2, lwd=1.5, col=colors[1], type="l")
lines(pred.loess("enet3", "auc", results5), lwd=1.5, col=colors[1], lty=2)
lines(pred.loess("sglasso3", "auc", results5), lwd=1.5, col=colors[2], lty=1)
lines(pred.loess("cmcp3", "auc", results5), lwd=1.5, col=colors[3], lty=1)
lines(pred.loess("gelasso3", "auc", results5), lwd=1.5, col=colors[4], lty=1)
lines(pred.loess("grridge", "auc", results5), lwd=1.5, col=colors[5], lty=1)
lines(xlim, rep(mean(results5$auc$ridge), 2), lwd=1.5,
      col=colors[5], lty=2)

# legend
leglabels <- c("enet", "sglasso", "cMCP", "gelasso", "ridge",
               "group-regularized", "regular")
legend("bottomright", merge=TRUE, seg.len=1, cex=1.3, fill=c(colors, 0, 0),
       lty=c(rep(NA, length(colors)), 1, 2),
       lwd=c(rep(NA, length(colors)), 1.5, 1.5),
       border=c(rep(1, length(colors)), 0, 0), legend=leglabels)

# briers
ylim <- range(sapply(names(results5$briers)[c(2, 5, 8, 11, 14, 17)], pred.loess,
                     "briers", results5)[seq(2, 12, 2)],
              mean(results5$briers$ridge))
xlim <- c(0, 500)
plot(pred.loess("gren3", "briers", results5), ylim=ylim, xlim=xlim,
     main="d)", xlab="Number of selected features", ylab="Brier skill score",
     cex.axis=1.5, cex.lab=2, cex.main=2, lwd=1.5, col=colors[1], type="l")
lines(pred.loess("enet3", "briers", results5), lwd=1.5, col=colors[1], lty=2)
lines(pred.loess("sglasso3", "briers", results5), lwd=1.5, col=colors[2], lty=1)
lines(pred.loess("cmcp3", "briers", results5), lwd=1.5, col=colors[3], lty=1)
lines(pred.loess("gelasso3", "briers", results5), lwd=1.5, col=colors[4], lty=1)
lines(pred.loess("grridge", "briers", results5), lwd=1.5, col=colors[5], lty=1)
lines(xlim, rep(mean(results5$briers$ridge), 2), lwd=1.5,
      col=colors[5], lty=2)
dev.off()



