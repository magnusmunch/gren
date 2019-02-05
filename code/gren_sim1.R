path.res <- ifelse(as.character(Sys.info()[1])!="Darwin", "~/EBEN/results/",
                   "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/results/")
path.graph <- "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/graphs/"

library(mvtnorm)
library(gren)
library(pROC)
library(GRridge)
library(grpreg)
library(SGL)
library(irr)

################################# scenario 1
n <- 100
ntest <- 1000
p <- 1000
G <- 5
rho <- 0.5
rg <- 50
beta.mean <- 0.07
f <- 1.6

Sigma <- diag(p)
for(i in 1:p) {
  for(j in 1:p) {
    Sigma[i, j] <- rho^abs(i - j)
  }
}
b1 <- 2*p*beta.mean/sum(rg*f^(c(1:G) - 1))
beta <- numeric(p)

part1 <- rep(c(1:G), each=p/G)
csel <- c(seq(1, 10, 2), seq(15, 50, 5), seq(60, 140, 10), seq(160, 200, 20))
nreps <- 100
psel1 <- kappa1 <- list(grridge=matrix(NA, nrow=nreps, ncol=length(csel)),
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
auc1 <- mse1 <- briers1 <- list(ridge=numeric(nreps),
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
lambdagest1 <- list(grridge=matrix(NA, nrow=nreps, ncol=G),
                    gren1=matrix(NA, nrow=nreps, ncol=G),
                    gren2=matrix(NA, nrow=nreps, ncol=G),
                    gren3=matrix(NA, nrow=nreps, ncol=G))
for(r in 1:nreps) {
  set.seed(2018 + r - 1)
  print(paste("rep", r))
  for(g in 1:G) {
    idz <- c(((g-1)*p/G + 1):(g*p/G - rg))
    ida <- c(((g*p/G - rg) + 1):(g*p/G))
    beta[idz] <- 0
    beta[ida] <- runif(rg, 0, b1*f^(g-1))
  }
  x <- rmvnorm(n, rep(0, p), Sigma)
  y <- rbinom(n, 1, as.numeric(1/(1 + exp(-x %*% beta))))

  xtest <- rmvnorm(ntest, rep(0, p), Sigma)
  ytest <- rbinom(ntest, 1, as.numeric(1/(1 + exp(-xtest %*% beta))))

  fit1.ridge <- cv.glmnet(x, y, alpha=0, standardize=FALSE)

  fit1.grridge <- vector("list", length(csel))
  invisible(capture.output(
    fit1.grridge[[1]] <- grridge(t(x), y, partitions=list(
      groups=CreatePartition(as.factor(part1))), selection=TRUE, maxsel=csel[1],
      trace=FALSE, standardizeX=FALSE)))
  for(s in 2:length(csel)) {
    invisible(capture.output(
      fit1.grridge[[s]] <- grridge(t(x), y, partitions=list(
        groups=CreatePartition(as.factor(part1))), selection=TRUE,
        maxsel=csel[s], optl=fit1.grridge[[1]]$optl, trace=FALSE,
        standardizeX=FALSE)))
  }

  fit1.gren1 <- gren(x, y, partitions=list(groups=part1), alpha=0.05,
                     trace=FALSE)
  fit1.gren2 <- gren(x, y, partitions=list(groups=part1), alpha=0.5,
                     trace=FALSE)
  fit1.gren3 <- gren(x, y, partitions=list(groups=part1), alpha=0.95,
                     trace=FALSE)

  fit1.sglasso1 <- SGL(list(x=x, y=y), part1, type="logit", alpha=0.05,
                       standardize=FALSE, nlam=100)
  fit1.sglasso2 <- SGL(list(x=x, y=y), part1, type="logit", alpha=0.5,
                       standardize=FALSE, nlam=100)
  fit1.sglasso3 <- SGL(list(x=x, y=y), part1, type="logit", alpha=0.95,
                       standardize=FALSE, nlam=100)

  fit1.cmcp1 <- grpreg(x, y, part1, penalty="cMCP", alpha=0.05)
  fit1.cmcp2 <- grpreg(x, y, part1, penalty="cMCP", alpha=0.5)
  fit1.cmcp3 <- grpreg(x, y, part1, penalty="cMCP", alpha=0.95)

  fit1.gelasso1 <- grpreg(x, y, part1, penalty="gel", alpha=0.05)
  fit1.gelasso2 <- grpreg(x, y, part1, penalty="gel", alpha=0.5)
  fit1.gelasso3 <- grpreg(x, y, part1, penalty="gel", alpha=0.95)

  pred1.ridge <- as.numeric(predict(fit1.ridge, xtest, "lambda.min"))

  pred1.grridge <- sapply(fit1.grridge, function(s) {
    predict.grridge(s, t(xtest))[, 3]})

  pred1.gren1 <- predict(fit1.gren1, xtest, type="groupreg",
                         s=fit1.gren1$freq.model$groupreg$lambda)
  pred1.gren2 <- predict(fit1.gren2, xtest, type="groupreg",
                         s=fit1.gren2$freq.model$groupreg$lambda)
  pred1.gren3 <- predict(fit1.gren3, xtest, type="groupreg",
                         s=fit1.gren3$freq.model$groupreg$lambda)

  pred1.enet1 <- predict(fit1.gren1, xtest, type="regular",
                         s=fit1.gren1$freq.model$regular$lambda)
  pred1.enet2 <- predict(fit1.gren2, xtest, type="regular",
                         s=fit1.gren2$freq.model$regular$lambda)
  pred1.enet3 <- predict(fit1.gren3, xtest, type="regular",
                         s=fit1.gren3$freq.model$regular$lambda)

  pred1.sglasso1 <- predictSGL(fit1.sglasso1, xtest)
  pred1.sglasso2 <- predictSGL(fit1.sglasso2, xtest)
  pred1.sglasso3 <- predictSGL(fit1.sglasso3, xtest)

  pred1.cmcp1 <- predict(fit1.cmcp1, xtest)
  pred1.cmcp2 <- predict(fit1.cmcp2, xtest)
  pred1.cmcp3 <- predict(fit1.cmcp3, xtest)

  pred1.gelasso1 <- predict(fit1.gelasso1, xtest)
  pred1.gelasso2 <- predict(fit1.gelasso2, xtest)
  pred1.gelasso3 <- predict(fit1.gelasso3, xtest)

  auc1$ridge[r] <- pROC::auc(ytest, pred1.ridge)

  auc1$grridge[r, ] <- apply(pred1.grridge, 2, function(pred) {
    pROC::auc(ytest, pred)})

  auc1$gren1[r, ] <- apply(pred1.gren1, 2, function(pred) {
    pROC::auc(ytest, pred)})
  auc1$gren2[r, ] <- apply(pred1.gren2, 2, function(pred) {
    pROC::auc(ytest, pred)})
  auc1$gren3[r, ] <- apply(pred1.gren3, 2, function(pred) {
    pROC::auc(ytest, pred)})

  auc1$enet1[r, ] <- apply(pred1.enet1, 2, function(pred) {
    pROC::auc(ytest, pred)})
  auc1$enet2[r, ] <- apply(pred1.enet2, 2, function(pred) {
    pROC::auc(ytest, pred)})
  auc1$enet3[r, ] <- apply(pred1.enet3, 2, function(pred) {
    pROC::auc(ytest, pred)})

  auc1$sglasso1[r, ] <- apply(pred1.sglasso1, 2, function(pred) {
    pROC::auc(ytest, pred)})
  auc1$sglasso2[r, ] <- apply(pred1.sglasso2, 2, function(pred) {
    pROC::auc(ytest, pred)})
  auc1$sglasso3[r, ] <- apply(pred1.sglasso3, 2, function(pred) {
    pROC::auc(ytest, pred)})

  auc1$cmcp1[r, ] <- apply(pred1.cmcp1, 2, function(pred) {
    pROC::auc(ytest, pred)})
  auc1$cmcp2[r, ] <- apply(pred1.cmcp2, 2, function(pred) {
    pROC::auc(ytest, pred)})
  auc1$cmcp3[r, ] <- apply(pred1.cmcp3, 2, function(pred) {
    pROC::auc(ytest, pred)})

  auc1$gelasso1[r, ] <- apply(pred1.gelasso1, 2, function(pred) {
    pROC::auc(ytest, pred)})
  auc1$gelasso2[r, ] <- apply(pred1.gelasso2, 2, function(pred) {
    pROC::auc(ytest, pred)})
  auc1$gelasso3[r, ] <- apply(pred1.gelasso3, 2, function(pred) {
    pROC::auc(ytest, pred)})

  const <- sum((ytest - mean(ytest))^2)
  briers1$ridge[r] <- 1 - sum((ytest - pred1.ridge)^2)/const

  briers1$grridge[r, ] <- apply(pred1.grridge, 2, function(pred) {
    1 - sum((ytest - pred)^2)/const})

  briers1$gren1[r, ] <- apply(pred1.gren1, 2, function(pred) {
    1 - sum((ytest - pred)^2)/const})
  briers1$gren2[r, ] <- apply(pred1.gren2, 2, function(pred) {
    1 - sum((ytest - pred)^2)/const})
  briers1$gren3[r, ] <- apply(pred1.gren3, 2, function(pred) {
    1 - sum((ytest - pred)^2)/const})

  briers1$enet1[r, ] <- apply(pred1.enet1, 2, function(pred) {
    1 - sum((ytest - pred)^2)/const})
  briers1$enet2[r, ] <- apply(pred1.enet2, 2, function(pred) {
    1 - sum((ytest - pred)^2)/const})
  briers1$enet3[r, ] <- apply(pred1.enet3, 2, function(pred) {
    1 - sum((ytest - pred)^2)/const})

  briers1$sglasso1[r, ] <- apply(pred1.sglasso1, 2, function(pred) {
    1 - sum((ytest - pred)^2)/const})
  briers1$sglasso2[r, ] <- apply(pred1.sglasso2, 2, function(pred) {
    1 - sum((ytest - pred)^2)/const})
  briers1$sglasso3[r, ] <- apply(pred1.sglasso3, 2, function(pred) {
    1 - sum((ytest - pred)^2)/const})

  briers1$cmcp1[r, ] <- apply(pred1.cmcp1, 2, function(pred) {
    1 - sum((ytest - pred)^2)/const})
  briers1$cmcp2[r, ] <- apply(pred1.cmcp2, 2, function(pred) {
    1 - sum((ytest - pred)^2)/const})
  briers1$cmcp3[r, ] <- apply(pred1.cmcp3, 2, function(pred) {
    1 - sum((ytest - pred)^2)/const})

  briers1$gelasso1[r, ] <- apply(pred1.gelasso1, 2, function(pred) {
    1 - sum((ytest - pred)^2)/const})
  briers1$gelasso2[r, ] <- apply(pred1.gelasso2, 2, function(pred) {
    1 - sum((ytest - pred)^2)/const})
  briers1$gelasso3[r, ] <- apply(pred1.gelasso3, 2, function(pred) {
    1 - sum((ytest - pred)^2)/const})

  mse1$ridge[r] <- mean((coef(fit1.ridge)[-1] - beta)^2)

  mse1$grridge[r, ] <- sapply(fit1.grridge, function(s) {
    mean((replace(rep(0, p), s$resEN$whichEN, s$resEN$betasEN) - beta)^2)})

  mse1$gren1[r, ] <- apply(fit1.gren1$freq.model$groupreg$beta, 2, function(b) {
    mean((b - beta)^2)})
  mse1$gren2[r, ] <- apply(fit1.gren2$freq.model$groupreg$beta, 2, function(b) {
    mean((b - beta)^2)})
  mse1$gren3[r, ] <- apply(fit1.gren3$freq.model$groupreg$beta, 2, function(b) {
    mean((b - beta)^2)})

  mse1$enet1[r, ] <- apply(fit1.gren1$freq.model$regular$beta, 2, function(b) {
    mean((b - beta)^2)})
  mse1$enet2[r, ] <- apply(fit1.gren2$freq.model$regular$beta, 2, function(b) {
    mean((b - beta)^2)})
  mse1$enet3[r, ] <- apply(fit1.gren3$freq.model$regular$beta, 2, function(b) {
    mean((b - beta)^2)})

  mse1$sglasso1[r, ] <- apply(fit1.sglasso1$beta, 2, function(b) {
    mean((b - beta)^2)})
  mse1$sglasso2[r, ] <- apply(fit1.sglasso2$beta, 2, function(b) {
    mean((b - beta)^2)})
  mse1$sglasso3[r, ] <- apply(fit1.sglasso3$beta, 2, function(b) {
    mean((b - beta)^2)})

  mse1$cmcp1[r, ] <- apply(fit1.cmcp1$beta[-1, ], 2, function(b) {
    mean((b - beta)^2)})
  mse1$cmcp2[r, ] <- apply(fit1.cmcp2$beta[-1, ], 2, function(b) {
    mean((b - beta)^2)})
  mse1$cmcp3[r, ] <- apply(fit1.cmcp3$beta[-1, ], 2, function(b) {
    mean((b - beta)^2)})

  mse1$gelasso1[r, ] <- apply(fit1.gelasso1$beta[-1, ], 2, function(b) {
    mean((b - beta)^2)})
  mse1$gelasso2[r, ] <- apply(fit1.gelasso2$beta[-1, ], 2, function(b) {
    mean((b - beta)^2)})
  mse1$gelasso3[r, ] <- apply(fit1.gelasso3$beta[-1, ], 2, function(b) {
    mean((b - beta)^2)})

  kappa1$grridge[r, ] <- sapply(fit1.grridge, function(s) {
    kappa2(cbind(beta!=0, replace(rep(FALSE, p), s$resEN$whichEN,
                                  TRUE)))$value})

  kappa1$gren1[r, ] <- apply(fit1.gren1$freq.model$groupreg$beta, 2,
                             function(b) {kappa2(cbind(beta!=0, b!=0))$value})
  kappa1$gren2[r, ] <- apply(fit1.gren2$freq.model$groupreg$beta, 2,
                             function(b) {kappa2(cbind(beta!=0, b!=0))$value})
  kappa1$gren3[r, ] <- apply(fit1.gren3$freq.model$groupreg$beta, 2,
                             function(b) {kappa2(cbind(beta!=0, b!=0))$value})

  kappa1$enet1[r, ] <- apply(fit1.gren1$freq.model$regular$beta, 2,
                             function(b) {kappa2(cbind(beta!=0, b!=0))$value})
  kappa1$enet2[r, ] <- apply(fit1.gren2$freq.model$regular$beta, 2,
                             function(b) {kappa2(cbind(beta!=0, b!=0))$value})
  kappa1$enet3[r, ] <- apply(fit1.gren3$freq.model$regular$beta, 2,
                             function(b) {kappa2(cbind(beta!=0, b!=0))$value})

  kappa1$sglasso1[r, ] <- apply(fit1.sglasso1$beta, 2, function(b) {
    kappa2(cbind(beta!=0, b!=0))$value})
  kappa1$sglasso2[r, ] <- apply(fit1.sglasso2$beta, 2, function(b) {
    kappa2(cbind(beta!=0, b!=0))$value})
  kappa1$sglasso3[r, ] <- apply(fit1.sglasso3$beta, 2, function(b) {
    kappa2(cbind(beta!=0, b!=0))$value})

  kappa1$cmcp1[r, ] <- apply(fit1.cmcp1$beta[-1, ], 2, function(b) {
    kappa2(cbind(beta!=0, b!=0))$value})
  kappa1$cmcp2[r, ] <- apply(fit1.cmcp2$beta[-1, ], 2, function(b) {
    kappa2(cbind(beta!=0, b!=0))$value})
  kappa1$cmcp3[r, ] <- apply(fit1.cmcp3$beta[-1, ], 2, function(b) {
    kappa2(cbind(beta!=0, b!=0))$value})

  kappa1$gelasso1[r, ] <- apply(fit1.gelasso1$beta[-1, ], 2, function(b) {
    kappa2(cbind(beta!=0, b!=0))$value})
  kappa1$gelasso2[r, ] <- apply(fit1.gelasso2$beta[-1, ], 2, function(b) {
    kappa2(cbind(beta!=0, b!=0))$value})
  kappa1$gelasso3[r, ] <- apply(fit1.gelasso3$beta[-1, ], 2, function(b) {
    kappa2(cbind(beta!=0, b!=0))$value})

  psel1$grridge[r, ] <- sapply(fit1.grridge, function(s) {
    length(s$resEN$whichEN)})

  psel1$gren1[r, ] <- fit1.gren1$freq.model$groupreg$df
  psel1$gren2[r, ] <- fit1.gren2$freq.model$groupreg$df
  psel1$gren3[r, ] <- fit1.gren3$freq.model$groupreg$df

  psel1$enet1[r, ] <- fit1.gren1$freq.model$regular$df
  psel1$enet2[r, ] <- fit1.gren2$freq.model$regular$df
  psel1$enet3[r, ] <- fit1.gren3$freq.model$regular$df

  psel1$sglasso1[r, ] <- apply(fit1.sglasso1$beta, 2, function(b) {sum(b!=0)})
  psel1$sglasso2[r, ] <- apply(fit1.sglasso2$beta, 2, function(b) {sum(b!=0)})
  psel1$sglasso3[r, ] <- apply(fit1.sglasso3$beta, 2, function(b) {sum(b!=0)})

  psel1$cmcp1[r, ] <- apply(fit1.cmcp1$beta, 2, function(b) {sum(b!=0)})
  psel1$cmcp2[r, ] <- apply(fit1.cmcp2$beta, 2, function(b) {sum(b!=0)})
  psel1$cmcp3[r, ] <- apply(fit1.cmcp3$beta, 2, function(b) {sum(b!=0)})

  psel1$gelasso1[r, ] <- apply(fit1.gelasso1$beta, 2, function(b) {sum(b!=0)})
  psel1$gelasso2[r, ] <- apply(fit1.gelasso2$beta, 2, function(b) {sum(b!=0)})
  psel1$gelasso3[r, ] <- apply(fit1.gelasso3$beta, 2, function(b) {sum(b!=0)})

  lambdagest1$grridge[r, ] <- fit1.grridge[[1]]$lambdamults$groups
  lambdagest1$gren1[r, ] <- fit1.gren1$lambdag$groups
  lambdagest1$gren2[r, ] <- fit1.gren2$lambdag$groups
  lambdagest1$gren3[r, ] <- fit1.gren3$lambdag$groups

  results1 <- list(auc=auc1, briers=briers1, mse=mse1, kappa=kappa1, psel=psel1,
                   lambdag=lambdagest1)
  save(results1, file=paste(path.res, "gren_sim1_res1.Rdata", sep=""))

}







path.res <- ifelse(as.character(Sys.info()[1])!="Darwin", "~/EBEN/results/",
                   "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/results/")
load(paste(path.res, "gren_sim1_res1.Rdata", sep=""))
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
png(paste(path.graph, "gren_sim1_res1_performance1.png", sep=""),
    units="in", width=14, height=10, res=120)
par(mfrow=c(2, 2), mar=c(5.1, 5.1, 4.1, 2.1))

# kappa
ylim <- range(sapply(names(results1$kappa)[-1], pred.loess, "kappa", results1)[
  seq(2, 2*(length(results1$kappa) - 1), 2)])
# xlim <- range(results1$psel)
xlim <- c(0, 250)
plot(pred.loess("gren1", "kappa", results1), ylim=ylim, xlim=xlim, 
     main="a)", xlab="Number of selected features", ylab="Cohen's kappa", 
     cex.axis=1.5, cex.lab=2, cex.main=2, lwd=1.5, col=colors[1], type="l")
lines(pred.loess("gren2", "kappa", results1), lwd=1.5, col=colors[1], lty=2)
lines(pred.loess("gren3", "kappa", results1), lwd=1.5, col=colors[1], lty=3)

lines(pred.loess("enet1", "kappa", results1), lwd=1.5, col=colors[2], lty=1)
lines(pred.loess("enet2", "kappa", results1), lwd=1.5, col=colors[2], lty=2)
lines(pred.loess("enet3", "kappa", results1), lwd=1.5, col=colors[2], lty=3)

lines(pred.loess("sglasso1", "kappa", results1), lwd=1.5, col=colors[3], lty=1)
lines(pred.loess("sglasso2", "kappa", results1), lwd=1.5, col=colors[3], lty=2)
lines(pred.loess("sglasso3", "kappa", results1), lwd=1.5, col=colors[3], lty=3)

lines(pred.loess("cmcp1", "kappa", results1), lwd=1.5, col=colors[4], lty=1)
lines(pred.loess("cmcp2", "kappa", results1), lwd=1.5, col=colors[4], lty=2)
lines(pred.loess("cmcp3", "kappa", results1), lwd=1.5, col=colors[4], lty=3)

lines(pred.loess("gelasso1", "kappa", results1), lwd=1.5, col=colors[5], lty=1)
lines(pred.loess("gelasso2", "kappa", results1), lwd=1.5, col=colors[5], lty=2)
lines(pred.loess("gelasso3", "kappa", results1), lwd=1.5, col=colors[5], lty=3)

lines(pred.loess("grridge", "kappa", results1), lwd=1.5, col=colors[6], lty=1)

# mse
ylim <- range(sapply(names(results1$mse)[-1], pred.loess, "mse", results1)[
  seq(2, 2*(length(results1$mse) - 1), 2)], mean(results1$mse$ridge))
xlim <- c(0, 250)
plot(pred.loess("gren1", "mse", results1), ylim=ylim, xlim=xlim, 
     main="b)", xlab="Number of selected features", ylab="MSE", cex.axis=1.5, 
     cex.lab=2, cex.main=2, lwd=1.5, col=colors[1], type="l")
lines(pred.loess("gren2", "mse", results1), lwd=1.5, col=colors[1], lty=2)
lines(pred.loess("gren3", "mse", results1), lwd=1.5, col=colors[1], lty=3)

lines(pred.loess("enet1", "mse", results1), lwd=1.5, col=colors[2], lty=1)
lines(pred.loess("enet2", "mse", results1), lwd=1.5, col=colors[2], lty=2)
lines(pred.loess("enet3", "mse", results1), lwd=1.5, col=colors[2], lty=3)

lines(pred.loess("sglasso1", "mse", results1), lwd=1.5, col=colors[3], lty=1)
lines(pred.loess("sglasso2", "mse", results1), lwd=1.5, col=colors[3], lty=2)
lines(pred.loess("sglasso3", "mse", results1), lwd=1.5, col=colors[3], lty=3)

lines(pred.loess("cmcp1", "mse", results1), lwd=1.5, col=colors[4], lty=1)
lines(pred.loess("cmcp2", "mse", results1), lwd=1.5, col=colors[4], lty=2)
lines(pred.loess("cmcp3", "mse", results1), lwd=1.5, col=colors[4], lty=3)

lines(pred.loess("gelasso1", "mse", results1), lwd=1.5, col=colors[5], lty=1)
lines(pred.loess("gelasso2", "mse", results1), lwd=1.5, col=colors[5], lty=2)
lines(pred.loess("gelasso3", "mse", results1), lwd=1.5, col=colors[5], lty=3)

lines(pred.loess("grridge", "mse", results1), lwd=1.5, col=colors[6], lty=1)
lines(xlim, rep(mean(results1$mse$ridge), 2), lwd=1.5, 
      col=colors[6], lty=2)

# auc
ylim <- range(sapply(names(results1$auc)[-1], pred.loess, "auc", results1)[
  seq(2, 2*(length(results1$auc) - 1), 2)], mean(results1$auc$ridge))
xlim <- c(0, 250)
plot(pred.loess("gren1", "auc", results1), ylim=ylim, xlim=xlim, 
     main="c)", xlab="Number of selected features", ylab="AUC", cex.axis=1.5, 
     cex.lab=2, cex.main=2, lwd=1.5, col=colors[1], type="l")
lines(pred.loess("gren2", "auc", results1), lwd=1.5, col=colors[1], lty=2)
lines(pred.loess("gren3", "auc", results1), lwd=1.5, col=colors[1], lty=3)

lines(pred.loess("enet1", "auc", results1), lwd=1.5, col=colors[2], lty=1)
lines(pred.loess("enet2", "auc", results1), lwd=1.5, col=colors[2], lty=2)
lines(pred.loess("enet3", "auc", results1), lwd=1.5, col=colors[2], lty=3)

lines(pred.loess("sglasso1", "auc", results1), lwd=1.5, col=colors[3], lty=1)
lines(pred.loess("sglasso2", "auc", results1), lwd=1.5, col=colors[3], lty=2)
lines(pred.loess("sglasso3", "auc", results1), lwd=1.5, col=colors[3], lty=3)

lines(pred.loess("cmcp1", "auc", results1), lwd=1.5, col=colors[4], lty=1)
lines(pred.loess("cmcp2", "auc", results1), lwd=1.5, col=colors[4], lty=2)
lines(pred.loess("cmcp3", "auc", results1), lwd=1.5, col=colors[4], lty=3)

lines(pred.loess("gelasso1", "auc", results1), lwd=1.5, col=colors[5], lty=1)
lines(pred.loess("gelasso2", "auc", results1), lwd=1.5, col=colors[5], lty=2)
lines(pred.loess("gelasso3", "auc", results1), lwd=1.5, col=colors[5], lty=3)

lines(pred.loess("grridge", "auc", results1), lwd=1.5, col=colors[6], lty=1)
lines(xlim, rep(mean(results1$auc$ridge), 2), lwd=1.5, 
      col=colors[6], lty=2)

# briers
ylim <- range(sapply(names(results1$briers)[-1], pred.loess, "briers", results1)[
  seq(2, 2*(length(results1$briers) - 1), 2)], mean(results1$briers$ridge))
xlim <- c(0, 250)
plot(pred.loess("gren1", "briers", results1), ylim=ylim, xlim=xlim, 
     main="d)", xlab="Number of selected features", ylab="Brier skill score", 
     cex.axis=1.5, cex.lab=2, cex.main=2, lwd=1.5, col=colors[1], type="l")
lines(pred.loess("gren2", "briers", results1), lwd=1.5, col=colors[1], lty=2)
lines(pred.loess("gren3", "briers", results1), lwd=1.5, col=colors[1], lty=3)

lines(pred.loess("enet1", "briers", results1), lwd=1.5, col=colors[2], lty=1)
lines(pred.loess("enet2", "briers", results1), lwd=1.5, col=colors[2], lty=2)
lines(pred.loess("enet3", "briers", results1), lwd=1.5, col=colors[2], lty=3)

lines(pred.loess("sglasso1", "briers", results1), lwd=1.5, col=colors[3], lty=1)
lines(pred.loess("sglasso2", "briers", results1), lwd=1.5, col=colors[3], lty=2)
lines(pred.loess("sglasso3", "briers", results1), lwd=1.5, col=colors[3], lty=3)

lines(pred.loess("cmcp1", "briers", results1), lwd=1.5, col=colors[4], lty=1)
lines(pred.loess("cmcp2", "briers", results1), lwd=1.5, col=colors[4], lty=2)
lines(pred.loess("cmcp3", "briers", results1), lwd=1.5, col=colors[4], lty=3)

lines(pred.loess("gelasso1", "briers", results1), lwd=1.5, col=colors[5], lty=1)
lines(pred.loess("gelasso2", "briers", results1), lwd=1.5, col=colors[5], lty=2)
lines(pred.loess("gelasso3", "briers", results1), lwd=1.5, col=colors[5], lty=3)

lines(pred.loess("grridge", "briers", results1), lwd=1.5, col=colors[6], lty=1)
lines(xlim, rep(mean(results1$briers$ridge), 2), lwd=1.5, 
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
png(paste(path.graph, "gren_sim1_res1_performance1_set1.png", sep=""),
    units="in", width=14, height=10, res=120)
par(mfrow=c(2, 2), mar=c(5.1, 5.1, 4.1, 2.1))

# kappa
ylim <- range(sapply(names(results1$kappa)[c(1, 2, 5, 8, 11, 14)], pred.loess, 
                     "kappa", results1)[seq(2, 12, 2)])
# xlim <- range(results1$psel)
xlim <- c(0, 250)
plot(pred.loess("gren1", "kappa", results1), ylim=ylim, xlim=xlim, 
     main="a)", xlab="Number of selected features", ylab="Cohen's kappa", 
     cex.axis=1.5, cex.lab=2, cex.main=2, lwd=1.5, col=colors[1], type="l")
lines(pred.loess("enet1", "kappa", results1), lwd=1.5, col=colors[1], lty=2)
lines(pred.loess("sglasso1", "kappa", results1), lwd=1.5, col=colors[2], lty=1)
lines(pred.loess("cmcp1", "kappa", results1), lwd=1.5, col=colors[3], lty=1)
lines(pred.loess("gelasso1", "kappa", results1), lwd=1.5, col=colors[4], lty=1)
lines(pred.loess("grridge", "kappa", results1), lwd=1.5, col=colors[5], lty=1)

# mse
ylim <- range(sapply(names(results1$mse)[c(2, 3, 6, 9, 12, 15)], pred.loess, 
                     "mse", results1)[seq(2, 12, 2)], mean(results1$mse$ridge))
xlim <- c(0, 250)
plot(pred.loess("gren1", "mse", results1), ylim=ylim, xlim=xlim, 
     main="b)", xlab="Number of selected features", ylab="MSE", cex.axis=1.5, 
     cex.lab=2, cex.main=2, lwd=1.5, col=colors[1], type="l")
lines(pred.loess("enet1", "mse", results1), lwd=1.5, col=colors[1], lty=2)
lines(pred.loess("sglasso1", "mse", results1), lwd=1.5, col=colors[2], lty=1)
lines(pred.loess("cmcp1", "mse", results1), lwd=1.5, col=colors[3], lty=1)
lines(pred.loess("gelasso1", "mse", results1), lwd=1.5, col=colors[4], lty=1)
lines(pred.loess("grridge", "mse", results1), lwd=1.5, col=colors[5], lty=1)
lines(xlim, rep(mean(results1$mse$ridge), 2), lwd=1.5, 
      col=colors[5], lty=2)

# auc
ylim <- range(sapply(names(results1$auc)[c(2, 3, 6, 9, 12, 15)], pred.loess, 
                     "auc", results1)[seq(2, 12, 2)], mean(results1$auc$ridge))
xlim <- c(0, 250)
plot(pred.loess("gren1", "auc", results1), ylim=ylim, xlim=xlim, 
     main="c)", xlab="Number of selected features", ylab="AUC", cex.axis=1.5, 
     cex.lab=2, cex.main=2, lwd=1.5, col=colors[1], type="l")
lines(pred.loess("enet1", "auc", results1), lwd=1.5, col=colors[1], lty=2)
lines(pred.loess("sglasso1", "auc", results1), lwd=1.5, col=colors[2], lty=1)
lines(pred.loess("cmcp1", "auc", results1), lwd=1.5, col=colors[3], lty=1)
lines(pred.loess("gelasso1", "auc", results1), lwd=1.5, col=colors[4], lty=1)
lines(pred.loess("grridge", "auc", results1), lwd=1.5, col=colors[5], lty=1)
lines(xlim, rep(mean(results1$auc$ridge), 2), lwd=1.5, 
      col=colors[5], lty=2)

# legend
leglabels <- c("enet", "sglasso", "cMCP", "gelasso", "ridge",
               "group-regularized", "regular")
legend("bottomright", merge=TRUE, seg.len=1, cex=1.3, fill=c(colors, 0, 0), 
       lty=c(rep(NA, length(colors)), 1, 2), 
       lwd=c(rep(NA, length(colors)), 1.5, 1.5),
       border=c(rep(1, length(colors)), 0, 0), legend=leglabels)

# briers
ylim <- range(sapply(names(results1$briers)[c(2, 3, 6, 9, 12, 15)], pred.loess, 
                     "briers", results1)[seq(2, 12, 2)], 
              mean(results1$briers$ridge))
xlim <- c(0, 250)
plot(pred.loess("gren1", "briers", results1), ylim=ylim, xlim=xlim, 
     main="d)", xlab="Number of selected features", ylab="Brier skill score", 
     cex.axis=1.5, cex.lab=2, cex.main=2, lwd=1.5, col=colors[1], type="l")
lines(pred.loess("enet1", "briers", results1), lwd=1.5, col=colors[1], lty=2)
lines(pred.loess("sglasso1", "briers", results1), lwd=1.5, col=colors[2], lty=1)
lines(pred.loess("cmcp1", "briers", results1), lwd=1.5, col=colors[3], lty=1)
lines(pred.loess("gelasso1", "briers", results1), lwd=1.5, col=colors[4], lty=1)
lines(pred.loess("grridge", "briers", results1), lwd=1.5, col=colors[5], lty=1)
lines(xlim, rep(mean(results1$briers$ridge), 2), lwd=1.5, 
      col=colors[5], lty=2)
dev.off()



### performance for alpha=0.5
colors <- bpy.colors(7)[-c(1, 7)]
png(paste(path.graph, "gren_sim1_res1_performance1_set2.png", sep=""),
    units="in", width=14, height=10, res=120)
par(mfrow=c(2, 2), mar=c(5.1, 5.1, 4.1, 2.1))

# kappa
ylim <- range(sapply(names(results1$kappa)[c(1, 3, 6, 9, 12, 15)], pred.loess, 
                     "kappa", results1)[seq(2, 12, 2)])
# xlim <- range(results1$psel)
xlim <- c(0, 250)
plot(pred.loess("gren2", "kappa", results1), ylim=ylim, xlim=xlim, 
     main="a)", xlab="Number of selected features", ylab="Cohen's kappa", 
     cex.axis=1.5, cex.lab=2, cex.main=2, lwd=1.5, col=colors[1], type="l")
lines(pred.loess("enet2", "kappa", results1), lwd=1.5, col=colors[1], lty=2)
lines(pred.loess("sglasso2", "kappa", results1), lwd=1.5, col=colors[2], lty=1)
lines(pred.loess("cmcp2", "kappa", results1), lwd=1.5, col=colors[3], lty=1)
lines(pred.loess("gelasso2", "kappa", results1), lwd=1.5, col=colors[4], lty=1)
lines(pred.loess("grridge", "kappa", results1), lwd=1.5, col=colors[5], lty=1)

# mse
ylim <- range(sapply(names(results1$mse)[c(2, 4, 7, 10, 13, 16)], pred.loess, 
                     "mse", results1)[seq(2, 12, 2)], mean(results1$mse$ridge))
xlim <- c(0, 250)
plot(pred.loess("gren2", "mse", results1), ylim=ylim, xlim=xlim, 
     main="b)", xlab="Number of selected features", ylab="MSE", cex.axis=1.5, 
     cex.lab=2, cex.main=2, lwd=1.5, col=colors[1], type="l")
lines(pred.loess("enet2", "mse", results1), lwd=1.5, col=colors[1], lty=2)
lines(pred.loess("sglasso2", "mse", results1), lwd=1.5, col=colors[2], lty=1)
lines(pred.loess("cmcp2", "mse", results1), lwd=1.5, col=colors[3], lty=1)
lines(pred.loess("gelasso2", "mse", results1), lwd=1.5, col=colors[4], lty=1)
lines(pred.loess("grridge", "mse", results1), lwd=1.5, col=colors[5], lty=1)
lines(xlim, rep(mean(results1$mse$ridge), 2), lwd=1.5, 
      col=colors[5], lty=2)

# auc
ylim <- range(sapply(names(results1$auc)[c(2, 4, 7, 10, 13, 16)], pred.loess, 
                     "auc", results1)[seq(2, 12, 2)], mean(results1$auc$ridge))
xlim <- c(0, 250)
plot(pred.loess("gren2", "auc", results1), ylim=ylim, xlim=xlim, 
     main="c)", xlab="Number of selected features", ylab="AUC", cex.axis=1.5, 
     cex.lab=2, cex.main=2, lwd=1.5, col=colors[1], type="l")
lines(pred.loess("enet2", "auc", results1), lwd=1.5, col=colors[1], lty=2)
lines(pred.loess("sglasso2", "auc", results1), lwd=1.5, col=colors[2], lty=1)
lines(pred.loess("cmcp2", "auc", results1), lwd=1.5, col=colors[3], lty=1)
lines(pred.loess("gelasso2", "auc", results1), lwd=1.5, col=colors[4], lty=1)
lines(pred.loess("grridge", "auc", results1), lwd=1.5, col=colors[5], lty=1)
lines(xlim, rep(mean(results1$auc$ridge), 2), lwd=1.5, 
      col=colors[5], lty=2)

# legend
leglabels <- c("enet", "sglasso", "cMCP", "gelasso", "ridge",
               "group-regularized", "regular")
legend("bottomright", merge=TRUE, seg.len=1, cex=1.3, fill=c(colors, 0, 0), 
       lty=c(rep(NA, length(colors)), 1, 2), 
       lwd=c(rep(NA, length(colors)), 1.5, 1.5),
       border=c(rep(1, length(colors)), 0, 0), legend=leglabels)

# briers
ylim <- range(sapply(names(results1$briers)[c(2, 4, 7, 10, 13, 16)], pred.loess, 
                     "briers", results1)[seq(2, 12, 2)], 
              mean(results1$briers$ridge))
xlim <- c(0, 250)
plot(pred.loess("gren2", "briers", results1), ylim=ylim, xlim=xlim, 
     main="d)", xlab="Number of selected features", ylab="Brier skill score", 
     cex.axis=1.5, cex.lab=2, cex.main=2, lwd=1.5, col=colors[1], type="l")
lines(pred.loess("enet2", "briers", results1), lwd=1.5, col=colors[1], lty=2)
lines(pred.loess("sglasso2", "briers", results1), lwd=1.5, col=colors[2], lty=1)
lines(pred.loess("cmcp2", "briers", results1), lwd=1.5, col=colors[3], lty=1)
lines(pred.loess("gelasso2", "briers", results1), lwd=1.5, col=colors[4], lty=1)
lines(pred.loess("grridge", "briers", results1), lwd=1.5, col=colors[5], lty=1)
lines(xlim, rep(mean(results1$briers$ridge), 2), lwd=1.5, 
      col=colors[5], lty=2)
dev.off()



### performance for alpha=0.95
colors <- bpy.colors(7)[-c(1, 7)]
png(paste(path.graph, "gren_sim1_res1_performance1_set3.png", sep=""),
    units="in", width=14, height=10, res=120)
par(mfrow=c(2, 2), mar=c(5.1, 5.1, 4.1, 2.1))

# kappa
ylim <- range(sapply(names(results1$kappa)[c(1, 4, 7, 10, 13, 16)], pred.loess, 
                     "kappa", results1)[seq(2, 12, 2)])
# xlim <- range(results1$psel)
xlim <- c(0, 250)
plot(pred.loess("gren3", "kappa", results1), ylim=ylim, xlim=xlim, 
     main="a)", xlab="Number of selected features", ylab="Cohen's kappa", 
     cex.axis=1.5, cex.lab=2, cex.main=2, lwd=1.5, col=colors[1], type="l")
lines(pred.loess("enet3", "kappa", results1), lwd=1.5, col=colors[1], lty=2)
lines(pred.loess("sglasso3", "kappa", results1), lwd=1.5, col=colors[2], lty=1)
lines(pred.loess("cmcp3", "kappa", results1), lwd=1.5, col=colors[3], lty=1)
lines(pred.loess("gelasso3", "kappa", results1), lwd=1.5, col=colors[4], lty=1)
lines(pred.loess("grridge", "kappa", results1), lwd=1.5, col=colors[5], lty=1)

# mse
ylim <- range(sapply(names(results1$mse)[c(2, 5, 8, 11, 14, 17)], pred.loess, 
                     "mse", results1)[seq(2, 12, 2)], mean(results1$mse$ridge))
xlim <- c(0, 250)
plot(pred.loess("gren3", "mse", results1), ylim=ylim, xlim=xlim, 
     main="b)", xlab="Number of selected features", ylab="MSE", cex.axis=1.5, 
     cex.lab=2, cex.main=2, lwd=1.5, col=colors[1], type="l")
lines(pred.loess("enet3", "mse", results1), lwd=1.5, col=colors[1], lty=2)
lines(pred.loess("sglasso3", "mse", results1), lwd=1.5, col=colors[2], lty=1)
lines(pred.loess("cmcp3", "mse", results1), lwd=1.5, col=colors[3], lty=1)
lines(pred.loess("gelasso3", "mse", results1), lwd=1.5, col=colors[4], lty=1)
lines(pred.loess("grridge", "mse", results1), lwd=1.5, col=colors[5], lty=1)
lines(xlim, rep(mean(results1$mse$ridge), 2), lwd=1.5, 
      col=colors[5], lty=2)

# auc
ylim <- range(sapply(names(results1$auc)[c(2, 5, 8, 11, 14, 17)], pred.loess, 
                     "auc", results1)[seq(2, 12, 2)], mean(results1$auc$ridge))
xlim <- c(0, 250)
plot(pred.loess("gren3", "auc", results1), ylim=ylim, xlim=xlim, 
     main="c)", xlab="Number of selected features", ylab="AUC", cex.axis=1.5, 
     cex.lab=2, cex.main=2, lwd=1.5, col=colors[1], type="l")
lines(pred.loess("enet3", "auc", results1), lwd=1.5, col=colors[1], lty=2)
lines(pred.loess("sglasso3", "auc", results1), lwd=1.5, col=colors[2], lty=1)
lines(pred.loess("cmcp3", "auc", results1), lwd=1.5, col=colors[3], lty=1)
lines(pred.loess("gelasso3", "auc", results1), lwd=1.5, col=colors[4], lty=1)
lines(pred.loess("grridge", "auc", results1), lwd=1.5, col=colors[5], lty=1)
lines(xlim, rep(mean(results1$auc$ridge), 2), lwd=1.5, 
      col=colors[5], lty=2)

# legend
leglabels <- c("enet", "sglasso", "cMCP", "gelasso", "ridge",
               "group-regularized", "regular")
legend("bottomright", merge=TRUE, seg.len=1, cex=1.3, fill=c(colors, 0, 0), 
       lty=c(rep(NA, length(colors)), 1, 2), 
       lwd=c(rep(NA, length(colors)), 1.5, 1.5),
       border=c(rep(1, length(colors)), 0, 0), legend=leglabels)

# briers
ylim <- range(sapply(names(results1$briers)[c(2, 5, 8, 11, 14, 17)], pred.loess, 
                     "briers", results1)[seq(2, 12, 2)], 
              mean(results1$briers$ridge))
xlim <- c(0, 250)
plot(pred.loess("gren3", "briers", results1), ylim=ylim, xlim=xlim, 
     main="d)", xlab="Number of selected features", ylab="Brier skill score", 
     cex.axis=1.5, cex.lab=2, cex.main=2, lwd=1.5, col=colors[1], type="l")
lines(pred.loess("enet3", "briers", results1), lwd=1.5, col=colors[1], lty=2)
lines(pred.loess("sglasso3", "briers", results1), lwd=1.5, col=colors[2], lty=1)
lines(pred.loess("cmcp3", "briers", results1), lwd=1.5, col=colors[3], lty=1)
lines(pred.loess("gelasso3", "briers", results1), lwd=1.5, col=colors[4], lty=1)
lines(pred.loess("grridge", "briers", results1), lwd=1.5, col=colors[5], lty=1)
lines(xlim, rep(mean(results1$briers$ridge), 2), lwd=1.5, 
      col=colors[5], lty=2)
dev.off()