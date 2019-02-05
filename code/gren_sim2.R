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

############################# scenario 2
# similar as model 1 in yuan and lin (2006)
n <- 100
ntest <- 1000
p <- 900
G <- 300
rho <- 0.5
Gactive <- 50

Sigma <- diag(G)
for(i in 1:G) {
  for(j in 1:G) {
    Sigma[i, j] <- rho^abs(i - j)
  }
}
beta.active.mean <- 0.3

nclass <- 3
q <- c(-Inf, qnorm((c(1:nclass)/nclass)[-nclass]), Inf)
beta <- numeric(p)
beta.active <- seq(beta.active.mean*Gactive/(Gactive/2 + 0.5),
                   beta.active.mean/(Gactive/2 + 0.5), length.out=Gactive)
for(g in 1:Gactive) {
  id <- ((g - 1)*nclass + 1):((g - 1)*nclass + nclass - 1)
  beta[id] <- beta.active[g]
}

part1 <- rep(c(1:G), each=p/G)
csel <- c(seq(1, 10, 2), seq(15, 50, 5), seq(60, 140, 10), seq(160, 200, 20))
nreps <- 100
psel2 <- kappa2 <- list(grridge=matrix(NA, nrow=nreps, ncol=length(csel)),
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
auc2 <- mse2 <- briers2 <- list(ridge=numeric(nreps),
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
lambdagest2 <- list(grridge=matrix(NA, nrow=nreps, ncol=G),
                    gren1=matrix(NA, nrow=nreps, ncol=G),
                    gren2=matrix(NA, nrow=nreps, ncol=G),
                    gren3=matrix(NA, nrow=nreps, ncol=G))
for(r in 1:nreps) {
  set.seed(2018 + r - 1)
  print(paste("rep", r))
  z <- rmvnorm(n, rep(0, G), Sigma)
  x <- matrix(0, nrow=n, ncol=p)
  for(g in 1:G) {
    id <- c(((g-1)*p/G + 1):(g*p/G))
    for(class in 1:nclass) {
      x[, id[class]] <- (z[, g] >= q[class] & z[, g] < q[class + 1])
    }
  }
  
  beta0 <- as.numeric(rep(1/3, p) %*% beta)
  y <- rbinom(n, 1, as.numeric(1/(1 + exp(-x %*% beta + beta0))))
  
  ztest <- rmvnorm(n, rep(0, G), Sigma)
  xtest <- matrix(0, nrow=ntest, ncol=p)
  for(g in 1:G) {
    id <- c(((g-1)*p/G + 1):(g*p/G))
    for(class in 1:nclass) {
      xtest[, id[class]] <- (ztest[, g] >= q[class] & ztest[, g] < q[class + 1])
    }
  }
  ytest <- rbinom(ntest, 1, as.numeric(1/(1 + exp(-xtest %*% beta + beta0))))
  
  fit2.ridge <- cv.glmnet(x, y, alpha=0, standardize=FALSE)
  
  fit2.grridge <- vector("list", length(csel))
  invisible(capture.output(
    fit2.grridge[[1]] <- grridge(t(x), y, partitions=list(
      groups=CreatePartition(as.factor(part1))), selection=TRUE, maxsel=csel[1], 
      trace=FALSE, standardizeX=FALSE)))
  for(s in 2:length(csel)) {
    invisible(capture.output(
      fit2.grridge[[s]] <- grridge(t(x), y, partitions=list(
        groups=CreatePartition(as.factor(part1))), selection=TRUE, 
        maxsel=csel[s], optl=fit2.grridge[[1]]$optl, trace=FALSE, 
        standardizeX=FALSE)))
  }
  
  fit2.gren1 <- gren(x, y, partitions=list(groups=part1), alpha=0.05, 
                     trace=FALSE, control=list(epsilon=0.001, maxit=100,
                                               maxit.opt=1000, maxit.vb=100))
  fit2.gren2 <- gren(x, y, partitions=list(groups=part1), alpha=0.5, 
                     trace=FALSE, control=list(epsilon=0.001, maxit=100,
                                               maxit.opt=1000, maxit.vb=100))
  fit2.gren3 <- gren(x, y, partitions=list(groups=part1), alpha=0.95,
                     trace=FALSE, control=list(epsilon=0.001, maxit=100,
                                               maxit.opt=1000, maxit.vb=100))
  
  fit2.sglasso1 <- SGL(list(x=x, y=y), part1, type="logit", alpha=0.05, 
                       standardize=FALSE, nlam=100)
  fit2.sglasso2 <- SGL(list(x=x, y=y), part1, type="logit", alpha=0.5, 
                       standardize=FALSE, nlam=100)
  fit2.sglasso3 <- SGL(list(x=x, y=y), part1, type="logit", alpha=0.95, 
                       standardize=FALSE, nlam=100)
  
  fit2.cmcp1 <- grpreg(x, y, part1, penalty="cMCP", alpha=0.05)
  fit2.cmcp2 <- grpreg(x, y, part1, penalty="cMCP", alpha=0.5)
  fit2.cmcp3 <- grpreg(x, y, part1, penalty="cMCP", alpha=0.95)
  
  fit2.gelasso1 <- grpreg(x, y, part1, penalty="gel", alpha=0.05)
  fit2.gelasso2 <- grpreg(x, y, part1, penalty="gel", alpha=0.5)
  fit2.gelasso3 <- grpreg(x, y, part1, penalty="gel", alpha=0.95)
  
  pred2.ridge <- as.numeric(predict(fit2.ridge, xtest, "lambda.min"))
  
  pred2.grridge <- sapply(fit2.grridge, function(s) {
    predict.grridge(s, t(xtest))[, 3]})
  
  pred2.gren1 <- predict(fit2.gren1, xtest, type="groupreg",
                         s=fit2.gren1$freq.model$groupreg$lambda)
  pred2.gren2 <- predict(fit2.gren2, xtest, type="groupreg",
                         s=fit2.gren2$freq.model$groupreg$lambda)
  pred2.gren3 <- predict(fit2.gren3, xtest, type="groupreg",
                         s=fit2.gren3$freq.model$groupreg$lambda)
  
  pred2.enet1 <- predict(fit2.gren1, xtest, type="regular",
                         s=fit2.gren1$freq.model$regular$lambda)
  pred2.enet2 <- predict(fit2.gren2, xtest, type="regular",
                         s=fit2.gren2$freq.model$regular$lambda)
  pred2.enet3 <- predict(fit2.gren3, xtest, type="regular",
                         s=fit2.gren3$freq.model$regular$lambda)
  
  pred2.sglasso1 <- predictSGL(fit2.sglasso1, xtest)
  pred2.sglasso2 <- predictSGL(fit2.sglasso2, xtest)
  pred2.sglasso3 <- predictSGL(fit2.sglasso3, xtest)
  
  pred2.cmcp1 <- predict(fit2.cmcp1, xtest)
  pred2.cmcp2 <- predict(fit2.cmcp2, xtest)
  pred2.cmcp3 <- predict(fit2.cmcp3, xtest)
  
  pred2.gelasso1 <- predict(fit2.gelasso1, xtest)
  pred2.gelasso2 <- predict(fit2.gelasso2, xtest)
  pred2.gelasso3 <- predict(fit2.gelasso3, xtest)
  
  auc2$ridge[r] <- pROC::auc(ytest, pred2.ridge)
  
  auc2$grridge[r, ] <- apply(pred2.grridge, 2, function(pred) {
    pROC::auc(ytest, pred)})
  
  auc2$gren1[r, ] <- apply(pred2.gren1, 2, function(pred) {
    pROC::auc(ytest, pred)})
  auc2$gren2[r, ] <- apply(pred2.gren2, 2, function(pred) {
    pROC::auc(ytest, pred)})
  auc2$gren3[r, ] <- apply(pred2.gren3, 2, function(pred) {
    pROC::auc(ytest, pred)})
  
  auc2$enet1[r, ] <- apply(pred2.enet1, 2, function(pred) {
    pROC::auc(ytest, pred)})
  auc2$enet2[r, ] <- apply(pred2.enet2, 2, function(pred) {
    pROC::auc(ytest, pred)})
  auc2$enet3[r, ] <- apply(pred2.enet3, 2, function(pred) {
    pROC::auc(ytest, pred)})
  
  auc2$sglasso1[r, ] <- apply(pred2.sglasso1, 2, function(pred) {
    pROC::auc(ytest, pred)})
  auc2$sglasso2[r, ] <- apply(pred2.sglasso2, 2, function(pred) {
    pROC::auc(ytest, pred)})
  auc2$sglasso3[r, ] <- apply(pred2.sglasso3, 2, function(pred) {
    pROC::auc(ytest, pred)})
  
  auc2$cmcp1[r, ] <- apply(pred2.cmcp1, 2, function(pred) {
    pROC::auc(ytest, pred)})
  auc2$cmcp2[r, ] <- apply(pred2.cmcp2, 2, function(pred) {
    pROC::auc(ytest, pred)})
  auc2$cmcp3[r, ] <- apply(pred2.cmcp3, 2, function(pred) {
    pROC::auc(ytest, pred)})
  
  auc2$gelasso1[r, ] <- apply(pred2.gelasso1, 2, function(pred) {
    pROC::auc(ytest, pred)})
  auc2$gelasso2[r, ] <- apply(pred2.gelasso2, 2, function(pred) {
    pROC::auc(ytest, pred)})
  auc2$gelasso3[r, ] <- apply(pred2.gelasso3, 2, function(pred) {
    pROC::auc(ytest, pred)})
  
  const <- sum((ytest - mean(ytest))^2)
  briers2$ridge[r] <- 1 - sum((ytest - pred2.ridge)^2)/const
  
  briers2$grridge[r, ] <- apply(pred2.grridge, 2, function(pred) {
    1 - sum((ytest - pred)^2)/const})
  
  briers2$gren1[r, ] <- apply(pred2.gren1, 2, function(pred) {
    1 - sum((ytest - pred)^2)/const})
  briers2$gren2[r, ] <- apply(pred2.gren2, 2, function(pred) {
    1 - sum((ytest - pred)^2)/const})
  briers2$gren3[r, ] <- apply(pred2.gren3, 2, function(pred) {
    1 - sum((ytest - pred)^2)/const})
  
  briers2$enet1[r, ] <- apply(pred2.enet1, 2, function(pred) {
    1 - sum((ytest - pred)^2)/const})
  briers2$enet2[r, ] <- apply(pred2.enet2, 2, function(pred) {
    1 - sum((ytest - pred)^2)/const})
  briers2$enet3[r, ] <- apply(pred2.enet3, 2, function(pred) {
    1 - sum((ytest - pred)^2)/const})
  
  briers2$sglasso1[r, ] <- apply(pred2.sglasso1, 2, function(pred) {
    1 - sum((ytest - pred)^2)/const})
  briers2$sglasso2[r, ] <- apply(pred2.sglasso2, 2, function(pred) {
    1 - sum((ytest - pred)^2)/const})
  briers2$sglasso3[r, ] <- apply(pred2.sglasso3, 2, function(pred) {
    1 - sum((ytest - pred)^2)/const})
  
  briers2$cmcp1[r, ] <- apply(pred2.cmcp1, 2, function(pred) {
    1 - sum((ytest - pred)^2)/const})
  briers2$cmcp2[r, ] <- apply(pred2.cmcp2, 2, function(pred) {
    1 - sum((ytest - pred)^2)/const})
  briers2$cmcp3[r, ] <- apply(pred2.cmcp3, 2, function(pred) {
    1 - sum((ytest - pred)^2)/const})
  
  briers2$gelasso1[r, ] <- apply(pred2.gelasso1, 2, function(pred) {
    1 - sum((ytest - pred)^2)/const})
  briers2$gelasso2[r, ] <- apply(pred2.gelasso2, 2, function(pred) {
    1 - sum((ytest - pred)^2)/const})
  briers2$gelasso3[r, ] <- apply(pred2.gelasso3, 2, function(pred) {
    1 - sum((ytest - pred)^2)/const})
  
  mse2$ridge[r] <- mean((coef(fit2.ridge)[-1] - beta)^2)
  
  mse2$grridge[r, ] <- sapply(fit2.grridge, function(s) {
    mean((replace(rep(0, p), s$resEN$whichEN, s$resEN$betasEN) - beta)^2)})
  
  mse2$gren1[r, ] <- apply(fit2.gren1$freq.model$groupreg$beta, 2, function(b) {
    mean((b - beta)^2)})
  mse2$gren2[r, ] <- apply(fit2.gren2$freq.model$groupreg$beta, 2, function(b) {
    mean((b - beta)^2)})
  mse2$gren3[r, ] <- apply(fit2.gren3$freq.model$groupreg$beta, 2, function(b) {
    mean((b - beta)^2)})
  
  mse2$enet1[r, ] <- apply(fit2.gren1$freq.model$regular$beta, 2, function(b) {
    mean((b - beta)^2)})
  mse2$enet2[r, ] <- apply(fit2.gren2$freq.model$regular$beta, 2, function(b) {
    mean((b - beta)^2)})
  mse2$enet3[r, ] <- apply(fit2.gren3$freq.model$regular$beta, 2, function(b) {
    mean((b - beta)^2)})
  
  mse2$sglasso1[r, ] <- apply(fit2.sglasso1$beta, 2, function(b) {
    mean((b - beta)^2)})
  mse2$sglasso2[r, ] <- apply(fit2.sglasso2$beta, 2, function(b) {
    mean((b - beta)^2)})
  mse2$sglasso3[r, ] <- apply(fit2.sglasso3$beta, 2, function(b) {
    mean((b - beta)^2)})
  
  mse2$cmcp1[r, ] <- apply(fit2.cmcp1$beta[-1, ], 2, function(b) {
    mean((b - beta)^2)})
  mse2$cmcp2[r, ] <- apply(fit2.cmcp2$beta[-1, ], 2, function(b) {
    mean((b - beta)^2)})
  mse2$cmcp3[r, ] <- apply(fit2.cmcp3$beta[-1, ], 2, function(b) {
    mean((b - beta)^2)})
  
  mse2$gelasso1[r, ] <- apply(fit2.gelasso1$beta[-1, ], 2, function(b) {
    mean((b - beta)^2)})
  mse2$gelasso2[r, ] <- apply(fit2.gelasso2$beta[-1, ], 2, function(b) {
    mean((b - beta)^2)})
  mse2$gelasso3[r, ] <- apply(fit2.gelasso3$beta[-1, ], 2, function(b) {
    mean((b - beta)^2)})
  
  kappa2$grridge[r, ] <- sapply(fit2.grridge, function(s) {
    kappa2(cbind(beta!=0, replace(rep(FALSE, p), s$resEN$whichEN, 
                                  TRUE)))$value})
  
  kappa2$gren1[r, ] <- apply(fit2.gren1$freq.model$groupreg$beta, 2, 
                             function(b) {kappa2(cbind(beta!=0, b!=0))$value})
  kappa2$gren2[r, ] <- apply(fit2.gren2$freq.model$groupreg$beta, 2, 
                             function(b) {kappa2(cbind(beta!=0, b!=0))$value})
  kappa2$gren3[r, ] <- apply(fit2.gren3$freq.model$groupreg$beta, 2, 
                             function(b) {kappa2(cbind(beta!=0, b!=0))$value})
  
  kappa2$enet1[r, ] <- apply(fit2.gren1$freq.model$regular$beta, 2, 
                             function(b) {kappa2(cbind(beta!=0, b!=0))$value})
  kappa2$enet2[r, ] <- apply(fit2.gren2$freq.model$regular$beta, 2, 
                             function(b) {kappa2(cbind(beta!=0, b!=0))$value})
  kappa2$enet3[r, ] <- apply(fit2.gren3$freq.model$regular$beta, 2, 
                             function(b) {kappa2(cbind(beta!=0, b!=0))$value})
  
  kappa2$sglasso1[r, ] <- apply(fit2.sglasso1$beta, 2, function(b) {
    kappa2(cbind(beta!=0, b!=0))$value})
  kappa2$sglasso2[r, ] <- apply(fit2.sglasso2$beta, 2, function(b) {
    kappa2(cbind(beta!=0, b!=0))$value})
  kappa2$sglasso3[r, ] <- apply(fit2.sglasso3$beta, 2, function(b) {
    kappa2(cbind(beta!=0, b!=0))$value})
  
  kappa2$cmcp1[r, ] <- apply(fit2.cmcp1$beta[-1, ], 2, function(b) {
    kappa2(cbind(beta!=0, b!=0))$value})
  kappa2$cmcp2[r, ] <- apply(fit2.cmcp2$beta[-1, ], 2, function(b) {
    kappa2(cbind(beta!=0, b!=0))$value})
  kappa2$cmcp3[r, ] <- apply(fit2.cmcp3$beta[-1, ], 2, function(b) {
    kappa2(cbind(beta!=0, b!=0))$value})
  
  kappa2$gelasso1[r, ] <- apply(fit2.gelasso1$beta[-1, ], 2, function(b) {
    kappa2(cbind(beta!=0, b!=0))$value})
  kappa2$gelasso2[r, ] <- apply(fit2.gelasso2$beta[-1, ], 2, function(b) {
    kappa2(cbind(beta!=0, b!=0))$value})
  kappa2$gelasso3[r, ] <- apply(fit2.gelasso3$beta[-1, ], 2, function(b) {
    kappa2(cbind(beta!=0, b!=0))$value})
  
  psel2$grridge[r, ] <- sapply(fit2.grridge, function(s) {
    length(s$resEN$whichEN)})
  
  psel2$gren1[r, ] <- fit2.gren1$freq.model$groupreg$df
  psel2$gren2[r, ] <- fit2.gren2$freq.model$groupreg$df
  psel2$gren3[r, ] <- fit2.gren3$freq.model$groupreg$df
  
  psel2$enet1[r, ] <- fit2.gren1$freq.model$regular$df
  psel2$enet2[r, ] <- fit2.gren2$freq.model$regular$df
  psel2$enet3[r, ] <- fit2.gren3$freq.model$regular$df
  
  psel2$sglasso1[r, ] <- apply(fit2.sglasso1$beta, 2, function(b) {sum(b!=0)})
  psel2$sglasso2[r, ] <- apply(fit2.sglasso2$beta, 2, function(b) {sum(b!=0)})
  psel2$sglasso3[r, ] <- apply(fit2.sglasso3$beta, 2, function(b) {sum(b!=0)})
  
  psel2$cmcp1[r, ] <- apply(fit2.cmcp1$beta, 2, function(b) {sum(b!=0)})
  psel2$cmcp2[r, ] <- apply(fit2.cmcp2$beta, 2, function(b) {sum(b!=0)})
  psel2$cmcp3[r, ] <- apply(fit2.cmcp3$beta, 2, function(b) {sum(b!=0)})
  
  psel2$gelasso1[r, ] <- apply(fit2.gelasso1$beta, 2, function(b) {sum(b!=0)})
  psel2$gelasso2[r, ] <- apply(fit2.gelasso2$beta, 2, function(b) {sum(b!=0)})
  psel2$gelasso3[r, ] <- apply(fit2.gelasso3$beta, 2, function(b) {sum(b!=0)})
  
  lambdagest2$grridge[r, ] <- fit2.grridge[[1]]$lambdamults$groups
  lambdagest2$gren1[r, ] <- fit2.gren1$lambdag$groups
  lambdagest2$gren2[r, ] <- fit2.gren2$lambdag$groups
  lambdagest2$gren3[r, ] <- fit2.gren3$lambdag$groups
  
  results2 <- list(auc=auc2, briers=briers2, mse=mse2, kappa=kappa2, psel=psel2,
                   lambdag=lambdagest2)
  save(results2, file=paste(path.res, "gren_sim2_res1.Rdata", sep=""))
  
}



path.res <- ifelse(as.character(Sys.info()[1])!="Darwin", "~/EBEN/results/",
                   "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/results/")

path.graph <- "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/graphs/"
load(paste(path.res, "gren_sim2_res1.Rdata", sep=""))
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

colors <- bpy.colors(8)[-c(1, 8)]
png(paste(path.graph, "gren_sim2_res1_performance1.png", sep=""),
    units="in", width=14, height=10, res=120)
par(mfrow=c(2, 2), mar=c(5.1, 5.1, 4.1, 2.1))

# kappa
ylim <- range(sapply(names(results2$kappa)[-1], pred.loess, "kappa", results2)[
  seq(2, 2*(length(results2$kappa) - 1), 2)])
plot(pred.loess("gren1", "kappa", results2), ylim=ylim, xlim=range(results2$psel), 
     main="a)", xlab="Number of selected features", ylab="Cohen's kappa", 
     cex.axis=1.5, cex.lab=2, cex.main=2, lwd=1.5, col=colors[1], type="l")
lines(pred.loess("gren2", "kappa", results2), lwd=1.5, col=colors[1], lty=2)
lines(pred.loess("gren3", "kappa", results2), lwd=1.5, col=colors[1], lty=3)

lines(pred.loess("enet1", "kappa", results2), lwd=1.5, col=colors[2], lty=1)
lines(pred.loess("enet2", "kappa", results2), lwd=1.5, col=colors[2], lty=2)
lines(pred.loess("enet3", "kappa", results2), lwd=1.5, col=colors[2], lty=3)

lines(pred.loess("sglasso1", "kappa", results2), lwd=1.5, col=colors[3], lty=1)
lines(pred.loess("sglasso2", "kappa", results2), lwd=1.5, col=colors[3], lty=2)
lines(pred.loess("sglasso3", "kappa", results2), lwd=1.5, col=colors[3], lty=3)

lines(pred.loess("cmcp1", "kappa", results2), lwd=1.5, col=colors[4], lty=1)
lines(pred.loess("cmcp2", "kappa", results2), lwd=1.5, col=colors[4], lty=2)
lines(pred.loess("cmcp3", "kappa", results2), lwd=1.5, col=colors[4], lty=3)

lines(pred.loess("gelasso1", "kappa", results2), lwd=1.5, col=colors[5], lty=1)
lines(pred.loess("gelasso2", "kappa", results2), lwd=1.5, col=colors[5], lty=2)
lines(pred.loess("gelasso3", "kappa", results2), lwd=1.5, col=colors[5], lty=3)

lines(pred.loess("grridge", "kappa", results2), lwd=1.5, col=colors[6], lty=1)

# mse
ylim <- range(sapply(names(results2$mse)[-1], pred.loess, "mse", results2)[
  seq(2, 2*(length(results2$mse) - 1), 2)], mean(results2$mse$ridge))
plot(pred.loess("gren1", "mse", results2), ylim=ylim, xlim=range(results2$psel), 
     main="b)", xlab="Number of selected features", ylab="MSE", cex.axis=1.5, 
     cex.lab=2, cex.main=2, lwd=1.5, col=colors[1], type="l")
lines(pred.loess("gren2", "mse", results2), lwd=1.5, col=colors[1], lty=2)
lines(pred.loess("gren3", "mse", results2), lwd=1.5, col=colors[1], lty=3)

lines(pred.loess("enet1", "mse", results2), lwd=1.5, col=colors[2], lty=1)
lines(pred.loess("enet2", "mse", results2), lwd=1.5, col=colors[2], lty=2)
lines(pred.loess("enet3", "mse", results2), lwd=1.5, col=colors[2], lty=3)

lines(pred.loess("sglasso1", "mse", results2), lwd=1.5, col=colors[3], lty=1)
lines(pred.loess("sglasso2", "mse", results2), lwd=1.5, col=colors[3], lty=2)
lines(pred.loess("sglasso3", "mse", results2), lwd=1.5, col=colors[3], lty=3)

lines(pred.loess("cmcp1", "mse", results2), lwd=1.5, col=colors[4], lty=1)
lines(pred.loess("cmcp2", "mse", results2), lwd=1.5, col=colors[4], lty=2)
lines(pred.loess("cmcp3", "mse", results2), lwd=1.5, col=colors[4], lty=3)

lines(pred.loess("gelasso1", "mse", results2), lwd=1.5, col=colors[5], lty=1)
lines(pred.loess("gelasso2", "mse", results2), lwd=1.5, col=colors[5], lty=2)
lines(pred.loess("gelasso3", "mse", results2), lwd=1.5, col=colors[5], lty=3)

lines(pred.loess("grridge", "mse", results2), lwd=1.5, col=colors[6], lty=1)
lines(range(results2$psel), rep(mean(results2$mse$ridge), 2), lwd=1.5, 
      col=colors[6], lty=2)

# auc
ylim <- range(sapply(names(results2$auc)[-1], pred.loess, "auc", results2)[
  seq(2, 2*(length(results2$auc) - 1), 2)], mean(results2$auc$ridge))
plot(pred.loess("gren1", "auc", results2), ylim=ylim, xlim=range(results2$psel), 
     main="c)", xlab="Number of selected features", ylab="AUC", cex.axis=1.5, 
     cex.lab=2, cex.main=2, lwd=1.5, col=colors[1], type="l")
lines(pred.loess("gren2", "auc", results2), lwd=1.5, col=colors[1], lty=2)
lines(pred.loess("gren3", "auc", results2), lwd=1.5, col=colors[1], lty=3)

lines(pred.loess("enet1", "auc", results2), lwd=1.5, col=colors[2], lty=1)
lines(pred.loess("enet2", "auc", results2), lwd=1.5, col=colors[2], lty=2)
lines(pred.loess("enet3", "auc", results2), lwd=1.5, col=colors[2], lty=3)

lines(pred.loess("sglasso1", "auc", results2), lwd=1.5, col=colors[3], lty=1)
lines(pred.loess("sglasso2", "auc", results2), lwd=1.5, col=colors[3], lty=2)
lines(pred.loess("sglasso3", "auc", results2), lwd=1.5, col=colors[3], lty=3)

lines(pred.loess("cmcp1", "auc", results2), lwd=1.5, col=colors[4], lty=1)
lines(pred.loess("cmcp2", "auc", results2), lwd=1.5, col=colors[4], lty=2)
lines(pred.loess("cmcp3", "auc", results2), lwd=1.5, col=colors[4], lty=3)

lines(pred.loess("gelasso1", "auc", results2), lwd=1.5, col=colors[5], lty=1)
lines(pred.loess("gelasso2", "auc", results2), lwd=1.5, col=colors[5], lty=2)
lines(pred.loess("gelasso3", "auc", results2), lwd=1.5, col=colors[5], lty=3)

lines(pred.loess("grridge", "auc", results2), lwd=1.5, col=colors[6], lty=1)
lines(range(results2$psel), rep(mean(results2$auc$ridge), 2), lwd=1.5, 
      col=colors[6], lty=2)

# briers
ylim <- range(sapply(names(results2$briers)[-1], pred.loess, "briers", results2)[
  seq(2, 2*(length(results2$briers) - 1), 2)], mean(results2$briers$ridge))
plot(pred.loess("gren1", "briers", results2), ylim=ylim, xlim=range(results2$psel), 
     main="d)", xlab="Number of selected features", ylab="Brier skill score", 
     cex.axis=1.5, cex.lab=2, cex.main=2, lwd=1.5, col=colors[1], type="l")
lines(pred.loess("gren2", "briers", results2), lwd=1.5, col=colors[1], lty=2)
lines(pred.loess("gren3", "briers", results2), lwd=1.5, col=colors[1], lty=3)

lines(pred.loess("enet1", "briers", results2), lwd=1.5, col=colors[2], lty=1)
lines(pred.loess("enet2", "briers", results2), lwd=1.5, col=colors[2], lty=2)
lines(pred.loess("enet3", "briers", results2), lwd=1.5, col=colors[2], lty=3)

lines(pred.loess("sglasso1", "briers", results2), lwd=1.5, col=colors[3], lty=1)
lines(pred.loess("sglasso2", "briers", results2), lwd=1.5, col=colors[3], lty=2)
lines(pred.loess("sglasso3", "briers", results2), lwd=1.5, col=colors[3], lty=3)

lines(pred.loess("cmcp1", "briers", results2), lwd=1.5, col=colors[4], lty=1)
lines(pred.loess("cmcp2", "briers", results2), lwd=1.5, col=colors[4], lty=2)
lines(pred.loess("cmcp3", "briers", results2), lwd=1.5, col=colors[4], lty=3)

lines(pred.loess("gelasso1", "briers", results2), lwd=1.5, col=colors[5], lty=1)
lines(pred.loess("gelasso2", "briers", results2), lwd=1.5, col=colors[5], lty=2)
lines(pred.loess("gelasso3", "briers", results2), lwd=1.5, col=colors[5], lty=3)

lines(pred.loess("grridge", "briers", results2), lwd=1.5, col=colors[6], lty=1)
lines(range(results2$psel), rep(mean(results2$briers$ridge), 2), lwd=1.5, 
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
png(paste(path.graph, "gren_sim2_res1_performance1_set1.png", sep=""),
    units="in", width=14, height=10, res=120)
par(mfrow=c(2, 2), mar=c(5.1, 5.1, 4.1, 2.1))

# kappa
ylim <- range(sapply(names(results2$kappa)[c(1, 2, 5, 8, 11, 14)], pred.loess,
                     "kappa", results2)[seq(2, 12, 2)])
# xlim <- range(results2$psel)
xlim <- c(0, 500)
plot(pred.loess("gren1", "kappa", results2), ylim=ylim, xlim=xlim,
     main="a)", xlab="Number of selected features", ylab="Cohen's kappa",
     cex.axis=1.5, cex.lab=2, cex.main=2, lwd=1.5, col=colors[1], type="l")
lines(pred.loess("enet1", "kappa", results2), lwd=1.5, col=colors[1], lty=2)
lines(pred.loess("sglasso1", "kappa", results2), lwd=1.5, col=colors[2], lty=1)
lines(pred.loess("cmcp1", "kappa", results2), lwd=1.5, col=colors[3], lty=1)
lines(pred.loess("gelasso1", "kappa", results2), lwd=1.5, col=colors[4], lty=1)
lines(pred.loess("grridge", "kappa", results2), lwd=1.5, col=colors[5], lty=1)

# mse
ylim <- range(sapply(names(results2$mse)[c(2, 3, 6, 9, 12, 15)], pred.loess,
                     "mse", results2)[seq(2, 12, 2)], mean(results2$mse$ridge))
xlim <- c(0, 500)
plot(pred.loess("gren1", "mse", results2), ylim=ylim, xlim=xlim,
     main="b)", xlab="Number of selected features", ylab="MSE", cex.axis=1.5,
     cex.lab=2, cex.main=2, lwd=1.5, col=colors[1], type="l")
lines(pred.loess("enet1", "mse", results2), lwd=1.5, col=colors[1], lty=2)
lines(pred.loess("sglasso1", "mse", results2), lwd=1.5, col=colors[2], lty=1)
lines(pred.loess("cmcp1", "mse", results2), lwd=1.5, col=colors[3], lty=1)
lines(pred.loess("gelasso1", "mse", results2), lwd=1.5, col=colors[4], lty=1)
lines(pred.loess("grridge", "mse", results2), lwd=1.5, col=colors[5], lty=1)
lines(xlim, rep(mean(results2$mse$ridge), 2), lwd=1.5,
      col=colors[5], lty=2)

# auc
ylim <- range(sapply(names(results2$auc)[c(2, 3, 6, 9, 12, 15)], pred.loess,
                     "auc", results2)[seq(2, 12, 2)], mean(results2$auc$ridge))
xlim <- c(0, 500)
plot(pred.loess("gren1", "auc", results2), ylim=ylim, xlim=xlim,
     main="c)", xlab="Number of selected features", ylab="AUC", cex.axis=1.5,
     cex.lab=2, cex.main=2, lwd=1.5, col=colors[1], type="l")
lines(pred.loess("enet1", "auc", results2), lwd=1.5, col=colors[1], lty=2)
lines(pred.loess("sglasso1", "auc", results2), lwd=1.5, col=colors[2], lty=1)
lines(pred.loess("cmcp1", "auc", results2), lwd=1.5, col=colors[3], lty=1)
lines(pred.loess("gelasso1", "auc", results2), lwd=1.5, col=colors[4], lty=1)
lines(pred.loess("grridge", "auc", results2), lwd=1.5, col=colors[5], lty=1)
lines(xlim, rep(mean(results2$auc$ridge), 2), lwd=1.5,
      col=colors[5], lty=2)

# legend
leglabels <- c("enet", "sglasso", "cMCP", "gelasso", "ridge",
               "group-regularized", "regular")
legend("bottomright", merge=TRUE, seg.len=1, cex=1.3, fill=c(colors, 0, 0),
       lty=c(rep(NA, length(colors)), 1, 2),
       lwd=c(rep(NA, length(colors)), 1.5, 1.5),
       border=c(rep(1, length(colors)), 0, 0), legend=leglabels)

# briers
ylim <- range(sapply(names(results2$briers)[c(2, 3, 6, 9, 12, 15)], pred.loess,
                     "briers", results2)[seq(2, 12, 2)],
              mean(results2$briers$ridge))
xlim <- c(0, 500)
plot(pred.loess("gren1", "briers", results2), ylim=ylim, xlim=xlim,
     main="d)", xlab="Number of selected features", ylab="Brier skill score",
     cex.axis=1.5, cex.lab=2, cex.main=2, lwd=1.5, col=colors[1], type="l")
lines(pred.loess("enet1", "briers", results2), lwd=1.5, col=colors[1], lty=2)
lines(pred.loess("sglasso1", "briers", results2), lwd=1.5, col=colors[2], lty=1)
lines(pred.loess("cmcp1", "briers", results2), lwd=1.5, col=colors[3], lty=1)
lines(pred.loess("gelasso1", "briers", results2), lwd=1.5, col=colors[4], lty=1)
lines(pred.loess("grridge", "briers", results2), lwd=1.5, col=colors[5], lty=1)
lines(xlim, rep(mean(results2$briers$ridge), 2), lwd=1.5,
      col=colors[5], lty=2)
dev.off()



### performance for alpha=0.5
colors <- bpy.colors(7)[-c(1, 7)]
png(paste(path.graph, "gren_sim2_res1_performance1_set2.png", sep=""),
    units="in", width=14, height=10, res=120)
par(mfrow=c(2, 2), mar=c(5.1, 5.1, 4.1, 2.1))

# kappa
ylim <- range(sapply(names(results2$kappa)[c(1, 3, 6, 9, 12, 15)], pred.loess,
                     "kappa", results2)[seq(2, 12, 2)])
# xlim <- range(results2$psel)
xlim <- c(0, 500)
plot(pred.loess("gren2", "kappa", results2), ylim=ylim, xlim=xlim,
     main="a)", xlab="Number of selected features", ylab="Cohen's kappa",
     cex.axis=1.5, cex.lab=2, cex.main=2, lwd=1.5, col=colors[1], type="l")
lines(pred.loess("enet2", "kappa", results2), lwd=1.5, col=colors[1], lty=2)
lines(pred.loess("sglasso2", "kappa", results2), lwd=1.5, col=colors[2], lty=1)
lines(pred.loess("cmcp2", "kappa", results2), lwd=1.5, col=colors[3], lty=1)
lines(pred.loess("gelasso2", "kappa", results2), lwd=1.5, col=colors[4], lty=1)
lines(pred.loess("grridge", "kappa", results2), lwd=1.5, col=colors[5], lty=1)

# mse
ylim <- range(sapply(names(results2$mse)[c(2, 4, 7, 10, 13, 16)], pred.loess,
                     "mse", results2)[seq(2, 12, 2)], mean(results2$mse$ridge))
xlim <- c(0, 500)
plot(pred.loess("gren2", "mse", results2), ylim=ylim, xlim=xlim,
     main="b)", xlab="Number of selected features", ylab="MSE", cex.axis=1.5,
     cex.lab=2, cex.main=2, lwd=1.5, col=colors[1], type="l")
lines(pred.loess("enet2", "mse", results2), lwd=1.5, col=colors[1], lty=2)
lines(pred.loess("sglasso2", "mse", results2), lwd=1.5, col=colors[2], lty=1)
lines(pred.loess("cmcp2", "mse", results2), lwd=1.5, col=colors[3], lty=1)
lines(pred.loess("gelasso2", "mse", results2), lwd=1.5, col=colors[4], lty=1)
lines(pred.loess("grridge", "mse", results2), lwd=1.5, col=colors[5], lty=1)
lines(xlim, rep(mean(results2$mse$ridge), 2), lwd=1.5,
      col=colors[5], lty=2)

# auc
ylim <- range(sapply(names(results2$auc)[c(2, 4, 7, 10, 13, 16)], pred.loess,
                     "auc", results2)[seq(2, 12, 2)], mean(results2$auc$ridge))
xlim <- c(0, 500)
plot(pred.loess("gren2", "auc", results2), ylim=ylim, xlim=xlim,
     main="c)", xlab="Number of selected features", ylab="AUC", cex.axis=1.5,
     cex.lab=2, cex.main=2, lwd=1.5, col=colors[1], type="l")
lines(pred.loess("enet2", "auc", results2), lwd=1.5, col=colors[1], lty=2)
lines(pred.loess("sglasso2", "auc", results2), lwd=1.5, col=colors[2], lty=1)
lines(pred.loess("cmcp2", "auc", results2), lwd=1.5, col=colors[3], lty=1)
lines(pred.loess("gelasso2", "auc", results2), lwd=1.5, col=colors[4], lty=1)
lines(pred.loess("grridge", "auc", results2), lwd=1.5, col=colors[5], lty=1)
lines(xlim, rep(mean(results2$auc$ridge), 2), lwd=1.5,
      col=colors[5], lty=2)

# legend
leglabels <- c("enet", "sglasso", "cMCP", "gelasso", "ridge",
               "group-regularized", "regular")
legend("bottomright", merge=TRUE, seg.len=1, cex=1.3, fill=c(colors, 0, 0),
       lty=c(rep(NA, length(colors)), 1, 2),
       lwd=c(rep(NA, length(colors)), 1.5, 1.5),
       border=c(rep(1, length(colors)), 0, 0), legend=leglabels)

# briers
ylim <- range(sapply(names(results2$briers)[c(2, 4, 7, 10, 13, 16)], pred.loess,
                     "briers", results2)[seq(2, 12, 2)],
              mean(results2$briers$ridge))
xlim <- c(0, 500)
plot(pred.loess("gren2", "briers", results2), ylim=ylim, xlim=xlim,
     main="d)", xlab="Number of selected features", ylab="Brier skill score",
     cex.axis=1.5, cex.lab=2, cex.main=2, lwd=1.5, col=colors[1], type="l")
lines(pred.loess("enet2", "briers", results2), lwd=1.5, col=colors[1], lty=2)
lines(pred.loess("sglasso2", "briers", results2), lwd=1.5, col=colors[2], lty=1)
lines(pred.loess("cmcp2", "briers", results2), lwd=1.5, col=colors[3], lty=1)
lines(pred.loess("gelasso2", "briers", results2), lwd=1.5, col=colors[4], lty=1)
lines(pred.loess("grridge", "briers", results2), lwd=1.5, col=colors[5], lty=1)
lines(xlim, rep(mean(results2$briers$ridge), 2), lwd=1.5,
      col=colors[5], lty=2)
dev.off()



### performance for alpha=0.95
colors <- bpy.colors(7)[-c(1, 7)]
png(paste(path.graph, "gren_sim2_res1_performance1_set3.png", sep=""),
    units="in", width=14, height=10, res=120)
par(mfrow=c(2, 2), mar=c(5.1, 5.1, 4.1, 2.1))

# kappa
ylim <- range(sapply(names(results2$kappa)[c(1, 4, 7, 10, 13, 16)], pred.loess,
                     "kappa", results2)[seq(2, 12, 2)])
# xlim <- range(results2$psel)
xlim <- c(0, 500)
plot(pred.loess("gren3", "kappa", results2), ylim=ylim, xlim=xlim,
     main="a)", xlab="Number of selected features", ylab="Cohen's kappa",
     cex.axis=1.5, cex.lab=2, cex.main=2, lwd=1.5, col=colors[1], type="l")
lines(pred.loess("enet3", "kappa", results2), lwd=1.5, col=colors[1], lty=2)
lines(pred.loess("sglasso3", "kappa", results2), lwd=1.5, col=colors[2], lty=1)
lines(pred.loess("cmcp3", "kappa", results2), lwd=1.5, col=colors[3], lty=1)
lines(pred.loess("gelasso3", "kappa", results2), lwd=1.5, col=colors[4], lty=1)
lines(pred.loess("grridge", "kappa", results2), lwd=1.5, col=colors[5], lty=1)

# mse
ylim <- range(sapply(names(results2$mse)[c(2, 5, 8, 11, 14, 17)], pred.loess,
                     "mse", results2)[seq(2, 12, 2)], mean(results2$mse$ridge))
xlim <- c(0, 500)
plot(pred.loess("gren3", "mse", results2), ylim=ylim, xlim=xlim,
     main="b)", xlab="Number of selected features", ylab="MSE", cex.axis=1.5,
     cex.lab=2, cex.main=2, lwd=1.5, col=colors[1], type="l")
lines(pred.loess("enet3", "mse", results2), lwd=1.5, col=colors[1], lty=2)
lines(pred.loess("sglasso3", "mse", results2), lwd=1.5, col=colors[2], lty=1)
lines(pred.loess("cmcp3", "mse", results2), lwd=1.5, col=colors[3], lty=1)
lines(pred.loess("gelasso3", "mse", results2), lwd=1.5, col=colors[4], lty=1)
lines(pred.loess("grridge", "mse", results2), lwd=1.5, col=colors[5], lty=1)
lines(xlim, rep(mean(results2$mse$ridge), 2), lwd=1.5,
      col=colors[5], lty=2)

# auc
ylim <- range(sapply(names(results2$auc)[c(2, 5, 8, 11, 14, 17)], pred.loess,
                     "auc", results2)[seq(2, 12, 2)], mean(results2$auc$ridge))
xlim <- c(0, 500)
plot(pred.loess("gren3", "auc", results2), ylim=ylim, xlim=xlim,
     main="c)", xlab="Number of selected features", ylab="AUC", cex.axis=1.5,
     cex.lab=2, cex.main=2, lwd=1.5, col=colors[1], type="l")
lines(pred.loess("enet3", "auc", results2), lwd=1.5, col=colors[1], lty=2)
lines(pred.loess("sglasso3", "auc", results2), lwd=1.5, col=colors[2], lty=1)
lines(pred.loess("cmcp3", "auc", results2), lwd=1.5, col=colors[3], lty=1)
lines(pred.loess("gelasso3", "auc", results2), lwd=1.5, col=colors[4], lty=1)
lines(pred.loess("grridge", "auc", results2), lwd=1.5, col=colors[5], lty=1)
lines(xlim, rep(mean(results2$auc$ridge), 2), lwd=1.5,
      col=colors[5], lty=2)

# legend
leglabels <- c("enet", "sglasso", "cMCP", "gelasso", "ridge",
               "group-regularized", "regular")
legend("bottomright", merge=TRUE, seg.len=1, cex=1.3, fill=c(colors, 0, 0),
       lty=c(rep(NA, length(colors)), 1, 2),
       lwd=c(rep(NA, length(colors)), 1.5, 1.5),
       border=c(rep(1, length(colors)), 0, 0), legend=leglabels)

# briers
ylim <- range(sapply(names(results2$briers)[c(2, 5, 8, 11, 14, 17)], pred.loess,
                     "briers", results2)[seq(2, 12, 2)],
              mean(results2$briers$ridge))
xlim <- c(0, 500)
plot(pred.loess("gren3", "briers", results2), ylim=ylim, xlim=xlim,
     main="d)", xlab="Number of selected features", ylab="Brier skill score",
     cex.axis=1.5, cex.lab=2, cex.main=2, lwd=1.5, col=colors[1], type="l")
lines(pred.loess("enet3", "briers", results2), lwd=1.5, col=colors[1], lty=2)
lines(pred.loess("sglasso3", "briers", results2), lwd=1.5, col=colors[2], lty=1)
lines(pred.loess("cmcp3", "briers", results2), lwd=1.5, col=colors[3], lty=1)
lines(pred.loess("gelasso3", "briers", results2), lwd=1.5, col=colors[4], lty=1)
lines(pred.loess("grridge", "briers", results2), lwd=1.5, col=colors[5], lty=1)
lines(xlim, rep(mean(results2$briers$ridge), 2), lwd=1.5,
      col=colors[5], lty=2)
dev.off()
