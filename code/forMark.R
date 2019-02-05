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
library(psych)

### load functions
source(paste(path.code, "grVBEM.R", sep=""))

### load data
# load(paste(path.data, "forMagnus.Rdata", sep=""))
load(paste(path.data, "forMagnusN88.Rdata", sep=""))

### data manipulation
n <- length(resp)
p <- nrow(mirnormcen_resp)
unpenal <- model.matrix(~ adjth + thscheme + age + pcrcdiff, data=datfr)[, -1]
u <- ncol(unpenal)

### create partitioning
part.grridge <- list(TS=CreatePartition(as.character(which(!is.na(whichin))),
                                        as.character(1:p)))
part.greben <- list(TS=rep(c(1, 2), times=unlist(lapply(partkeep$TS, length)))[
  order(unlist(partkeep$TS))])

### cross validation
set.seed(2018)
nfolds <- n
rest <- n %% nfolds
foldid <- sample(rep(1:nfolds, times=round(c(rep(
  n %/% nfolds + as.numeric(rest!=0), times=rest),
  rep(n %/% nfolds, times=nfolds - rest)))))

psel <- c(seq(1, 5, 1), seq(7, 15, 2), seq(20, 40, 5), seq(50, 90, 10))
methods <- c("ridge", "grridge", "enet1", "enet2", "enet3", "greben1",
             "greben2", "greben3")
pred.ridge <- numeric(n)
for(m in 2:length(methods)) {
  assign(paste("pred.", methods[m], sep=""), matrix(NA, nrow=n,
                                                    ncol=length(psel)))
  assign(paste("psel.", methods[m], sep=""), matrix(NA, nrow=n,
                                                    ncol=length(psel)))
}
for(k in sort(unique(foldid))) {
  cat(paste("Fold ", k, "\n"))

  xtrain <- t(mirnormcen_resp)[foldid!=k, ]
  xtest <- matrix(t(mirnormcen_resp)[foldid==k, ], ncol=p, byrow=TRUE)
  ytrain <- (as.numeric(resp) - 1)[foldid!=k]
  utrain1 <- datfr[foldid!=k, ]
  utest1 <- datfr[foldid==k, ]
  utrain2 <- unpenal[foldid!=k, ]
  utest2 <- matrix(unpenal[foldid==k, ], ncol=u, byrow=TRUE)

  cv.ridge <- cv.glmnet(cbind(utrain2, xtrain), ytrain, alpha=0,
                        standardize=FALSE, family="binomial", intercept=TRUE,
                        penalty.factor=c(rep(0, u), rep(1, p)))
  cv.grridge <- vector("list", length(psel))
  cv.grridge[[1]] <- grridge(mirnormcen_resp, resp, part.grridge3, optl=NULL,
                             unpenal = ~1 + adjth + thscheme + age + pcrcdiff,
                             niter=1, method="exact", dataunpen=datfr,
                             innfold=10, savepredobj="all", comparelasso=FALSE,
                             optllasso=NULL, compareunpenal=TRUE,
                             selection=TRUE, maxsel=psel[1])
  for(s in 2:length(psel)) {
    cv.grridge[[s]] <- grridge(mirnormcen_resp, resp, part.grridge3,
                               optl=cv.grridge[[1]]$optl,
                               unpenal = ~1 + adjth + thscheme + age + pcrcdiff,
                               niter=1, method="exact", dataunpen=datfr,
                               innfold=10, savepredobj="all",
                               comparelasso=FALSE, optllasso=NULL,
                               compareunpenal=TRUE, selection=TRUE,
                               maxsel=psel[s])
  }
  cv1.greben <- grEBEN3(xtrain, ytrain, rep(1, length(ytrain)),
                        unpenalized=utrain2, partitions=part.greben3, alpha=0.05,
                        psel=psel, nfolds=10, trace=TRUE)
  cv2.greben <- grEBEN3(xtrain, ytrain, rep(1, length(ytrain)),
                        unpenalized=utrain2, partitions=part.greben3, alpha=0.5,
                        psel=psel, nfolds=10, trace=TRUE) # 0.251227 1.175005
  cv3.greben <- grEBEN3(xtrain, ytrain, rep(1, length(ytrain)),
                        unpenalized=utrain2, partitions=part.greben3, alpha=0.95,
                        psel=psel, nfolds=10, trace=TRUE)

  pred.enet1[foldid==k, ] <- 1/(1 + exp(-cbind(1, utest2, xtest) %*%
                                          cv1.greben$beta.nogroups))
  pred.greben1[foldid==k, ] <- 1/(1 + exp(-cbind(1, utest2, xtest) %*%
                                            cv1.greben$beta))
  pred.enet2[foldid==k, ] <- 1/(1 + exp(-cbind(1, utest2, xtest) %*%
                                          cv2.greben$beta.nogroups))
  pred.greben2[foldid==k, ] <- 1/(1 + exp(-cbind(1, utest2, xtest) %*%
                                            cv2.greben$beta))
  pred.enet3[foldid==k, ] <- 1/(1 + exp(-cbind(1, utest2, xtest) %*%
                                          cv3.greben$beta.nogroups))
  pred.greben3[foldid==k, ] <- 1/(1 + exp(-cbind(1, utest2, xtest) %*%
                                            cv3.greben$beta))
  pred.grridge[foldid==k, ] <- sapply(cv.grridge, function(s) {
    predict.grridge(s, t(xtest), dataunpennew=utest1)[, 3]})
  pred.ridge[foldid==k] <- predict(cv.ridge, cbind(utest2, xtest),
                                   type="response", s="lambda.min")

  psel.grridge[foldid==k, ] <- sapply(cv.grridge, function(s) {
    length(s$resEN$whichEN)})
  psel.enet1[foldid==k, ] <- colSums(cv1.greben$beta.nogroups!=0) - u - 1
  psel.greben1[foldid==k, ] <- colSums(cv1.greben$beta!=0) - u - 1
  psel.enet2[foldid==k, ] <- colSums(cv2.greben$beta.nogroups!=0) - u - 1
  psel.greben2[foldid==k, ] <- colSums(cv2.greben$beta!=0) - u - 1
  psel.enet3[foldid==k, ] <- colSums(cv3.greben$beta.nogroups!=0) - u - 1
  psel.greben3[foldid==k, ] <- colSums(cv3.greben$beta!=0) - u - 1

}

results5 <- list(pred=list(pred.ridge, pred.grridge, pred.enet1, pred.enet2,
                           pred.enet3, pred.greben1, pred.greben2,
                           pred.greben3),
                 psel=list(psel.grridge, psel.enet1, psel.enet2, psel.enet3,
                           psel.greben1, psel.greben2, psel.greben3))
save(results5, file=paste(path.res, "grEBEN_mirseq_Maarten_res5.Rdata", sep=""))

### cross-validation results
load(paste(path.res, "grEBEN_mirseq_Maarten_res4.Rdata", sep=""))

auc <- lapply(results4$pred[-1], function(l) {
  apply(l, 2, function(preds) {pROC::roc(as.numeric(resp) - 1, preds)$auc})})
auc <- c(pROC::roc(as.numeric(resp) - 1, results4$pred[[1]])$auc, auc)

briers <- lapply(results4$pred[-1], function(l) {
  apply(l, 2, function(preds) {
    1 - sum((as.numeric(resp) - 1 - preds)^2)/
      sum((as.numeric(resp) - 1 - mean(as.numeric(resp) - 1))^2)})})
briers <- c(1 - sum((as.numeric(resp) - 1 - results4$pred[[1]])^2)/
              sum((as.numeric(resp) - 1 - mean(as.numeric(resp) - 1))^2),
            briers)

psel <- lapply(results4$psel, function(l) {colMeans(l)})

leglabels <- c("ridge", expression(paste("enet, ", alpha==0.05)),
               expression(paste("enet, ", alpha==0.5)),
               expression(paste("enet, ", alpha==0.95)),
               "group-regularized", "not group-regularized")

png(paste(path.graph, "grEBEN_mirseq_Maarten_res4_performance.png", sep=""),
    units="in", width=12, height=6, res=120)
par(mfrow=c(1, 2))
plot(psel[[1]], auc[[2]], type="l", xlim=range(psel), ylim=range(auc), col=2,
     xlab="Number of selected variables", ylab="AUC", main="a)")
lines(range(psel), rep(auc[[1]], 2), col=2, lty=2)
lines(psel[[2]], auc[[3]], col=3, lty=2)
lines(psel[[3]], auc[[4]], col=4, lty=2)
lines(psel[[4]], auc[[5]], col=5, lty=2)
lines(psel[[5]], auc[[6]], col=3)
lines(psel[[6]], auc[[7]], col=4)
lines(psel[[7]], auc[[8]], col=5)

plot(psel[[1]], briers[[2]], type="l", xlim=range(psel), ylim=range(briers),
     col=2, xlab="Number of selected variables", ylab="Brier skill score",
     main="b)")
lines(range(psel), rep(briers[[1]], 2), col=2, lty=2)
lines(psel[[2]], briers[[3]], col=3, lty=2)
lines(psel[[3]], briers[[4]], col=4, lty=2)
lines(psel[[4]], briers[[5]], col=5, lty=2)
lines(psel[[5]], briers[[6]], col=3)
lines(psel[[6]], briers[[7]], col=4)
lines(psel[[7]], briers[[8]], col=5)
legend("bottomleft", legend=leglabels, fill=c(2:5, 0, 0),
       lty=c(rep(NA, 4), 1, 2), border=c(rep(1, 4), 0, 0), merge=TRUE,
       seg.len=1)
dev.off()