#!/usr/bin/env Rscript

### installation of packages
if(substr(system('git log -n 1 --format="%h %aN %s %ad"', intern=TRUE), 1, 7)!=
   substr(packageDescription("gren")$GithubSHA1, 1, 7)) {
  if(!("devtools" %in% installed.packages())) {
    install.packages("devtools")
  }
  library(devtools)
  install_github("magnusmunch/gren/rpackage", local=FALSE,
                 auth_token=Sys.getenv("GITHUB_PAT"))
}

### libraries
library(gren)
library(GRridge)
library(grpreg)
library(Biobase)
library(SGL)
library(grpregOverlap)

### load data
load("data/forMagnusN88.Rdata")

### create model matrix for unpenalized covariates
unpenal <- model.matrix(~ adjth + thscheme + age + pcrcdiff, data=datfr)[, -1]

### create partitioning based on FDR <= 0.05
diff.twogroup <- rep(c(2, 1), times=unlist(lapply(partkeep$TS, length)))[
  order(unlist(partkeep$TS))]

### create partitioning based on FDR <= 0.05 and FDR <= 0.001
miRNA.BFDR <- as.character(TumMirs$miRNA)[TumMirs$BFDR_PNminP < 0.001]
miRNA.TumMirs <- as.character(TumMirs$miRNA)
miRNA <- as.character(sapply(rownames(mirnormcen_resp), function(s) {
  strsplit(s, split=" ")[[1]][1]}))
diff.threegroup <- miRNA %in% miRNA.BFDR + miRNA %in% miRNA.TumMirs + 1

### target vector
benefit <- as.numeric(resp) - 1

### mirna data
micrornas <- scale(t(mirnormcen_resp))
colnames(micrornas) <- miRNA

### randomly split in test and train data (30% and 70%, respectively)
set.seed(2019)
id.train <- c(sample(which(benefit==0), size=floor(sum(benefit==0)*0.7)),
              sample(which(benefit==1), size=floor(sum(benefit==1)*0.7)))

################################### model 1 ####################################
### fitting the models
set.seed(2019)
ytrain <- benefit[id.train]
xtrain <- micrornas[id.train, ]
utrain <- unpenal[id.train, ]
ytest <- benefit[-id.train]
xtest <- scale(micrornas[-id.train, ])
utest <- utrain[-id.train, ]
part <- diff.threegroup
p <- ncol(xtrain)

sum(apply(xtrain, 2, sd)==0)
test <- glmnet(xtrain, ytrain, family="binomial")
fit1.gren1 <- gren(xtrain, ytrain, partitions=list(part=part), alpha=0.05, 
                   trace=FALSE)
fit1.gren2 <- gren(xtrain, ytrain, partitions=list(part=part), alpha=0.5, 
                   trace=FALSE)
fit1.gren3 <- gren(xtrain, ytrain, partitions=list(part=part), alpha=0.95, 
                   trace=FALSE)

fit1.grridge <- grridge(t(xtrain), ytrain, list(part=split(1:p, part)))

fit1.sgl1 <- cvSGL(list(x=xtrain, y=ytrain), part, type="logit", alpha=0.05)
fit1.sgl1$fit$type <- "logit"
fit1.sgl2 <- cvSGL(list(x=xtrain, y=ytrain), part, type="logit", alpha=0.5)
fit1.sgl2$fit$type <- "logit"
fit1.sgl3 <- cvSGL(list(x=xtrain, y=ytrain), part, type="logit", alpha=0.95)
fit1.sgl3$fit$type <- "logit"

fit1.cmcp1 <- cv.grpreg(xtrain, ytrain, part, penalty="cMCP", 
                        family="binomial", alpha=0.05)
fit1.cmcp2 <- cv.grpreg(xtrain, ytrain, part, penalty="cMCP", 
                        family="binomial", alpha=0.5)
fit1.cmcp3 <- cv.grpreg(xtrain, ytrain, part, penalty="cMCP", 
                        family="binomial", alpha=0.95)

fit1.gel1 <- cv.grpreg(xtrain, ytrain, part, penalty="gel", family="binomial", 
                       alpha=0.05)
fit1.gel2 <- cv.grpreg(xtrain, ytrain, part, penalty="gel", family="binomial", 
                       alpha=0.5)
fit1.gel3 <- cv.grpreg(xtrain, ytrain, part, penalty="gel", family="binomial", 
                       alpha=0.95)

save(fit1.grridge, fit1.gren1, fit1.gren2, fit1.gren3, fit1.sgl1, fit1.sgl2,
     fit1.sgl3, fit1.cmcp1, fit1.cmcp2, fit1.cmcp3, fit1.gel1, fit1.gel2,
     fit1.gel3, file="results/metabolomics_alzheimer_fit1.Rdata")

### cross validation
set.seed(2019)
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