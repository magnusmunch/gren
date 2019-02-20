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
micrornas <- t(mirnormcen_resp)
colnames(micrornas) <- miRNA

### randomly split in test and train data (30% and 70%, respectively)
set.seed(2019)
id.train <- c(sample(which(benefit==0), size=floor(sum(benefit==0)*0.7)),
              sample(which(benefit==1), size=floor(sum(benefit==1)*0.7)))

################################### model 1 ####################################
### fitting the models
set.seed(2019)
ytrain <- benefit[id.train]
# we remove the constant micrornas
xtrain <- scale(micrornas[id.train, ])[, apply(micrornas[id.train, ], 2, sd)!=0]
utrain <- unpenal[id.train, ]
ytest <- benefit[-id.train]
xtest <- scale(micrornas[-id.train, ])[, apply(micrornas[id.train, ], 2, sd)!=0]
# set the constant micrornas to 0
xtest[is.nan(xtest)] <- 0
utest <- utrain[-id.train, ]
part <- diff.threegroup[apply(micrornas[id.train, ], 2, sd)!=0]
p <- ncol(xtrain)

fit1.gren1 <- gren(xtrain, ytrain, partitions=list(part=part), alpha=0.05, 
                   standardize=TRUE, trace=FALSE)
fit1.gren2 <- gren(xtrain, ytrain, partitions=list(part=part), alpha=0.5, 
                   standardize=TRUE, trace=FALSE)
fit1.gren3 <- gren(xtrain, ytrain, partitions=list(part=part), alpha=0.95, 
                   standardize=TRUE, trace=FALSE)

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
     fit1.gel3, file="results/micrornaseq_colorectal_cancer_fit1.Rdata")

### prediction on test set
pred1 <- data.frame(ridge=predict.grridge(fit1.grridge, t(xtest))[, 1], 
                    grridge=predict.grridge(fit1.grridge, t(xtest))[, 2],
                    gren1=predict(fit1.gren1$freq.model$groupreg, xtest, 
                                  type="response"),
                    
                    enet1=predict(fit1.gren1$freq.model$regular, xtest, 
                                  type="response"),
                    
                    sgl1=predictSGL(fit1.sgl1$fit, xtest),
                    
                    cmcp1=predict(fit1.cmcp1$fit, xtest, type="response"), 
                    
                    gel1=predict(fit1.gel1$fit, xtest, type="response")
                    )
psel1 <- c(ridge=p, grridge=p,
           gren1=fit1.gren1$freq.model$groupreg$df,
           
           enet1=fit1.gren1$freq.model$regular$df,
           
           sgl1=colSums(fit1.sgl1$fit$beta!=0),
           
           cmcp1=colSums(fit1.cmcp1$fit$beta[-1, ]!=0), 
           
           gel1=colSums(fit1.gel1$fit$beta[-1, ]!=0))
auc1 <- apply(pred1, 2, function(m) {pROC::auc(ytest, m)})
auc1.gren1 <- auc1[substr(names(auc1), 1, 4)=="gren"]
auc1.enet1 <- auc1[substr(names(auc1), 1, 4)=="enet"]
auc1.sgl1 <- auc1[substr(names(auc1), 1, 3)=="sgl"]
auc1.cmcp1 <- auc1[substr(names(auc1), 1, 4)=="cmcp"]
auc1.gel1 <- auc1[substr(names(auc1), 1, 3)=="gel"]
auc1.ridge <- auc1[substr(names(auc1), 1, 5)=="ridge"]
auc1.grridge <- auc1[substr(names(auc1), 1, 7)=="grridge"]

psel1.gren1 <- psel1[substr(names(psel1), 1, 4)=="gren"]
psel1.enet1 <- psel1[substr(names(psel1), 1, 4)=="enet"]
psel1.sgl1 <- psel1[substr(names(psel1), 1, 3)=="sgl"]
psel1.cmcp1 <- psel1[substr(names(psel1), 1, 4)=="cmcp"]
psel1.gel1 <- psel1[substr(names(psel1), 1, 3)=="gel"]
psel1.ridge <- psel1[substr(names(psel1), 1, 5)=="ridge"]
psel1.grridge <- psel1[substr(names(psel1), 1, 7)=="grridge"]

plot(psel1.gren1, auc1.gren1, col=1, type="l", xlim=range(c(psel1.gren1 ,
                                                            psel1.enet1 ,
                                                            psel1.sgl1 ,
                                                            psel1.cmcp1 ,
                                                            psel1.gel1 ,
                                                            psel1.ridge ,
                                                            psel1.grridge)),
     ylim=range(c(auc1.gren1 ,
                  auc1.enet1 ,
                  auc1.sgl1 , 
                  auc1.cmcp1 ,
                  auc1.gel1 , 
                  auc1.ridge ,
                  auc1.grridge)))
lines(psel1.enet1, auc1.enet1, col=2)
lines(psel1.sgl1, auc1.sgl1, col=3)
lines(psel1.cmcp1, auc1.cmcp1, col=4)
lines(psel1.gel1, auc1.gel1, col=5)
abline(h=auc1.ridge, lty=2, col=6)
abline(h=auc1.grridge, lty=2, col=7)


pred1 <- data.frame(ridge=predict.grridge(fit1.grridge, t(xtest))[, 1], 
                    grridge=predict.grridge(fit1.grridge, t(xtest))[, 2],
                    gren1=predict(fit1.gren1$freq.model$groupreg, xtest, 
                                  type="response"),
                    gren2=predict(fit1.gren2$freq.model$groupreg, xtest, 
                                  type="response"), 
                    gren3=predict(fit1.gren3$freq.model$groupreg, xtest, 
                                  type="response"), 
                    enet1=predict(fit1.gren1$freq.model$regular, xtest, 
                                  type="response"),
                    enet2=predict(fit1.gren2$freq.model$regular, xtest, 
                                  type="response"), 
                    enet3=predict(fit1.gren3$freq.model$regular, xtest, 
                                  type="response"), 
                    sgl1=predictSGL(fit1.sgl1$fit, xtest),
                    sgl2=predictSGL(fit1.sgl2$fit, xtest),
                    sgl3=predictSGL(fit1.sgl3$fit, xtest),
                    cmcp1=predict(fit1.cmcp1$fit, xtest, type="response"), 
                    cmcp2=predict(fit1.cmcp2$fit, xtest, type="response"), 
                    cmcp3=predict(fit1.cmcp2$fit, xtest, type="response"),
                    gel1=predict(fit1.gel1$fit, xtest, type="response"), 
                    gel2=predict(fit1.gel2$fit, xtest, type="response"), 
                    gel3=predict(fit1.gel3$fit, xtest, type="response"))
psel1 <- c(ridge=p, grridge=p,
           gren1=fit1.gren1$freq.model$groupreg$df,
           gren2=fit1.gren2$freq.model$groupreg$df,
           gren3=fit1.gren3$freq.model$groupreg$df,
           enet1=fit1.gren1$freq.model$regular$df,
           enet2=fit1.gren2$freq.model$regular$df,
           enet3=fit1.gren3$freq.model$regular$df,
           sgl1=colSums(fit1.sgl1$fit$beta!=0),
           sgl2=colSums(fit1.sgl2$fit$beta!=0),
           sgl3=colSums(fit1.sgl3$fit$beta!=0),
           cmcp1=colSums(fit1.cmcp1$fit$beta[-1, ]!=0), 
           cmcp2=colSums(fit1.cmcp2$fit$beta[-1, ]!=0), 
           cmcp3=colSums(fit1.cmcp3$fit$beta[-1, ]!=0),
           gel1=colSums(fit1.gel1$fit$beta[-1, ]!=0), 
           gel2=colSums(fit1.gel1$fit$beta[-1, ]!=0), 
           gel3=colSums(fit1.gel1$fit$beta[-1, ]!=0))
auc1 <- apply(pred1, 2, function(m) {pROC::auc(ytest, m)})
res1 <- rbind(pred1, psel1, auc1)
rownames(res1) <- c(paste0("pred", c(1:length(ytest))), paste0("psel", 1),
                    paste0("auc", 1))
write.table(res1, file="results/metabolomics_alzheimer_res1.csv")

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