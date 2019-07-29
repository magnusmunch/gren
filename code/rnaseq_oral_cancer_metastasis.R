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
library(CoRF)
library(gren)
library(GRridge)
library(grpregOverlap)
library(SGL)
library(randomForestSRC)
library(freeknotsplines)
library(microbenchmark)

### load data
data(LNM_Example)

### data preparation
pv <- CoDataTrain$pvalsVUmc
selectgenes <- which(pv<=0.1)
p <- length(selectgenes)
CDT <- CoDataTrain[selectgenes, ]
Tr <- TrainData[, selectgenes]
n <- nrow(Tr)

grridge2gren <- function(partit){
  vec <- array(dim=length(unlist(partit)))
  for(i in 1:length(partit)){
    vec[partit[[i]]] <- i
  }
  return(vec)
}
pvpart0 <- grridge2gren(CreatePartition(CDT$pvalsVUmc, ngroup=5,
                                        decreasing=FALSE))
corpart0 <- grridge2gren(CreatePartition(CDT$Corrs, ngroup=5, decreasing=TRUE))
partLNM <- list(corr=corpart0, pv=pvpart0)

kn1 <- knots(ecdf(CDT$Corrs))
kn2 <- knots(ecdf(CDT$pvalsVUmc))
spl1 <- freelsgen(kn1, 1:length(kn1), degree=1, numknot=4, seed=2019, stream=0)
spl2 <- freelsgen(kn2, 1:length(kn2), degree=1, numknot=4, seed=2019, stream=0)

corr2 <- (CDT$Corrs <= spl1@optknot[1]) + (CDT$Corrs <= spl1@optknot[2]) +
  (CDT$Corrs <= spl1@optknot[3]) + (CDT$Corrs <= spl1@optknot[4]) + 1
pv2 <- (CDT$pvalsVUmc <= spl2@optknot[1]) + (CDT$pvalsVUmc <= spl2@optknot[2]) +
  (CDT$pvalsVUmc <= spl2@optknot[3]) + (CDT$pvalsVUmc <= spl2@optknot[4]) + 1

################################### model 1 ####################################
### fitting the models
set.seed(2019)
ytrain <- RespTrain
xtrain <- Tr
ytest <- as.numeric(levels(RespValidation))[RespValidation]
xtest <- ValidationData[, selectgenes]
part <- partLNM
part2 <- list(corr=corr2, pv=pv2)
part.ogl <- unlist(lapply(list(corr=split(1:p, part$corr), 
                               pv=split(1:p, part$pv)), function(part) {
                                 lapply(part, function(s) {
                                   colnames(xtrain)[s]})}), recursive=FALSE)
n <- nrow(xtrain)
p <- ncol(xtrain)

bench <- microbenchmark(
fit.gren1 <- gren(xtrain, ytrain, partitions=part, alpha=0.05, trace=FALSE),
fit.gren2 <- gren(xtrain, ytrain, partitions=part, alpha=0.5, trace=FALSE),
fit.gren3 <- gren(xtrain, ytrain, partitions=part, alpha=0.95, trace=FALSE),

fit.grridge <- grridge(t(xtrain), ytrain, list(corr=split(1:p, part$corr), 
                                                pv=split(1:p, part$pv))), {

fit.sgl1 <- cvSGL(list(x=xtrain, y=ytrain), as.numeric(as.factor(paste(
  corpart0, pvpart0))), type="logit", alpha=0.05)
fit.sgl1$fit$type <- "logit"}, {
fit.sgl2 <- cvSGL(list(x=xtrain, y=ytrain), as.numeric(as.factor(paste(
  corpart0, pvpart0))), type="logit", alpha=0.5)
fit.sgl2$fit$type <- "logit"}, {
fit.sgl3 <- cvSGL(list(x=xtrain, y=ytrain), as.numeric(as.factor(paste(
  corpart0, pvpart0))), type="logit", alpha=0.95)
fit.sgl3$fit$type <- "logit"},

fit.cmcp1 <- cv.grpreg(xtrain, ytrain, as.numeric(as.factor(paste(
  corpart0, pvpart0))), penalty="cMCP", family="binomial", alpha=0.05),
fit.cmcp2 <- cv.grpreg(xtrain, ytrain, as.numeric(as.factor(paste(
  corpart0, pvpart0))), penalty="cMCP", family="binomial", alpha=0.5),
fit.cmcp3 <- cv.grpreg(xtrain, ytrain, as.numeric(as.factor(paste(
  corpart0, pvpart0))), penalty="cMCP", family="binomial", alpha=0.95),

fit.gel1 <- cv.grpreg(xtrain, ytrain, as.numeric(as.factor(paste(
  corpart0, pvpart0))), penalty="gel", family="binomial", alpha=0.05),
fit.gel2 <- cv.grpreg(xtrain, ytrain, as.numeric(as.factor(paste(
  corpart0, pvpart0))), penalty="gel", family="binomial", alpha=0.5),
fit.gel3 <- cv.grpreg(xtrain, ytrain, as.numeric(as.factor(paste(
  corpart0, pvpart0))), penalty="gel", family="binomial", alpha=0.95),

fit.ocmcp1 <- cv.grpregOverlap(xtrain, ytrain, part.ogl, penalty="cMCP", 
                                family="binomial", alpha=0.05),
fit.ocmcp2 <- cv.grpregOverlap(xtrain, ytrain, part.ogl, penalty="cMCP", 
                                family="binomial", alpha=0.5),
fit.ocmcp3 <- cv.grpregOverlap(xtrain, ytrain, part.ogl, penalty="cMCP", 
                                family="binomial", alpha=0.95),

fit.ogel1 <- cv.grpregOverlap(xtrain, ytrain, part.ogl, penalty="gel", 
                               family="binomial", alpha=0.05),
fit.ogel2 <- cv.grpregOverlap(xtrain, ytrain, part.ogl, penalty="gel", 
                               family="binomial", alpha=0.5), {
# ogel3 gives just the saturated model, so we create cv.grpregOverlap ourselves
fit.ogel3 <- grpregOverlap(xtrain, ytrain, part.ogl, penalty="gel", 
                            family="binomial", alpha=0.95)
fit.ogel3 <- list(cve=NA, cvse=NA, lambda=fit.ogel3$lambda, fit=fit.ogel3, 
                   min=1, lambda.min=fit.ogel3$lambda, null.dev=NA, pe=NA)
class(fit.ogel3) <- c("cv.grpregOverlap", "cv.grpreg")}, 

fit.rf <- rfsrc(y ~ ., data=data.frame(y=ytrain, x=xtrain), 
                var.used="all.trees", ntree=5000, importance="none"),
times=1, control=list(order="inorder"))

fit.gren4 <- gren(xtrain, ytrain, partitions=part2, alpha=0.05, trace=FALSE)
fit.gren5 <- gren(xtrain, ytrain, partitions=part2, alpha=0.5, trace=FALSE)
fit.gren6 <- gren(xtrain, ytrain, partitions=part2, alpha=0.95, trace=FALSE)

fit.grridge2 <- grridge(t(xtrain), ytrain, list(corr=split(1:p, part2$corr), 
                                                pv=split(1:p, part2$pv)))

rownames(bench) <- c(paste0("gren", 1:3), "grridge", paste0("sgl", 1:3),
                     paste0("cmcp", 1:3), paste0("gel", 1:3),
                     paste0("ocmcp", 1:3), paste0("ogel", 1:3), "rf")
write.table(bench, file="results/rnaseq_oral_cancer_metastasis_bench1.csv")

save(fit.grridge, fit.gren1, fit.gren2, fit.gren3, fit.sgl1, fit.sgl2,
     fit.sgl3, fit.cmcp1, fit.cmcp2, fit.cmcp3, fit.gel1, fit.gel2,
     fit.gel3, fit.rf, fit.gren4, fit.gren5, fit.gren6, fit.grridge2,
     file="results/rnaseq_oral_cancer_metastasis_fit1.Rdata")

### obtaining predictions for test data
pred <- data.frame(ridge=predict.grridge(fit.grridge, t(xtest))[, 1], 
                   grridge=predict.grridge(fit.grridge, t(xtest))[, 2],
                   gren1=predict(fit.gren1, xtest, type="groupreg"),
                   gren2=predict(fit.gren2, xtest, type="groupreg"), 
                   gren3=predict(fit.gren3, xtest, type="groupreg"), 
                   enet1=predict(fit.gren1, xtest, type="regular"),
                   enet2=predict(fit.gren2, xtest, type="regular"), 
                   enet3=predict(fit.gren3, xtest, type="regular"), 
                   sgl1=predictSGL(fit.sgl1$fit, xtest),
                   sgl2=predictSGL(fit.sgl2$fit, xtest),
                   sgl3=predictSGL(fit.sgl3$fit, xtest),
                   cmcp1=predict(fit.cmcp1$fit, xtest, type="response"), 
                   cmcp2=predict(fit.cmcp2$fit, xtest, type="response"), 
                   cmcp3=predict(fit.cmcp2$fit, xtest, type="response"),
                   gel1=predict(fit.gel1$fit, xtest, type="response"), 
                   gel2=predict(fit.gel2$fit, xtest, type="response"), 
                   gel3=predict(fit.gel3$fit, xtest, type="response"),
                   ocmcp1=predict(fit.ocmcp1$fit, xtest, type="response"), 
                   ocmcp2=predict(fit.ocmcp2$fit, xtest, type="response"), 
                   ocmcp3=predict(fit.ocmcp2$fit, xtest, type="response"),
                   ogel1=predict(fit.ogel1$fit, xtest, type="response"), 
                   ogel2=predict(fit.ogel2$fit, xtest, type="response"), 
                   ogel3=predict(fit.ogel3$fit, xtest, type="response"),
                   rf=predict(fit.rf, data.frame(x=xtest))$predicted,
                   gren4=predict(fit.gren4, xtest, type="groupreg"),
                   gren5=predict(fit.gren5, xtest, type="groupreg"),
                   gren6=predict(fit.gren6, xtest, type="groupreg"),
                   grridge2=predict.grridge(fit.grridge2, t(xtest))[, 2])
psel <- c(ridge=p, grridge=p,
          gren1=fit.gren1$freq.model$groupreg$df,
          gren2=fit.gren2$freq.model$groupreg$df,
          gren3=fit.gren3$freq.model$groupreg$df,
          enet1=fit.gren1$freq.model$regular$df,
          enet2=fit.gren2$freq.model$regular$df,
          enet3=fit.gren3$freq.model$regular$df,
          sgl1=colSums(fit.sgl1$fit$beta!=0),
          sgl2=colSums(fit.sgl2$fit$beta!=0),
          sgl3=colSums(fit.sgl3$fit$beta!=0),
          cmcp1=colSums(fit.cmcp1$fit$beta[-1, ]!=0), 
          cmcp2=colSums(fit.cmcp2$fit$beta[-1, ]!=0), 
          cmcp3=colSums(fit.cmcp3$fit$beta[-1, ]!=0),
          gel1=colSums(fit.gel1$fit$beta[-1, ]!=0), 
          gel2=colSums(fit.gel1$fit$beta[-1, ]!=0), 
          gel3=colSums(fit.gel1$fit$beta[-1, ]!=0),
          ocmcp1=colSums(fit.ocmcp1$fit$beta[-1, ]!=0), 
          ocmcp2=colSums(fit.ocmcp2$fit$beta[-1, ]!=0), 
          ocmcp3=colSums(fit.ocmcp3$fit$beta[-1, ]!=0),
          ogel1=colSums(fit.ogel1$fit$beta[-1, ]!=0), 
          ogel2=colSums(fit.ogel1$fit$beta[-1, ]!=0), 
          ogel3=colSums(fit.ogel1$fit$beta[-1, ]!=0),
          rf=p,
          gren4=fit.gren4$freq.model$groupreg$df,
          gren5=fit.gren5$freq.model$groupreg$df,
          gren6=fit.gren6$freq.model$groupreg$df,
          grridge2=p)
auc <- apply(pred, 2, function(m) {pROC::auc(ytest, m)})
briers <- apply(pred, 2, function(m) {
  1 - mean((m - ytest)^2)/mean((mean(ytest) - ytest)^2)})
res <- rbind(pred, psel, auc, briers)
rownames(res) <- c(paste0("pred", c(1:length(ytest))), paste0("psel", 1),
                    paste0("auc", 1), paste0("briers", 1))
write.table(res, file="results/rnaseq_oral_cancer_metastasis_res1.csv")
