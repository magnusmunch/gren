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
library(GRridge)
library(gren)
library(grpregOverlap)
library(SGL)
library(pROC)


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
  #partit <- TS
  vec <- array(dim=length(unlist(partit)))
  for(i in 1:length(partit)){
    vec[partit[[i]]] <- i
  }
  return(vec)
}
pvpart0 <- grridge2gren(CreatePartition(CDT$pvalsVUmc, ngroup=5,
                                        decreasing=FALSE))
corpart0 <- grridge2gren(CreatePartition(CDT$Corrs, ngroup=5, decreasing=TRUE))
#sds
#roepmanpart <- CDT$RoepmanGenes+1
partLNM <- list(corr=corpart0, pv=pvpart0)

################################### model 1 ####################################
### fitting the models
set.seed(2019)
ytrain <- RespTrain
xtrain <- Tr
ytest <- as.numeric(levels(RespValidation))[RespValidation]
xtest <- ValidationData[, selectgenes]
part <- partLNM
part.ogl <- unlist(lapply(list(corr=split(1:p, part$corr), 
                               pv=split(1:p, part$pv)), function(part) {
                                 lapply(part, function(s) {
                                   colnames(xtrain)[s]})}), recursive=FALSE)
n <- nrow(xtrain)
p <- ncol(xtrain)

bench1 <- microbenchmark(
fit1.gren1 <- gren(xtrain, ytrain, partitions=part, alpha=0.05, trace=FALSE),
fit1.gren2 <- gren(xtrain, ytrain, partitions=part, alpha=0.5, trace=FALSE),
fit1.gren3 <- gren(xtrain, ytrain, partitions=part, alpha=0.95, trace=FALSE),

fit1.grridge <- grridge(t(xtrain), ytrain, list(corr=split(1:p, part$corr), 
                                                pv=split(1:p, part$pv))), {

fit1.sgl1 <- cvSGL(list(x=xtrain, y=ytrain), as.numeric(as.factor(paste(
  corpart0, pvpart0))), type="logit", alpha=0.05)
fit1.sgl1$fit$type <- "logit"}, {
fit1.sgl2 <- cvSGL(list(x=xtrain, y=ytrain), as.numeric(as.factor(paste(
  corpart0, pvpart0))), type="logit", alpha=0.5)
fit1.sgl2$fit$type <- "logit"}, {
fit1.sgl3 <- cvSGL(list(x=xtrain, y=ytrain), as.numeric(as.factor(paste(
  corpart0, pvpart0))), type="logit", alpha=0.95)
fit1.sgl3$fit$type <- "logit"},

fit1.cmcp1 <- cv.grpreg(xtrain, ytrain, as.numeric(as.factor(paste(
  corpart0, pvpart0))), penalty="cMCP", family="binomial", alpha=0.05),
fit1.cmcp2 <- cv.grpreg(xtrain, ytrain, as.numeric(as.factor(paste(
  corpart0, pvpart0))), penalty="cMCP", family="binomial", alpha=0.5),
fit1.cmcp3 <- cv.grpreg(xtrain, ytrain, as.numeric(as.factor(paste(
  corpart0, pvpart0))), penalty="cMCP", family="binomial", alpha=0.95),

fit1.gel1 <- cv.grpreg(xtrain, ytrain, as.numeric(as.factor(paste(
  corpart0, pvpart0))), penalty="gel", family="binomial", alpha=0.05),
fit1.gel2 <- cv.grpreg(xtrain, ytrain, as.numeric(as.factor(paste(
  corpart0, pvpart0))), penalty="gel", family="binomial", alpha=0.5),
fit1.gel3 <- cv.grpreg(xtrain, ytrain, as.numeric(as.factor(paste(
  corpart0, pvpart0))), penalty="gel", family="binomial", alpha=0.95),

fit1.ocmcp1 <- cv.grpregOverlap(xtrain, ytrain, part.ogl, penalty="cMCP", 
                                family="binomial", alpha=0.05),
fit1.ocmcp2 <- cv.grpregOverlap(xtrain, ytrain, part.ogl, penalty="cMCP", 
                                family="binomial", alpha=0.5),
fit1.ocmcp3 <- cv.grpregOverlap(xtrain, ytrain, part.ogl, penalty="cMCP", 
                                family="binomial", alpha=0.95),

fit1.ogel1 <- cv.grpregOverlap(xtrain, ytrain, part.ogl, penalty="gel", 
                               family="binomial", alpha=0.05),
fit1.ogel2 <- cv.grpregOverlap(xtrain, ytrain, part.ogl, penalty="gel", 
                               family="binomial", alpha=0.5), {
# ogel3 gives just the saturated model, so we create cv.grpregOverlap ourselves
fit1.ogel3 <- grpregOverlap(xtrain, ytrain, part.ogl, penalty="gel", 
                            family="binomial", alpha=0.95)
fit1.ogel3 <- list(cve=NA, cvse=NA, lambda=fit1.ogel3$lambda, fit=fit1.ogel3, 
                   min=1, lambda.min=fit1.ogel3$lambda, null.dev=NA, pe=NA)
class(fit1.ogel3) <- c("cv.grpregOverlap", "cv.grpreg")}, times=1, 
control=list(order="inorder"))

rownames(bench1) <- c(paste0("gren", 1:3), "grridge", paste0("sgl", 1:3),
                      paste0("cmcp", 1:3), paste0("gel", 1:3),
                      paste0("ocmcp", 1:3), paste0("ogel", 1:3))
write.table(bench1, file="results/rnaseq_lymph_node_metastasis_bench1.csv")

save(fit1.grridge, fit1.gren1, fit1.gren2, fit1.gren3, fit1.sgl1, fit1.sgl2,
     fit1.sgl3, fit1.cmcp1, fit1.cmcp2, fit1.cmcp3, fit1.gel1, fit1.gel2,
     fit1.gel3, file="results/rnaseq_lymph_node_metastasis_fit1.Rdata")

### obtaining predictions for test data
pred1 <- data.frame(ridge=predict.grridge(fit1.grridge, t(xtest))[, 1], 
                    grridge=predict.grridge(fit1.grridge, t(xtest))[, 2],
                    gren1=predict(fit1.gren1, xtest, type="groupreg"),
                    gren2=predict(fit1.gren2, xtest, type="groupreg"), 
                    gren3=predict(fit1.gren3, xtest, type="groupreg"), 
                    enet1=predict(fit1.gren1, xtest, type="regular"),
                    enet2=predict(fit1.gren2, xtest, type="regular"), 
                    enet3=predict(fit1.gren3, xtest, type="regular"), 
                    sgl1=predictSGL(fit1.sgl1$fit, xtest),
                    sgl2=predictSGL(fit1.sgl2$fit, xtest),
                    sgl3=predictSGL(fit1.sgl3$fit, xtest),
                    cmcp1=predict(fit1.cmcp1$fit, xtest, type="response"), 
                    cmcp2=predict(fit1.cmcp2$fit, xtest, type="response"), 
                    cmcp3=predict(fit1.cmcp2$fit, xtest, type="response"),
                    gel1=predict(fit1.gel1$fit, xtest, type="response"), 
                    gel2=predict(fit1.gel2$fit, xtest, type="response"), 
                    gel3=predict(fit1.gel3$fit, xtest, type="response"),
                    ocmcp1=predict(fit1.ocmcp1$fit, xtest, type="response"), 
                    ocmcp2=predict(fit1.ocmcp2$fit, xtest, type="response"), 
                    ocmcp3=predict(fit1.ocmcp2$fit, xtest, type="response"),
                    ogel1=predict(fit1.ogel1$fit, xtest, type="response"), 
                    ogel2=predict(fit1.ogel2$fit, xtest, type="response"), 
                    ogel3=predict(fit1.ogel3$fit, xtest, type="response"))
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
           gel3=colSums(fit1.gel1$fit$beta[-1, ]!=0),
           ocmcp1=colSums(fit1.ocmcp1$fit$beta[-1, ]!=0), 
           ocmcp2=colSums(fit1.ocmcp2$fit$beta[-1, ]!=0), 
           ocmcp3=colSums(fit1.ocmcp3$fit$beta[-1, ]!=0),
           ogel1=colSums(fit1.ogel1$fit$beta[-1, ]!=0), 
           ogel2=colSums(fit1.ogel1$fit$beta[-1, ]!=0), 
           ogel3=colSums(fit1.ogel1$fit$beta[-1, ]!=0))
auc1 <- apply(pred1, 2, function(m) {pROC::auc(ytest, m)})
res1 <- rbind(pred1, psel1, auc1)
rownames(res1) <- c(paste0("pred", c(1:length(ytest))), paste0("psel", 1),
                    paste0("auc", 1))
write.table(res1, file="results/rnaseq_lymph_node_metastasis_res1.csv")
