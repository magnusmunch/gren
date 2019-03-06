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
library(Biobase)
library(gren)
library(GRridge)
library(grpregOverlap)
library(SGL)

### load data
load("data/ESetMbolCSFPR2.Rdata")
load("data/NwkDegreePartition.Rdata")

### select the patients with APOE E4 allel and controls without APOE E4
pheno <- pData(ESetMbolCSFPR2)
pheno.apoe <- pheno[(pheno$D_diag_name=="Probable AD" & pheno$APOE=="E4YES") |
                      (pheno$D_diag_name=="Subjectieve klachten" & 
                         pheno$APOE=="E4NO"), ]

metabol <- t(exprs(ESetMbolCSFPR2))
metabol.apoe <- metabol[(pheno$D_diag_name=="Probable AD" & 
                           pheno$APOE=="E4YES") |
                          (pheno$D_diag_name=="Subjectieve klachten" & 
                             pheno$APOE=="E4NO"), ]

### randomly split in test and train data (30% and 70%, respectively)
set.seed(2019)
id.train <- c(sample(which(pheno.apoe$D_diag_name=="Probable AD"), 
                     size=floor(sum(pheno.apoe$D_diag_name=="Probable AD")*
                                  0.7)),
              sample(which(pheno.apoe$D_diag_name=="Subjectieve klachten"), 
                     size=floor(sum(
                       pheno.apoe$D_diag_name=="Subjectieve klachten")*0.7)))

### create co-data
feat <- fData(ESetMbolCSFPR2)

# platformcode gives the type of metabolites
platformcode <- feat$PlatformCode2

# quality gives groups based on a quality measure of the assayed features
part.RSDqc <- CreatePartition(feat$RSDqc, ngroup=3)
quality <- rep(c(1:3), times=sapply(part.RSDqc, length))[
  order(unlist(part.RSDqc))]

# sds gives groups based on the standard deviations of the features
part.sds.apoe.train <- CreatePartition(apply(metabol.apoe[id.train, ], 1, sd), 
                                       ngroup=3)
sds.apoe.train <- rep(c(1:3), sapply(part.sds.apoe.train, length))[
  order(unlist(part.sds.apoe.train))]

# degree gives the degree classes of the differential network
degree <- NetworkDegreeClass[, 2]

### transform data
alzheim.apoe <- as.numeric(pheno.apoe$D_diag_name) - 1

################################### model 1 ####################################
### fitting the models
ytrain <- alzheim.apoe[id.train]
xtrain <- scale(metabol.apoe[id.train, ])
ytest <- alzheim.apoe[-id.train]
xtest <- scale(metabol.apoe[-id.train, ])
p <- ncol(xtrain)
part <- list(quality=quality, degree=degree)
part.gl <- as.numeric(as.factor(paste(quality, degree)))
part.ogl <- unlist(lapply(list(quality=split(1:p, part$quality), 
                               degree=split(1:p, part$degree)), function(part) {
                                 lapply(part, function(s) {
                                   colnames(xtrain)[s]})}), recursive=FALSE)

fit.gren1 <- gren(xtrain, ytrain, partitions=part, alpha=0.05, trace=FALSE)
fit.gren2 <- gren(xtrain, ytrain, partitions=part, alpha=0.5, trace=FALSE)
fit.gren3 <- gren(xtrain, ytrain, partitions=part, alpha=0.95, trace=FALSE)

fit.grridge <- grridge(t(xtrain), ytrain, list(quality=split(1:p, quality), 
                                               degree=split(1:p, degree)))

fit.sgl1 <- cvSGL(list(x=xtrain, y=ytrain), part.gl, type="logit", alpha=0.05)
fit.sgl1$fit$type <- "logit"
fit.sgl2 <- cvSGL(list(x=xtrain, y=ytrain), part.gl, type="logit", alpha=0.5)
fit.sgl2$fit$type <- "logit"
fit.sgl3 <- cvSGL(list(x=xtrain, y=ytrain), part.gl, type="logit", alpha=0.95)
fit.sgl3$fit$type <- "logit"

fit.cmcp1 <- cv.grpreg(xtrain, ytrain, part.gl, penalty="cMCP", 
                       family="binomial", alpha=0.05)
fit.cmcp2 <- cv.grpreg(xtrain, ytrain, part.gl, penalty="cMCP", 
                       family="binomial", alpha=0.5)
fit.cmcp3 <- cv.grpreg(xtrain, ytrain, part.gl, penalty="cMCP", 
                       family="binomial", alpha=0.95)

fit.gel1 <- cv.grpreg(xtrain, ytrain, part.gl, penalty="gel", 
                      family="binomial", alpha=0.05)
fit.gel2 <- cv.grpreg(xtrain, ytrain, part.gl, penalty="gel", 
                      family="binomial", alpha=0.5)
fit.gel3 <- cv.grpreg(xtrain, ytrain, part.gl, penalty="gel", 
                      family="binomial", alpha=0.95)

fit.ocmcp1 <- cv.grpregOverlap(xtrain, ytrain, part.ogl, penalty="cMCP", 
                               family="binomial", alpha=0.05)
fit.ocmcp2 <- cv.grpregOverlap(xtrain, ytrain, part.ogl, penalty="cMCP", 
                               family="binomial", alpha=0.5)
fit.ocmcp3 <- cv.grpregOverlap(xtrain, ytrain, part.ogl, penalty="cMCP", 
                               family="binomial", alpha=0.95)

fit.ogel1 <- cv.grpregOverlap(xtrain, ytrain, part.ogl, penalty="gel", 
                              family="binomial", alpha=0.05)
fit.ogel2 <- cv.grpregOverlap(xtrain, ytrain, part.ogl, penalty="gel", 
                              family="binomial", alpha=0.5)
fit.ogel3 <- cv.grpregOverlap(xtrain, ytrain, part.ogl, penalty="gel", 
                               family="binomial", alpha=0.95)

save(fit.grridge, fit.gren1, fit.gren2, fit.gren3, fit.sgl1, fit.sgl2,
     fit.sgl3, fit.cmcp1, fit.cmcp2, fit.cmcp3, fit.gel1, fit.gel2,
     fit.gel3, fit.ocmcp1, fit.ocmcp2, fit.ocmcp3, fit.ogel1, fit.ogel2, 
     fit.ogel3, file="results/metabolomics_alzheimer_fit1.Rdata")

### prediction on test set
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
                   ogel3=predict(fit.ogel3$fit, xtest, type="response"))
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
          ogel3=colSums(fit.ogel1$fit$beta[-1, ]!=0))
auc <- apply(pred, 2, function(m) {pROC::auc(ytest, m)})
briers <- apply(pred, 2, function(m) {
  1 - mean((m - ytest)^2)/mean((mean(ytest) - ytest)^2)})
res <- rbind(pred, psel, auc, briers)
rownames(res) <- c(paste0("pred", c(1:length(ytest))), paste0("psel", 1),
                    paste0("auc", 1), paste0("briers", 1))
write.table(res, file="results/metabolomics_alzheimer_res1.csv")
