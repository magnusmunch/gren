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

# quality gives groups based on a quality measure of the assayd features
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
set.seed(2019)
ytrain <- alzheim.apoe[id.train]
xtrain <- scale(metabol.apoe[id.train, ])
ytest <- alzheim.apoe[-id.train]
xtest <- scale(metabol.apoe[-id.train, ])
part <- platformcode
p <- ncol(xtrain)

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

### prediction on test set
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


################################### model 2 ####################################
### fitting the models
set.seed(2018)
ytrain <- alzheim.apoe[id.train]
xtrain <- scale(metabol.apoe[id.train, ])
ytest <- alzheim.apoe[-id.train]
xtest <- scale(metabol.apoe[-id.train, ])
part <- quality
p <- ncol(xtrain)

fit2.gren1 <- gren(xtrain, ytrain, partitions=list(part=part), alpha=0.05, 
                   trace=FALSE)
fit2.gren2 <- gren(xtrain, ytrain, partitions=list(part=part), alpha=0.5, 
                   trace=FALSE)
fit2.gren3 <- gren(xtrain, ytrain, partitions=list(part=part), alpha=0.95, 
                   trace=FALSE)

fit2.grridge <- grridge(t(xtrain), ytrain, list(part=split(1:p, part)))

fit2.sgl1 <- cvSGL(list(x=xtrain, y=ytrain), part, type="logit", alpha=0.05)
fit2.sgl1$fit$type <- "logit"
fit2.sgl2 <- cvSGL(list(x=xtrain, y=ytrain), part, type="logit", alpha=0.5)
fit2.sgl2$fit$type <- "logit"
fit2.sgl3 <- cvSGL(list(x=xtrain, y=ytrain), part, type="logit", alpha=0.95)
fit2.sgl3$fit$type <- "logit"

fit2.cmcp1 <- cv.grpreg(xtrain, ytrain, part, penalty="cMCP", 
                        family="binomial", alpha=0.05)
fit2.cmcp2 <- cv.grpreg(xtrain, ytrain, part, penalty="cMCP", 
                        family="binomial", alpha=0.5)
fit2.cmcp3 <- cv.grpreg(xtrain, ytrain, part, penalty="cMCP", 
                        family="binomial", alpha=0.95)

fit2.gel1 <- cv.grpreg(xtrain, ytrain, part, penalty="gel", family="binomial", 
                       alpha=0.05)
fit2.gel2 <- cv.grpreg(xtrain, ytrain, part, penalty="gel", family="binomial", 
                       alpha=0.5)
fit2.gel3 <- cv.grpreg(xtrain, ytrain, part, penalty="gel", family="binomial", 
                       alpha=0.95)

save(fit2.grridge, fit2.gren1, fit2.gren2, fit2.gren3, fit2.sgl1, fit2.sgl2,
     fit2.sgl3, fit2.cmcp1, fit2.cmcp2, fit2.cmcp3, fit2.gel1, fit2.gel2,
     fit2.gel3, file="results/metabolomics_alzheimer_fit2.Rdata")

### prediction on test set
pred2 <- data.frame(ridge=predict.grridge(fit2.grridge, t(xtest))[, 1], 
                    grridge=predict.grridge(fit2.grridge, t(xtest))[, 2],
                    gren1=predict(fit2.gren1$freq.model$groupreg, xtest, 
                                  type="response"),
                    gren2=predict(fit2.gren2$freq.model$groupreg, xtest, 
                                  type="response"), 
                    gren3=predict(fit2.gren3$freq.model$groupreg, xtest, 
                                  type="response"), 
                    enet1=predict(fit2.gren1$freq.model$regular, xtest, 
                                  type="response"),
                    enet2=predict(fit2.gren2$freq.model$regular, xtest, 
                                  type="response"), 
                    enet3=predict(fit2.gren3$freq.model$regular, xtest, 
                                  type="response"), 
                    sgl1=predictSGL(fit2.sgl1$fit, xtest),
                    sgl2=predictSGL(fit2.sgl2$fit, xtest),
                    sgl3=predictSGL(fit2.sgl3$fit, xtest),
                    cmcp1=predict(fit2.cmcp1$fit, xtest, type="response"), 
                    cmcp2=predict(fit2.cmcp2$fit, xtest, type="response"), 
                    cmcp3=predict(fit2.cmcp2$fit, xtest, type="response"),
                    gel1=predict(fit2.gel1$fit, xtest, type="response"), 
                    gel2=predict(fit2.gel2$fit, xtest, type="response"), 
                    gel3=predict(fit2.gel3$fit, xtest, type="response"))
psel2 <- c(ridge=p, grridge=p,
           gren1=fit2.gren1$freq.model$groupreg$df,
           gren2=fit2.gren2$freq.model$groupreg$df,
           gren3=fit2.gren3$freq.model$groupreg$df,
           enet1=fit2.gren1$freq.model$regular$df,
           enet2=fit2.gren2$freq.model$regular$df,
           enet3=fit2.gren3$freq.model$regular$df,
           sgl1=colSums(fit2.sgl1$fit$beta!=0),
           sgl2=colSums(fit2.sgl2$fit$beta!=0),
           sgl3=colSums(fit2.sgl3$fit$beta!=0),
           cmcp1=colSums(fit2.cmcp1$fit$beta[-1, ]!=0), 
           cmcp2=colSums(fit2.cmcp2$fit$beta[-1, ]!=0), 
           cmcp3=colSums(fit2.cmcp3$fit$beta[-1, ]!=0),
           gel1=colSums(fit2.gel1$fit$beta[-1, ]!=0), 
           gel2=colSums(fit2.gel1$fit$beta[-1, ]!=0), 
           gel3=colSums(fit2.gel1$fit$beta[-1, ]!=0))
auc2 <- apply(pred2, 2, function(m) {pROC::auc(ytest, m)})
res2 <- rbind(pred2, psel2, auc2)
rownames(res2) <- c(paste0("pred", c(1:length(ytest))), paste0("psel", 1),
                    paste0("auc", 1))
write.table(res2, file="results/metabolomics_alzheimer_res2.csv")


################################### model 3 ####################################
### fitting the models
set.seed(2018)
ytrain <- alzheim.apoe[id.train]
xtrain <- scale(metabol.apoe[id.train, ])
ytest <- alzheim.apoe[-id.train]
xtest <- scale(metabol.apoe[-id.train, ])
part <- degree
p <- ncol(xtrain)

fit3.gren1 <- gren(xtrain, ytrain, partitions=list(part=part), alpha=0.05, 
                   trace=FALSE)
fit3.gren2 <- gren(xtrain, ytrain, partitions=list(part=part), alpha=0.5, 
                   trace=FALSE)
fit3.gren3 <- gren(xtrain, ytrain, partitions=list(part=part), alpha=0.95, 
                   trace=FALSE)

fit3.grridge <- grridge(t(xtrain), ytrain, list(part=split(1:p, part)))

fit3.sgl1 <- cvSGL(list(x=xtrain, y=ytrain), part, type="logit", alpha=0.05)
fit3.sgl1$fit$type <- "logit"
fit3.sgl2 <- cvSGL(list(x=xtrain, y=ytrain), part, type="logit", alpha=0.5)
fit3.sgl2$fit$type <- "logit"
fit3.sgl3 <- cvSGL(list(x=xtrain, y=ytrain), part, type="logit", alpha=0.95)
fit3.sgl3$fit$type <- "logit"

fit3.cmcp1 <- cv.grpreg(xtrain, ytrain, part, penalty="cMCP", 
                        family="binomial", alpha=0.05)
fit3.cmcp2 <- cv.grpreg(xtrain, ytrain, part, penalty="cMCP", 
                        family="binomial", alpha=0.5)
fit3.cmcp3 <- cv.grpreg(xtrain, ytrain, part, penalty="cMCP", 
                        family="binomial", alpha=0.95)

fit3.gel1 <- cv.grpreg(xtrain, ytrain, part, penalty="gel", family="binomial", 
                       alpha=0.05)
fit3.gel2 <- cv.grpreg(xtrain, ytrain, part, penalty="gel", family="binomial", 
                       alpha=0.5)
fit3.gel3 <- cv.grpreg(xtrain, ytrain, part, penalty="gel", family="binomial", 
                       alpha=0.95)

save(fit3.grridge, fit3.gren1, fit3.gren2, fit3.gren3, fit3.sgl1, fit3.sgl2,
     fit3.sgl3, fit3.cmcp1, fit3.cmcp2, fit3.cmcp3, fit3.gel1, fit3.gel2,
     fit3.gel3, file="results/metabolomics_alzheimer_fit3.Rdata")

### prediction on test set
pred3 <- data.frame(ridge=predict.grridge(fit3.grridge, t(xtest))[, 1], 
                    grridge=predict.grridge(fit3.grridge, t(xtest))[, 2],
                    gren1=predict(fit3.gren1$freq.model$groupreg, xtest, 
                                  type="response"),
                    gren2=predict(fit3.gren2$freq.model$groupreg, xtest, 
                                  type="response"), 
                    gren3=predict(fit3.gren3$freq.model$groupreg, xtest, 
                                  type="response"), 
                    enet1=predict(fit3.gren1$freq.model$regular, xtest, 
                                  type="response"),
                    enet2=predict(fit3.gren2$freq.model$regular, xtest, 
                                  type="response"), 
                    enet3=predict(fit3.gren3$freq.model$regular, xtest, 
                                  type="response"), 
                    sgl1=predictSGL(fit3.sgl1$fit, xtest),
                    sgl2=predictSGL(fit3.sgl2$fit, xtest),
                    sgl3=predictSGL(fit3.sgl3$fit, xtest),
                    cmcp1=predict(fit3.cmcp1$fit, xtest, type="response"), 
                    cmcp2=predict(fit3.cmcp2$fit, xtest, type="response"), 
                    cmcp3=predict(fit3.cmcp2$fit, xtest, type="response"),
                    gel1=predict(fit3.gel1$fit, xtest, type="response"), 
                    gel2=predict(fit3.gel2$fit, xtest, type="response"), 
                    gel3=predict(fit3.gel3$fit, xtest, type="response"))
psel3 <- c(ridge=p, grridge=p,
           gren1=fit3.gren1$freq.model$groupreg$df,
           gren2=fit3.gren2$freq.model$groupreg$df,
           gren3=fit3.gren3$freq.model$groupreg$df,
           enet1=fit3.gren1$freq.model$regular$df,
           enet2=fit3.gren2$freq.model$regular$df,
           enet3=fit3.gren3$freq.model$regular$df,
           sgl1=colSums(fit3.sgl1$fit$beta!=0),
           sgl2=colSums(fit3.sgl2$fit$beta!=0),
           sgl3=colSums(fit3.sgl3$fit$beta!=0),
           cmcp1=colSums(fit3.cmcp1$fit$beta[-1, ]!=0), 
           cmcp2=colSums(fit3.cmcp2$fit$beta[-1, ]!=0), 
           cmcp3=colSums(fit3.cmcp3$fit$beta[-1, ]!=0),
           gel1=colSums(fit3.gel1$fit$beta[-1, ]!=0), 
           gel2=colSums(fit3.gel1$fit$beta[-1, ]!=0), 
           gel3=colSums(fit3.gel1$fit$beta[-1, ]!=0))
auc3 <- apply(pred3, 2, function(m) {pROC::auc(ytest, m)})
res3 <- rbind(pred3, psel3, auc3)
rownames(res3) <- c(paste0("pred", c(1:length(ytest))), paste0("psel", 1),
                    paste0("auc", 1))
write.table(res3, file="results/metabolomics_alzheimer_res3.csv")


################################### model 4 ####################################
### fitting the models
set.seed(2018)
ytrain <- alzheim.apoe[id.train]
xtrain <- scale(metabol.apoe[id.train, ])
ytest <- alzheim.apoe[-id.train]
xtest <- scale(metabol.apoe[-id.train, ])
p <- ncol(xtrain)
part.gren <- list(platformcode=platformcode, quality=quality, degree=degree)
part.grridge <- list(platformcode=split(1:p, platformcode), 
                     quality=split(1:p, quality),
                     degree=split(1:p, degree))
part.gl <- as.numeric(factor(paste(platformcode, quality, degree)))
part.ogl <- unlist(lapply(part.grridge, function(part) {
  lapply(part, function(s) {colnames(xtrain)[s]})}), recursive=FALSE)

fit4.gren1 <- gren(xtrain, ytrain, partitions=part.gren, alpha=0.05, 
                   trace=FALSE)
fit4.gren2 <- gren(xtrain, ytrain, partitions=part.gren, alpha=0.5, 
                   trace=FALSE)
fit4.gren3 <- gren(xtrain, ytrain, partitions=part.gren, alpha=0.95, 
                   trace=FALSE)

fit4.grridge <- grridge(t(xtrain), ytrain, part.grridge)

fit4.sgl1 <- cvSGL(list(x=xtrain, y=ytrain), part.gl, type="logit", alpha=0.05)
fit4.sgl1$fit$type <- "logit"
fit4.sgl2 <- cvSGL(list(x=xtrain, y=ytrain), part.gl, type="logit", alpha=0.5)
fit4.sgl2$fit$type <- "logit"
fit4.sgl3 <- cvSGL(list(x=xtrain, y=ytrain), part.gl, type="logit", alpha=0.95)
fit4.sgl3$fit$type <- "logit"

fit4.cmcp1 <- cv.grpreg(xtrain, ytrain, part.gl, penalty="cMCP", 
                        family="binomial", alpha=0.05)
fit4.cmcp2 <- cv.grpreg(xtrain, ytrain, part.gl, penalty="cMCP", 
                        family="binomial", alpha=0.5)
fit4.cmcp3 <- cv.grpreg(xtrain, ytrain, part.gl, penalty="cMCP", 
                        family="binomial", alpha=0.95)

fit4.gel1 <- cv.grpreg(xtrain, ytrain, part.gl, penalty="gel", 
                       family="binomial", alpha=0.05)
fit4.gel2 <- cv.grpreg(xtrain, ytrain, part.gl, penalty="gel", 
                       family="binomial", alpha=0.5)
fit4.gel3 <- cv.grpreg(xtrain, ytrain, part.gl, penalty="gel", 
                       family="binomial", alpha=0.95)

fit4.ocmcp1 <- cv.grpregOverlap(xtrain, ytrain, part.ogl, penalty="cMCP", 
                                family="binomial", alpha=0.05)
fit4.ocmcp2 <- cv.grpregOverlap(xtrain, ytrain, part.ogl, penalty="cMCP", 
                                family="binomial", alpha=0.5)
fit4.ocmcp3 <- cv.grpregOverlap(xtrain, ytrain, part.ogl, penalty="cMCP", 
                                family="binomial", alpha=0.95)

fit4.ogel1 <- cv.grpregOverlap(xtrain, ytrain, part.ogl, penalty="gel", 
                                family="binomial", alpha=0.05)
fit4.ogel2 <- cv.grpregOverlap(xtrain, ytrain, part.ogl, penalty="gel", 
                                family="binomial", alpha=0.5)
fit4.ogel3 <- cv.grpregOverlap(xtrain, ytrain, part.ogl, penalty="gel", 
                                family="binomial", alpha=0.95)

save(fit4.grridge, fit4.gren1, fit4.gren2, fit4.gren3, fit4.sgl1, fit4.sgl2,
     fit4.sgl3, fit4.cmcp1, fit4.cmcp2, fit4.cmcp3, fit4.gel1, fit4.gel2,
     fit4.gel3, fit4.ocmcp1, fit4.ocmcp2, fit4.ocmcp3, fit4.ogel1, fit4.ogel2, 
     fit4.ogel3, file="results/metabolomics_alzheimer_fit4.Rdata")

### prediction on test set
pred4 <- data.frame(ridge=predict.grridge(fit4.grridge, t(xtest))[, 1], 
                    grridge=predict.grridge(fit4.grridge, t(xtest))[, 2],
                    gren1=predict(fit4.gren1$freq.model$groupreg, xtest, 
                                  type="response"),
                    gren2=predict(fit4.gren2$freq.model$groupreg, xtest, 
                                  type="response"), 
                    gren3=predict(fit4.gren3$freq.model$groupreg, xtest, 
                                  type="response"), 
                    enet1=predict(fit4.gren1$freq.model$regular, xtest, 
                                  type="response"),
                    enet2=predict(fit4.gren2$freq.model$regular, xtest, 
                                  type="response"), 
                    enet3=predict(fit4.gren3$freq.model$regular, xtest, 
                                  type="response"), 
                    sgl1=predictSGL(fit4.sgl1$fit, xtest),
                    sgl2=predictSGL(fit4.sgl2$fit, xtest),
                    sgl3=predictSGL(fit4.sgl3$fit, xtest),
                    cmcp1=predict(fit4.cmcp1$fit, xtest, type="response"), 
                    cmcp2=predict(fit4.cmcp2$fit, xtest, type="response"), 
                    cmcp3=predict(fit4.cmcp2$fit, xtest, type="response"),
                    gel1=predict(fit4.gel1$fit, xtest, type="response"), 
                    gel2=predict(fit4.gel2$fit, xtest, type="response"), 
                    gel3=predict(fit4.gel3$fit, xtest, type="response"),
                    ocmcp1=predict(fit4.ocmcp1$fit, xtest, type="response"), 
                    ocmcp2=predict(fit4.ocmcp2$fit, xtest, type="response"), 
                    ocmcp3=predict(fit4.ocmcp2$fit, xtest, type="response"),
                    ogel1=predict(fit4.ogel1$fit, xtest, type="response"), 
                    ogel2=predict(fit4.ogel2$fit, xtest, type="response"), 
                    ogel3=predict(fit4.ogel3$fit, xtest, type="response"))
psel4 <- c(ridge=p, grridge=p,
           gren1=fit4.gren1$freq.model$groupreg$df,
           gren2=fit4.gren2$freq.model$groupreg$df,
           gren3=fit4.gren3$freq.model$groupreg$df,
           enet1=fit4.gren1$freq.model$regular$df,
           enet2=fit4.gren2$freq.model$regular$df,
           enet3=fit4.gren3$freq.model$regular$df,
           sgl1=colSums(fit4.sgl1$fit$beta!=0),
           sgl2=colSums(fit4.sgl2$fit$beta!=0),
           sgl3=colSums(fit4.sgl3$fit$beta!=0),
           cmcp1=colSums(fit4.cmcp1$fit$beta[-1, ]!=0), 
           cmcp2=colSums(fit4.cmcp2$fit$beta[-1, ]!=0), 
           cmcp3=colSums(fit4.cmcp3$fit$beta[-1, ]!=0),
           gel1=colSums(fit4.gel1$fit$beta[-1, ]!=0), 
           gel2=colSums(fit4.gel1$fit$beta[-1, ]!=0), 
           gel3=colSums(fit4.gel1$fit$beta[-1, ]!=0),
           ocmcp1=colSums(fit4.ocmcp1$fit$beta[-1, ]!=0), 
           ocmcp2=colSums(fit4.ocmcp2$fit$beta[-1, ]!=0), 
           ocmcp3=colSums(fit4.ocmcp3$fit$beta[-1, ]!=0),
           ogel1=colSums(fit4.ogel1$fit$beta[-1, ]!=0), 
           ogel2=colSums(fit4.ogel1$fit$beta[-1, ]!=0), 
           ogel3=colSums(fit4.ogel1$fit$beta[-1, ]!=0))
auc4 <- apply(pred4, 2, function(m) {pROC::auc(ytest, m)})
res4 <- rbind(pred4, psel4, auc4)
rownames(res4) <- c(paste0("pred", c(1:length(ytest))), paste0("psel", 1),
                    paste0("auc", 1))
write.table(res4, file="results/metabolomics_alzheimer_res4.csv")