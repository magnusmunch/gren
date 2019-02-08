#!/usr/bin/env Rscript

### installation of packages
# if(!("devtools" %in% installed.packages())) {
#   install.packages("devtools")
# }
# library(devtools)
# install_github("magnusmunch/gren/rpackage", local=FALSE,
#                auth_token=Sys.getenv("GITHUB_PAT"))

### libraries
library(gren)
library(GRridge)
library(grpreg)
library(SGL)
library(Biobase)

### load data
load("data/ESetMbolCSFPR2.Rdata")

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

### create co-data
feat <- fData(ESetMbolCSFPR2)

# platformcode gives the type of metabolites
platformcode <- feat$PlatformCode2

# quality gives groups based on a quality measure of the assayd features
part.RSDqc <- CreatePartition(feat$RSDqc, ngroup=3)
quality <- rep(c(1:3), sapply(part.RSDqc, length))[unlist(part.RSDqc)]

# sds gives groups based on the standard deviations of the features
part.sds.apoe <- CreatePartition(apply(metabol.apoe, 1, sd), ngroup=3)
sds.apoe <- rep(c(1:3), sapply(part.sds.apoe, length))[unlist(part.sds.apoe)]

### transform data
alzheim.apoe <- as.numeric(pheno.apoe$D_diag_name) - 1
metabol.apoe.scaled <- scale(metabol.apoe)

################################### model 1 ####################################
### fitting the models
set.seed(2018)
y <- alzheim.apoe
x <- metabol.apoe.scaled
part <- platformcode
n <- nrow(x)
p <- ncol(x)

fit1.gren1 <- gren(x, y, partitions=list(part=part), alpha=0.05)
fit1.gren2 <- gren(x, y, partitions=list(part=part), alpha=0.5)
fit1.gren3 <- gren(x, y, partitions=list(part=part), alpha=0.95)

fit1.grridge <- grridge(t(x), y, list(part=split(1:p, part)))
fit1.sgl1 <- cvSGL(list(x=x, y=y), part, type="logit", alpha=0.05)
fit1.sgl2 <- cvSGL(list(x=x, y=y), part, type="logit", alpha=0.5)
fit1.sgl3 <- cvSGL(list(x=x, y=y), part, type="logit", alpha=0.95)

fit1.cmcp1 <- cv.grpreg(x, y, part, penalty="cMCP", family="binomial", 
                        alpha=0.05)
fit1.cmcp2 <- cv.grpreg(x, y, part, penalty="cMCP", family="binomial", 
                        alpha=0.5)
fit1.cmcp3 <- cv.grpreg(x, y, part, penalty="cMCP", family="binomial", 
                        alpha=0.95)

fit1.gel1 <- cv.grpreg(x, y, part, penalty="gel", family="binomial", alpha=0.05)
fit1.gel2 <- cv.grpreg(x, y, part, penalty="gel", family="binomial", alpha=0.5)
fit1.gel3 <- cv.grpreg(x, y, part, penalty="gel", family="binomial", alpha=0.95)

save(fit1.grridge, fit1.gren1, fit1.gren2, fit1.gren3, fit1.sgl1, fit1.sgl2,
     fit1.sgl3, fit1.cmcp1, fit1.cmcp2, fit1.cmcp3, fit1.gel1, fit1.gel2,
     fit1.gel3, file="results/metabolomics_alzheimer_fit1.Rdata")

### cross-validation of performance
set.seed(2018)
nfolds <- 10
rest <- n %% nfolds
foldid <- sample(rep(1:nfolds, times=round(c(rep(
  n %/% nfolds + as.numeric(rest!=0), times=rest),
  rep(n %/% nfolds, times=nfolds - rest)))))
pred1 <- data.frame(ridge=rep(NA, n), grridge=rep(NA, n),
                    gren1=rep(NA, n), gren2=rep(NA, n), gren3=rep(NA, n),
                    enet1=rep(NA, n), enet2=rep(NA, n), enet3=rep(NA, n),
                    sgl1=rep(NA, n), sgl2=rep(NA, n), sgl3=rep(NA, n),
                    cmcp1=rep(NA, n), cmcp2=rep(NA, n), cmcp3=rep(NA, n),
                    gel1=rep(NA, n), gel2=rep(NA, n), gel3=rep(NA, n))
psel1 <- data.frame(ridge=rep(p, nfolds), grridge=rep(p, nfolds),
                    gren1=rep(NA, nfolds), gren2=rep(NA, nfolds), 
                    gren3=rep(NA, nfolds), enet1=rep(NA, nfolds), 
                    enet2=rep(NA, nfolds), enet3=rep(NA, nfolds),
                    sgl1=rep(NA, nfolds), sgl2=rep(NA, nfolds), 
                    sgl3=rep(NA, nfolds), cmcp1=rep(NA, nfolds), 
                    cmcp2=rep(NA, nfolds), cmcp3=rep(NA, nfolds),
                    gel1=rep(NA, nfolds), gel2=rep(NA, nfolds), 
                    gel3=rep(NA, nfolds))
for(k in 1:nfolds) {
  xtrain <- matrix(x[foldid!=k, ], ncol=p)
  xtest <- matrix(x[foldid==k, ], ncol=p)
  ytrain <- y[foldid!=k]
  ytest <- y[foldid==k]
  
  ### fit methods
  cv1.gren1 <- gren(xtrain, ytrain, partitions=list(part=part), alpha=0.05)
  cv1.gren2 <- gren(xtrain, ytrain, partitions=list(part=part), alpha=0.5)
  cv1.gren3 <- gren(xtrain, ytrain, partitions=list(part=part), alpha=0.95)
  
  cv1.grridge <- grridge(t(xtrain), ytrain, list(part=split(1:p, part)))
  cv1.sgl1 <- cvSGL(list(x=xtrain, y=ytrain), part, type="logit", alpha=0.05)
  cv1.sgl2 <- cvSGL(list(x=xtrain, y=ytrain), part, type="logit", alpha=0.5)
  cv1.sgl3 <- cvSGL(list(x=xtrain, y=ytrain), part, type="logit", alpha=0.95)
  
  cv1.cmcp1 <- cv.grpreg(xtrain, ytrain, part, penalty="cMCP", 
                         family="binomial", alpha=0.05)
  cv1.cmcp2 <- cv.grpreg(xtrain, ytrain, part, penalty="cMCP", 
                         family="binomial", alpha=0.5)
  cv1.cmcp3 <- cv.grpreg(xtrain, ytrain, part, penalty="cMCP", 
                         family="binomial", alpha=0.95)
  
  cv1.gel1 <- cv.grpreg(xtrain, ytrain, part, penalty="gel", family="binomial",
                        alpha=0.05)
  cv1.gel2 <- cv.grpreg(xtrain, ytrain, part, penalty="gel", family="binomial",
                        alpha=0.5)
  cv1.gel3 <- cv.grpreg(xtrain, ytrain, part, penalty="gel", family="binomial",
                        alpha=0.95)
  
  pred1$ridge[foldid==k] <- predict.grridge(cv1.grridge, t(xtest))[, 1]
  pred1$grridge[foldid==k] <- predict.grridge(cv1.grridge, t(xtest))[, 2]
  
  pred1$enet1[foldid==k] <- predict(cv1.gren1, xtest, type="regular")
  pred1$enet2[foldid==k] <- predict(cv1.gren2, xtest, type="regular")
  pred1$enet3[foldid==k] <- predict(cv1.gren3, xtest, type="regular")
  
  pred1$gren1[foldid==k] <- predict(cv1.gren1, xtest, type="groupreg")
  pred1$gren2[foldid==k] <- predict(cv1.gren2, xtest, type="groupreg")
  pred1$gren3[foldid==k] <- predict(cv1.gren3, xtest, type="groupreg")
  
  pred1$sgl1[foldid==k] <- 1/(1 + exp(-xtest %*% cv1.sgl1$fit$beta[, which.min(
    cv1.sgl1$lldiff)]))
  pred1$sgl2[foldid==k] <- 1/(1 + exp(-xtest %*% cv1.sgl2$fit$beta[, which.min(
    cv1.sgl2$lldiff)]))
  pred1$sgl3[foldid==k] <- 1/(1 + exp(-xtest %*% cv1.sgl3$fit$beta[, which.min(
    cv1.sgl3$lldiff)]))
  
  pred1$cmcp1[foldid==k] <- predict(cv1.cmcp1, xtest, type="response")
  pred1$cmcp2[foldid==k] <- predict(cv1.cmcp2, xtest, type="response")
  pred1$cmcp3[foldid==k] <- predict(cv1.cmcp3, xtest, type="response")
  
  pred1$gel1[foldid==k] <- predict(cv1.gel1, xtest, type="response")
  pred1$gel2[foldid==k] <- predict(cv1.gel2, xtest, type="response")
  pred1$gel3[foldid==k] <- predict(cv1.gel3, xtest, type="response")
  
  psel1$enet1[k] <- sum(coef(cv1.gren1, type="regular")!=0) - 1
  psel1$enet2[k] <- sum(coef(cv1.gren2, type="regular")!=0) - 1
  psel1$enet3[k] <- sum(coef(cv1.gren3, type="regular")!=0) - 1
  
  psel1$gren1[k] <- sum(coef(cv1.gren1, type="groupreg")!=0) - 1
  psel1$gren2[k] <- sum(coef(cv1.gren2, type="groupreg")!=0) - 1
  psel1$gren3[k] <- sum(coef(cv1.gren3, type="groupreg")!=0) - 1
  
  psel1$sgl1[k] <- sum(cv1.sgl1$fit$beta[, which.min(cv1.sgl1$lldiff)]!=0)
  psel1$sgl2[k] <- sum(cv1.sgl2$fit$beta[, which.min(cv1.sgl2$lldiff)]!=0)
  psel1$sgl3[k] <- sum(cv1.sgl3$fit$beta[, which.min(cv1.sgl3$lldiff)]!=0)
  
  psel1$cmcp1[k] <- sum(cv1.cmcp1$fit$beta[, cv1.cmcp1$min]!=0) - 1
  psel1$cmcp2[k] <- sum(cv1.cmcp2$fit$beta[, cv1.cmcp2$min]!=0) - 1
  psel1$cmcp3[k] <- sum(cv1.cmcp3$fit$beta[, cv1.cmcp3$min]!=0) - 1
  
  psel1$gel1[k] <- sum(cv1.gel1$fit$beta[, cv1.gel1$min]!=0) - 1
  psel1$gel2[k] <- sum(cv1.gel2$fit$beta[, cv1.gel2$min]!=0) - 1
  psel1$gel3[k] <- sum(cv1.gel3$fit$beta[, cv1.gel3$min]!=0) - 1
  
  res1 <- rbind(pred1, psel1)
  rownames(res1) <- c(paste0("pred", c(1:n)), paste0("psel", 1:nfolds))
  write.table(res1, file="results/metabolomics_alzheimer_res1.csv")
  
}


################################### model 2 ####################################
### fitting the models
set.seed(2018)
y <- alzheim.apoe
x <- metabol.apoe.scaled
part <- quality
n <- nrow(x)
p <- ncol(x)

fit2.gren1 <- gren(x, y, partitions=list(part=part), alpha=0.05)
fit2.gren2 <- gren(x, y, partitions=list(part=part), alpha=0.5)
fit2.gren3 <- gren(x, y, partitions=list(part=part), alpha=0.95)

fit2.grridge <- grridge(t(x), y, list(part=split(1:p, part)))
fit2.sgl1 <- cvSGL(list(x=x, y=y), part, type="logit", alpha=0.05)
fit2.sgl2 <- cvSGL(list(x=x, y=y), part, type="logit", alpha=0.5)
fit2.sgl3 <- cvSGL(list(x=x, y=y), part, type="logit", alpha=0.95)

fit2.cmcp1 <- cv.grpreg(x, y, part, penalty="cMCP", family="binomial", 
                        alpha=0.05)
fit2.cmcp2 <- cv.grpreg(x, y, part, penalty="cMCP", family="binomial", 
                        alpha=0.5)
fit2.cmcp3 <- cv.grpreg(x, y, part, penalty="cMCP", family="binomial", 
                        alpha=0.95)

fit2.gel1 <- cv.grpreg(x, y, part, penalty="gel", family="binomial", alpha=0.05)
fit2.gel2 <- cv.grpreg(x, y, part, penalty="gel", family="binomial", alpha=0.5)
fit2.gel3 <- cv.grpreg(x, y, part, penalty="gel", family="binomial", alpha=0.95)

save(fit2.grridge, fit2.gren1, fit2.gren2, fit2.gren3, fit2.sgl1, fit2.sgl2,
     fit2.sgl3, fit2.cmcp1, fit2.cmcp2, fit2.cmcp3, fit2.gel1, fit2.gel2,
     fit2.gel3, file="results/metabolomics_alzheimer_fit2.Rdata")

### cross-validation of performance
set.seed(2018)
nfolds <- 10
rest <- n %% nfolds
foldid <- sample(rep(1:nfolds, times=round(c(rep(
  n %/% nfolds + as.numeric(rest!=0), times=rest),
  rep(n %/% nfolds, times=nfolds - rest)))))
pred2 <- data.frame(ridge=rep(NA, n), grridge=rep(NA, n),
                    gren1=rep(NA, n), gren2=rep(NA, n), gren3=rep(NA, n),
                    enet1=rep(NA, n), enet2=rep(NA, n), enet3=rep(NA, n),
                    sgl1=rep(NA, n), sgl2=rep(NA, n), sgl3=rep(NA, n),
                    cmcp1=rep(NA, n), cmcp2=rep(NA, n), cmcp3=rep(NA, n),
                    gel1=rep(NA, n), gel2=rep(NA, n), gel3=rep(NA, n))
psel2 <- data.frame(ridge=rep(p, nfolds), grridge=rep(p, nfolds),
                    gren1=rep(NA, nfolds), gren2=rep(NA, nfolds), 
                    gren3=rep(NA, nfolds), enet1=rep(NA, nfolds), 
                    enet2=rep(NA, nfolds), enet3=rep(NA, nfolds),
                    sgl1=rep(NA, nfolds), sgl2=rep(NA, nfolds), 
                    sgl3=rep(NA, nfolds), cmcp1=rep(NA, nfolds), 
                    cmcp2=rep(NA, nfolds), cmcp3=rep(NA, nfolds),
                    gel1=rep(NA, nfolds), gel2=rep(NA, nfolds), 
                    gel3=rep(NA, nfolds))
for(k in 1:nfolds) {
  xtrain <- matrix(x[foldid!=k, ], ncol=p)
  xtest <- matrix(x[foldid==k, ], ncol=p)
  ytrain <- y[foldid!=k]
  ytest <- y[foldid==k]
  
  ### fit methods
  cv2.gren1 <- gren(xtrain, ytrain, partitions=list(part=part), alpha=0.05)
  cv2.gren2 <- gren(xtrain, ytrain, partitions=list(part=part), alpha=0.5)
  cv2.gren3 <- gren(xtrain, ytrain, partitions=list(part=part), alpha=0.95)
  
  cv2.grridge <- grridge(t(xtrain), ytrain, list(part=split(1:p, part)))
  cv2.sgl1 <- cvSGL(list(x=xtrain, y=ytrain), part, type="logit", alpha=0.05)
  cv2.sgl2 <- cvSGL(list(x=xtrain, y=ytrain), part, type="logit", alpha=0.5)
  cv2.sgl3 <- cvSGL(list(x=xtrain, y=ytrain), part, type="logit", alpha=0.95)
  
  cv2.cmcp1 <- cv.grpreg(xtrain, ytrain, part, penalty="cMCP", 
                         family="binomial", alpha=0.05)
  cv2.cmcp2 <- cv.grpreg(xtrain, ytrain, part, penalty="cMCP", 
                         family="binomial", alpha=0.5)
  cv2.cmcp3 <- cv.grpreg(xtrain, ytrain, part, penalty="cMCP", 
                         family="binomial", alpha=0.95)
  
  cv2.gel1 <- cv.grpreg(xtrain, ytrain, part, penalty="gel", family="binomial",
                        alpha=0.05)
  cv2.gel2 <- cv.grpreg(xtrain, ytrain, part, penalty="gel", family="binomial",
                        alpha=0.5)
  cv2.gel3 <- cv.grpreg(xtrain, ytrain, part, penalty="gel", family="binomial",
                        alpha=0.95)
  
  pred2$ridge[foldid==k] <- predict.grridge(cv2.grridge, t(xtest))[, 1]
  pred2$grridge[foldid==k] <- predict.grridge(cv2.grridge, t(xtest))[, 2]
  
  pred2$enet1[foldid==k] <- predict(cv2.gren1, xtest, type="regular")
  pred2$enet2[foldid==k] <- predict(cv2.gren2, xtest, type="regular")
  pred2$enet3[foldid==k] <- predict(cv2.gren3, xtest, type="regular")
  
  pred2$gren1[foldid==k] <- predict(cv2.gren1, xtest, type="groupreg")
  pred2$gren2[foldid==k] <- predict(cv2.gren2, xtest, type="groupreg")
  pred2$gren3[foldid==k] <- predict(cv2.gren3, xtest, type="groupreg")
  
  pred2$sgl1[foldid==k] <- 1/(1 + exp(-xtest %*% cv2.sgl1$fit$beta[, which.min(
    cv2.sgl1$lldiff)]))
  pred2$sgl2[foldid==k] <- 1/(1 + exp(-xtest %*% cv2.sgl2$fit$beta[, which.min(
    cv2.sgl2$lldiff)]))
  pred2$sgl3[foldid==k] <- 1/(1 + exp(-xtest %*% cv2.sgl3$fit$beta[, which.min(
    cv2.sgl3$lldiff)]))
  
  pred2$cmcp1[foldid==k] <- predict(cv2.cmcp1, xtest, type="response")
  pred2$cmcp2[foldid==k] <- predict(cv2.cmcp2, xtest, type="response")
  pred2$cmcp3[foldid==k] <- predict(cv2.cmcp3, xtest, type="response")
  
  pred2$gel1[foldid==k] <- predict(cv2.gel1, xtest, type="response")
  pred2$gel2[foldid==k] <- predict(cv2.gel2, xtest, type="response")
  pred2$gel3[foldid==k] <- predict(cv2.gel3, xtest, type="response")
  
  psel2$enet1[k] <- sum(coef(cv2.gren1, type="regular")!=0) - 1
  psel2$enet2[k] <- sum(coef(cv2.gren2, type="regular")!=0) - 1
  psel2$enet3[k] <- sum(coef(cv2.gren3, type="regular")!=0) - 1
  
  psel2$gren1[k] <- sum(coef(cv2.gren1, type="groupreg")!=0) - 1
  psel2$gren2[k] <- sum(coef(cv2.gren2, type="groupreg")!=0) - 1
  psel2$gren3[k] <- sum(coef(cv2.gren3, type="groupreg")!=0) - 1
  
  psel2$sgl1[k] <- sum(cv2.sgl1$fit$beta[, which.min(cv2.sgl1$lldiff)]!=0)
  psel2$sgl2[k] <- sum(cv2.sgl2$fit$beta[, which.min(cv2.sgl2$lldiff)]!=0)
  psel2$sgl3[k] <- sum(cv2.sgl3$fit$beta[, which.min(cv2.sgl3$lldiff)]!=0)
  
  psel2$cmcp1[k] <- sum(cv2.cmcp1$fit$beta[, cv2.cmcp1$min]!=0) - 1
  psel2$cmcp2[k] <- sum(cv2.cmcp2$fit$beta[, cv2.cmcp2$min]!=0) - 1
  psel2$cmcp3[k] <- sum(cv2.cmcp3$fit$beta[, cv2.cmcp3$min]!=0) - 1
  
  psel2$gel1[k] <- sum(cv2.gel1$fit$beta[, cv2.gel1$min]!=0) - 1
  psel2$gel2[k] <- sum(cv2.gel2$fit$beta[, cv2.gel2$min]!=0) - 1
  psel2$gel3[k] <- sum(cv2.gel3$fit$beta[, cv2.gel3$min]!=0) - 1
  
  res2 <- rbind(pred2, psel2)
  rownames(res2) <- c(paste0("pred", c(1:n)), paste0("psel", 1:nfolds))
  write.table(res2, file="results/metabolomics_alzheimer_res2.csv")
  
}