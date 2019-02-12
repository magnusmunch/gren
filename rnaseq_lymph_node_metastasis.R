#!/usr/bin/env Rscript

### installation of packages
if(!("CoRF" %in% installed.packages())) {
  if(!("devtools" %in% installed.packages())) {
    install.packages("devtools")
  }
  library(devtools)
  install_github("DennisBeest/CoRF", local=FALSE,
                 auth_token=Sys.getenv("GITHUB_PAT"))
}

### libraries
library(CoRF)
library(GRridge)
library(gren)
library(grpreg)
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

### fitting the models
nfolds <- 1
set.seed(2018)
y <- RespTrain
x <- Tr
part <- partLNM
n <- nrow(x)
p <- ncol(x)

fit1.gren1 <- gren(x, y, partitions=part, alpha=0.05, trace=FALSE)
fit1.gren2 <- gren(x, y, partitions=part, alpha=0.5, trace=FALSE)
fit1.gren3 <- gren(x, y, partitions=part, alpha=0.95, trace=FALSE)

fit1.grridge <- grridge(t(x), y, list(corr=split(1:p, part$corr), 
                                      pv=split(1:p, part$pv)))

fit1.cmcp1 <- cv.grpreg(x, y, as.numeric(as.factor(paste(corpart0, pvpart0))), 
                        penalty="cMCP", family="binomial", alpha=0.05)
fit1.cmcp2 <- cv.grpreg(x, y, as.numeric(as.factor(paste(corpart0, pvpart0))), 
                        penalty="cMCP", family="binomial", alpha=0.5)
fit1.cmcp3 <- cv.grpreg(x, y, as.numeric(as.factor(paste(corpart0, pvpart0))), 
                        penalty="cMCP", family="binomial", alpha=0.95)

fit1.gel1 <- cv.grpreg(x, y, as.numeric(as.factor(paste(corpart0, pvpart0))), 
                       penalty="gel", family="binomial", alpha=0.05)
fit1.gel2 <- cv.grpreg(x, y, as.numeric(as.factor(paste(corpart0, pvpart0))), 
                       penalty="gel", family="binomial", alpha=0.5)
fit1.gel3 <- cv.grpreg(x, y, as.numeric(as.factor(paste(corpart0, pvpart0))), 
                       penalty="gel", family="binomial", alpha=0.95)

### obtaining predictions for test data
ytest <- as.numeric(levels(RespValidation))[RespValidation]
xtest <- ValidationData[, selectgenes]
pred1 <- data.frame(ridge=predict.grridge(fit1.grridge, t(xtest))[, 1], 
                    grridge=predict.grridge(fit1.grridge, t(xtest))[, 2],
                    gren1=as.numeric(predict(fit1.gren1, xtest, 
                                             type="groupreg")),
                    gren2=as.numeric(predict(fit1.gren2, xtest, 
                                             type="groupreg")), 
                    gren3=as.numeric(predict(fit1.gren3, xtest, 
                                             type="groupreg")), 
                    enet1=as.numeric(predict(fit1.gren1, xtest, 
                                             type="regular")),
                    enet2=as.numeric(predict(fit1.gren2, xtest, 
                                             type="regular")), 
                    enet3=as.numeric(predict(fit1.gren3, xtest, 
                                             type="regular")), 
                    cmcp1=predict(fit1.cmcp1, xtest, type="response"), 
                    cmcp2=predict(fit1.cmcp2, xtest, type="response"), 
                    cmcp3=predict(fit1.cmcp2, xtest, type="response"),
                    gel1=predict(fit1.gel1, xtest, type="response"), 
                    gel2=predict(fit1.gel2, xtest, type="response"), 
                    gel3=predict(fit1.gel3, xtest, type="response"))

psel1 <- data.frame(ridge=p, grridge=p,
                    gren1=sum(coef(fit1.gren1$freq.model$groupreg, 
                                   s=fit1.gren1$lambda)[-1, ]!=0),
                    gren2=sum(coef(fit1.gren2$freq.model$groupreg, 
                                   s=fit1.gren2$lambda)[-1, ]!=0),
                    gren3=sum(coef(fit1.gren3$freq.model$groupreg, 
                                   s=fit1.gren3$lambda)[-1, ]!=0),
                    enet1=sum(coef(fit1.gren1$freq.model$regular, 
                                   s=fit1.gren1$lambda)[-1, ]!=0),
                    enet2=sum(coef(fit1.gren2$freq.model$regular, 
                                   s=fit1.gren2$lambda)[-1, ]!=0),
                    enet3=sum(coef(fit1.gren3$freq.model$regular, 
                                   s=fit1.gren3$lambda)[-1, ]!=0),
                    cmcp1=sum(fit1.cmcp1$fit$beta[-1, fit1.cmcp1$min]!=0), 
                    cmcp2=sum(fit1.cmcp2$fit$beta[-1, fit1.cmcp2$min]!=0), 
                    cmcp3=sum(fit1.cmcp3$fit$beta[-1, fit1.cmcp3$min]!=0),
                    gel1=sum(fit1.gel1$fit$beta[-1, fit1.gel1$min]!=0), 
                    gel2=sum(fit1.gel1$fit$beta[-1, fit1.gel1$min]!=0), 
                    gel3=sum(fit1.gel1$fit$beta[-1, fit1.gel1$min]!=0))

res1 <- rbind(pred1, psel1)
rownames(res1) <- c(paste0("pred", c(1:length(ytest))), 
                    paste0("psel", 1:nfolds))
write.table(res1, file="results/rnaseq_lymph_node_metastasis_res1.csv")


