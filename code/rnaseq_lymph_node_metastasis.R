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

################################### model 1 ####################################
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

save(fit1.grridge, fit1.gren1, fit1.gren2, fit1.gren3, 
     fit1.cmcp1, fit1.cmcp2, fit1.cmcp3, fit1.gel1, fit1.gel2,
     fit1.gel3, file="results/rnaseq_lymph_node_metastasis_fit1.Rdata")

### obtaining predictions for test data
ytest <- as.numeric(levels(RespValidation))[RespValidation]
xtest <- ValidationData[, selectgenes]
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
           cmcp1=colSums(fit1.cmcp1$fit$beta[-1, ]!=0), 
           cmcp2=colSums(fit1.cmcp2$fit$beta[-1, ]!=0), 
           cmcp3=colSums(fit1.cmcp3$fit$beta[-1, ]!=0),
           gel1=colSums(fit1.gel1$fit$beta[-1, ]!=0), 
           gel2=colSums(fit1.gel1$fit$beta[-1, ]!=0), 
           gel3=colSums(fit1.gel1$fit$beta[-1, ]!=0))

res1 <- rbind(pred1, psel1)
rownames(res1) <- c(paste0("pred", c(1:length(ytest))), 
                    paste0("psel", 1:nfolds))
write.table(res1, file="results/rnaseq_lymph_node_metastasis_res1.csv")
