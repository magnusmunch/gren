#!/usr/bin/env Rscript

### installation of gren if updated on github
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
utest <- unpenal[-id.train, ]
part <- diff.twogroup[apply(micrornas[id.train, ], 2, sd)!=0]
p <- ncol(xtrain)
u <- ncol(utrain)

fit1.gren1 <- gren(xtrain, ytrain, unpenalized=utrain, 
                   partitions=list(part=part), alpha=0.05, standardize=TRUE, 
                   trace=FALSE)
fit1.gren2 <- gren(xtrain, ytrain, unpenalized=utrain, 
                   partitions=list(part=part), alpha=0.5, standardize=TRUE, 
                   trace=FALSE)
fit1.gren3 <- gren(xtrain, ytrain, unpenalized=utrain, 
                   partitions=list(part=part), alpha=0.95, standardize=TRUE, 
                   trace=FALSE)

init1.grridge <- grridge(t(scale(micrornas)), benefit, 
                         list(part=split(1:ncol(micrornas), diff.twogroup)),
                         unpenal= ~ 1 + adjth2 + thscheme2 + thscheme3 + age + 
                           pcrcdiff3, 
                         optl=NULL, dataunpen=as.data.frame(unpenal))

fit1.grridge <- grridge(t(xtrain), ytrain, list(part=split(1:p, part)),
                        unpenal= ~ 1 + adjth2 + thscheme2 + thscheme3 + age + 
                          pcrcdiff3, optl=init1.grridge$optl, 
                        dataunpen=as.data.frame(utrain))

fit1.sgl1 <- cvSGL(list(x=xtrain, y=ytrain), part, type="logit", alpha=0.05)
fit1.sgl1$fit$type <- "logit"
fit1.sgl2 <- cvSGL(list(x=xtrain, y=ytrain), part, type="logit", alpha=0.5)
fit1.sgl2$fit$type <- "logit"
fit1.sgl3 <- cvSGL(list(x=xtrain, y=ytrain), part, type="logit", alpha=0.95)
fit1.sgl3$fit$type <- "logit"

fit1.cmcp1 <- cv.grpreg(cbind(utrain, xtrain), ytrain, 
                        c(rep(0, ncol(utrain)), part), penalty="cMCP", 
                        family="binomial", alpha=0.05)
fit1.cmcp2 <- cv.grpreg(cbind(utrain, xtrain), ytrain, 
                        c(rep(0, ncol(utrain)), part), penalty="cMCP", 
                        family="binomial", alpha=0.5)
fit1.cmcp3 <- cv.grpreg(cbind(utrain, xtrain), ytrain, 
                        c(rep(0, ncol(utrain)), part), penalty="cMCP", 
                        family="binomial", alpha=0.95)

fit1.gel1 <- cv.grpreg(cbind(utrain, xtrain), ytrain, 
                       c(rep(0, ncol(utrain)), part), penalty="gel", 
                       family="binomial", alpha=0.05)
# gel2 and gel3 give just the saturated model, so we create cv.grpreg ourselves
fit1.gel2 <- grpreg(cbind(utrain, xtrain), ytrain, 
                    c(rep(0, ncol(utrain)), part), penalty="gel", 
                    family="binomial", alpha=0.5)
fit1.gel2 <- list(cve=NA, cvse=NA, lambda=fit1.gel2$lambda, fit=fit1.gel2, 
                   min=1, lambda.min=fit1.gel2$lambda, null.dev=NA, pe=NA)
class(fit1.gel2) <- "cv.grpreg"
fit1.gel3 <- grpreg(cbind(utrain, xtrain), ytrain, 
                    c(rep(0, ncol(utrain)), part), penalty="gel", 
                    family="binomial", alpha=0.95)
fit1.gel3 <- list(cve=NA, cvse=NA, lambda=fit1.gel3$lambda, fit=fit1.gel3, 
                  min=1, lambda.min=fit1.gel3$lambda, null.dev=NA, pe=NA)
class(fit1.gel3) <- "cv.grpreg"

save(fit1.grridge, fit1.gren1, fit1.gren2, fit1.gren3, fit1.sgl1, fit1.sgl2,
     fit1.sgl3, fit1.cmcp1, fit1.cmcp2, fit1.cmcp3, fit1.gel1, fit1.gel2,
     fit1.gel3, file="results/micrornaseq_colorectal_cancer_fit1.Rdata")

### prediction on test set
pred1 <- data.frame(ridge=predict.grridge(fit1.grridge, t(xtest), 
                                          dataunpennew=as.data.frame(utest))[
                                            , 1], 
                    grridge=predict.grridge(fit1.grridge, t(xtest),
                                            dataunpennew=as.data.frame(utest))[
                                              , 2],
                    gren1=predict(fit1.gren1, xtest, utest, type="groupreg"),
                    gren2=predict(fit1.gren2, xtest, utest, type="groupreg"), 
                    gren3=predict(fit1.gren3, xtest, utest, type="groupreg"), 
                    enet1=predict(fit1.gren1, xtest, utest, type="regular"),
                    enet2=predict(fit1.gren2, xtest, utest, type="regular"), 
                    enet3=predict(fit1.gren3, xtest, utest, type="regular"), 
                    sgl1=predictSGL(fit1.sgl1$fit, xtest),
                    sgl2=predictSGL(fit1.sgl2$fit, xtest),
                    sgl3=predictSGL(fit1.sgl3$fit, xtest),
                    cmcp1=predict(fit1.cmcp1$fit, cbind(utest, xtest), 
                                  type="response"), 
                    cmcp2=predict(fit1.cmcp2$fit, cbind(utest, xtest), 
                                  type="response"), 
                    cmcp3=predict(fit1.cmcp3$fit, cbind(utest, xtest), 
                                  type="response"),
                    gel1=predict(fit1.gel1$fit, cbind(utest, xtest), 
                                 type="response"), 
                    gel2=predict(fit1.gel2$fit, cbind(utest, xtest), 
                                 type="response"), 
                    gel3=predict(fit1.gel3$fit, cbind(utest, xtest), 
                                 type="response"))
psel1 <- c(ridge=p, grridge=p,
           gren1=fit1.gren1$freq.model$groupreg$df - u,
           gren2=fit1.gren2$freq.model$groupreg$df - u,
           gren3=fit1.gren3$freq.model$groupreg$df - u,
           enet1=fit1.gren1$freq.model$regular$df - u,
           enet2=fit1.gren2$freq.model$regular$df - u,
           enet3=fit1.gren3$freq.model$regular$df - u,
           sgl1=colSums(fit1.sgl1$fit$beta!=0),
           sgl2=colSums(fit1.sgl2$fit$beta!=0),
           sgl3=colSums(fit1.sgl3$fit$beta!=0),
           cmcp1=colSums(fit1.cmcp1$fit$beta[-c(1:(u + 1)), ]!=0), 
           cmcp2=colSums(fit1.cmcp2$fit$beta[-c(1:(u + 1)), ]!=0), 
           cmcp3=colSums(fit1.cmcp3$fit$beta[-c(1:(u + 1)), ]!=0),
           gel1=colSums(fit1.gel1$fit$beta[-c(1:(u + 1)), ]!=0), 
           gel2=colSums(as.matrix(fit1.gel2$fit$beta[-c(1:(u + 1)), ]!=0)), 
           gel3=colSums(as.matrix(fit1.gel3$fit$beta[-c(1:(u + 1)), ]!=0)))
auc1 <- apply(pred1, 2, function(m) {pROC::auc(ytest, m)})
res1 <- rbind(pred1, psel1, auc1)
rownames(res1) <- c(paste0("pred", c(1:length(ytest))), paste0("psel", 1),
                    paste0("auc", 1))
write.table(res1, file="results/micrornaseq_colorectal_cancer_res1.csv")


################################### model 2 ####################################
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
utest <- unpenal[-id.train, ]
part <- diff.threegroup[apply(micrornas[id.train, ], 2, sd)!=0]
p <- ncol(xtrain)

fit2.gren1 <- gren(xtrain, ytrain, unpenalized=utrain, 
                   partitions=list(part=part), alpha=0.05, standardize=TRUE, 
                   trace=FALSE)
fit2.gren2 <- gren(xtrain, ytrain, unpenalized=utrain, 
                   partitions=list(part=part), alpha=0.5, standardize=TRUE, 
                   trace=FALSE)
fit2.gren3 <- gren(xtrain, ytrain, unpenalized=utrain, 
                   partitions=list(part=part), alpha=0.95, standardize=TRUE, 
                   trace=FALSE)

init2.grridge <- grridge(t(scale(micrornas)), benefit, 
                         list(part=split(1:ncol(micrornas), diff.threegroup)),
                         unpenal= ~ 1 + adjth2 + thscheme2 + thscheme3 + age + 
                           pcrcdiff3, optl=NULL, 
                         dataunpen=as.data.frame(unpenal))

fit2.grridge <- grridge(t(xtrain), ytrain, list(part=split(1:p, part)),
                        unpenal= ~ 1 + adjth2 + thscheme2 + thscheme3 + age + 
                          pcrcdiff3, optl=init2.grridge$optl, 
                        dataunpen=as.data.frame(utrain))

fit2.sgl1 <- cvSGL(list(x=xtrain, y=ytrain), part, type="logit", alpha=0.05)
fit2.sgl1$fit$type <- "logit"
fit2.sgl2 <- cvSGL(list(x=xtrain, y=ytrain), part, type="logit", alpha=0.5)
fit2.sgl2$fit$type <- "logit"
fit2.sgl3 <- cvSGL(list(x=xtrain, y=ytrain), part, type="logit", alpha=0.95)
fit2.sgl3$fit$type <- "logit"

fit2.cmcp1 <- cv.grpreg(cbind(utrain, xtrain), ytrain, 
                        c(rep(0, ncol(utrain)), part), penalty="cMCP", 
                        family="binomial", alpha=0.05)
fit2.cmcp2 <- cv.grpreg(cbind(utrain, xtrain), ytrain, 
                        c(rep(0, ncol(utrain)), part), penalty="cMCP", 
                        family="binomial", alpha=0.5)
fit2.cmcp3 <- cv.grpreg(cbind(utrain, xtrain), ytrain, 
                        c(rep(0, ncol(utrain)), part), penalty="cMCP", 
                        family="binomial", alpha=0.95)

fit2.gel1 <- cv.grpreg(cbind(utrain, xtrain), ytrain, 
                       c(rep(0, ncol(utrain)), part), penalty="gel", 
                       family="binomial", alpha=0.05)
# gel2 and gel3 give just the saturated model, so we create cv.grpreg ourselves
fit2.gel2 <- grpreg(cbind(utrain, xtrain), ytrain, 
                    c(rep(0, ncol(utrain)), part), penalty="gel", 
                    family="binomial", alpha=0.5)
fit2.gel2 <- list(cve=NA, cvse=NA, lambda=fit2.gel2$lambda, fit=fit2.gel2, 
                  min=1, lambda.min=fit2.gel2$lambda, null.dev=NA, pe=NA)
class(fit2.gel2) <- "cv.grpreg"
fit2.gel3 <- grpreg(cbind(utrain, xtrain), ytrain, 
                    c(rep(0, ncol(utrain)), part), penalty="gel", 
                    family="binomial", alpha=0.95)
fit2.gel3 <- list(cve=NA, cvse=NA, lambda=fit2.gel3$lambda, fit=fit2.gel3, 
                  min=1, lambda.min=fit2.gel3$lambda, null.dev=NA, pe=NA)
class(fit2.gel3) <- "cv.grpreg"

save(fit2.grridge, fit2.gren1, fit2.gren2, fit2.gren3, fit2.sgl1, fit2.sgl2,
     fit2.sgl3, fit2.cmcp1, fit2.cmcp2, fit2.cmcp3, fit2.gel1, fit2.gel2,
     fit2.gel3, file="results/micrornaseq_colorectal_cancer_fit2.Rdata")

### prediction on test set
pred2 <- data.frame(ridge=predict.grridge(fit2.grridge, t(xtest), 
                                          dataunpennew=as.data.frame(utest))[
                                            , 1], 
                    grridge=predict.grridge(fit2.grridge, t(xtest),
                                            dataunpennew=as.data.frame(utest))[
                                              , 2],
                    gren1=predict(fit2.gren1, xtest, utest, type="groupreg"),
                    gren2=predict(fit2.gren2, xtest, utest, type="groupreg"), 
                    gren3=predict(fit2.gren3, xtest, utest, type="groupreg"), 
                    enet1=predict(fit2.gren1, xtest, utest, type="regular"),
                    enet2=predict(fit2.gren2, xtest, utest, type="regular"), 
                    enet3=predict(fit2.gren3, xtest, utest, type="regular"), 
                    sgl1=predictSGL(fit2.sgl1$fit, xtest),
                    sgl2=predictSGL(fit2.sgl2$fit, xtest),
                    sgl3=predictSGL(fit2.sgl3$fit, xtest),
                    cmcp1=predict(fit2.cmcp1$fit, cbind(utest, xtest), 
                                  type="response"), 
                    cmcp2=predict(fit2.cmcp2$fit, cbind(utest, xtest), 
                                  type="response"), 
                    cmcp3=predict(fit2.cmcp3$fit, cbind(utest, xtest), 
                                  type="response"),
                    gel1=predict(fit2.gel1$fit, cbind(utest, xtest), 
                                 type="response"), 
                    gel2=predict(fit2.gel2$fit, cbind(utest, xtest), 
                                 type="response"), 
                    gel3=predict(fit2.gel3$fit, cbind(utest, xtest), 
                                 type="response"))
psel2 <- c(ridge=p, grridge=p,
           gren1=fit2.gren1$freq.model$groupreg$df - u,
           gren2=fit2.gren2$freq.model$groupreg$df - u,
           gren3=fit2.gren3$freq.model$groupreg$df - u,
           enet1=fit2.gren1$freq.model$regular$df - u,
           enet2=fit2.gren2$freq.model$regular$df - u,
           enet3=fit2.gren3$freq.model$regular$df - u,
           sgl1=colSums(fit2.sgl1$fit$beta!=0),
           sgl2=colSums(fit2.sgl2$fit$beta!=0),
           sgl3=colSums(fit2.sgl3$fit$beta!=0),
           cmcp1=colSums(fit2.cmcp1$fit$beta[-c(1:(u + 1)), ]!=0), 
           cmcp2=colSums(fit2.cmcp2$fit$beta[-c(1:(u + 1)), ]!=0), 
           cmcp3=colSums(fit2.cmcp3$fit$beta[-c(1:(u + 1)), ]!=0),
           gel1=colSums(fit2.gel1$fit$beta[-c(1:(u + 1)), ]!=0), 
           gel2=colSums(as.matrix(fit2.gel2$fit$beta[-c(1:(u + 1)), ]!=0)), 
           gel3=colSums(as.matrix(fit2.gel3$fit$beta[-c(1:(u + 1)), ]!=0)))
auc2 <- apply(pred2, 2, function(m) {pROC::auc(ytest, m)})
res2 <- rbind(pred2, psel2, auc2)
rownames(res2) <- c(paste0("pred", c(1:length(ytest))), paste0("psel", 1),
                    paste0("auc", 1))
write.table(res2, file="results/micrornaseq_colorectal_cancer_res2.csv")


################################### model 3 ####################################
###  create 100 random splits of features into groups
set.seed(2019)
nsplits <- 100
y <- benefit
x <- scale(micrornas)
unpenal <- unpenal
p <- ncol(x)
u <- ncol(unpenal)
part <- diff.threegroup

multipliers <- data.frame(grridge=matrix(NA, nrow=nsplits, 
                                         ncol=length(unique(part))),
                          gren1=matrix(NA, nrow=nsplits, 
                                       ncol=length(unique(part))),
                          gren2=matrix(NA, nrow=nsplits, 
                                       ncol=length(unique(part))),
                          gren3=matrix(NA, nrow=nsplits, 
                                       ncol=length(unique(part))))

for(k in 1:nsplits) {
  cat(paste("split ", k, "\n"))
  set.seed(2019 + k)
  
  part.train <- sample(part)

  fit3.gren1 <- gren(x, y, unpenalized=unpenal, 
                     partitions=list(part=part.train), alpha=0.05, 
                     standardize=TRUE, trace=FALSE)
  fit3.gren2 <- gren(x, y, unpenalized=unpenal, 
                     partitions=list(part=part.train), alpha=0.5, 
                     standardize=TRUE, trace=FALSE)
  fit3.gren3 <- gren(x, y, unpenalized=unpenal, 
                     partitions=list(part=part.train), alpha=0.95, 
                     standardize=TRUE, trace=FALSE)
  
  fit3.grridge <- grridge(t(x), y, list(part=split(1:p, part.train)),
                          unpenal= ~ 1 + adjth2 + thscheme2 + thscheme3 + age + 
                            pcrcdiff3, dataunpen=as.data.frame(unpenal))

  multipliers[k, c(1:length(unique(part)))] <- fit3.grridge$lambdamults$part
  multipliers[k, c((length(unique(part)) + 1):(2*length(unique(part))))] <- 
    fit3.gren1$lambdag$part
  multipliers[k, c((2*length(unique(part)) + 1):(3*length(unique(part))))] <- 
    fit3.gren2$lambdag$part
  multipliers[k, c((3*length(unique(part)) + 1):(4*length(unique(part))))] <- 
    fit3.gren3$lambdag$part

}
res3 <- multipliers
rownames(res3) <- paste0("split", c(1:nsplits))
write.table(res3, file="results/micrornaseq_colorectal_cancer_res3.csv")


################################### model 4 ####################################
###  create 100 random splits of features into groups
set.seed(2019)
nsplits <- 100
y <- benefit
x <- scale(micrornas)
unpenal <- unpenal
p <- ncol(x)
u <- ncol(unpenal)
ngroups <- 10
part <- rep(1:ngroups, times=round(c(rep(
  p %/% ngroup + as.numeric((p %% ngroups)!=0), times=p %% ngroups),
  rep(p %/% ngroups, times=ngroups - p %% ngroups))))

multipliers <- data.frame(grridge=matrix(NA, nrow=nsplits, 
                                         ncol=length(unique(part))),
                          gren1=matrix(NA, nrow=nsplits, 
                                       ncol=length(unique(part))),
                          gren2=matrix(NA, nrow=nsplits, 
                                       ncol=length(unique(part))),
                          gren3=matrix(NA, nrow=nsplits, 
                                       ncol=length(unique(part))))

for(k in 1:nsplits) {
  cat(paste("split ", k, "\n"))
  set.seed(2019 + k)
  
  part.train <- sample(part)
  
  fit3.gren1 <- gren(x, y, unpenalized=unpenal, 
                     partitions=list(part=part.train), alpha=0.05, 
                     standardize=TRUE, trace=FALSE)
  fit3.gren2 <- gren(x, y, unpenalized=unpenal, 
                     partitions=list(part=part.train), alpha=0.5, 
                     standardize=TRUE, trace=FALSE)
  fit3.gren3 <- gren(x, y, unpenalized=unpenal, 
                     partitions=list(part=part.train), alpha=0.95, 
                     standardize=TRUE, trace=FALSE)
  
  fit3.grridge <- grridge(t(x), y, list(part=split(1:p, part.train)),
                          unpenal= ~ 1 + adjth2 + thscheme2 + thscheme3 + age + 
                            pcrcdiff3, dataunpen=as.data.frame(unpenal))
  
  multipliers[k, c(1:length(unique(part)))] <- fit3.grridge$lambdamults$part
  multipliers[k, c((length(unique(part)) + 1):(2*length(unique(part))))] <- 
    fit3.gren1$lambdag$part
  multipliers[k, c((2*length(unique(part)) + 1):(3*length(unique(part))))] <- 
    fit3.gren2$lambdag$part
  multipliers[k, c((3*length(unique(part)) + 1):(4*length(unique(part))))] <- 
    fit3.gren3$lambdag$part
  
}
res3 <- multipliers
rownames(res3) <- paste0("split", c(1:nsplits))
write.table(res3, file="results/micrornaseq_colorectal_cancer_res3.csv")

################################### model 3 ####################################
###  create 100 random splits of features into 3 groups, same size as FDR groups
set.seed(2019)
nsplits <- 100
y <- benefit
x <- scale(micrornas)
unpenal <- unpenal
p <- ncol(x)
u <- ncol(unpenal)
part <- diff.threegroup

multipliers <- data.frame(grridge=matrix(NA, nrow=nsplits, 
                                         ncol=length(unique(part))),
                          gren1=matrix(NA, nrow=nsplits, 
                                       ncol=length(unique(part))),
                          gren2=matrix(NA, nrow=nsplits, 
                                       ncol=length(unique(part))),
                          gren3=matrix(NA, nrow=nsplits, 
                                       ncol=length(unique(part))))

for(k in 1:nsplits) {
  cat(paste("split ", k, "\n"))
  set.seed(2019 + k)
  
  part.train <- sample(part)
  
  fit3.gren1 <- gren(x, y, unpenalized=unpenal, 
                     partitions=list(part=part.train), alpha=0.05, 
                     standardize=TRUE, trace=FALSE)
  fit3.gren2 <- gren(x, y, unpenalized=unpenal, 
                     partitions=list(part=part.train), alpha=0.5, 
                     standardize=TRUE, trace=FALSE)
  fit3.gren3 <- gren(x, y, unpenalized=unpenal, 
                     partitions=list(part=part.train), alpha=0.95, 
                     standardize=TRUE, trace=FALSE)
  
  fit3.grridge <- grridge(t(x), y, list(part=split(1:p, part.train)),
                          unpenal= ~ 1 + adjth2 + thscheme2 + thscheme3 + age + 
                            pcrcdiff3, dataunpen=as.data.frame(unpenal))
  
  multipliers[k, c(1:length(unique(part)))] <- fit3.grridge$lambdamults$part
  multipliers[k, c((length(unique(part)) + 1):(2*length(unique(part))))] <- 
    fit3.gren1$lambdag$part
  multipliers[k, c((2*length(unique(part)) + 1):(3*length(unique(part))))] <- 
    fit3.gren2$lambdag$part
  multipliers[k, c((3*length(unique(part)) + 1):(4*length(unique(part))))] <- 
    fit3.gren3$lambdag$part
  
}
res3 <- multipliers
rownames(res3) <- paste0("split", c(1:nsplits))
write.table(res3, file="results/micrornaseq_colorectal_cancer_res3.csv")


################################### model 4 ####################################
###  create 100 random splits of features into 10 evenly sized groups
set.seed(2019)
nsplits <- 100
y <- benefit
x <- scale(micrornas)
unpenal <- unpenal
p <- ncol(x)
u <- ncol(unpenal)
ngroups <- 10
part <- rep(1:ngroups, times=round(c(rep(
  p %/% ngroups + as.numeric((p %% ngroups)!=0), times=p %% ngroups),
  rep(p %/% ngroups, times=ngroups - p %% ngroups))))

multipliers <- data.frame(grridge=matrix(NA, nrow=nsplits, 
                                         ncol=length(unique(part))),
                          gren1=matrix(NA, nrow=nsplits, 
                                       ncol=length(unique(part))),
                          gren2=matrix(NA, nrow=nsplits, 
                                       ncol=length(unique(part))),
                          gren3=matrix(NA, nrow=nsplits, 
                                       ncol=length(unique(part))))

for(k in 1:nsplits) {
  cat(paste("split ", k, "\n"))
  set.seed(2019 + k)
  
  part.train <- sample(part)
  
  fit4.gren1 <- gren(x, y, unpenalized=unpenal, 
                     partitions=list(part=part.train), alpha=0.05, 
                     standardize=TRUE, trace=FALSE)
  fit4.gren2 <- gren(x, y, unpenalized=unpenal, 
                     partitions=list(part=part.train), alpha=0.5, 
                     standardize=TRUE, trace=FALSE)
  fit4.gren3 <- gren(x, y, unpenalized=unpenal, 
                     partitions=list(part=part.train), alpha=0.95, 
                     standardize=TRUE, trace=FALSE)
  
  fit4.grridge <- grridge(t(x), y, list(part=split(1:p, part.train)),
                          unpenal= ~ 1 + adjth2 + thscheme2 + thscheme3 + age + 
                            pcrcdiff3, dataunpen=as.data.frame(unpenal))
  
  multipliers[k, c(1:length(unique(part)))] <- fit4.grridge$lambdamults$part
  multipliers[k, c((length(unique(part)) + 1):(2*length(unique(part))))] <- 
    fit4.gren1$lambdag$part
  multipliers[k, c((2*length(unique(part)) + 1):(3*length(unique(part))))] <- 
    fit4.gren2$lambdag$part
  multipliers[k, c((3*length(unique(part)) + 1):(4*length(unique(part))))] <- 
    fit4.gren3$lambdag$part
  
}
res4 <- multipliers
rownames(res4) <- paste0("split", c(1:nsplits))
write.table(res4, file="results/micrornaseq_colorectal_cancer_res4.csv")

################################### model 5 ####################################
###  create 100 random splits of features into 10 evenly sized groups
set.seed(2019)
nsamples <- 100
psel <- 25
y <- benefit
x <- scale(micrornas)
unpenal <- unpenal
p <- ncol(x)
u <- ncol(unpenal)
ngroups <- 10
part <- diff.threegroup
p <- ncol(x)

for(k in 1:nsamples) {
  cat(paste("sample ", k, "\n"))
  
  set.seed(2019 + k)
  
  id.boot <- c(sample(which(y==0), replace=TRUE), 
               sample(which(y==1), replace=TRUE))
  # remove the constant micrornas
  is.const <- apply(x[id.boot, ], 2, sd)==0
  xboot <- scale(x[id.boot, !is.const])
  yboot <- y[id.boot]
  uboot <- unpenal[id.boot, ]
  part.boot <- part[!is.const]
  
  fit5.gren1 <- gren(xboot, yboot, unpenalized=uboot, 
                     partitions=list(part=part.boot), alpha=0.05, 
                     standardize=TRUE, trace=FALSE)
  fit5.gren2 <- gren(xboot, yboot, unpenalized=uboot, 
                     partitions=list(part=part.boot), alpha=0.5, 
                     standardize=TRUE, trace=FALSE, psel=psel)
  fit5.gren3 <- gren(xboot, yboot, unpenalized=uboot, 
                     partitions=list(part=part.boot), alpha=0.95, 
                     standardize=TRUE, trace=FALSE, psel=psel)
  
  # sel.grridge[[k]] <- names(boot.grridge$resEN$whichEN)
  sel.gren1[[k]] <- rownames(coef(boot.gren1$freq.model$groupreg))[
    as.numeric(coef(boot.gren1$freq.model$groupreg))!=0][-1]
  sel.gren2[[k]] <- rownames(coef(boot.gren2$freq.model$groupreg))[
    as.numeric(coef(boot.gren2$freq.model$groupreg))!=0][-1]
  sel.gren3[[k]] <- rownames(coef(boot.gren3$freq.model$groupreg))[
    as.numeric(coef(boot.gren3$freq.model$groupreg))!=0][-1]
  sel.enet1[[k]] <- rownames(coef(boot.gren1$freq.model$regular))[
    as.numeric(coef(boot.gren1$freq.model$regular))!=0][-1]
  sel.enet2[[k]] <- rownames(coef(boot.gren2$freq.model$regular))[
    as.numeric(coef(boot.gren2$freq.model$regular))!=0][-1]
  sel.enet3[[k]] <- rownames(coef(boot.gren3$freq.model$regular))[
    as.numeric(coef(boot.gren3$freq.model$regular))!=0][-1]
  
}




res4 <- multipliers
rownames(res4) <- paste0("split", c(1:nsplits))
write.table(res4, file="results/micrornaseq_colorectal_cancer_res4.csv")
