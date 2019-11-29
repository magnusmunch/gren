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
library(SGL)
<<<<<<< HEAD
=======
library(randomForestSRC)
library(harmonicmeanp)
library(freeknotsplines)
>>>>>>> 31db1da7fc892df40e05ce4287d82f38e942befc
library(foreach)
library(doParallel)
library(microbenchmark)

### parallelisation
parallel <- TRUE

### load data
load("data/forMagnusN88.Rdata")
<<<<<<< HEAD
=======
codata <- read.table("data/results_all.txt", header=TRUE)
>>>>>>> 31db1da7fc892df40e05ce4287d82f38e942befc

### create model matrix for unpenalized covariates
unpenal <- model.matrix(~ adjth + thscheme + age + pcrcdiff, data=datfr)[, -1]

<<<<<<< HEAD
### create partitioning based on FDR <= 0.05 and FDR <= 0.001
miRNA.BFDR <- as.character(TumMirs$miRNA)[TumMirs$BFDR_PNminP < 0.001]
miRNA.TumMirs <- as.character(TumMirs$miRNA)
miRNA <- as.character(sapply(rownames(mirnormcen_resp), function(s) {
  strsplit(s, split=" ")[[1]][1]}))
diff.expr <- miRNA %in% miRNA.BFDR + miRNA %in% miRNA.TumMirs + 1

=======
>>>>>>> 31db1da7fc892df40e05ce4287d82f38e942befc
### target vector
benefit <- as.numeric(resp) - 1

### mirna data
micrornas <- t(mirnormcen_resp)
<<<<<<< HEAD
colnames(micrornas) <- miRNA
=======
colnames(micrornas) <- sapply(strsplit(colnames(micrornas), " "), "[[", 1)


### create partitioning based on FDR <= 0.05 and FDR <= 0.001
# miRNA.BFDR <- as.character(TumMirs$miRNA)[TumMirs$BFDR_PNminP < 0.001]
# miRNA.TumMirs <- as.character(TumMirs$miRNA)
# miRNA <- as.character(sapply(rownames(mirnormcen_resp), function(s) {
#   strsplit(s, split=" ")[[1]][1]}))
# diff.expr <- miRNA %in% miRNA.BFDR + miRNA %in% miRNA.TumMirs + 1
# 
# sum(codata$BFDR_MNminM <= 0.001 & codata$BFDR_PNminP <= 0.001 & 
#       codata$miRNA %in% colnames(micrornas))
# 
# sum(codata$BFDR_MNminM > 0.001 & codata$BFDR_MNminM <= 0.05 &
#       codata$BFDR_PNminP > 0.001 & codata$BFDR_PNminP <= 0.05 &
#       codata$miRNA %in% colnames(micrornas))

hmfdr <- apply(cbind(codata$BFDR_MNminM, codata$BFDR_PNminP), 1, p.hmp)
diff.expr <- rep(4, ncol(micrornas))
diff.expr[colnames(micrornas) %in% codata$miRNA[hmfdr <= 0.001]] <- 1
diff.expr[colnames(micrornas) %in% 
            codata$miRNA[hmfdr > 0.001 & hmfdr <= 0.05]] <- 2
diff.expr[colnames(micrornas) %in% codata$miRNA[hmfdr > 0.05]] <- 3

kn <- knots(ecdf(hmfdr))
spl <- freelsgen(kn, 1:length(kn), degree=1, numknot=2, seed=2019, stream=0)
diff.expr2 <- rep(4, ncol(micrornas))
diff.expr2[colnames(micrornas) %in% codata$miRNA[
  hmfdr <= spl@optknot[1]]] <- 1
diff.expr2[colnames(micrornas) %in% codata$miRNA[
  hmfdr > spl@optknot[1] & hmfdr <= spl@optknot[2]]] <- 2
diff.expr2[colnames(micrornas) %in% codata$miRNA[hmfdr > spl@optknot[2]]] <- 3
>>>>>>> 31db1da7fc892df40e05ce4287d82f38e942befc

### randomly split in test and train data (30% and 70%, respectively)
set.seed(2019)
id.train <- c(sample(which(benefit==0), size=floor(sum(benefit==0)*0.7)),
              sample(which(benefit==1), size=floor(sum(benefit==1)*0.7)))

################################### model 1 ####################################
### preparing data
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
part <- diff.expr[apply(micrornas[id.train, ], 2, sd)!=0]
<<<<<<< HEAD
=======
part2 <- diff.expr2[apply(micrornas[id.train, ], 2, sd)!=0]
>>>>>>> 31db1da7fc892df40e05ce4287d82f38e942befc
p <- ncol(xtrain)
u <- ncol(utrain)

### running and benchmarking model
bench <- microbenchmark(
fit.gren1 <- gren(xtrain, ytrain, unpenalized=utrain,
                  partitions=list(part=part), alpha=0.05, standardize=TRUE,
                  trace=FALSE),
fit.gren2 <- gren(xtrain, ytrain, unpenalized=utrain,
                  partitions=list(part=part), alpha=0.5, standardize=TRUE,
                  trace=FALSE),
fit.gren3 <- gren(xtrain, ytrain, unpenalized=utrain,
                  partitions=list(part=part), alpha=0.95, standardize=TRUE,
                  trace=FALSE), {

init.grridge <- grridge(t(scale(micrornas)), benefit,
                        list(part=split(1:ncol(micrornas), diff.expr)),
                        unpenal= ~ 1 + adjth2 + thscheme2 + thscheme3 + age +
                          pcrcdiff3, optl=NULL,
                        dataunpen=as.data.frame(unpenal));

fit.grridge <- grridge(t(xtrain), ytrain, list(part=split(1:p, part)),
                       unpenal= ~ 1 + adjth2 + thscheme2 + thscheme3 + age +
                         pcrcdiff3, optl=init.grridge$optl,
                       dataunpen=as.data.frame(utrain))}, {

fit.sgl1 <- cvSGL(list(x=xtrain, y=ytrain), part, type="logit", alpha=0.05);
fit.sgl1$fit$type <- "logit"}, {
fit.sgl2 <- cvSGL(list(x=xtrain, y=ytrain), part, type="logit", alpha=0.5);
fit.sgl2$fit$type <- "logit"}, {
fit.sgl3 <- cvSGL(list(x=xtrain, y=ytrain), part, type="logit", alpha=0.95);
fit.sgl3$fit$type <- "logit"},

fit.cmcp1 <- cv.grpreg(cbind(utrain, xtrain), ytrain,
                       c(rep(0, ncol(utrain)), part), penalty="cMCP",
                       family="binomial", alpha=0.05),
fit.cmcp2 <- cv.grpreg(cbind(utrain, xtrain), ytrain,
                       c(rep(0, ncol(utrain)), part), penalty="cMCP",
                       family="binomial", alpha=0.5),
fit.cmcp3 <- cv.grpreg(cbind(utrain, xtrain), ytrain,
                       c(rep(0, ncol(utrain)), part), penalty="cMCP",
                       family="binomial", alpha=0.95),

fit.gel1 <- cv.grpreg(cbind(utrain, xtrain), ytrain,
                      c(rep(0, ncol(utrain)), part), penalty="gel",
                      family="binomial", alpha=0.05), {
# gel2 and gel3 give just the saturated model, so we create cv.grpreg ourselves
fit.gel2 <- grpreg(cbind(utrain, xtrain), ytrain,
                   c(rep(0, ncol(utrain)), part), penalty="gel",
                   family="binomial", alpha=0.5);
fit.gel2 <- list(cve=NA, cvse=NA, lambda=fit.gel2$lambda, fit=fit.gel2,
                 min=1, lambda.min=fit.gel2$lambda, null.dev=NA, pe=NA);
class(fit.gel2) <- "cv.grpreg"}, {
fit.gel3 <- grpreg(cbind(utrain, xtrain), ytrain,
                   c(rep(0, ncol(utrain)), part), penalty="gel",
                   family="binomial", alpha=0.95);
fit.gel3 <- list(cve=NA, cvse=NA, lambda=fit.gel3$lambda, fit=fit.gel3,
                 min=1, lambda.min=fit.gel3$lambda, null.dev=NA, pe=NA);
<<<<<<< HEAD
class(fit.gel3) <- "cv.grpreg"}, times=1, control=list(order="inorder"))

rownames(bench) <- c(paste0("gren", 1:3), "grridge", paste0("sgl", 1:3),
                      paste0("cmcp", 1:3), paste0("gel", 1:3))
=======
class(fit.gel3) <- "cv.grpreg"}, times=1, control=list(order="inorder"),
fit.rf <- rfsrc(y ~ ., data=data.frame(y=ytrain, u=utrain, x=xtrain), 
                var.used="all.trees", ntree=5000, importance="none"))

fit.gren4 <- gren(xtrain, ytrain, unpenalized=utrain,
                  partitions=list(part=part2), alpha=0.05, standardize=TRUE,
                  trace=FALSE)
fit.gren5 <- gren(xtrain, ytrain, unpenalized=utrain,
                  partitions=list(part=part2), alpha=0.5, standardize=TRUE,
                  trace=FALSE)
fit.gren6 <- gren(xtrain, ytrain, unpenalized=utrain,
                  partitions=list(part=part2), alpha=0.95, standardize=TRUE,
                  trace=FALSE)

init.grridge2 <- grridge(t(scale(micrornas)), benefit,
                         list(part=split(1:ncol(micrornas), diff.expr2)),
                         unpenal= ~ 1 + adjth2 + thscheme2 + thscheme3 + age +
                           pcrcdiff3, optl=NULL,
                         dataunpen=as.data.frame(unpenal))

fit.grridge2 <- grridge(t(xtrain), ytrain, list(part=split(1:p, part2)),
                        unpenal= ~ 1 + adjth2 + thscheme2 + thscheme3 + age +
                          pcrcdiff3, optl=init.grridge2$optl,
                        dataunpen=as.data.frame(utrain))

rownames(bench) <- c(paste0("gren", 1:3), "grridge", paste0("sgl", 1:3),
                      paste0("cmcp", 1:3), paste0("gel", 1:3), "rf")
>>>>>>> 31db1da7fc892df40e05ce4287d82f38e942befc
write.table(bench, file="results/micrornaseq_colorectal_cancer_bench1.csv")

save(fit.grridge, fit.gren1, fit.gren2, fit.gren3, fit.sgl1, fit.sgl2,
     fit.sgl3, fit.cmcp1, fit.cmcp2, fit.cmcp3, fit.gel1, fit.gel2,
<<<<<<< HEAD
     fit.gel3, file="results/micrornaseq_colorectal_cancer_fit1.Rdata")
=======
     fit.gel3, fit.rf, fit.gren4, fit.gren5, fit.gren6, fit.grridge2,
     file="results/micrornaseq_colorectal_cancer_fit1.Rdata")
>>>>>>> 31db1da7fc892df40e05ce4287d82f38e942befc

### prediction on test set
pred <- data.frame(ridge=predict.grridge(fit.grridge, t(xtest), FALSE, 
                                         as.data.frame(utest))[, 1],
                   grridge=predict.grridge(fit.grridge, t(xtest), FALSE,
                                           as.data.frame(utest))[, 2],
                   gren1=predict(fit.gren1, xtest, utest, type="groupreg"),
                   gren2=predict(fit.gren2, xtest, utest, type="groupreg"),
                   gren3=predict(fit.gren3, xtest, utest, type="groupreg"),
                   enet1=predict(fit.gren1, xtest, utest, type="regular"),
                   enet2=predict(fit.gren2, xtest, utest, type="regular"),
                   enet3=predict(fit.gren3, xtest, utest, type="regular"),
                   sgl1=predictSGL(fit.sgl1$fit, xtest),
                   sgl2=predictSGL(fit.sgl2$fit, xtest),
                   sgl3=predictSGL(fit.sgl3$fit, xtest),
                   cmcp1=predict(fit.cmcp1$fit, cbind(utest, xtest),
                                 type="response"),
                   cmcp2=predict(fit.cmcp2$fit, cbind(utest, xtest),
                                 type="response"),
                   cmcp3=predict(fit.cmcp3$fit, cbind(utest, xtest),
                                 type="response"),
                   gel1=predict(fit.gel1$fit, cbind(utest, xtest),
                                type="response"),
                   gel2=predict(fit.gel2$fit, cbind(utest, xtest),
                                type="response"),
                   gel3=predict(fit.gel3$fit, cbind(utest, xtest),
<<<<<<< HEAD
                                type="response"))
=======
                                type="response"),
                   rf=predict(fit.rf, data.frame(u=utest, x=xtest))$predicted,
                   gren4=predict(fit.gren4, xtest, utest, type="groupreg"),
                   gren5=predict(fit.gren5, xtest, utest, type="groupreg"),
                   gren6=predict(fit.gren6, xtest, utest, type="groupreg"),
                   grridge2=predict.grridge(fit.grridge2, t(xtest), FALSE,
                                            as.data.frame(utest))[, 2])

>>>>>>> 31db1da7fc892df40e05ce4287d82f38e942befc
psel <- c(ridge=p, grridge=p,
          gren1=fit.gren1$freq.model$groupreg$df - u,
          gren2=fit.gren2$freq.model$groupreg$df - u,
          gren3=fit.gren3$freq.model$groupreg$df - u,
          enet1=fit.gren1$freq.model$regular$df - u,
          enet2=fit.gren2$freq.model$regular$df - u,
          enet3=fit.gren3$freq.model$regular$df - u,
          sgl1=colSums(fit.sgl1$fit$beta!=0),
          sgl2=colSums(fit.sgl2$fit$beta!=0),
          sgl3=colSums(fit.sgl3$fit$beta!=0),
          cmcp1=colSums(fit.cmcp1$fit$beta[-c(1:(u + 1)), ]!=0),
          cmcp2=colSums(fit.cmcp2$fit$beta[-c(1:(u + 1)), ]!=0),
          cmcp3=colSums(fit.cmcp3$fit$beta[-c(1:(u + 1)), ]!=0),
          gel1=colSums(fit.gel1$fit$beta[-c(1:(u + 1)), ]!=0),
          gel2=colSums(as.matrix(fit.gel2$fit$beta[-c(1:(u + 1)), ]!=0)),
<<<<<<< HEAD
          gel3=colSums(as.matrix(fit.gel3$fit$beta[-c(1:(u + 1)), ]!=0)))
=======
          gel3=colSums(as.matrix(fit.gel3$fit$beta[-c(1:(u + 1)), ]!=0)),
          rf=p,
          gren4=fit.gren4$freq.model$groupreg$df - u,
          gren5=fit.gren5$freq.model$groupreg$df - u,
          gren6=fit.gren6$freq.model$groupreg$df - u,
          grridge2=p)
>>>>>>> 31db1da7fc892df40e05ce4287d82f38e942befc
auc <- apply(pred, 2, function(m) {pROC::auc(ytest, m)})
briers <- apply(pred, 2, function(m) {
  1 - mean((m - ytest)^2)/mean((mean(ytest) - ytest)^2)})
res <- rbind(pred, psel, auc, briers)
rownames(res) <- c(paste0("pred", c(1:length(ytest))), paste0("psel", 1),
                   paste0("auc", 1), paste0("briers", 1))
<<<<<<< HEAD
=======

>>>>>>> 31db1da7fc892df40e05ce4287d82f38e942befc
write.table(res, file="results/micrornaseq_colorectal_cancer_res1.csv")


################################### model 2 ####################################
###  create random splits of features into 3 groups of the original sizes
nsplits <- 100
y <- benefit
x <- scale(micrornas)
unpenal <- unpenal
p <- ncol(x)
u <- ncol(unpenal)
part <- diff.expr

mults <- data.frame(grridge=matrix(NA, nrow=nsplits,
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

  fit.gren1 <- gren(x, y, unpenalized=unpenal,
                    partitions=list(part=part.train), alpha=0.05,
                    standardize=TRUE, trace=FALSE)
  fit.gren2 <- gren(x, y, unpenalized=unpenal, 
                    partitions=list(part=part.train), alpha=0.5,
                    standardize=TRUE, trace=FALSE)
  fit.gren3 <- gren(x, y, unpenalized=unpenal,
                    partitions=list(part=part.train), alpha=0.95,
                    standardize=TRUE, trace=FALSE)

  fit.grridge <- grridge(t(x), y, list(part=split(1:p, part.train)),
                         unpenal= ~ 1 + adjth2 + thscheme2 + thscheme3 + age +
                           pcrcdiff3, dataunpen=as.data.frame(unpenal))

  mults[k, c(1:length(unique(part)))] <- fit.grridge$lambdamults$part
  mults[k, c((length(unique(part)) + 1):(2*length(unique(part))))] <-
    fit.gren1$lambdag$part
  mults[k, c((2*length(unique(part)) + 1):(3*length(unique(part))))] <-
    fit.gren2$lambdag$part
  mults[k, c((3*length(unique(part)) + 1):(4*length(unique(part))))] <-
    fit.gren3$lambdag$part

}
res <- mults
rownames(res) <- paste0("split", c(1:nsplits))
write.table(res, file="results/micrornaseq_colorectal_cancer_res2.csv")


################################### model 3 ####################################
###  create random splits of features into 10 evenly sized groups
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

mults <- data.frame(grridge=matrix(NA, nrow=nsplits,
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

  fit.gren1 <- gren(x, y, unpenalized=unpenal,
                    partitions=list(part=part.train), alpha=0.05,
                    standardize=TRUE, trace=FALSE)
  fit.gren2 <- gren(x, y, unpenalized=unpenal,
                    partitions=list(part=part.train), alpha=0.5,
                    standardize=TRUE, trace=FALSE)
  fit.gren3 <- gren(x, y, unpenalized=unpenal,
                    partitions=list(part=part.train), alpha=0.95,
                    standardize=TRUE, trace=FALSE)

  fit.grridge <- grridge(t(x), y, list(part=split(1:p, part.train)),
                         unpenal= ~ 1 + adjth2 + thscheme2 + thscheme3 + age +
                           pcrcdiff3, dataunpen=as.data.frame(unpenal))

  mults[k, c(1:length(unique(part)))] <- fit.grridge$lambdamults$part
  mults[k, c((length(unique(part)) + 1):(2*length(unique(part))))] <-
    fit.gren1$lambdag$part
  mults[k, c((2*length(unique(part)) + 1):(3*length(unique(part))))] <-
    fit.gren2$lambdag$part
  mults[k, c((3*length(unique(part)) + 1):(4*length(unique(part))))] <-
    fit.gren3$lambdag$part

}
res <- mults
rownames(res) <- paste0("split", c(1:nsplits))
write.table(res, file="results/micrornaseq_colorectal_cancer_res3.csv")


################################### model 4 ####################################
###  create random bootstrap samples of the data
nsamples <- 50
psel <- 25
y <- benefit
x <- scale(micrornas)
unpenal <- unpenal
p <- ncol(x)
u <- ncol(unpenal)
ngroups <- 10
part <- diff.expr

################################################################################
######### Error in { : task 46 failed - "inv(): matrix seems singular" #########
################################################################################
# RcppArmadillo in .i() error, thrown by op_inv::apply_noalias() in 
# op_inv_meat.hpp. Maybe increase lambda by hand?

### analyse samples in parallel
ncores <- min(detectCores() - 1, nsamples)
cluster <- makeForkCluster(ncores)
if(parallel) {
  registerDoParallel(cluster)
} else {
  registerDoSEQ()
}

selnames <- foreach(k=c(1:nsamples), .errorhandling="pass") %dopar% {
  set.seed(2019 + k)
  
  id.boot <- c(sample(which(y==0), replace=TRUE), 
               sample(which(y==1), replace=TRUE))
  # remove the constant micrornas
  is.const <- apply(x[id.boot, ], 2, sd)==0
  xboot <- scale(x[id.boot, !is.const])
  yboot <- y[id.boot]
  uboot <- unpenal[id.boot, ]
  part.boot <- part[!is.const]
  
  fit.gren1 <- gren(xboot, yboot, unpenalized=uboot, 
                    partitions=list(part=part.boot), alpha=0.05, 
                    standardize=TRUE, trace=FALSE, psel=psel)
  fit.gren2 <- gren(xboot, yboot, unpenalized=uboot, 
                    partitions=list(part=part.boot), alpha=0.5, 
                    standardize=TRUE, trace=FALSE, psel=psel)
  fit.gren3 <- gren(xboot, yboot, unpenalized=uboot, 
                    partitions=list(part=part.boot), alpha=0.95, 
                    standardize=TRUE, trace=FALSE, psel=psel)

  list(gren1=colnames(xboot)[fit.gren1$freq.model$groupreg$beta[-c(1:u), ]!=0],
       gren2=colnames(xboot)[fit.gren2$freq.model$groupreg$beta[-c(1:u), ]!=0],
       gren3=colnames(xboot)[fit.gren3$freq.model$groupreg$beta[-c(1:u), ]!=0],
       enet1=colnames(xboot)[fit.gren1$freq.model$regular$beta[-c(1:u), ]!=0],
       enet2=colnames(xboot)[fit.gren2$freq.model$regular$beta[-c(1:u), ]!=0],
       enet3=colnames(xboot)[fit.gren3$freq.model$regular$beta[-c(1:u), ]!=0])
    
}
if(parallel) {stopCluster(cluster)}

################################## DEBGUGGING ##################################
test1 <- selnames[-46]
test2 <- sapply(c(paste0("gren", 1:3), paste0("enet", 1:3)), function(m) {
  sapply(test1, function(s) {s[[grep(m, names(s))]]}, simplify=FALSE)},
  simplify=FALSE)
intersect <- sapply(test2, function(m) {
  combn(1:length(test1), 2, function(s) {
    length(intersect(m[[s[1]]], m[[s[2]]]))})})
res <- intersect
rownames(res) <- paste0("combn", c(1:choose(length(test1), 2)))
write.table(res, file="results/micrornaseq_colorectal_cancer_res4.csv")
################################################################################

selnames <- sapply(c(paste0("gren", 1:3), paste0("enet", 1:3)), function(m) {
  sapply(selnames, function(s) {s[[grep(m, names(s))]]}, simplify=FALSE)},
  simplify=FALSE)

### calculate the size of all intersections of selected features
intersect <- sapply(selnames, function(m) {
  combn(1:length(m), 2, function(s) {
    length(intersect(m[[s[1]]], m[[s[2]]]))})})

res <- intersect
rownames(res) <- paste0("combn", c(1:choose(length(selnames), 2)))
write.table(res, file="results/micrornaseq_colorectal_cancer_res4.csv")
