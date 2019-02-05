##################################### notes ####################################
# -results3 and fitted1 are with the 103 subjects                              #
# -results4 and fitted2 are with the 88 subjects                               #
# -results5 is with the 88 subjects and partition in three groups              #
# -results6 is with the 88 subjects and partition in three random groups       #
################################################################################

# save(pInc, file="~/results.Rdata")

### paths
path.data <- ifelse(as.character(Sys.info()[1])!="Darwin", "~/EBEN/data/",
                    "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/data/")
path.code <- ifelse(as.character(Sys.info()[1])!="Darwin", "~/EBEN/code/",
                    "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/code/")
path.res <- ifelse(as.character(Sys.info()[1])!="Darwin", "~/EBEN/results/",
                   "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/results/")
path.graph <- "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/graphs/"

### loading libraries
library(gren)
library(GRridge)
library(sp)
library(grpreg)
library(SGL)

### load functions
# source(paste(path.code, "grVBEM.R", sep=""))

### load data
# load(paste(path.data, "forMagnus.Rdata", sep=""))
load(paste(path.data, "forMagnusN88.Rdata", sep=""))

### set color scheme for graphs
# colors <- bpy.colors(6)[-c(1, 6)]
colors <- bpy.colors(8)[-c(1, 2, 8)]

### data manipulation
n <- length(resp)
p <- nrow(mirnormcen_resp)
unpenal <- model.matrix(~ adjth + thscheme + age + pcrcdiff, data=datfr)[, -1]
u <- ncol(unpenal)

### create partitioning
part.grridge <- list(TS=CreatePartition(as.character(which(!is.na(whichin))),
                                        as.character(1:p)))
part.greben <- list(TS=rep(c(1, 2), times=unlist(lapply(partkeep$TS, length)))[
  order(unlist(partkeep$TS))])

miRNA.BFDR <- as.character(TumMirs$miRNA)[TumMirs$BFDR_PNminP < 0.001]
miRNA.TumMirs <- as.character(TumMirs$miRNA)
miRNA <- as.character(sapply(rownames(mirnormcen_resp), function(s) {
  strsplit(s, split=" ")[[1]][1]}))

TS <- miRNA %in% miRNA.BFDR + miRNA %in% miRNA.TumMirs + 1
part.greben2 <- list(TS=miRNA %in% miRNA.BFDR + miRNA %in% miRNA.TumMirs + 1)
part.grridge2 <- list(TS=CreatePartition(as.factor(part.greben2$TS)))

## fitting models
set.seed(2018)
# fit1.enet <- cv.glmnet(cbind(unpenal, t(mirnormcen_resp)), resp, alpha=0.05,
#                        standardize=FALSE, family="binomial", intercept=TRUE,
#                        penalty.factor=c(rep(0, u), rep(1, p)))
# fit2.enet <- cv.glmnet(cbind(unpenal, t(mirnormcen_resp)), resp, alpha=0.5,
#                        standardize=FALSE, family="binomial", intercept=TRUE,
#                        penalty.factor=c(rep(0, u), rep(1, p)))
# fit3.enet <- cv.glmnet(cbind(unpenal, t(mirnormcen_resp)), resp, alpha=0.95,
#                        standardize=FALSE, family="binomial", intercept=TRUE,
#                        penalty.factor=c(rep(0, u), rep(1, p)))
# fit1.ridge <- cv.glmnet(cbind(unpenal, t(mirnormcen_resp)), resp, alpha=0,
#                         standardize=FALSE, family="binomial", intercept=TRUE,
#                         penalty.factor=c(rep(0, u), rep(1, p)))
# psel <- c(seq(1, 5, 1), seq(7, 15, 2), seq(20, 40, 5), seq(50, 90, 10))
# fit1.grridge <- vector("list", length(psel))
# fit1.grridge[[1]] <- grridge(mirnormcen_resp, resp, part.grridge, optl=NULL,
#                              unpenal = ~1 + adjth + thscheme + age + pcrcdiff,
#                              niter=1, method="exact", dataunpen = datfr,
#                              innfold=10, savepredobj="all", comparelasso=FALSE,
#                              optllasso=NULL, compareunpenal=TRUE,
#                              selection=TRUE, maxsel=psel[1])
# for(s in 2:length(psel)) {
#   fit1.grridge[[s]] <- grridge(mirnormcen_resp, resp, part.grridge,
#                                optl=fit1.grridge[[1]]$optl,
#                                unpenal = ~1 + adjth + thscheme + age + pcrcdiff,
#                                niter=1, method="exact", dataunpen = datfr,
#                                innfold=10, savepredobj="all",
#                                comparelasso=FALSE, optllasso=NULL,
#                                compareunpenal=TRUE, selection=TRUE,
#                                maxsel=psel[s])
# }
# fit1.greben <- grEBEN3(t(mirnormcen_resp), as.numeric(resp) - 1,
#                        rep(1, length(resp)), unpenalized=unpenal,
#                        partitions=part.greben, alpha=0.05, psel=psel,
#                        nfolds=NULL, posterior=FALSE, ELBO=FALSE, eps=0.001,
#                        maxiter=500, trace=TRUE)
# fit2.greben <- grEBEN3(t(mirnormcen_resp), as.numeric(resp) - 1,
#                        rep(1, length(resp)), unpenalized=unpenal,
#                        partitions=part.greben, alpha=0.5, psel=psel,
#                        nfolds=NULL, posterior=FALSE, ELBO=FALSE, eps=0.001,
#                        maxiter=500, trace=TRUE)
# fit3.greben <- grEBEN3(t(mirnormcen_resp), as.numeric(resp) - 1,
#                        rep(1, length(resp)), unpenalized=unpenal,
#                        partitions=part.greben, alpha=0.95, psel=psel,
#                        nfolds=NULL, posterior=FALSE, ELBO=FALSE, eps=0.001,
#                        maxiter=500, trace=TRUE)
part.grplasso <- c(rep(1, u), part.greben$TS + 1)
seq.grplasso <- grpreg(t(mirnormcen_resp), as.numeric(resp) - 1,
                       group=part.greben$TS, family="binomial")$lambda
lambda.grplasso <- exp(seq(log(max(seq.grplasso)), log(10e-4*max(seq.grplasso)), 
                    length.out=100))
fit.grplasso <- cv.grpreg(cbind(unpenal, t(mirnormcen_resp)), 
                          as.numeric(resp) - 1, group=part.grplasso, 
                          family="binomial", lambda=lambda.grplasso, 
                          group.multiplier=c(0, sqrt(
                            rle(sort(part.greben$TS))$lengths)))
(grplasso.nzero <- apply(coef(fit.grplasso, lambda=fit.grplasso$lambda), 2, 
                         function(b) {sum(b!=0)}))

library(oem)
fit.sglasso <- oem(cbind(unpenal, t(mirnormcen_resp)), as.numeric(resp) - 1,
                   family="binomial", penalty="sparse.grp.lasso", tau=0.5,
                   groups=part.grplasso, penalty.factor=c(rep(0, u), rep(1, p)))
install_version("sglOptim", version="1.3.6", repos="http://cran.us.r-project.org")
install_version("msgl", version="2.3.6", repos="http://cran.us.r-project.org")
library(msgl)
fit.sglasso <- msgl::fit(cbind(unpenal, t(mirnormcen_resp)), resp,
                    grouping=part.grplasso, lambda=1)

coef(fit.sglasso)[[1]]@x


install_github("jlaria/sglfast")
library(sglfast)




fit1.grridge2 <- vector("list", length(psel))
fit1.grridge2[[1]] <- grridge(mirnormcen_resp, resp, part.grridge2, optl=NULL,
                              unpenal = ~1 + adjth + thscheme + age + pcrcdiff,
                              niter=1, method="exact", dataunpen = datfr,
                              innfold=10, savepredobj="all", comparelasso=FALSE,
                              optllasso=NULL, compareunpenal=TRUE,
                              selection=TRUE, maxsel=psel[1])
for(s in 2:length(psel)) {
  fit1.grridge2[[s]] <- grridge(mirnormcen_resp, resp, part.grridge2,
                                optl=fit1.grridge2[[1]]$optl, unpenal = ~1 +
                                  adjth + thscheme + age + pcrcdiff, niter=1,
                                method="exact", dataunpen = datfr, innfold=10,
                                savepredobj="all", comparelasso=FALSE,
                                optllasso=NULL, compareunpenal=TRUE,
                                selection=TRUE, maxsel=psel[s])
}

fit1.greben2 <- grEBEN3(t(mirnormcen_resp), as.numeric(resp) - 1,
                        rep(1, length(resp)), unpenalized=unpenal,
                        partitions=part.greben2, alpha=0.05, psel=psel,
                        nfolds=NULL, posterior=FALSE, ELBO=FALSE, eps=0.001,
                        maxiter=500, trace=TRUE)
fit2.greben2 <- grEBEN3(t(mirnormcen_resp), as.numeric(resp) - 1,
                        rep(1, length(resp)), unpenalized=unpenal,
                        partitions=part.greben2, alpha=0.5, psel=psel,
                        nfolds=NULL, posterior=FALSE, ELBO=FALSE, eps=0.001,
                        maxiter=500, trace=TRUE)
fit3.greben2 <- grEBEN3(t(mirnormcen_resp), as.numeric(resp) - 1,
                        rep(1, length(resp)), unpenalized=unpenal,
                        partitions=part.greben2, alpha=0.95, psel=psel,
                        nfolds=NULL, posterior=FALSE, ELBO=FALSE, eps=0.001,
                        maxiter=500, trace=TRUE)

save(fit1.ridge, fit1.grridge, fit1.greben, fit2.greben, fit3.greben,
     fit1.grridge2, fit1.greben2, fit2.greben2, fit3.greben2,
     file=paste(path.res, "grEBEN_mirseq_Maarten_fitted2.Rdata", sep=""))

mirsel <- t(mirnormcen_resp)[, which(miRNA %in% miRNA.TumMirs)]
pmirsel <- ncol(mirsel)

fit1.ridge2 <- cv.glmnet(cbind(unpenal, mirsel), resp, alpha=0,
                         standardize=FALSE, family="binomial", intercept=TRUE,
                         penalty.factor=c(rep(0, u), rep(1, pmirsel)))
fit1.enet <- cv.glmnet(cbind(unpenal, mirsel), resp, alpha=0.05,
                       standardize=FALSE, family="binomial", intercept=TRUE,
                       penalty.factor=c(rep(0, u), rep(1, pmirsel)))
fit2.enet <- cv.glmnet(cbind(unpenal, mirsel), resp, alpha=0.5,
                       standardize=FALSE, family="binomial", intercept=TRUE,
                       penalty.factor=c(rep(0, u), rep(1, pmirsel)))
fit3.enet <- cv.glmnet(cbind(unpenal, mirsel), resp, alpha=0.95,
                       standardize=FALSE, family="binomial", intercept=TRUE,
                       penalty.factor=c(rep(0, u), rep(1, pmirsel)))



# ### cross validation
# set.seed(2018)
# nfolds <- n
# rest <- n %% nfolds
# foldid <- sample(rep(1:nfolds, times=round(c(rep(
#   n %/% nfolds + as.numeric(rest!=0), times=rest),
#   rep(n %/% nfolds, times=nfolds - rest)))))
# 
# psel <- c(seq(1, 5, 1), seq(7, 15, 2), seq(20, 40, 5), seq(50, 90, 10))
# methods <- c("ridge", "grridge", "enet1", "enet2", "enet3", "greben1",
#              "greben2", "greben3")
# pred.ridge <- numeric(n)
# for(m in 2:length(methods)) {
#   assign(paste("pred.", methods[m], sep=""), matrix(NA, nrow=n,
#                                                     ncol=length(psel)))
#   assign(paste("psel.", methods[m], sep=""), matrix(NA, nrow=n,
#                                                     ncol=length(psel)))
# }
# 
# fitinit.grridge <- grridge(mirnormcen_resp, resp, part.grridge, optl=NULL,
#                            unpenal = ~1 + adjth + thscheme + age + pcrcdiff,
#                            niter=1, method="exact", dataunpen=datfr,
#                            innfold=10, savepredobj="all", comparelasso=FALSE,
#                            optllasso=NULL, compareunpenal=TRUE,
#                            selection=TRUE, maxsel=psel[1])
# fitinit.ridge <- cv.glmnet(cbind(unpenal, t(mirnormcen_resp)), resp, alpha=0,
#                            standardize=FALSE, family="binomial", intercept=TRUE,
#                            penalty.factor=c(rep(0, u), rep(1, p)))
# fitinit1.enet <- cv.glmnet(cbind(unpenal, t(mirnormcen_resp)), resp, alpha=0.05,
#                            standardize=FALSE, family="binomial", intercept=TRUE,
#                            penalty.factor=c(rep(0, u), rep(1, p)))
# fitinit2.enet <- cv.glmnet(cbind(unpenal, t(mirnormcen_resp)), resp, alpha=0.5,
#                            standardize=FALSE, family="binomial", intercept=TRUE,
#                            penalty.factor=c(rep(0, u), rep(1, p)))
# fitinit3.enet <- cv.glmnet(cbind(unpenal, t(mirnormcen_resp)), resp, alpha=0.95,
#                            standardize=FALSE, family="binomial", intercept=TRUE,
#                            penalty.factor=c(rep(0, u), rep(1, p)))
#
# for(k in sort(unique(foldid))) {
#   cat(paste("Fold ", k, "\n"))
# 
#   xtrain <- t(mirnormcen_resp)[foldid!=k, ]
#   xtest <- matrix(t(mirnormcen_resp)[foldid==k, ], ncol=p, byrow=TRUE)
#   ytrain <- (as.numeric(resp) - 1)[foldid!=k]
#   utrain1 <- datfr[foldid!=k, ]
#   utest1 <- datfr[foldid==k, ]
#   utrain2 <- unpenal[foldid!=k, ]
#   utest2 <- matrix(unpenal[foldid==k, ], ncol=u, byrow=TRUE)
# 
#   cv.ridge <- glmnet(cbind(utrain2, xtrain), ytrain, alpha=0, 
#                      standardize=FALSE, family="binomial", intercept=TRUE,
#                      penalty.factor=c(rep(0, u), rep(1, p)), 
#                      lambda=fitinit.ridge$lambda.min)
#   cv.grridge <- vector("list", length(psel))
#   cv.grridge[[1]] <- grridge(t(xtrain), ytrain, part.grridge, optl=fitinit.grridge$optl,
#                              unpenal = ~1 + adjth + thscheme + age + pcrcdiff,
#                              niter=1, method="exact", dataunpen=utrain1,
#                              innfold=10, savepredobj="all", comparelasso=FALSE,
#                              optllasso=NULL, compareunpenal=TRUE,
#                              selection=TRUE, maxsel=psel[1])
#   for(s in 2:length(psel)) {
#     cv.grridge[[s]] <- grridge(t(xtrain), ytrain, part.grridge,
#                                optl=fitinit.grridge$optl,
#                                unpenal = ~1 + adjth + thscheme + age + pcrcdiff,
#                                niter=1, method="exact", dataunpen=utrain1,
#                                innfold=10, savepredobj="all",
#                                comparelasso=FALSE, optllasso=NULL,
#                                compareunpenal=TRUE, selection=TRUE,
#                                maxsel=psel[s])
#   }
#   cv1.greben <- grEBEN3(xtrain, ytrain, rep(1, length(ytrain)),
#                         unpenalized=utrain2, partitions=part.greben, alpha=0.05,
#                         psel=psel, trace=TRUE,
#                         lambda=fitinit1.enet$lambda.min)
#   cv2.greben <- grEBEN3(xtrain, ytrain, rep(1, length(ytrain)),
#                         unpenalized=utrain2, partitions=part.greben, alpha=0.5,
#                         psel=psel, trace=TRUE,
#                         lambda=fitinit2.enet$lambda.min) # 0.251227 1.175005
#   cv3.greben <- grEBEN3(xtrain, ytrain, rep(1, length(ytrain)),
#                         unpenalized=utrain2, partitions=part.greben, alpha=0.95,
#                         psel=psel, trace=TRUE,
#                         lambda=fitinit3.enet$lambda.min)
# 
#   pred.enet1[foldid==k, ] <- 1/(1 + exp(-cbind(1, utest2, xtest) %*%
#                                           cv1.greben$beta.nogroups))
#   pred.greben1[foldid==k, ] <- 1/(1 + exp(-cbind(1, utest2, xtest) %*%
#                                             cv1.greben$beta))
#   pred.enet2[foldid==k, ] <- 1/(1 + exp(-cbind(1, utest2, xtest) %*%
#                                           cv2.greben$beta.nogroups))
#   pred.greben2[foldid==k, ] <- 1/(1 + exp(-cbind(1, utest2, xtest) %*%
#                                             cv2.greben$beta))
#   pred.enet3[foldid==k, ] <- 1/(1 + exp(-cbind(1, utest2, xtest) %*%
#                                           cv3.greben$beta.nogroups))
#   pred.greben3[foldid==k, ] <- 1/(1 + exp(-cbind(1, utest2, xtest) %*%
#                                             cv3.greben$beta))
#   pred.grridge[foldid==k, ] <- sapply(cv.grridge, function(s) {
#     predict.grridge(s, t(xtest), dataunpennew=utest1)[, 3]})
#   pred.ridge[foldid==k] <- predict(cv.ridge, cbind(utest2, xtest),
#                                    type="response")
# 
#   psel.grridge[foldid==k, ] <- sapply(cv.grridge, function(s) {
#     length(s$resEN$whichEN)})
#   psel.enet1[foldid==k, ] <- colSums(cv1.greben$beta.nogroups!=0) - u - 1
#   psel.greben1[foldid==k, ] <- colSums(cv1.greben$beta!=0) - u - 1
#   psel.enet2[foldid==k, ] <- colSums(cv2.greben$beta.nogroups!=0) - u - 1
#   psel.greben2[foldid==k, ] <- colSums(cv2.greben$beta!=0) - u - 1
#   psel.enet3[foldid==k, ] <- colSums(cv3.greben$beta.nogroups!=0) - u - 1
#   psel.greben3[foldid==k, ] <- colSums(cv3.greben$beta!=0) - u - 1
# 
# }
# 
# results4 <- list(pred=list(pred.ridge, pred.grridge, pred.enet1, pred.enet2,
#                            pred.enet3, pred.greben1, pred.greben2,
#                            pred.greben3),
#                  psel=list(psel.grridge, psel.enet1, psel.enet2, psel.enet3,
#                            psel.greben1, psel.greben2, psel.greben3))
# save(results4, file=paste(path.res, "grEBEN_mirseq_Maarten_res4.Rdata", sep=""))

# # random groups
# nsplits <- 100
# mult.grridge <- mult.greben1 <- mult.greben2 <- mult.greben3 <- 
#   matrix(NA, ncol=3, nrow=nsplits)
# for(k in 1:nsplits) {
#   cat(paste("Split ", k, "\n"))
#   
#   set.seed(2018 + k)
#   part.greben3 <- list(random=sample(part.greben2$TS))
#   part.grridge3 <- list(random=CreatePartition(as.factor(part.greben3$random)))
#   
#   cv1.greben <- grEBEN3(t(mirnormcen_resp), as.numeric(resp) - 1,
#                         rep(1, length(resp)), unpenalized=unpenal,
#                         partitions=part.greben3, alpha=0.05, psel=NULL,
#                         nfolds=NULL, posterior=FALSE, ELBO=FALSE, eps=0.001,
#                         maxiter=500, trace=TRUE)
#   cv2.greben <- grEBEN3(t(mirnormcen_resp), as.numeric(resp) - 1,
#                         rep(1, length(resp)), unpenalized=unpenal,
#                         partitions=part.greben3, alpha=0.5, psel=NULL,
#                         nfolds=NULL, posterior=FALSE, ELBO=FALSE, eps=0.001,
#                         maxiter=500, trace=TRUE)
#   cv3.greben <- grEBEN3(t(mirnormcen_resp), as.numeric(resp) - 1,
#                         rep(1, length(resp)), unpenalized=unpenal,
#                         partitions=part.greben3, alpha=0.95, psel=NULL,
#                         nfolds=NULL, posterior=FALSE, ELBO=FALSE, eps=0.001,
#                         maxiter=500, trace=TRUE)
# 
# 
#   cv.grridge <- grridge(mirnormcen_resp, resp, part.grridge3, optl=NULL,
#                         unpenal = ~1 + adjth + thscheme + age + pcrcdiff,
#                         niter=1, method="exact", dataunpen = datfr,
#                         innfold=10, savepredobj="all", comparelasso=FALSE,
#                         optllasso=NULL, compareunpenal=TRUE)
#   
#   mult.grridge[k, ] <- cv.grridge$lambdamults$random
#   mult.greben1[k, ] <- cv1.greben$lambdag$random[, cv1.greben$nouteriter + 1]
#   mult.greben2[k, ] <- cv2.greben$lambdag$random[, cv2.greben$nouteriter + 1]
#   mult.greben3[k, ] <- cv3.greben$lambdag$random[, cv3.greben$nouteriter + 1]
#   
# }

# results6 <- list(mult=list(mult.grridge, mult.greben1, mult.greben2,
#                            mult.greben3))
# # save(results6, file=paste(path.res, "grEBEN_mirseq_Maarten_res6.Rdata", sep=""))


### stability selection
K <- 50
psel <- 25
methods <- c("grridge", "enet1", "enet2", "enet3", "gren1",
             "gren2", "gren3")
for(s in 1:length(methods)) {
  assign(paste("sel.", methods[s], sep=""), vector("list", K))
}

for(k in 1:K) {
  cat(paste("Sample ", k, "\n"))
  
  set.seed(2018 + k)

  id <- c(sample(which(resp==levels(resp)[1]), replace=TRUE),
          sample(which(resp==levels(resp)[2]), replace=TRUE))
  xtrain <- t(mirnormcen_resp)[id, ]
  ytrain <- (as.numeric(resp) - 1)[id]
  mtrain <- rep(1, length(ytrain))
  utrain1 <- datfr[id, ]
  utrain2 <- unpenal[id, ]

  which.const1 <- apply(xtrain, 2, sd)==0
  if(any(which.const1)) {
    xtrain <- xtrain[, !which.const1]
    boot.TS <- TS[!which.const1]
  } else {
    boot.TS <- TS
  }

  which.const2.1 <- apply(utrain1, 2, sd)==0
  which.const2.2 <- apply(utrain2, 2, sd)==0
  if(any(which.const2.1)) {
    utrain1 <- utrain[, !which.const2.1]
    utrain2 <- utrain[, !which.const2.2]
  }
  
  boot.parTS <- CreatePartition(as.factor(boot.TS))

  # boot.grridge <- grridge(t(xtrain), ytrain, list(TS=boot.parTS),
  #                         ~1 + adjth + thscheme + age + pcrcdiff, niter=1,
  #                         method="exact", dataunpen=utrain1, innfold=10,
  #                         savepredobj="all", comparelasso=FALSE, optllasso=NULL,
  #                         compareunpenal=TRUE, selectionEN=TRUE, maxsel=psel)
  boot.gren1 <- tryCatch(gren(xtrain, ytrain, mtrain, utrain2, list(TS=boot.TS),
                              alpha=0.05, psel=psel, compare=TRUE), 
                         error=function(e) {
                           gren(xtrain, ytrain, mtrain, utrain2, lambda=0.01,
                                list(TS=boot.TS), alpha=0.05, psel=psel, 
                                compare=TRUE)}
                         )
  boot.gren2 <- tryCatch(gren(xtrain, ytrain, mtrain, utrain2, list(TS=boot.TS),
                              alpha=0.5, psel=psel, compare=TRUE), 
                         error=function(e) {
                           gren(xtrain, ytrain, mtrain, utrain2, lambda=0.01,
                                list(TS=boot.TS), alpha=0.5, psel=psel, 
                                compare=TRUE)}
                         )
  boot.gren3 <- tryCatch(gren(xtrain, ytrain, mtrain, utrain2, list(TS=boot.TS),
                              alpha=0.95, psel=psel, compare=TRUE), 
                         error=function(e) {
                           gren(xtrain, ytrain, mtrain, utrain2, lambda=0.01,
                                list(TS=boot.TS), alpha=0.95, psel=psel, 
                                compare=TRUE)})
  # boot.gren2 <- gren(xtrain, ytrain, mtrain, utrain2, list(TS=boot.TS),
  #                    alpha=0.5, psel=psel, compare=TRUE)
  # boot.gren3 <- gren(xtrain, ytrain, mtrain, utrain2, list(TS=boot.TS),
  #                    alpha=0.95, psel=psel, compare=TRUE)

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

# results7 <- list(sel=list(sel.grridge, sel.enet1, sel.enet2, sel.enet3,
#                           sel.gren1, sel.gren2, sel.gren3))
results7 <- list(sel=list(sel.enet1, sel.enet2, sel.enet3,
                          sel.gren1, sel.gren2, sel.gren3))
save(results7, file=paste(path.res, "grEBEN_mirseq_Maarten_res7.Rdata", sep=""))

# # random groups with 10 groups
# nsplits <- 100
# ngroups <- 10
# rest <- p %% ngroups
# mult.grridge <- mult.greben1 <- mult.greben2 <- mult.greben3 <-
#   matrix(NA, ncol=ngroups, nrow=nsplits)
# for(k in 1:nsplits) {
#   cat(paste("Split ", k, "\n"))
# 
#   set.seed(2018 + k)
#   part.greben3 <- list(random=sample(rep(c(1:ngroups), times=round(c(
#     rep(p %/% ngroups + as.numeric(rest!=0), times=rest), 
#     rep(p %/% ngroups, times=ngroups - rest))))))
#   part.grridge3 <- list(random=CreatePartition(as.factor(part.greben3$random)))
# 
#   cv1.greben <- grEBEN3(t(mirnormcen_resp), as.numeric(resp) - 1,
#                         rep(1, length(resp)), unpenalized=unpenal,
#                         partitions=part.greben3, alpha=0.05, psel=NULL,
#                         nfolds=NULL, posterior=FALSE, ELBO=FALSE, eps=0.001,
#                         maxiter=500, trace=TRUE)
#   cv2.greben <- grEBEN3(t(mirnormcen_resp), as.numeric(resp) - 1,
#                         rep(1, length(resp)), unpenalized=unpenal,
#                         partitions=part.greben3, alpha=0.5, psel=NULL,
#                         nfolds=NULL, posterior=FALSE, ELBO=FALSE, eps=0.001,
#                         maxiter=500, trace=TRUE)
#   cv3.greben <- grEBEN3(t(mirnormcen_resp), as.numeric(resp) - 1,
#                         rep(1, length(resp)), unpenalized=unpenal,
#                         partitions=part.greben3, alpha=0.95, psel=NULL,
#                         nfolds=NULL, posterior=FALSE, ELBO=FALSE, eps=0.001,
#                         maxiter=500, trace=TRUE)
# 
#   cv.grridge <- grridge(mirnormcen_resp, resp, part.grridge3, optl=NULL,
#                         unpenal = ~1 + adjth + thscheme + age + pcrcdiff,
#                         niter=1, method="exact", dataunpen = datfr,
#                         innfold=10, savepredobj="all", comparelasso=FALSE,
#                         optllasso=NULL, compareunpenal=TRUE)
# 
#   mult.grridge[k, ] <- cv.grridge$lambdamults$random
#   mult.greben1[k, ] <- cv1.greben$lambdag$random[, cv1.greben$nouteriter + 1]
#   mult.greben2[k, ] <- cv2.greben$lambdag$random[, cv2.greben$nouteriter + 1]
#   mult.greben3[k, ] <- cv3.greben$lambdag$random[, cv3.greben$nouteriter + 1]
# 
# }
# 
# results8 <- list(mult=list(mult.grridge, mult.greben1, mult.greben2,
#                            mult.greben3))
# save(results8, file=paste(path.res, "grEBEN_mirseq_Maarten_res8.Rdata", sep=""))

### random groups results
# load(paste(path.res, "grEBEN_mirseq_Maarten_res6.Rdata", sep=""))
# 
# boxdata <- data.frame(mults=c(as.numeric(results6$mult[[1]]), 
#                               as.numeric(results6$mult[[2]]), 
#                               as.numeric(results6$mult[[3]]), 
#                               as.numeric(results6$mult[[4]])),
#                       group=rep(rep(c(1:3), each=100), 4),
#                       method=rep(c(1:4), each=300))
# 
# leglabels <- c("GRridge", expression(paste("gren, ", alpha==0.05)),
#                expression(paste("gren, ", alpha==0.5)),
#                expression(paste("gren, ", alpha==0.95)))
# 
# png(paste(path.graph, "grEBEN_mirseq_Maarten_res6_boxplot.png", sep=""),
#     units="in", width=8, height=6, res=120)
# par(mar=c(5.1, 5.6, 4.1, 2.1))
# boxplot(mults ~ interaction(method, group), data=boxdata, 
#         at=c(c(1:4), c(6:9), c(11:14)), xaxt="n", col=rep(colors, 3), 
#         outline=FALSE, cex.lab=2, cex.names=1.5, 
#         ylab=expression({lambda^{"'"}}[g]), 
#         cex.axis=1.5, boxlwd=0.5)#,  border=rep(colors, 3))
# axis(1, c(2.5, 7.5, 12.5), c("Group 1", "Group 2", "Group 3"), tick=FALSE,
#      cex.axis=1.5)
# abline(h=1, lty=2, lwd=1.5)
# legend("topleft", legend=leglabels, fill=colors, border=rep(1, 4),
#        seg.len=1, cex=1.3)
# dev.off()
# 
# tabmed <- rbind(apply(results6$mult[[1]], 2, function(x) {round(median(x), 2)}),
#                 apply(results6$mult[[2]], 2, function(x) {round(median(x), 2)}),
#                 apply(results6$mult[[3]], 2, function(x) {round(median(x), 2)}),
#                 apply(results6$mult[[4]], 2, function(x) {round(median(x), 2)}))

# # 10 random groups results
# load(paste(path.res, "grEBEN_mirseq_Maarten_res8.Rdata", sep=""))
# 
# boxdata <- data.frame(mults=c(as.numeric(results8$mult[[1]]),
#                               as.numeric(results8$mult[[2]]),
#                               as.numeric(results8$mult[[3]]),
#                               as.numeric(results8$mult[[4]])),
#                       group=rep(rep(c(1:10), each=100), 4),
#                       method=rep(c(1:4), each=1000))
# 
# leglabels <- c("GRridge", expression(paste("gren, ", alpha==0.05)),
#                expression(paste("gren, ", alpha==0.5)),
#                expression(paste("gren, ", alpha==0.95)))
# 
# png(paste(path.graph, "grEBEN_mirseq_Maarten_res8_boxplot.png", sep=""),
#     units="in", width=12, height=6, res=120)
# par(mar=c(5.1, 5.6, 4.1, 2.1))
# boxplot(mults ~ interaction(method, group), data=boxdata,
#         at=as.numeric(sapply(c(1:10), function(s) {
#           (s - 1)*4 + (s - 1) + c(1:4)})), xaxt="n", col=rep(colors, 10),
#         outline=FALSE, cex.lab=2, cex.names=1.5,
#         ylab=expression({lambda^{"'"}}[g]), xlab="Group",
#         cex.axis=1.5, boxlwd=0.5)#,  border=rep(colors, 3))
# axis(1, 2.5 + c(0:9)*5, c(1:10), tick=FALSE,
#      cex.axis=1.5)
# abline(h=1, lty=2, lwd=1.5)
# legend("bottomright", legend=leglabels, fill=colors, border=rep(1, 4),
#        seg.len=1, cex=1.3)
# dev.off()
# 
# tabmed <- rbind(apply(results8$mult[[1]], 2, function(x) {round(median(x), 2)}),
#                 apply(results8$mult[[2]], 2, function(x) {round(median(x), 2)}),
#                 apply(results8$mult[[3]], 2, function(x) {round(median(x), 2)}),
#                 apply(results8$mult[[4]], 2, function(x) {round(median(x), 2)}))

### fitted models results
# load(paste(path.res, "grEBEN_mirseq_Maarten_fitted2.Rdata", sep=""))
# bsel1.greben <- apply(fit1.greben$beta, 2, function(b) {
#   which(b[-c(1:(u + 1))]!=0)})
# bsel2.greben <- apply(fit2.greben$beta, 2, function(b) {
#   which(b[-c(1:(u + 1))]!=0)})
# bsel3.greben <- apply(fit3.greben$beta, 2, function(b) {
#   which(b[-c(1:(u + 1))]!=0)})
# bsel1.enet <- apply(fit1.greben$beta.nogroups, 2, function(b) {
#   which(b[-c(1:(u + 1))]!=0)})
# bsel2.enet <- apply(fit2.greben$beta.nogroups, 2, function(b) {
#   which(b[-c(1:(u + 1))]!=0)})
# bsel3.enet <- apply(fit3.greben$beta.nogroups, 2, function(b) {
#   which(b[-c(1:(u + 1))]!=0)})
# bsel1.grridge <- sapply(fit1.grridge, function(s) {s$resEN$whichEN})
# 
# psel1.greben <- colSums(fit1.greben$beta!=0) - u - 1
# psel2.greben <- colSums(fit2.greben$beta!=0) - u - 1
# psel3.greben <- colSums(fit3.greben$beta!=0) - u - 1
# psel1.enet <- colSums(fit1.greben$beta.nogroups!=0) - u - 1
# psel2.enet <- colSums(fit2.greben$beta.nogroups!=0) - u - 1
# psel3.enet <- colSums(fit3.greben$beta.nogroups!=0) - u - 1
# psel1.grridge <- sapply(fit1.grridge, function(s) {length(s$resEN$whichEN)})
# 
# kap1.greben <- sapply(bsel1.greben, function(s) {
#   cohen.kappa(table(!is.na(whichin), replace(rep(FALSE, p), s, TRUE)))$kappa})
# kap2.greben <- sapply(bsel2.greben, function(s) {
#   cohen.kappa(table(!is.na(whichin), replace(rep(FALSE, p), s, TRUE)))$kappa})
# kap3.greben <- sapply(bsel3.greben, function(s) {
#   cohen.kappa(table(!is.na(whichin), replace(rep(FALSE, p), s, TRUE)))$kappa})
# kap1.enet <- sapply(bsel1.enet, function(s) {
#   cohen.kappa(table(!is.na(whichin), replace(rep(FALSE, p), s, TRUE)))$kappa})
# kap2.enet <- sapply(bsel2.enet, function(s) {
#   cohen.kappa(table(!is.na(whichin), replace(rep(FALSE, p), s, TRUE)))$kappa})
# kap3.enet <- sapply(bsel3.enet, function(s) {
#   cohen.kappa(table(!is.na(whichin), replace(rep(FALSE, p), s, TRUE)))$kappa})
# kap1.grridge <- sapply(bsel1.grridge, function(s) {
#   cohen.kappa(table(!is.na(whichin), replace(rep(FALSE, p), s, TRUE)))$kappa})
# 
# per1.greben <- sapply(bsel1.greben, function(s) {
#   mean(s %in% part.grridge$TS$VarIn)})
# per2.greben <- sapply(bsel2.greben, function(s) {
#   mean(s %in% part.grridge$TS$VarIn)})
# per3.greben <- sapply(bsel3.greben, function(s) {
#   mean(s %in% part.grridge$TS$VarIn)})
# per1.enet <- sapply(bsel1.enet, function(s) {
#   mean(s %in% part.grridge$TS$VarIn)})
# per2.enet <- sapply(bsel2.enet, function(s) {
#   mean(s %in% part.grridge$TS$VarIn)})
# per3.enet <- sapply(bsel3.enet, function(s) {
#   mean(s %in% part.grridge$TS$VarIn)})
# per1.grridge <- sapply(bsel1.grridge, function(s) {
#   mean(s %in% part.grridge$TS$VarIn)})
# 
# # cohens kappa against selected variables 2 groups
# plot(sort(psel1.greben), kap1.greben[order(psel1.greben)], type="l",
#      xlim=range(psel), col=2,
#      ylim=range(kap1.greben, kap2.greben, kap3.greben, kap1.enet, kap2.enet,
#                 kap3.enet, kap1.grridge), ylab="Cohen's kappa",
#      xlab="Number of selected variables")
# lines(sort(psel2.greben), kap2.greben[order(psel2.greben)], type="l", col=3)
# lines(sort(psel3.greben), kap3.greben[order(psel3.greben)], type="l", col=4)
# lines(psel1.enet, kap1.enet, type="l", col=2, lty=2)
# lines(psel2.enet, kap2.enet, type="l", col=3, lty=2)
# lines(psel3.enet, kap3.enet, type="l", col=4, lty=2)
# lines(psel1.grridge, kap1.grridge, type="l", col=5)
# legend("topleft", fill=c(c(2:5), 0, 0), lty=c(rep(NA, 4), 1, 2),
#        border=c(rep(1, 4), 0, 0), merge=TRUE, seg.len=1,
#        legend=c("EN, alpha=0.05", "EN, alpha=0.5", "EN, alpha=0.95", "GRridge",
#                 "Group-regularized", "Not group-regularized"))
# 
# # percentage in TS group, two groups
# plot(psel1.greben, per1.greben, type="l", xlim=range(psel), col=2,
#      ylim=range(per1.greben, per2.greben, per3.greben, per1.enet, per2.enet,
#                 per3.enet, per1.grridge), ylab="Percentage variables in TS",
#      xlab="Number of selected variables")
# lines(psel2.greben, per2.greben, type="l", col=3)
# lines(sort(psel3.greben), per3.greben[order(psel3.greben)], type="l", col=4)
# lines(psel1.enet, per1.enet, type="l", col=2, lty=2)
# lines(psel2.enet, per2.enet, type="l", col=3, lty=2)
# lines(psel3.enet, per3.enet, type="l", col=4, lty=2)
# lines(psel1.grridge, per1.grridge, type="l", col=5)
# legend("bottomleft", fill=c(c(2:5), 0, 0), lty=c(rep(NA, 4), 1, 2),
#        border=c(rep(1, 4), 0, 0), merge=TRUE, seg.len=1,
#        legend=c("EN, alpha=0.05", "EN, alpha=0.5", "EN, alpha=0.95", "GRridge",
#                 "Group-regularized", "Not group-regularized"))
# 
# # jaccard indexes of two groups
# line1 <- sapply(1:length(bsel1.greben), function(s) {
#   length(intersect(bsel1.greben[[s]], bsel2.greben[[s]]))/
#     length(union(bsel1.greben[[s]], bsel2.greben[[s]]))})
# line2 <- sapply(1:length(bsel1.greben), function(s) {
#   length(intersect(bsel1.greben[[s]], bsel3.greben[[s]]))/
#     length(union(bsel1.greben[[s]], bsel3.greben[[s]]))})
# line3 <- sapply(1:length(bsel2.greben), function(s) {
#   length(intersect(bsel2.greben[[s]], bsel3.greben[[s]]))/
#     length(union(bsel2.greben[[s]], bsel3.greben[[s]]))})
# line4 <- sapply(1:length(bsel1.greben), function(s) {
#   length(intersect(bsel1.greben[[s]], bsel1.grridge[[s]]))/
#     length(union(bsel1.greben[[s]], bsel1.grridge[[s]]))})
# line5 <- sapply(1:length(bsel2.greben), function(s) {
#   length(intersect(bsel2.greben[[s]], bsel1.grridge[[s]]))/
#     length(union(bsel2.greben[[s]], bsel1.grridge[[s]]))})
# line6 <- sapply(1:length(bsel3.greben), function(s) {
#   length(intersect(bsel3.greben[[s]], bsel1.grridge[[s]]))/
#     length(union(bsel3.greben[[s]], bsel1.grridge[[s]]))})
# 
# xlim <- range(psel)
# ylim <- c(min(line1, line2, line3, line4, line5, line6), 1)
# quartz()
# par(mfrow=c(4, 4), mar=c(2.1,4.2,2.1,2.1))
# plot(1, type="n", axes=F, xlab="", ylab="")
# text(1, 1, labels="EN, alpha=0.05")
# plot(1, type="n", axes=F, xlab="", ylab="")
# plot(1, type="n", axes=F, xlab="", ylab="")
# plot(1, type="n", axes=F, xlab="", ylab="")
# plot(psel, line1, type="l", xlim=xlim, ylim=ylim, xlab="", ylab="Jaccard index")
# plot(1, type="n", axes=F, xlab="", ylab="")
# text(1, 1, labels="EN, alpha=0.5")
# plot(1, type="n", axes=F, xlab="", ylab="")
# plot(1, type="n", axes=F, xlab="", ylab="")
# plot(psel, line2, type="l", xlim=xlim, ylim=ylim, xlab="", ylab="Jaccard index")
# plot(psel, line3, type="l", xlim=xlim, ylim=ylim, xlab="", ylab="")
# plot(1, type="n", axes=F, xlab="", ylab="")
# text(1, 1, labels="EN, alpha=0.95")
# par(mar=c(4.1,4.2,2.1,2.1))
# plot(1, type="n", axes=F, xlab="", ylab="")
# plot(psel, line4, type="l", xlim=xlim, ylim=ylim,
#      xlab="Number of selected variables", ylab="Jaccard index")
# plot(psel, line5, type="l", xlim=xlim, ylim=ylim,
#      xlab="Number of selected variables", ylab="")
# plot(psel, line6, type="l", xlim=xlim,
#      xlab="Number of selected variables", ylim=ylim, ylab="")
# plot(1, type="n", axes=F, xlab="", ylab="")
# text(1, 1, labels="GRridge")
# par(mfrow=c(1, 1), mar=c(5.1,4.1,4.1,2.1))
# 
# # comparing two to three groups settings
# bsel1.grridge2 <- sapply(fit1.grridge2, function(s) {s$resEN$whichEN})
# bsel1.greben2 <- apply(fit1.greben2$beta[-c(1:(u + 1)), ], 2, function(b) {
#   which(b!=0)})
# bsel2.greben2 <- apply(fit2.greben2$beta[-c(1:(u + 1)), ], 2, function(b) {
#   which(b!=0)})
# bsel3.greben2 <- apply(fit3.greben2$beta[-c(1:(u + 1)), ], 2, function(b) {
#   which(b!=0)})
# # 
# jaccard <- function(set1, set2) {
#   length(intersect(set1, set2))/length(union(set1, set2))
# }
# jac1.grridge <- sapply(1:length(fit1.grridge), function(s) {
#   jaccard(bsel1.grridge2[[s]], bsel1.grridge[[s]])})
# jac1.greben <- sapply(1:ncol(fit1.greben$beta), function(s) {
#   jaccard(bsel1.greben2[[s]], bsel1.greben[[s]])})
# jac2.greben <- sapply(1:ncol(fit2.greben$beta), function(s) {
#   jaccard(bsel2.greben2[[s]], bsel2.greben[[s]])})
# jac3.greben <- sapply(1:ncol(fit3.greben$beta), function(s) {
#   jaccard(bsel3.greben2[[s]], bsel3.greben[[s]])})
# 
# png(paste(path.graph, "grEBEN_mirseq_Maarten_jac1.png", sep=""), units="in",
#     width=5, height=5, res=200)
# plot(psel, jac1.greben, type="l", xlim=range(psel), col=2,
#      ylim=range(jac1.greben, jac2.greben, jac3.greben, jac1.grridge),
#      ylab="Jaccard index", xlab="Number of selected variables")
# lines(psel, jac2.greben, type="l", col=3)
# lines(psel, jac3.greben, type="l", col=4)
# lines(psel, jac1.grridge, type="l", col=5)
# legend("topright", lty=1, col=c(2:5),
#        legend=c("EN, alpha=0.05", "EN, alpha=0.5", "EN, alpha=0.95", "GRridge"))
# dev.off()




# # barplot of multipliers
# load(paste(path.res, "grEBEN_mirseq_Maarten_fitted2.Rdata", sep=""))
# leg.lab <- c("GRridge", expression(paste("gren, ", alpha==0.05)),
#              expression(paste("gren, ", alpha==0.5)),
#              expression(paste("gren, ", alpha==0.95)),
#              "not group-regularized")
# png(paste(path.graph, "grEBEN_mirseq_Maarten_bar1.png", sep=""), units="in",
#     width=14, height=5, res=200)
# par(mfrow=c(1, 2), mar=c(5.1, 5.6, 4.1, 2.1))
# barplot(rbind(fit1.grridge[[1]]$lambdamults$TS,
#               fit1.greben$lambdag$TS[, fit1.greben$nouteriter + 1],
#               fit2.greben$lambdag$TS[, fit2.greben$nouteriter + 1],
#               fit3.greben$lambdag$TS[, fit3.greben$nouteriter + 1]), beside=TRUE,
#         names.arg=c("high", "low"), col=colors, 
#         ylab=expression({lambda^{"'"}}[g]), xlab="microRNA expression", 
#         main="a)", legend.text=leg.lab, cex.axis=1.5, cex.names=1.5, cex.lab=2,
#         args.legend=list(x="topleft", fill=c(colors, 0),
#                          lty=c(rep(NA, 4), 2), lwd=c(rep(NA, 4), 1.5),
#                          border=c(rep(1, 4), 0), merge=TRUE, seg.len=1, 
#                          cex=1.3), cex.main=2)
# abline(h=1, lty=2, lwd=1.5)
# barplot(rbind(rev(fit1.grridge2[[1]]$lambdamults$TS),
#               rev(fit1.greben2$lambdag$TS[, fit1.greben2$nouteriter + 1]),
#               rev(fit2.greben2$lambdag$TS[, fit2.greben2$nouteriter + 1]),
#               rev(fit3.greben2$lambdag$TS[, fit3.greben2$nouteriter + 1])),
#         beside=TRUE, col=colors, 
#         names.arg=c("high", "medium", "low"),
#         ylab=expression({lambda^{"'"}}[g]), xlab="microRNA expression", 
#         main="b)", cex.axis=1.5, cex.names=1.5, cex.lab=2, cex.main=2)
# abline(h=1, lty=2, lwd=1.5)
# dev.off()





### cross-validation results
# leglabels <- c("ridge", expression(paste("enet, ", alpha==0.05)),
#                expression(paste("enet, ", alpha==0.5)),
#                expression(paste("enet, ", alpha==0.95)),
#                "group-regularized", "not group-regularized")
#
# # two groups
# load(paste(path.res, "grEBEN_mirseq_Maarten_res4.Rdata", sep=""))
# 
# auc4 <- lapply(results4$pred[-1], function(l) {
#   apply(l, 2, function(preds) {pROC::roc(as.numeric(resp) - 1, preds)$auc})})
# auc4 <- c(pROC::roc(as.numeric(resp) - 1, results4$pred[[1]])$auc, auc4)
# 
# briers4 <- lapply(results4$pred[-1], function(l) {
#   apply(l, 2, function(preds) {
#     1 - sum((as.numeric(resp) - 1 - preds)^2)/
#       sum((as.numeric(resp) - 1 - mean(as.numeric(resp) - 1))^2)})})
# briers4 <- c(1 - sum((as.numeric(resp) - 1 - results4$pred[[1]])^2)/
#                sum((as.numeric(resp) - 1 - mean(as.numeric(resp) - 1))^2),
#              briers4)
# 
# psel4 <- lapply(results4$psel, function(l) {colMeans(l)})
# 
# png(paste(path.graph, "grEBEN_mirseq_Maarten_res4_performance.png", sep=""),
#     units="in", width=14, height=6, res=120)
# par(mfrow=c(1, 2), mar=c(5.1, 5.1, 4.1, 2.1))
# plot(psel4[[1]], auc4[[2]], type="l", xlim=range(psel4), ylim=range(auc4), 
#      col=colors[1], xlab="Number of selected features", ylab="AUC", main="a)", 
#      lwd=1.5, cex.axis=1.5, cex.names=1.5, cex.lab=2, cex.main=2)
# lines(range(psel4), rep(auc4[[1]], 2), col=colors[1], lty=2, lwd=1.5)
# lines(psel4[[2]], auc4[[3]], col=colors[2], lty=2, lwd=1.5)
# lines(psel4[[3]], auc4[[4]], col=colors[3], lty=2, lwd=1.5)
# lines(psel4[[4]], auc4[[5]], col=colors[4], lty=2, lwd=1.5)
# lines(psel4[[5]], auc4[[6]], col=colors[2], lwd=1.5)
# lines(psel4[[6]], auc4[[7]], col=colors[3], lwd=1.5)
# lines(psel4[[7]], auc4[[8]], col=colors[4], lwd=1.5)
# 
# plot(psel4[[1]], briers4[[2]], type="l", xlim=range(psel4), ylim=range(briers4),
#      col=colors[1], xlab="Number of selected features", 
#      ylab="Brier skill score", main="b)", cex.axis=1.5, cex.names=1.5, 
#      cex.lab=2, lwd=1.5, cex.main=2)
# lines(range(psel4), rep(briers4[[1]], 2), col=colors[1], lty=2, lwd=1.5)
# lines(psel4[[2]], briers4[[3]], col=colors[2], lty=2, lwd=1.5)
# lines(psel4[[3]], briers4[[4]], col=colors[3], lty=2, lwd=1.5)
# lines(psel4[[4]], briers4[[5]], col=colors[4], lty=2, lwd=1.5)
# lines(psel4[[5]], briers4[[6]], col=colors[2], lwd=1.5)
# lines(psel4[[6]], briers4[[7]], col=colors[3], lwd=1.5)
# lines(psel4[[7]], briers4[[8]], col=colors[4], lwd=1.5)
# legend("bottomleft", legend=leglabels, fill=c(colors, 0, 0),
#        lty=c(rep(NA, 4), 1, 2), lwd=c(rep(NA, 4), 1.5, 1.5), 
#        border=c(rep(1, 4), 0, 0), merge=TRUE, seg.len=1, cex=1.3)
# dev.off()
# 
# rocs4 <- lapply(results4$pred[-1], function(l) {
#   apply(l, 2, function(preds) {pROC::roc(as.numeric(resp) - 1, preds)})})
# rocs4 <- c(list(pROC::roc(as.numeric(resp) - 1, results4$pred[[1]])), rocs4)
# 
# png(paste(path.graph, "grEBEN_mirseq_Maarten_res4_rocs.png", sep=""),
#     units="in", width=14, height=6, res=120)
# par(mfrow=c(1, 2), mar=c(5.1, 5.1, 4.1, 2.1))
# plot(rocs4[[1]], asp=1, col=colors[1], lty=2, main="a)", cex.axis=1.5, 
#      cex.names=1.5, cex.lab=2, lwd=1.5, cex.main=2)
# plot(rocs4[[3]][[11]], add=TRUE, col=colors[2], lty=2, lwd=1.5)
# plot(rocs4[[4]][[11]], add=TRUE, col=colors[3], lty=2, lwd=1.5)
# plot(rocs4[[5]][[11]], add=TRUE, col=colors[4], lty=2, lwd=1.5)
# 
# plot(rocs4[[2]][[11]], col=colors[1], main="b)", cex.axis=1.5, cex.names=1.5, 
#      cex.lab=2, lwd=1.5, cex.main=2)
# plot(rocs4[[6]][[11]], add=TRUE, col=colors[2], lwd=1.5)
# plot(rocs4[[7]][[11]], add=TRUE, col=colors[3], lwd=1.5)
# plot(rocs4[[8]][[11]], add=TRUE, col=colors[4], lwd=1.5)
# legend("bottomright", legend=leglabels, fill=c(colors, 0, 0),
#        lty=c(rep(NA, 4), 1, 2), lwd=c(rep(NA, 4), 1.5, 1.5), 
#        border=c(rep(1, 4), 0, 0), merge=TRUE, seg.len=1, cex=1.3)
# dev.off()
# 
# 
# # three groups
# load(paste(path.res, "grEBEN_mirseq_Maarten_res5.Rdata", sep=""))
# 
# auc5 <- lapply(results5$pred[-1], function(l) {
#   apply(l, 2, function(preds) {pROC::roc(as.numeric(resp) - 1, preds)$auc})})
# auc5 <- c(pROC::roc(as.numeric(resp) - 1, results5$pred[[1]])$auc, auc5)
# 
# briers5 <- lapply(results5$pred[-1], function(l) {
#   apply(l, 2, function(preds) {
#     1 - sum((as.numeric(resp) - 1 - preds)^2)/
#       sum((as.numeric(resp) - 1 - mean(as.numeric(resp) - 1))^2)})})
# briers5 <- c(1 - sum((as.numeric(resp) - 1 - results5$pred[[1]])^2)/
#                sum((as.numeric(resp) - 1 - mean(as.numeric(resp) - 1))^2),
#              briers5)
# 
# psel5 <- lapply(results5$psel, function(l) {colMeans(l)})
# 
# png(paste(path.graph, "grEBEN_mirseq_Maarten_res5_performance.png", sep=""),
#     units="in", width=14, height=6, res=120)
# par(mfrow=c(1, 2), mar=c(5.1, 5.1, 4.1, 2.1))
# plot(psel5[[1]], auc5[[2]], type="l", xlim=range(psel5), ylim=range(auc5), 
#      col=colors[1], xlab="Number of selected features", ylab="AUC", main="a)", 
#      lwd=1.5, cex.axis=1.5, cex.names=1.5, cex.lab=2, cex.main=2)
# lines(range(psel5), rep(auc5[[1]], 2), col=colors[1], lty=2, lwd=1.5)
# lines(psel5[[2]], auc5[[3]], col=colors[2], lty=2, lwd=1.5)
# lines(psel5[[3]], auc5[[4]], col=colors[3], lty=2, lwd=1.5)
# lines(psel5[[4]], auc5[[5]], col=colors[4], lty=2, lwd=1.5)
# lines(psel5[[5]], auc5[[6]], col=colors[2], lwd=1.5)
# lines(psel5[[6]], auc5[[7]], col=colors[3], lwd=1.5)
# lines(psel5[[7]], auc5[[8]], col=colors[4], lwd=1.5)
# 
# plot(psel5[[1]], briers5[[2]], type="l", xlim=range(psel5), ylim=range(briers5),
#      col=colors[1], xlab="Number of selected features", 
#      ylab="Brier skill score", main="b)", cex.axis=1.5, cex.names=1.5, 
#      cex.lab=2, lwd=1.5, cex.main=2)
# lines(range(psel5), rep(briers5[[1]], 2), col=colors[1], lty=2, lwd=1.5)
# lines(psel5[[2]], briers5[[3]], col=colors[2], lty=2, lwd=1.5)
# lines(psel5[[3]], briers5[[4]], col=colors[3], lty=2, lwd=1.5)
# lines(psel5[[4]], briers5[[5]], col=colors[4], lty=2, lwd=1.5)
# lines(psel5[[5]], briers5[[6]], col=colors[2], lwd=1.5)
# lines(psel5[[6]], briers5[[7]], col=colors[3], lwd=1.5)
# lines(psel5[[7]], briers5[[8]], col=colors[4], lwd=1.5)
# legend("bottomleft", legend=leglabels, fill=c(colors, 0, 0),
#        lty=c(rep(NA, 4), 1, 2), lwd=c(rep(NA, 4), 1.5, 1.5), 
#        border=c(rep(1, 4), 0, 0), merge=TRUE, seg.len=1, cex=1.3)
# dev.off()
# 
# rocs5 <- lapply(results5$pred[-1], function(l) {
#   apply(l, 2, function(preds) {pROC::roc(as.numeric(resp) - 1, preds)})})
# rocs5 <- c(list(pROC::roc(as.numeric(resp) - 1, results5$pred[[1]])), rocs5)
# 
# png(paste(path.graph, "grEBEN_mirseq_Maarten_res5_rocs.png", sep=""),
#     units="in", width=14, height=6, res=120)
# par(mfrow=c(1, 2), mar=c(5.1, 5.1, 4.1, 2.1))
# plot(rocs5[[1]], asp=1, col=colors[1], lty=2, main="a)", cex.axis=1.5, 
#      cex.names=1.5, cex.lab=2, lwd=1.5, cex.main=2)
# plot(rocs5[[3]][[11]], add=TRUE, col=colors[2], lty=2, lwd=1.5)
# plot(rocs5[[4]][[11]], add=TRUE, col=colors[3], lty=2, lwd=1.5)
# plot(rocs5[[5]][[11]], add=TRUE, col=colors[4], lty=2, lwd=1.5)
# 
# plot(rocs5[[2]][[11]], col=colors[1], main="b)", cex.axis=1.5, cex.names=1.5, 
#      cex.lab=2, lwd=1.5, cex.main=2)
# plot(rocs5[[6]][[11]], add=TRUE, col=colors[2], lwd=1.5)
# plot(rocs5[[7]][[11]], add=TRUE, col=colors[3], lwd=1.5)
# plot(rocs5[[8]][[11]], add=TRUE, col=colors[4], lwd=1.5)
# legend("bottomright", legend=leglabels, fill=c(colors, 0, 0),
#        lty=c(rep(NA, 4), 1, 2), lwd=c(rep(NA, 4), 1.5, 1.5), 
#        border=c(rep(1, 4), 0, 0), merge=TRUE, seg.len=1, cex=1.3)
# dev.off()



# combined
# leglabels <- c("ridge", expression(paste("enet, ", alpha==0.05)),
#                expression(paste("enet, ", alpha==0.5)),
#                expression(paste("enet, ", alpha==0.95)),
#                "group-regularized", "not group-regularized")
# 
# load(paste(path.res, "grEBEN_mirseq_Maarten_res4.Rdata", sep=""))
# 
# auc4 <- lapply(results4$pred[-1], function(l) {
#   apply(l, 2, function(preds) {pROC::roc(as.numeric(resp) - 1, preds)$auc})})
# auc4 <- c(pROC::roc(as.numeric(resp) - 1, results4$pred[[1]])$auc, auc4)
# 
# briers4 <- lapply(results4$pred[-1], function(l) {
#   apply(l, 2, function(preds) {
#     1 - sum((as.numeric(resp) - 1 - preds)^2)/
#       sum((as.numeric(resp) - 1 - mean(as.numeric(resp) - 1))^2)})})
# briers4 <- c(1 - sum((as.numeric(resp) - 1 - results4$pred[[1]])^2)/
#                sum((as.numeric(resp) - 1 - mean(as.numeric(resp) - 1))^2),
#              briers4)
# 
# psel4 <- lapply(results4$psel, function(l) {colMeans(l)})
# 
# load(paste(path.res, "grEBEN_mirseq_Maarten_res5.Rdata", sep=""))
# 
# auc5 <- lapply(results5$pred[-1], function(l) {
#   apply(l, 2, function(preds) {pROC::roc(as.numeric(resp) - 1, preds)$auc})})
# auc5 <- c(pROC::roc(as.numeric(resp) - 1, results5$pred[[1]])$auc, auc5)
# 
# briers5 <- lapply(results5$pred[-1], function(l) {
#   apply(l, 2, function(preds) {
#     1 - sum((as.numeric(resp) - 1 - preds)^2)/
#       sum((as.numeric(resp) - 1 - mean(as.numeric(resp) - 1))^2)})})
# briers5 <- c(1 - sum((as.numeric(resp) - 1 - results5$pred[[1]])^2)/
#                sum((as.numeric(resp) - 1 - mean(as.numeric(resp) - 1))^2),
#              briers5)
# 
# psel5 <- lapply(results5$psel, function(l) {colMeans(l)})
# 
# png(paste(path.graph, "grEBEN_mirseq_Maarten_combined_performance.png", sep=""),
#     units="in", width=14, height=12, res=120)
# par(mfrow=c(2, 2), mar=c(5.1, 5.1, 4.1, 2.1))
# plot(psel4[[1]], auc4[[2]], type="l", xlim=range(psel4), ylim=range(auc4), 
#      col=colors[1], xlab="Number of selected features", ylab="AUC", main="a)", 
#      lwd=1.5, cex.axis=1.5, cex.lab=2, cex.main=2)
# lines(range(psel4), rep(auc4[[1]], 2), col=colors[1], lty=2, lwd=1.5)
# lines(psel4[[2]], auc4[[3]], col=colors[2], lty=2, lwd=1.5)
# lines(psel4[[3]], auc4[[4]], col=colors[3], lty=2, lwd=1.5)
# lines(psel4[[4]], auc4[[5]], col=colors[4], lty=2, lwd=1.5)
# lines(psel4[[5]], auc4[[6]], col=colors[2], lwd=1.5)
# lines(psel4[[6]], auc4[[7]], col=colors[3], lwd=1.5)
# lines(psel4[[7]], auc4[[8]], col=colors[4], lwd=1.5)
# 
# plot(psel4[[1]], briers4[[2]], type="l", xlim=range(psel4), ylim=range(briers4),
#      col=colors[1], xlab="Number of selected features", 
#      ylab="Brier skill score", main="b)", cex.axis=1.5,
#      cex.lab=2, lwd=1.5, cex.main=2)
# lines(range(psel4), rep(briers4[[1]], 2), col=colors[1], lty=2, lwd=1.5)
# lines(psel4[[2]], briers4[[3]], col=colors[2], lty=2, lwd=1.5)
# lines(psel4[[3]], briers4[[4]], col=colors[3], lty=2, lwd=1.5)
# lines(psel4[[4]], briers4[[5]], col=colors[4], lty=2, lwd=1.5)
# lines(psel4[[5]], briers4[[6]], col=colors[2], lwd=1.5)
# lines(psel4[[6]], briers4[[7]], col=colors[3], lwd=1.5)
# lines(psel4[[7]], briers4[[8]], col=colors[4], lwd=1.5)
# 
# plot(psel5[[1]], auc5[[2]], type="l", xlim=range(psel5), ylim=range(auc5), 
#      col=colors[1], xlab="Number of selected features", ylab="AUC", main="c)", 
#      lwd=1.5, cex.axis=1.5, cex.lab=2, cex.main=2)
# lines(range(psel5), rep(auc5[[1]], 2), col=colors[1], lty=2, lwd=1.5)
# lines(psel5[[2]], auc5[[3]], col=colors[2], lty=2, lwd=1.5)
# lines(psel5[[3]], auc5[[4]], col=colors[3], lty=2, lwd=1.5)
# lines(psel5[[4]], auc5[[5]], col=colors[4], lty=2, lwd=1.5)
# lines(psel5[[5]], auc5[[6]], col=colors[2], lwd=1.5)
# lines(psel5[[6]], auc5[[7]], col=colors[3], lwd=1.5)
# lines(psel5[[7]], auc5[[8]], col=colors[4], lwd=1.5)
# 
# plot(psel5[[1]], briers5[[2]], type="l", xlim=range(psel5), ylim=range(briers5),
#      col=colors[1], xlab="Number of selected features", 
#      ylab="Brier skill score", main="d)", cex.axis=1.5,
#      cex.lab=2, lwd=1.5, cex.main=2)
# lines(range(psel5), rep(briers5[[1]], 2), col=colors[1], lty=2, lwd=1.5)
# lines(psel5[[2]], briers5[[3]], col=colors[2], lty=2, lwd=1.5)
# lines(psel5[[3]], briers5[[4]], col=colors[3], lty=2, lwd=1.5)
# lines(psel5[[4]], briers5[[5]], col=colors[4], lty=2, lwd=1.5)
# lines(psel5[[5]], briers5[[6]], col=colors[2], lwd=1.5)
# lines(psel5[[6]], briers5[[7]], col=colors[3], lwd=1.5)
# lines(psel5[[7]], briers5[[8]], col=colors[4], lwd=1.5)
# legend("bottomleft", legend=leglabels, fill=c(colors, 0, 0),
#        lty=c(rep(NA, 4), 1, 2), lwd=c(rep(NA, 4), 1.5, 1.5), 
#        border=c(rep(1, 4), 0, 0), merge=TRUE, seg.len=1, cex=1.3)
# dev.off()



### histogram with bootstrap variable selection
load(paste(path.res, "grEBEN_mirseq_Maarten_res7.Rdata", sep=""))
colors <- bpy.colors(6)[-c(1, 6)]

leg.lab <- c("not group-regularized", "group-regularized")

int <- lapply(results7$sel, function(l) {
  sapply(combn(l, 2, simplify=FALSE), function(x) {
    length(intersect(x[[1]], x[[2]])) - 5}, simplify=TRUE)})

png(paste(path.graph, "gren_mirseq_Maarten_overlap.png", sep=""),
    units="in", width=10, height=8, res=120)
par(mar=c(5.1, 5.1, 4.1, 2.1))
layout(matrix(c(1, 0, 1,  3, 2, 3, 2, 0), nrow = 2, ncol = 4))
barplot(rbind(sapply(c(1:max(unlist(int))), function(s) {sum(int[[1]]==s)}),
              sapply(c(1:max(unlist(int))), function(s) {sum(int[[4]]==s)}))/
          choose(50, 2),
        xlab="Number of overlapping features", ylab="Density", beside=TRUE,
        col=colors[c(2:3)], names.arg=c(1:max(unlist(int))), main="a)", 
        cex.axis=1.5, cex.names=1.5, cex.lab=2, cex.main=2, legend.text=leg.lab,
        args.legend=list(x="topright", fill=colors[c(2:3)],
                         border=colors[c(2:3)], cex=1.5))
barplot(rbind(sapply(c(1:max(unlist(int))), function(s) {sum(int[[2]]==s)}),
              sapply(c(1:max(unlist(int))), function(s) {sum(int[[5]]==s)}))/
          choose(50, 2),
        xlab="Number of overlapping features", ylab="Density", beside=TRUE,
        col=colors[c(2:3)], names.arg=c(1:max(unlist(int))), main="b)", 
        cex.axis=1.5, cex.names=1.5, cex.lab=2, cex.main=2)
barplot(rbind(sapply(c(1:max(unlist(int))), function(s) {sum(int[[3]]==s)}),
              sapply(c(1:max(unlist(int))), function(s) {sum(int[[6]]==s)}))/
          choose(50, 2),
        xlab="Number of overlapping features", ylab="Density", beside=TRUE,
        col=colors[c(2:3)], names.arg=c(1:max(unlist(int))), main="c)", 
        cex.axis=1.5, cex.names=1.5, cex.lab=2, cex.main=2)
dev.off()

png(paste(path.graph, "gren_mirseq_Maarten_overlap1.png", sep=""),
    units="in", width=7, height=6, res=120)
par(mar=c(5.1, 5.1, 4.1, 2.1))
barplot(rbind(sapply(c(1:max(unlist(int))), function(s) {sum(int[[1]]==s)}),
              sapply(c(1:max(unlist(int))), function(s) {sum(int[[4]]==s)}))/
          choose(50, 2),
        xlab="Number of overlapping features", ylab="Density", beside=TRUE,
        col=colors[c(2:3)], names.arg=c(1:max(unlist(int))), main="a)", 
        cex.axis=1.5, cex.names=1.5, cex.lab=2, cex.main=2, legend.text=leg.lab,
        args.legend=list(x="topright", fill=colors[c(2:3)],
                         border=colors[c(2:3)], cex=1.5))
dev.off()

png(paste(path.graph, "gren_mirseq_Maarten_overlap2.png", sep=""),
    units="in", width=7, height=6, res=120)
par(mar=c(5.1, 5.1, 4.1, 2.1))
barplot(rbind(sapply(c(1:max(unlist(int))), function(s) {sum(int[[2]]==s)}),
              sapply(c(1:max(unlist(int))), function(s) {sum(int[[5]]==s)}))/
          choose(50, 2),
        xlab="Number of overlapping features", ylab="Density", beside=TRUE,
        col=colors[c(2:3)], names.arg=c(1:max(unlist(int))), main="b)", 
        cex.axis=1.5, cex.names=1.5, cex.lab=2, cex.main=2, legend.text=leg.lab,
        args.legend=list(x="topright", fill=colors[c(2:3)],
                         border=colors[c(2:3)], cex=1.5))
dev.off()

png(paste(path.graph, "gren_mirseq_Maarten_overlap3.png", sep=""),
    units="in", width=7, height=6, res=120)
par(mar=c(5.1, 5.1, 4.1, 2.1))
barplot(rbind(sapply(c(1:max(unlist(int))), function(s) {sum(int[[3]]==s)}),
              sapply(c(1:max(unlist(int))), function(s) {sum(int[[6]]==s)}))/
          choose(50, 2),
        xlab="Number of overlapping features", ylab="Density", beside=TRUE,
        col=colors[c(2:3)], names.arg=c(1:max(unlist(int))), main="c)", 
        cex.axis=1.5, cex.names=1.5, cex.lab=2, cex.main=2, legend.text=leg.lab,
        args.legend=list(x="topright", fill=colors[c(2:3)],
                         border=colors[c(2:3)], cex=1.5))
dev.off()
