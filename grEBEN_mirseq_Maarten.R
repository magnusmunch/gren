##################################### notes ####################################
# -results3 and fitted1 are with the 103 subjects                              #
# -results4 and fitted2 are with the 88 subjects                               #
# -results5 is with the 88 subjects and partition in three groups              #
# -results6 is with the 88 subjects and partition in three random groups       #
################################################################################

### paths
path.data <- ifelse(as.character(Sys.info()[1])!="Darwin", "~/EBEN/data/",
                    "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/data/")
path.code <- ifelse(as.character(Sys.info()[1])!="Darwin", "~/EBEN/code/",
                    "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/code/")
path.res <- ifelse(as.character(Sys.info()[1])!="Darwin", "~/EBEN/results/",
                   "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/results/")
path.graph <- "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/graphs/"

### loading libraries
library(GRridge)
library(glmnet)
library(pROC)
library(irr)

### load functions
source(paste(path.code, "grVBEM.R", sep=""))

### load data
# load(paste(path.data, "forMagnus.Rdata", sep=""))
load(paste(path.data, "forMagnusN88.Rdata", sep=""))

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

part.greben2 <- list(TS=miRNA %in% miRNA.BFDR + miRNA %in% miRNA.TumMirs + 1)
part.grridge2 <- list(TS=CreatePartition(as.factor(part.greben2$TS)))

# ## fitting models
set.seed(2018)
fit1.enet <- cv.glmnet(cbind(unpenal, t(mirnormcen_resp)), resp, alpha=0.05,
                       standardize=FALSE, family="binomial", intercept=TRUE,
                       penalty.factor=c(rep(0, u), rep(1, p)))
fit2.enet <- cv.glmnet(cbind(unpenal, t(mirnormcen_resp)), resp, alpha=0.5,
                       standardize=FALSE, family="binomial", intercept=TRUE,
                       penalty.factor=c(rep(0, u), rep(1, p)))
fit3.enet <- cv.glmnet(cbind(unpenal, t(mirnormcen_resp)), resp, alpha=0.95,
                       standardize=FALSE, family="binomial", intercept=TRUE,
                       penalty.factor=c(rep(0, u), rep(1, p)))
fit1.ridge <- cv.glmnet(cbind(unpenal, t(mirnormcen_resp)), resp, alpha=0,
                        standardize=FALSE, family="binomial", intercept=TRUE,
                        penalty.factor=c(rep(0, u), rep(1, p)))
psel <- c(seq(1, 5, 1), seq(7, 15, 2), seq(20, 40, 5), seq(50, 90, 10))
fit1.grridge <- vector("list", length(psel))
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
# 
# fit1.grridge2 <- vector("list", length(psel))
# fit1.grridge2[[1]] <- grridge(mirnormcen_resp, resp, part.grridge2, optl=NULL,
#                               unpenal = ~1 + adjth + thscheme + age + pcrcdiff,
#                               niter=1, method="exact", dataunpen = datfr,
#                               innfold=10, savepredobj="all", comparelasso=FALSE,
#                               optllasso=NULL, compareunpenal=TRUE,
#                               selection=TRUE, maxsel=psel[1])
# for(s in 2:length(psel)) {
#   fit1.grridge2[[s]] <- grridge(mirnormcen_resp, resp, part.grridge2,
#                                 optl=fit1.grridge2[[1]]$optl, unpenal = ~1 +
#                                   adjth + thscheme + age + pcrcdiff, niter=1,
#                                 method="exact", dataunpen = datfr, innfold=10,
#                                 savepredobj="all", comparelasso=FALSE,
#                                 optllasso=NULL, compareunpenal=TRUE,
#                                 selection=TRUE, maxsel=psel[s])
# }
# 
# fit1.greben2 <- grEBEN3(t(mirnormcen_resp), as.numeric(resp) - 1,
#                         rep(1, length(resp)), unpenalized=unpenal,
#                         partitions=part.greben2, alpha=0.05, psel=psel,
#                         nfolds=NULL, posterior=FALSE, ELBO=FALSE, eps=0.001,
#                         maxiter=500, trace=TRUE)
# fit2.greben2 <- grEBEN3(t(mirnormcen_resp), as.numeric(resp) - 1,
#                         rep(1, length(resp)), unpenalized=unpenal,
#                         partitions=part.greben2, alpha=0.5, psel=psel,
#                         nfolds=NULL, posterior=FALSE, ELBO=FALSE, eps=0.001,
#                         maxiter=500, trace=TRUE)
# fit3.greben2 <- grEBEN3(t(mirnormcen_resp), as.numeric(resp) - 1,
#                         rep(1, length(resp)), unpenalized=unpenal,
#                         partitions=part.greben2, alpha=0.95, psel=psel,
#                         nfolds=NULL, posterior=FALSE, ELBO=FALSE, eps=0.001,
#                         maxiter=500, trace=TRUE)
# 
# save(fit1.ridge, fit1.grridge, fit1.greben, fit2.greben, fit3.greben,
#      fit1.grridge2, fit1.greben2, fit2.greben2, fit3.greben2,
#      file=paste(path.res, "grEBEN_mirseq_Maarten_fitted2.Rdata", sep=""))
# 
# mirsel <- t(mirnormcen_resp)[, which(miRNA %in% miRNA.TumMirs)]
# pmirsel <- ncol(mirsel)
# 
# fit1.ridge2 <- cv.glmnet(cbind(unpenal, mirsel), resp, alpha=0,
#                          standardize=FALSE, family="binomial", intercept=TRUE,
#                          penalty.factor=c(rep(0, u), rep(1, pmirsel)))
# fit1.enet <- cv.glmnet(cbind(unpenal, mirsel), resp, alpha=0.05,
#                        standardize=FALSE, family="binomial", intercept=TRUE,
#                        penalty.factor=c(rep(0, u), rep(1, pmirsel)))
# fit2.enet <- cv.glmnet(cbind(unpenal, mirsel), resp, alpha=0.5,
#                        standardize=FALSE, family="binomial", intercept=TRUE,
#                        penalty.factor=c(rep(0, u), rep(1, pmirsel)))
# fit3.enet <- cv.glmnet(cbind(unpenal, mirsel), resp, alpha=0.95,
#                        standardize=FALSE, family="binomial", intercept=TRUE,
#                        penalty.factor=c(rep(0, u), rep(1, pmirsel)))



### cross validation
set.seed(2018)
nfolds <- n
rest <- n %% nfolds
foldid <- sample(rep(1:nfolds, times=round(c(rep(
  n %/% nfolds + as.numeric(rest!=0), times=rest),
  rep(n %/% nfolds, times=nfolds - rest)))))

psel <- c(seq(1, 5, 1), seq(7, 15, 2), seq(20, 40, 5), seq(50, 90, 10))
methods <- c("ridge", "grridge", "enet1", "enet2", "enet3", "greben1",
             "greben2", "greben3")
pred.ridge <- numeric(n)
for(m in 2:length(methods)) {
  assign(paste("pred.", methods[m], sep=""), matrix(NA, nrow=n,
                                                    ncol=length(psel)))
  assign(paste("psel.", methods[m], sep=""), matrix(NA, nrow=n,
                                                    ncol=length(psel)))
}

fitinit.grridge <- grridge(mirnormcen_resp, resp, part.grridge, optl=NULL,
                           unpenal = ~1 + adjth + thscheme + age + pcrcdiff,
                           niter=1, method="exact", dataunpen=datfr,
                           innfold=10, savepredobj="all", comparelasso=FALSE,
                           optllasso=NULL, compareunpenal=TRUE,
                           selection=TRUE, maxsel=psel[1])
fitinit.ridge <- cv.glmnet(cbind(unpenal, t(mirnormcen_resp)), resp, alpha=0,
                           standardize=FALSE, family="binomial", intercept=TRUE,
                           penalty.factor=c(rep(0, u), rep(1, p)))
fitinit1.enet <- cv.glmnet(cbind(unpenal, t(mirnormcen_resp)), resp, alpha=0.05,
                           standardize=FALSE, family="binomial", intercept=TRUE,
                           penalty.factor=c(rep(0, u), rep(1, p)))
fitinit2.enet <- cv.glmnet(cbind(unpenal, t(mirnormcen_resp)), resp, alpha=0.5,
                           standardize=FALSE, family="binomial", intercept=TRUE,
                           penalty.factor=c(rep(0, u), rep(1, p)))
fitinit3.enet <- cv.glmnet(cbind(unpenal, t(mirnormcen_resp)), resp, alpha=0.95,
                           standardize=FALSE, family="binomial", intercept=TRUE,
                           penalty.factor=c(rep(0, u), rep(1, p)))

for(k in sort(unique(foldid))) {
  cat(paste("Fold ", k, "\n"))

  xtrain <- t(mirnormcen_resp)[foldid!=k, ]
  xtest <- matrix(t(mirnormcen_resp)[foldid==k, ], ncol=p, byrow=TRUE)
  ytrain <- (as.numeric(resp) - 1)[foldid!=k]
  utrain1 <- datfr[foldid!=k, ]
  utest1 <- datfr[foldid==k, ]
  utrain2 <- unpenal[foldid!=k, ]
  utest2 <- matrix(unpenal[foldid==k, ], ncol=u, byrow=TRUE)

  cv.ridge <- glmnet(cbind(utrain2, xtrain), ytrain, alpha=0, 
                     standardize=FALSE, family="binomial", intercept=TRUE,
                     penalty.factor=c(rep(0, u), rep(1, p)), 
                     lambda=fitinit.ridge$lambda.min)
  cv.grridge <- vector("list", length(psel))
  cv.grridge[[1]] <- grridge(t(xtrain), ytrain, part.grridge, optl=fitinit.grridge$optl,
                             unpenal = ~1 + adjth + thscheme + age + pcrcdiff,
                             niter=1, method="exact", dataunpen=utrain1,
                             innfold=10, savepredobj="all", comparelasso=FALSE,
                             optllasso=NULL, compareunpenal=TRUE,
                             selection=TRUE, maxsel=psel[1])
  for(s in 2:length(psel)) {
    cv.grridge[[s]] <- grridge(t(xtrain), ytrain, part.grridge,
                               optl=fitinit.grridge$optl,
                               unpenal = ~1 + adjth + thscheme + age + pcrcdiff,
                               niter=1, method="exact", dataunpen=utrain1,
                               innfold=10, savepredobj="all",
                               comparelasso=FALSE, optllasso=NULL,
                               compareunpenal=TRUE, selection=TRUE,
                               maxsel=psel[s])
  }
  cv1.greben <- grEBEN3(xtrain, ytrain, rep(1, length(ytrain)),
                        unpenalized=utrain2, partitions=part.greben, alpha=0.05,
                        psel=psel, trace=TRUE,
                        lambda=fitinit1.enet$lambda.min)
  cv2.greben <- grEBEN3(xtrain, ytrain, rep(1, length(ytrain)),
                        unpenalized=utrain2, partitions=part.greben, alpha=0.5,
                        psel=psel, trace=TRUE,
                        lambda=fitinit2.enet$lambda.min) # 0.251227 1.175005
  cv3.greben <- grEBEN3(xtrain, ytrain, rep(1, length(ytrain)),
                        unpenalized=utrain2, partitions=part.greben, alpha=0.95,
                        psel=psel, trace=TRUE,
                        lambda=fitinit3.enet$lambda.min)

  pred.enet1[foldid==k, ] <- 1/(1 + exp(-cbind(1, utest2, xtest) %*%
                                          cv1.greben$beta.nogroups))
  pred.greben1[foldid==k, ] <- 1/(1 + exp(-cbind(1, utest2, xtest) %*%
                                            cv1.greben$beta))
  pred.enet2[foldid==k, ] <- 1/(1 + exp(-cbind(1, utest2, xtest) %*%
                                          cv2.greben$beta.nogroups))
  pred.greben2[foldid==k, ] <- 1/(1 + exp(-cbind(1, utest2, xtest) %*%
                                            cv2.greben$beta))
  pred.enet3[foldid==k, ] <- 1/(1 + exp(-cbind(1, utest2, xtest) %*%
                                          cv3.greben$beta.nogroups))
  pred.greben3[foldid==k, ] <- 1/(1 + exp(-cbind(1, utest2, xtest) %*%
                                            cv3.greben$beta))
  pred.grridge[foldid==k, ] <- sapply(cv.grridge, function(s) {
    predict.grridge(s, t(xtest), dataunpennew=utest1)[, 3]})
  pred.ridge[foldid==k] <- predict(cv.ridge, cbind(utest2, xtest),
                                   type="response")

  psel.grridge[foldid==k, ] <- sapply(cv.grridge, function(s) {
    length(s$resEN$whichEN)})
  psel.enet1[foldid==k, ] <- colSums(cv1.greben$beta.nogroups!=0) - u - 1
  psel.greben1[foldid==k, ] <- colSums(cv1.greben$beta!=0) - u - 1
  psel.enet2[foldid==k, ] <- colSums(cv2.greben$beta.nogroups!=0) - u - 1
  psel.greben2[foldid==k, ] <- colSums(cv2.greben$beta!=0) - u - 1
  psel.enet3[foldid==k, ] <- colSums(cv3.greben$beta.nogroups!=0) - u - 1
  psel.greben3[foldid==k, ] <- colSums(cv3.greben$beta!=0) - u - 1

}

results4 <- list(pred=list(pred.ridge, pred.grridge, pred.enet1, pred.enet2,
                           pred.enet3, pred.greben1, pred.greben2,
                           pred.greben3),
                 psel=list(psel.grridge, psel.enet1, psel.enet2, psel.enet3,
                           psel.greben1, psel.greben2, psel.greben3))
save(results4, file=paste(path.res, "grEBEN_mirseq_Maarten_res4.Rdata", sep=""))

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
# load(paste(path.res, "grEBEN_mirseq_Maarten_res6.Rdata", sep=""))

# png(paste(path.graph, "grEBEN_mirseq_Maarten_res4_boxplot.png", sep=""),
#     units="in", width=6, height=6, res=120)
# par(mfrow=c(2, 2))
# boxplot(results6$mult[[1]], ylim=range(results6$mult), main="a)", 
#         ylab=expression(paste(lambda[g], "'")), xlab="Random groups")
# boxplot(results6$mult[[2]], ylim=range(results6$mult), main="b)", 
#         ylab=expression(paste(lambda[g], "'")), xlab="Random groups")
# boxplot(results6$mult[[3]], ylim=range(results6$mult), main="c)", 
#         ylab=expression(paste(lambda[g], "'")), xlab="Random groups")
# boxplot(results6$mult[[4]], ylim=range(results6$mult), main="d)", 
#         ylab=expression(paste(lambda[g], "'")), xlab="Random groups")
# dev.off()


# ### fitted models results
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
bsel1.grridge2 <- sapply(fit1.grridge2, function(s) {s$resEN$whichEN})
bsel1.greben2 <- apply(fit1.greben2$beta[-c(1:(u + 1)), ], 2, function(b) {
  which(b!=0)})
bsel2.greben2 <- apply(fit2.greben2$beta[-c(1:(u + 1)), ], 2, function(b) {
  which(b!=0)})
bsel3.greben2 <- apply(fit3.greben2$beta[-c(1:(u + 1)), ], 2, function(b) {
  which(b!=0)})
# 
jaccard <- function(set1, set2) {
  length(intersect(set1, set2))/length(union(set1, set2))
}
jac1.grridge <- sapply(1:length(fit1.grridge), function(s) {
  jaccard(bsel1.grridge2[[s]], bsel1.grridge[[s]])})
jac1.greben <- sapply(1:ncol(fit1.greben$beta), function(s) {
  jaccard(bsel1.greben2[[s]], bsel1.greben[[s]])})
jac2.greben <- sapply(1:ncol(fit2.greben$beta), function(s) {
  jaccard(bsel2.greben2[[s]], bsel2.greben[[s]])})
jac3.greben <- sapply(1:ncol(fit3.greben$beta), function(s) {
  jaccard(bsel3.greben2[[s]], bsel3.greben[[s]])})

png(paste(path.graph, "grEBEN_mirseq_Maarten_jac1.png", sep=""), units="in",
    width=5, height=5, res=200)
plot(psel, jac1.greben, type="l", xlim=range(psel), col=2,
     ylim=range(jac1.greben, jac2.greben, jac3.greben, jac1.grridge),
     ylab="Jaccard index", xlab="Number of selected variables")
lines(psel, jac2.greben, type="l", col=3)
lines(psel, jac3.greben, type="l", col=4)
lines(psel, jac1.grridge, type="l", col=5)
legend("topright", lty=1, col=c(2:5),
       legend=c("EN, alpha=0.05", "EN, alpha=0.5", "EN, alpha=0.95", "GRridge"))
dev.off()

# # barplot of multipliers
# leg.lab <- c("GRridge", expression(paste("greben, ", alpha==0.05)),
#              expression(paste("greben, ", alpha==0.5)),
#              expression(paste("greben, ", alpha==0.95)),
#              "not group-regularized")
# png(paste(path.graph, "grEBEN_mirseq_Maarten_bar1.png", sep=""), units="in",
#     width=12, height=5, res=200)
# par(bg=NA, mfrow=c(1, 2))
# barplot(rbind(fit1.grridge[[1]]$lambdamults$TS,
#               fit1.greben$lambdag$TS[, fit1.greben$nouteriter + 1],
#               fit2.greben$lambdag$TS[, fit2.greben$nouteriter + 1],
#               fit3.greben$lambdag$TS[, fit3.greben$nouteriter + 1]), beside=TRUE,
#         names.arg=c(expression(BFDR<0.05), expression(BFDR>=0.05)),
#         ylab=expression(paste(lambda[g], "'")), main="a)",
#         args.legend=list(x="topleft", fill=c(gray.colors(4), 0),
#                          lty=c(rep(NA, 4), 2), border=c(rep(1, 4), 0),
#                          merge=TRUE, seg.len=1), legend.text=leg.lab)
# abline(h=1, lty=2)
# barplot(rbind(rev(fit1.grridge2[[1]]$lambdamults$TS),
#               rev(fit1.greben2$lambdag$TS[, fit1.greben2$nouteriter + 1]),
#               rev(fit2.greben2$lambdag$TS[, fit2.greben2$nouteriter + 1]),
#               rev(fit3.greben2$lambdag$TS[, fit3.greben2$nouteriter + 1])),
#         beside=TRUE,
#         names.arg=c(expression(BFDR < 0.001),
#                     expression(0.001 <= BFDR ~ phantom(x) < 0.05),
#                     expression(BFDR>=0.05)),
#         ylab=expression(paste(lambda[g], "'")), main="b)")
# abline(h=1, lty=2)
# dev.off()
# 
# part.greben$TS
# vars.benet <- apply(fit1.greben$beta.nogroups, 2, function(b) {
#   c(var(b[part.greben$TS==1]), var(b[part.greben$TS==2]))})
# means.benet <- apply(fit1.greben$beta.nogroups, 2, function(b) {
#   c(mean(b[part.greben$TS==1]), mean(b[part.greben$TS==2]))})
# vars.benet[1, ]/vars.benet[2, ]
# means.benet[1, ]/means.benet[2, ]
# 
# logllen <- function(loglambda, x) {
#   n <- length(x)
#   loglambda1 <- loglambda[1]
#   loglambda2 <- loglambda[2]
#   logll <- 0.5*n*loglambda2 + 
#     n*dnorm(0.5*exp(loglambda1)/sqrt(exp(loglambda2)), log=TRUE) -
#     n*pnorm(-0.5*exp(loglambda1)/sqrt(exp(loglambda2)), log.p=TRUE) - 
#     0.5*exp(loglambda1)*sum(abs(x)) - 0.5*exp(loglambda2)*sum(x^2)
#   return(-logll)
# }
# 
# 
# plot(as.numeric(coef(fit1.ridge)), fit1.greben$beta.nogroups[, 1])
# opt1 <- optim(c(0, 0), logllen, x=as.matrix(coef(fit1.ridge))[part.greben$TS==1, ])
# opt2 <- optim(c(0, 0), logllen, x=as.matrix(coef(fit1.ridge))[part.greben$TS==2, ])
# alpha1 <- exp(opt1$par)[1]/(2*exp(opt1$par)[2] + exp(opt1$par)[1])
# alpha2 <- exp(opt2$par)[1]/(2*exp(opt2$par)[2] + exp(opt2$par)[1])
# lambda1 <- (0.5*exp(opt1$par)[1] + exp(opt1$par)[2])/length(resp)
# lambda2 <- (0.5*exp(opt2$par)[1] + exp(opt2$par)[2])/length(resp)
# 
# exp(opt1$par)
# exp(opt2$par)



# 
# 
### cross-validation results
load(paste(path.res, "grEBEN_mirseq_Maarten_res5.Rdata", sep=""))

auc <- lapply(results5$pred[-1], function(l) {
  apply(l, 2, function(preds) {pROC::roc(as.numeric(resp) - 1, preds)$auc})})
auc <- c(pROC::roc(as.numeric(resp) - 1, results5$pred[[1]])$auc, auc)

briers <- lapply(results5$pred[-1], function(l) {
  apply(l, 2, function(preds) {
    1 - sum((as.numeric(resp) - 1 - preds)^2)/
      sum((as.numeric(resp) - 1 - mean(as.numeric(resp) - 1))^2)})})
briers <- c(1 - sum((as.numeric(resp) - 1 - results5$pred[[1]])^2)/
              sum((as.numeric(resp) - 1 - mean(as.numeric(resp) - 1))^2),
            briers)

psel <- lapply(results5$psel, function(l) {colMeans(l)})

leglabels <- c("ridge", expression(paste("enet, ", alpha==0.05)),
               expression(paste("enet, ", alpha==0.5)),
               expression(paste("enet, ", alpha==0.95)),
               "group-regularized", "not group-regularized")

png(paste(path.graph, "grEBEN_mirseq_Maarten_res5_performance.png", sep=""),
    units="in", width=12, height=6, res=120)
par(mfrow=c(1, 2))
plot(psel[[1]], auc[[2]], type="l", xlim=range(psel), ylim=range(auc), col=2,
     xlab="Number of selected variables", ylab="AUC", main="a)")
lines(range(psel), rep(auc[[1]], 2), col=2, lty=2)
lines(psel[[2]], auc[[3]], col=3, lty=2)
lines(psel[[3]], auc[[4]], col=4, lty=2)
lines(psel[[4]], auc[[5]], col=5, lty=2)
lines(psel[[5]], auc[[6]], col=3)
lines(psel[[6]], auc[[7]], col=4)
lines(psel[[7]], auc[[8]], col=5)

plot(psel[[1]], briers[[2]], type="l", xlim=range(psel), ylim=range(briers),
     col=2, xlab="Number of selected variables", ylab="Brier skill score",
     main="b)")
lines(range(psel), rep(briers[[1]], 2), col=2, lty=2)
lines(psel[[2]], briers[[3]], col=3, lty=2)
lines(psel[[3]], briers[[4]], col=4, lty=2)
lines(psel[[4]], briers[[5]], col=5, lty=2)
lines(psel[[5]], briers[[6]], col=3)
lines(psel[[6]], briers[[7]], col=4)
lines(psel[[7]], briers[[8]], col=5)
legend("bottomleft", legend=leglabels, fill=c(2:5, 0, 0),
       lty=c(rep(NA, 4), 1, 2), border=c(rep(1, 4), 0, 0), merge=TRUE,
       seg.len=1)
dev.off()

rocs <- lapply(results5$pred[-1], function(l) {
  apply(l, 2, function(preds) {pROC::roc(as.numeric(resp) - 1, preds)})})
rocs <- c(list(pROC::roc(as.numeric(resp) - 1, results5$pred[[1]])), rocs)

png(paste(path.graph, "grEBEN_mirseq_Maarten_res5_rocs.png", sep=""),
    units="in", width=12, height=6, res=120)
par(mfrow=c(1, 2))
par(mfrow=c(1, 2))
plot(rocs[[1]], asp=1, col=2, lty=2, main="a)")
plot(rocs[[3]][[11]], add=TRUE, col=3, lty=2)
plot(rocs[[4]][[11]], add=TRUE, col=4, lty=2)
plot(rocs[[5]][[11]], add=TRUE, col=5, lty=2)

plot(rocs[[2]][[11]], col=2, main="b)")
plot(rocs[[6]][[11]], add=TRUE, col=3)
plot(rocs[[7]][[11]], add=TRUE, col=4)
plot(rocs[[8]][[11]], add=TRUE, col=5)
legend("bottomright", legend=leglabels, fill=c(2:5, 0, 0),
       lty=c(rep(NA, 4), 1, 2), border=c(rep(1, 4), 0, 0), merge=TRUE,
       seg.len=1)
dev.off()








load(paste(path.res, "grEBEN_mirseq_Maarten_fitted2.Rdata", sep=""))

best1.enet <- as.numeric(coef(fit1.enet, s="lambda.min"))[-c(1:(u + 1))]
best1.ridge <- as.numeric(coef(fit1.ridge, s="lambda.min"))[-c(1:(u + 1))]

lambda.ridge <- mlenet(best1.ridge)
lambda1.enet <- mlenet(best1.enet)


lambda3.enet <- c(mlenet(best1.enet[part.greben2$TS==1], alpha=0.5)[2],
                  mlenet(best1.enet[part.greben2$TS==2], alpha=0.5)[2],
                  mlenet(best1.enet[part.greben2$TS==3], alpha=0.5)[2])

lambdag.enet <- lambda3.enet/exp(sum(rle(sort(part.greben2$TS))$lengths*log(
  lambda3.enet))/length(part.greben$TS))

lambda2.ridge <- mlenet(best1.ridge, alpha=0.5)
lambda2.enet <- mlenet(best1.enet, alpha=0.5)

hist(best1.ridge, breaks=50, freq=FALSE, ylim=c(0, 130))
curve(denet(x, lambda1=lambda.ridge[1], lambda2=lambda.ridge[2], log=FALSE), 
      add=TRUE, col=2)
curve(denet(x, lambda1=lambda2.ridge[1]*lambda2.ridge[2], 
            lambda2=0.5*(1 - lambda2.ridge[1])*lambda2.ridge[2], log=FALSE), 
      add=TRUE, col=3)

hist(best1.enet, breaks=50, freq=FALSE)
curve(denet(x, lambda1=lambda1.enet[1], lambda2=lambda1.enet[2], log=FALSE), 
      add=TRUE, col=2)
curve(denet(x, lambda1=lambda2.enet[1]*lambda2.enet[2], 
            lambda2=0.5*(1 - lambda2.enet[1])*lambda2.enet[2], log=FALSE), 
      add=TRUE, col=3)



set.seed(2018)
n <- 200
p <- 1000
alpha <- 0.5
lambda <- 100
G <- 4
# lambdag <- exp(c(-1.0625, -0.0625, 0.4375, 0.6875))
# lambdag <- exp(log(seq(0.3, 1.8, length.out=G)) -
#                  mean(log(seq(0.3, 1.8, length.out=G))))
lambdag <- exp(seq(-2, 2, length.out=4))
part <- list(groups=rep(c(1:G), each=p/G))
g <- 4
q <- 0.5
# beta <- as.numeric(sapply(1:G, function(g) {
#   c(rep(0, q*p/G), renbeta(ceiling((1 - q)*p/G), lambda*lambdag[g]*alpha, 
#                                  0.5*(1 - alpha)*lambda*lambdag[g]))}))
beta <- as.numeric(sapply(1:G, function(g) {
  b <- renbeta(p/G, lambda*lambdag[g]*alpha, 0.5*(1 - alpha)*lambda*lambdag[g]);
  b[abs(b)<=quantile(abs(b), q)] <- 0
  return(b)}))
x <- matrix(rnorm(n*p), ncol=p, nrow=n)
pblock <- 25
rho <- 0.7
sigma <- matrix(rho, ncol=pblock, nrow=pblock); diag(sigma) <- 1
x <- do.call(cbind, replicate(p/pblock, rmvnorm(n, mean=rep(0, pblock),
                                                sigma=sigma), simplify=FALSE))
prob <- as.numeric(exp(x %*% beta)/(1 + exp(x %*% beta)))
y <- rbinom(n, 1, prob)

ntest <- 1000
xtest <- do.call(cbind, replicate(p/pblock, rmvnorm(ntest, mean=rep(
  0, pblock), sigma=sigma), simplify=FALSE))
probtest <- as.numeric(exp(xtest %*% beta)/(1 + exp(xtest %*% beta)))
ytest <- rbinom(ntest, 1, probtest)

test1.enet <- glmnet(x, y, alpha=0.5, standardize=FALSE)
test.ridge <- cv.glmnet(x, y, alpha=0, standardize=FALSE)
test1.greben <- grEBEN3(x, y, rep(1, length(y)), partitions=part, alpha=0.5,
                        trace=TRUE)
test1.greben2 <- glmnet(x, y, alpha=0.5, standardize=FALSE, 
                        penalty.factor=rep(test1.greben$lambdag$groups[
                          , test1.greben$nouteriter + 1], each=p/G))
psel <- c(seq(1, 10, 2), seq(15, 50, 5), seq(60, 140, 10), seq(160, 200, 20))
test.grridge <- vector("list", length(psel))
test.grridge[[1]] <- grridge(t(x), y, partitions=list(groups=CreatePartition(
  as.factor(part$groups))), selection=TRUE, maxsel=psel[1])
for(s in 2:length(psel)) {
  test.grridge[[s]] <- grridge(t(x), y, partitions=list(groups=CreatePartition(
    as.factor(part$groups))), selection=TRUE, maxsel=psel[s], 
    optl=test.grridge[[1]]$optl)
}

auc <- list(ridge=pROC::roc(ytest, as.numeric(predict(test.ridge, xtest, 
                                                      "lambda.min")))$auc,
            grridge=sapply(test.grridge, function(s) {
              pROC::roc(ytest, predict.grridge(s, t(xtest))[, 3])$auc}),
            enet1=apply(predict(test1.enet, xtest, type="response"), 2, 
                        function(s) {pROC::roc(ytest, s)$auc}),
            greben1=apply(predict(test1.greben2, xtest, type="response"), 2, 
                          function(s) {pROC::roc(ytest, s)$auc}))
psel <- list(ridge=NA, 
             grridge=sapply(test.grridge, function(s) {
               length(s$resEN$whichEN)}),
             enet1=test1.enet$df, 
             greben1=test1.greben2$df)
kappa <- list(ridge=NA, 
              grridge=sapply(test.grridge, function(s) {
                kappa2(cbind(beta!=0, replace(rep(FALSE, p), s$resEN$whichEN, 
                                              TRUE)))$value}),
              enet1=apply(test1.enet$beta, 2, function(b) {
                kappa2(cbind(beta!=0, b!=0))$value}),
              greben1=apply(test1.greben2$beta, 2, function(b) {
                kappa2(cbind(beta!=0, b!=0))$value}))

plot(sort(psel[[2]]), auc[[2]][order(psel[[2]])], type="l", 
     ylim=range(auc, na.rm=TRUE), xlim=range(psel, na.rm=TRUE), col=2)
lines(range(psel, na.rm=TRUE), rep(auc[[1]], 2), col=2, lty=2)
lines(sort(psel[[3]]), auc[[3]][order(psel[[3]])], col=3, lty=2)
lines(sort(psel[[4]]), auc[[4]][order(psel[[4]])], col=3)

plot(sort(psel[[2]]), kappa[[2]][order(psel[[2]])], type="l", 
     ylim=range(kappa, na.rm=TRUE), xlim=range(psel, na.rm=TRUE), col=2)
lines(sort(psel[[3]]), kappa[[3]][order(psel[[3]])], col=3, lty=2)
lines(sort(psel[[4]]), kappa[[4]][order(psel[[4]])], col=3)




