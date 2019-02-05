##############################  preamble  #############################
# grEBEN on real data                                                 #
# version: 01                                                         #
# author: Magnus M?nch                                                #
# created: 16-10-2017                                                 #
# last edited: 16-10-2017                                             #
#######################################################################

##########################  data description  #########################
# A deep sequencing analysis on small non-coding ribonucleic acid     #
# (miRNAseq) was performed on 56 samples (24 women with high-grade    #
# cervical intraepithelial neoplasia (CIN3) and 32 healthy women) for #
# the purpose of finding relevant screening markers for cervical      #
# cancer screening. The next generation sequencing analysis resulted  #
# in 2,576 transcripts. The data was normalized and pre-processed,    #
# rendering 772 transcripts. More detail description of the data sets #
# and the preprocessing step are available on the supplementary       #
# material of this following publication "Better diagnostic           #
# signatures from RNAseq data through use of auxiliary co-data",      #
# Bioinformatics (2017).                                              #
#######################################################################

### paths
path.code <- as.character(ifelse(Sys.info()[1]=="Darwin", "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/code/" ,
                                 "~/EBEN/code/"))
path.res <- as.character(ifelse(Sys.info()[1]=="Darwin", "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/results/" ,"~/EBEN/results/"))
path.data <- as.character(ifelse(Sys.info()[1]=="Darwin","/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/data/" ,"~/EBEN/data/"))
path.graph <- "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/graphs/"
path.pack <- as.character("/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/rpackages/")

### libraries
library(glmnet)
library(penalized)
library(GRridge)
library(pROC)

### functions
# source grENVB functions
source(paste(path.code, "grVBEM.R", sep=""))

### loading data
data(dataVerlaat)

### data manipulation
data.mirna <- t(datcenVerlaat)
norm.mirna <- apply(data.mirna, 2, function(x) {(x - mean(x))/sd(x)})
resp.cin3 <- respVerlaat
n <- nrow(norm.mirna)
p <- ncol(norm.mirna)
m <- rep(1, n)

### partitions
annotation <- as.numeric(CpGann)
partitions1 <- list(annotation=CreatePartition(CpGann))
partitions2 <- list(annotation=annotation)
partitions3 <- list(pvalues=CreatePartition(pvalFarkas, decreasing=FALSE, uniform=TRUE, ngroup=10))
partitions4 <- list(pvalues=rep(1:length(partitions3$pvalues), 
                                times=unlist(lapply(partitions3$pvalues, length)))[
                                  order(unlist(partitions3$pvalues))])
partitions5 <- list(annotation=partitions1$annotation,
                    pvalues=partitions3$pvalues)
partitions6 <- list(annotation=partitions2$annotation,
                    pvalues=partitions4$pvalues)

### using annotation as grouping information
# fitting models
enet.pen <- cv.pen(norm.mirna, resp.cin3, unpenalized=NULL, intercept=TRUE, psel=NULL)
enet0.05.pen <- cv.glmnet(norm.mirna, resp.cin3, family="binomial", alpha=0.05,
                          grouped=FALSE, standardize=FALSE, nfolds=length(resp.cin3))
fit.ridge <- cv.glmnet(norm.mirna, resp.cin3, family="binomial", alpha=0, grouped=FALSE,
                       standardize=FALSE, nfolds=length(resp.cin3))
fit.lasso <- cv.glmnet(norm.mirna, resp.cin3, family="binomial", alpha=1, grouped=FALSE,
                       standardize=FALSE, nfolds=length(resp.cin3))
fit.enet <- glmnet(norm.mirna, resp.cin3, family="binomial", standardize=FALSE,
                   alpha=enet.pen$alpha[which.min(enet.pen$cvll)],
                   lambda=enet.pen$lambda[which.min(enet.pen$cvll)])
fit1.grridge <- grridge(t(norm.mirna), resp.cin3, partitions1)
fit1.greben <- grEBEN(norm.mirna, resp.cin3, rep(1, length(resp.cin3)), partitions=partitions2,
                      lambda=0.001)

# using p-values
fit2.grridge <- grridge(t(norm.mirna), resp.cin3, partitions3, monotone=FALSE)
fit2.greben <- grEBEN(norm.mirna, resp.cin3, rep(1, length(resp.cin3)), partitions=partitions4,
                      lambda1=enet.pen$lambda1bayes[which.min(enet.pen$cvll)],
                      lambda2=enet.pen$lambda2bayes[which.min(enet.pen$cvll)],
                      monotone=list(monotone=FALSE, decreasing=FALSE))
fit2.2.greben <- grEBEN(norm.mirna, resp.cin3, rep(1, length(resp.cin3)), partitions=partitions4,
                        lambda1=2*n*0.05*enet0.05.pen$lambda.min,
                        lambda2=n*0.95*enet0.05.pen$lambda.min,
                        monotone=list(monotone=FALSE, decreasing=FALSE))

fit2.2.grridge <- grridge(t(norm.mirna), resp.cin3, partitions3, monotone=TRUE)
fit2.3.greben <- grEBEN(norm.mirna, resp.cin3, rep(1, length(resp.cin3)), partitions=partitions4,
                      lambda1=enet.pen$lambda1bayes[which.min(enet.pen$cvll)],
                      lambda2=enet.pen$lambda2bayes[which.min(enet.pen$cvll)],
                      monotone=list(monotone=TRUE, decreasing=FALSE))

fit2.4.greben <- grEBEN(norm.mirna, resp.cin3, rep(1, length(resp.cin3)), partitions=partitions4,
                        lambda1=2*n*0.05*enet0.05.pen$lambda.min,
                        lambda2=n*0.95*enet0.05.pen$lambda.min,
                        monotone=list(monotone=TRUE, decreasing=FALSE))
fit2.5.greben <- grEBEN(norm.mirna, resp.cin3, rep(1, length(resp.cin3)), partitions=partitions4,
                        lambda1=fit2.grridge$optl/n,
                        lambda2=fit2.grridge$optl/n,
                        monotone=list(monotone=FALSE, decreasing=FALSE))

fit2.grvbem2 <- grVBEM2(norm.mirna, resp.cin3, rep(1, length(resp.cin3)), 
                        partitions=partitions4, intercept=TRUE, posterior=FALSE,
                        lambda1=enet.pen$lambda1bayes[which.min(enet.pen$cvll)], 
                        lambda2=enet.pen$lambda2bayes[which.min(enet.pen$cvll)],   
                        eps=0.001, maxiter=500, trace=TRUE)
fit2.bfgs <- grEBEN(norm.mirna, resp.cin3, rep(1, length(resp.cin3)), partitions=partitions4,
                    lambda1=enet.pen$lambda1bayes[which.min(enet.pen$cvll)],
                    lambda2=enet.pen$lambda2bayes[which.min(enet.pen$cvll)],
                    monotone=list(monotone=FALSE, decreasing=FALSE))


fit2.greben$conv

plot(fit2.greben$lambdag$pvalues[, fit2.greben$nouteriter], fit2.grvbem2$lambdag$pvalues[, fit2.grvbem2$nouteriter + 1])
plot(fit2.greben$lambdag$pvalues[, fit2.greben$nouteriter], fit2.grridge$lambdamults$pvalues)
plot(fit2.grvbem2$lambdag$pvalues[, fit2.grvbem2$nouteriter + 1], fit2.grridge$lambdamults$pvalues)




fit2.2.grridge <- grridge(t(norm.mirna), resp.cin3, partitions3, monotone=TRUE)
fit2.2.greben <- grEBEN(norm.mirna, resp.cin3, rep(1, length(resp.cin3)), partitions=partitions4,
                      lambda1=enet.pen$lambda1bayes[which.min(enet.pen$cvll)],
                      lambda2=enet.pen$lambda2bayes[which.min(enet.pen$cvll)],
                      monotone=list(monotone=TRUE, decreasing=FALSE))


plot(fit2.grridge$lambdamults$pvalues, fit2.greben$lambdag$pvalues[, fit2.greben$nouteriter + 1])
barplot(rbind(fit2.grridge$lambdamults$pvalues, 
              fit2.greben$lambdag$pvalues[, fit2.greben$nouteriter + 1]),
        beside=TRUE)

exp(sum(rle(sort(partitions4$pvalues))$lengths*log(fit2.greben$lambdag$pvalues[, fit2.greben$nouteriter + 1])))
exp(sum(unlist(lapply(partitions3$pvalues, length))*log(fit2.grridge$lambdamults$pvalues)))

plot(partitions4$pvalues, pvalFarkas)
plot(rep(1:10, times=unlist(lapply(partitions3$pvalues, length)))[order(unlist(partitions3$pvalues))], pvalFarkas)





fit3.grridge <- grridge(t(norm.mirna), resp.cin3, partitions5, monotone=c(FALSE, TRUE))
fit3.greben <- grEBEN(norm.mirna, resp.cin3, rep(1, length(resp.cin3)), partitions=partitions6,
                      lambda1=enet.pen$lambda1bayes[which.min(enet.pen$cvll)],
                      lambda2=enet.pen$lambda2bayes[which.min(enet.pen$cvll)],
                      monotone=list(monotone=c(FALSE, TRUE), decreasing=c(FALSE, FALSE)))



grmagn1 <- function(par, lambda2, sizes, sum1) {
  s <- par[1]
  loglambdag <- par[-1]
  magn <- sum((0.5*lambda2*sum1*exp(loglambdag) - 0.5*sizes + s*sizes)^2) + sum(sizes*loglambdag)^2
  return(magn)
}

grmagn3 <- function(par, lambda2, nparts, partsind, partsmat, sizes, G, sum1) {
  
  # s <- rep(par[1:nparts], times=unlist(G))
  # loglambdag <- par[-c(1:nparts)]
  s <- par[1]
  loglambdag <- par[-1]
  
  loglambdamult <- rowSums(sapply(1:nparts, function(part) {
    loglambdag[partsind==part][partsmat[, part]]}))
  
  partsum <- sum((0.5*lambda2*unlist(sapply(1:nparts, function(part) {
    tapply(exp(loglambdamult)*sum1, partsmat[, part], sum)})) + (s - 0.5)*sizes)^2)
  constrsum <- sum(sapply(1:nparts, function(part) {
    sum(loglambdag[partsind==part]*sizes[partsind==part])^2}))
  magn <- partsum + constrsum
  return(magn)
  
}



x=norm.mirna
y=resp.cin3
m=rep(1, length(resp.cin3)) 
partitions=partitions4
intercept=TRUE
posterior=FALSE
lambda1=enet.pen$lambda1bayes[which.min(enet.pen$cvll)]
lambda2=enet.pen$lambda2bayes[which.min(enet.pen$cvll)]  
eps=0.001
maxiter=100
trace=TRUE

grVBEM2 <- function(x, y, m, partitions, lambda1=NULL, lambda2=NULL, intercept, posterior, eps, 
                    maxiter, trace=TRUE) {
  
  if(!is.list(partitions) & is.null(names(partitions))) {
    partitions <- list(partition1=partitions)
  } else if(is.null(names(partitions))) {
    names(partitions) <- paste("partition", 1:length(partitions), sep="")
  }
  
  partnames <- names(partitions)
  
  # assigning fixed (throughout algorithm) variables
  nparts <- length(partitions)
  sizes <- lapply(partitions, function(part) {rle(sort(part))$lengths})
  G <- lapply(partitions, function(part) {length(unique(part))})
  p <- ncol(x)
  n <- nrow(x)
  kappa <- y - m/2
  
  # if no penalty parameters are given we estimate them by cross-validation
  if(is.null(lambda1) | is.null(lambda2)) {
    if(trace) {cat("\r", "Estimating global lambda1 and lambda2 by cross-validation", sep="")}
    srt <- proc.time()[3]
    opt.glob <- cv.pen(x, y, intercept)
    cv.time <- proc.time()[3] - srt
    lambda1 <- opt.glob$lambda1[which.min(opt.glob$cvll)]
    lambda2 <- opt.glob$lambda2[which.min(opt.glob$cvll)]
    if(trace) {cat("\n", "Global lambda1 and lambda2 estimated at ", round(lambda1, 2), " and ", 
                   round(lambda2, 2), " in ", round(cv.time, 2), " seconds", sep="")}
  }
  
  # in the multiplier setting phi does not change
  phi <- 0.25*lambda1^2/lambda2
  
  # starting values for lambdag, lagrange multiplier s and possible approximate hessian B
  lambdagnew <- lambdag <- lambdagold <- lambdagseq <- lapply(G, function(gpart) {rep(1, gpart)})
  # s <- rep(0, times=nparts)
  s <- 0
  
  # starting values for mu and sigma
  fit.pen <- penalized(y, x, unpenalized=formula(ifelse(intercept, "~1", "~0")), 
                       model="logistic", lambda1=0, 
                       lambda2=2*(lambda1 + lambda2), trace=FALSE)
  #muold <- coef(fit.pen, which="all")
  if(intercept) {
    xadj <- cbind(1, x)
  } else {
    xadj <- x
  }
  b0 <- coef(fit.pen, which="all")
  pred0 <- as.numeric(exp(xadj %*% b0)/(1 + exp(xadj %*% b0)))
  w <- sqrt(pred0*(1 - pred0))
  xw <- xadj*w
  svdxw <- svd(xw)
  d <- svdxw$d
  v <- svdxw$v
  invmat <- d^2/(d^2 + 4*(lambda1 + lambda2))^2
  sigmaold <- t(t(v)*invmat) %*% t(v)
  
  # rest of the starting values follow from that
  muold <- as.numeric(sigmaold %*% (t(xadj) %*% as.matrix(kappa)))
  ci <- as.numeric(sqrt(colSums(t(xadj) * (sigmaold %*% t(xadj))) + (colSums(t(xadj)*muold))^2))
  chi <- as.numeric(0.5*(lambda1 + lambda2)*(diag(sigmaold) + muold^2))[(intercept + 1):(intercept + p)]
  dsigmaold <- diag(sigmaold)
  
  # sum is needed in optimisation routine
  intsec <- do.call("paste", c(partitions, sep=" "))
  uintsec <- unique(intsec)
  sum1 <- sapply(1:length(uintsec), function(cursec) {
    ind <- which(intsec==uintsec[[cursec]]) + intercept;
    sum((diag(sigmaold)[ind] + muold[ind]^2)*(1 + sqrt(phi/chi[ind - intercept])))})
  intsizes <- as.numeric(t(table(as.data.frame(matrix(unlist(partitions), ncol=nparts)))))
  partsmat <- unique(matrix(unlist(partitions), ncol=nparts))
  partsind <- rep(1:nparts, times=unlist(G))
  
  # remove intercept for iterations
  if(intercept) {
    muold <- muold[-1]
    dsigmaold <- dsigmaold[-1]
  }
  
  # keeping track of things:
  lambdamult <- exp(rowSums(sapply(1:nparts, function(part) {
    log(unlist(lambdag))[partsind==part][partsmat[, part]]})))
  # lambdagprod <- sapply(1:length(uintsec), function(cursec) {prod(sapply(1:nparts, function(part) {
  #   lambdag[[part]][as.numeric(strsplit(uintsec[cursec], split=" ")[[1]][part])]}))})
  lowermllseq <- 0.5*sum(sapply(1:nparts, function(part) {sum(sizes[[part]]*log(lambdag[[part]]))})) - 
    0.5*lambda2*sum(sum1*lambdamult)
  # lowermllseq <- sum(sizes*log(lambdagold)) - 0.5*lambda2*sum(lambdagold*sum1)
  niter2seq <- numeric(0)
  
  # outer loop of algorithm:
  conv1 <- FALSE
  iter1 <- 0
  if(trace) {cat("\n", "Estimating penalty multipliers by empirical Bayes", "\n", sep="")}
  srt <- proc.time()[3]
  while(!conv1 & (iter1 < maxiter)) {
    
    iter1 <- iter1 + 1
    
    # estimating new lambdag
    if(length(partitions)==1) {
      # local_opts <- list(algorithm="NLOPT_LD_MMA", xtol_rel= 1.0e-7)
      # opts <- list(algorithm="NLOPT_LD_AUGLAG", xtol_rel=1.0e-7, maxeval=1000,
      #              local_opts=local_opts)
      # opt <- nloptr(x0=log(lambdag), eval_f=grlowermll, eval_g_eq=grconstr,
      #               opts=opts, lambda2=lambda2, sizes=as.numeric(sizes), sum1=sum1)
      opt <- optim(par=c(s, log(unlist(lambdag, use.names=FALSE))), fn=grmagn1, lambda2=lambda2, 
                   sizes=sizes[[1]], sum1=sum1, method="BFGS", control=list(maxit=5000))
      s <- opt$par[1]
      lambdagnew <- list(exp(opt$par[-1]))
      lambdagseq <- sapply(1:nparts, function(part) {cbind(lambdagseq[[part]], lambdagnew[[part]])}, simplify=FALSE)
      names(lambdagseq) <- partnames  
    } else {
      # opt <- optim(par=c(s, log(unlist(lambdag, use.names=FALSE))), fn=grmagn2, lambda2=lambda2, nparts=nparts, 
      #              sizes=sizes, G=G, uintsec=uintsec, sum1=sum1, method="Nelder-Mead", control=list(maxit=5000))
      # s <- opt$par[1:nparts]
      # lambdagnew <- split(exp(opt$par[-c(1:nparts)]), factor(rep(partnames, unlist(G)), 
      #                                                        levels=partnames))
      # lambdagseq <- sapply(1:nparts, function(part) {cbind(lambdagseq[[part]], lambdagnew[[part]])}, simplify=FALSE)
      # names(lambdagseq) <- partnames
      
      opt <- optim(par=c(s, log(unlist(lambdag, use.names=FALSE))), fn=grmagn3, lambda2=lambda2, nparts=nparts, 
                   partsind=partsind, partsmat=partsmat, sizes=unlist(sizes), G=G, sum1=sum1, 
                   method="Nelder-Mead", control=list(maxit=5000))
      s <- opt$par[1]
      lambdagnew <- split(exp(opt$par[-1]), factor(rep(partnames, unlist(G)), levels=partnames))
      # s <- opt$par[1:nparts]
      # lambdagnew <- split(exp(opt$par[-c(1:nparts)]), factor(rep(partnames, unlist(G)), levels=partnames))
      lambdagseq <- sapply(1:nparts, function(part) {cbind(lambdagseq[[part]], lambdagnew[[part]])}, simplify=FALSE)
      names(lambdagseq) <- partnames
    }
    
    # inner loop of algorithm:
    conv2 <- 0
    iter2 <- 0
    while(!conv2 & (iter2 < maxiter)) {
      iter2 <- iter2 + 1
      
      # estimating new model parameters
      lambdamultvec <- apply(sapply(1:nparts, function(part) {lambdagnew[[part]][partitions[[part]]]}), 1, prod)
      newparam <- est_param2(x, kappa, m, n, p, ci, phi, chi, lambdamultvec*lambda2, intercept)
      
      dsigma <- as.numeric(newparam$dsigma)
      mu <- as.numeric(newparam$mu)
      ci <- as.numeric(newparam$ci)
      chi <- as.numeric(newparam$chi)
      
      # checking convergence of inner loop
      conv2 <- max(c(abs((mu - muold)/muold), abs((dsigma - dsigmaold)/dsigmaold))) < eps
      
      muold <- mu
      dsigmaold <- dsigma
    }
    
    niter2seq <- c(niter2seq, iter2)
    
    # sum is needed in optimisation routine
    sum1 <- sapply(1:length(uintsec), function(cursec) {
      ind <- which(intsec==uintsec[[cursec]]);
      sum((dsigmaold[ind] + muold[ind]^2)*(1 + sqrt(phi/chi[ind])))})
    
    # keeping track of lower bound on marginal log likelihood
    lambdamult <- exp(rowSums(sapply(1:nparts, function(part) {
      log(unlist(lambdagnew))[partsind==part][partsmat[, part]]})))
    lowermllseq <- c(lowermllseq, 0.5*sum(sapply(1:nparts, function(part) {sum(sizes[[part]]*log(lambdag[[part]]))})) - 
                       0.5*lambda2*sum(sum1*lambdamult))
    # lowermllseq <- c(lowermllseq, sum(sizes*log(lambdag)) - 0.5*lambda2*sum(lambdag*sum1))
    
    # checking convergence of outer loop:
    conv1 <- max(abs((unlist(lambdagnew) - unlist(lambdag))/unlist(lambdag))) < eps
    
    # updating lambdag for new iteration
    lambdagold <- lambdag
    lambdag <- lambdagnew
    
    # printing progress
    if(trace) {cat("\r", "Penalty multipliers estimated at ", paste(partnames, lapply(lambdag, function(part) {
      paste(round(part, 2), collapse=", ")}), sep=": ", collapse=" and "), "      ", sep="")}
    # cat("\r", "Penalty multipliers estimated at ", paste(round(lambdag[-G], 2), collapse=", "), 
    #     " and ", round(lambdag[G], 2), "      ", sep="")
    
  }
  if(posterior) {
    newparam <- est_param(x, kappa, m, n, p, ci, rep(phi, p), chi, lambdamultvec*lambda2, intercept)
    dsigma <- newparam$sigma
    mu <- newparam$mu
    ci <- newparam$ci
    chi <- newparam$chi
  }
  eb.time <- proc.time()[3] - srt
  if(trace) {cat("\r", "Penalty multipliers estimated at ", paste(partnames, lapply(lambdag, function(part) {
    paste(round(part, 2), collapse=", ")}), sep=": ", collapse=" and "), " in ", 
    round(eb.time, 2), " seconds ", sep="")}
  
  out <- list(mu=mu, sigma=dsigma, ci=ci, chi=chi, lambda1=lambda1, lambda2=lambda2, 
              lambdag=lambdagseq, lowermll=lowermllseq, nouteriter=iter1, ninneriter=niter2seq, 
              conv=conv1)
  
  return(out)
  
}





