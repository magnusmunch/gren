##############################  preamble  #############################
# grMCEM implemented in R and C++                                     #
# version: 01                                                         #
# author: Magnus M?nch                                                #
# created: 14-03-2017                                                 #
# last edited: 17-10-2017                                             #
#######################################################################

###############################  notes  ###############################
# 28-03-2017: Added multiple partitions                               #
# 15-03-2017: First implementation of group-regularized elastic net   #
#             with lambda1 penalty multiplier, the square root of     #
#             lambda2 penalty multiplier                              #
#######################################################################

# paths
path.code <- as.character(ifelse(Sys.info()[1]=="Darwin", "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/code/" ,
                                 "~/EBEN/code/"))
path.graph <- "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/graphs/"
path.res <- as.character(ifelse(Sys.info()[1]=="Darwin", "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/results/" ,
                                "~/EBEN/results/"))

### libraries
library(Rcpp)
library(glmnet)

### functions
# source function for variational Bayes
sourceCpp(paste(path.code, "ENVB2.cpp", sep=""))

# various functions of the elastic net prior
dtau <- function(x, lambda1, lambda2, log=FALSE) {
  scale <- 8*lambda2/lambda1^2
  if(log) {
    dens <- -(log(0.5) - 0.5*log(x*(x > 1)) - 0.5*log(scale) - 0.5*log(pi) - x/scale - 
                pnorm(-sqrt(2/scale), log.p=TRUE))
  } else {
    dens <- 0.5*sqrt((x > 1)/(ifelse(x==0, 1, x)*scale*pi))*exp(-x/scale)/pnorm(-sqrt(2/scale))
  }
  return(dens)
}

ptau <- function(q, lambda1, lambda2, log.p=FALSE) {
  scale <- 8*lambda2/lambda1^2
  denom <- pnorm(-sqrt(2/scale))
  if(log.p) {
    prob <- log(1 - pnorm(-sqrt(2*ifelse(q <= 1, 1, q)/scale)/denom))
  } else {
    prob <- 1 - pnorm(-sqrt(2*ifelse(q <= 1, 1, q)/scale)/denom)
  }
  return(prob)
}

qtau <- function(p, lambda1, lambda2, log.p=FALSE) {
  scale <- 8*lambda2/lambda1^2
  if(log.p) {
    aux <- pnorm(-sqrt(2/scale), log.p=log.p)
    tau <- 0.5*scale*qnorm(exp(aux) - exp(aux + p))^2
  } else {
    tau <- 0.5*scale*qnorm((1 - p)*pnorm(-sqrt(2/scale)))^2
  }
  return(tau)
}

rtau <- function(n, lambda1, lambda2, log.p=FALSE) {
  u <- runif(n)
  if(log.p) {
    tau <- qtau(log(u), lambda1, lambda2, log.p=log.p)
  } else {
    tau <- qtau(u, lambda1, lambda2, log.p=log.p)
  }
  return(tau)
}

denet <- function(x, lambda1=1, lambda2=1, log=FALSE) {
  if(log) {
    dens <- 0.5*log(lambda2) - log(2) - 0.5*lambda1*abs(x) - 0.5*lambda2*x^2 +
      dnorm(0.5*lambda1/sqrt(lambda2), log=TRUE) - pnorm(-0.5*lambda1/sqrt(lambda2), log.p=TRUE)
  } else {
    dens <- 0.5*sqrt(lambda2)*exp(-0.5*(lambda1*abs(x) + lambda2*x^2))*
      exp(dnorm(0.5*lambda1/sqrt(lambda2), log=TRUE) - 
            pnorm(-0.5*lambda1/sqrt(lambda2), log.p=TRUE))
  }
  return(dens)
}

renet <- function(n, lambda1=1, lambda2=1, log.p=FALSE) {
  if(lambda1==0) {
    beta <- rnorm(n, 0, sqrt(1/lambda2))
  } else if(lambda2==0) {
    u <- runif(n) - 0.5
    beta <- -2*sign(u)*log(1 - 2*abs(u))/lambda1
  } else {
    tau <- rtau(n, lambda1, lambda2, log.p=log.p)
    sigma <- (tau - 1)/(lambda2*tau)
    beta <- rnorm(n, 0, sqrt(sigma))
  }
  return(beta)
}

varenet <- function(lambda1, lambda2) {
  return(lambda1^2/(4*lambda2^2) + 1/lambda2 - lambda1/(2*lambda2^(3/2))*
           exp(dnorm(lambda1/(2*sqrt(lambda2)), log=TRUE) - 
                 pnorm(-lambda1/(2*sqrt(lambda2)), log.p=TRUE)))
}

renbeta <- function(n, lambda1, lambda2, log.p=FALSE) {
  return(renet(n, lambda1, lambda2, log.p))
}

logllenet <- function(lambda, alpha=NULL, x) {
  n <- length(x)
  if(is.null(alpha)) {
    lambda1 <- lambda[1]
    lambda2 <- lambda[2]
  } else {
    lambda1 <- alpha*lambda
    lambda2 <- 0.5*(1 - alpha)*lambda
  }
  logll <- 0.5*n*log(lambda2) + n*dnorm(0.5*lambda1/sqrt(lambda2), log=TRUE) -
    n*pnorm(-0.5*lambda1/sqrt(lambda2), log.p=TRUE) - 0.5*lambda1*sum(abs(x)) -
    0.5*lambda2*sum(x^2)
  return(-logll)
}

mlenet <- function(x, alpha=NULL) {
  n <- length(x)
  if(is.null(alpha)) {
    lambda.start <- rep((n - 1)/(2*sum(x^2)), 2)
    opt <- constrOptim(lambda.start, logllenet, alpha=alpha, x=x, grad=NULL, 
                       ui=diag(2), ci=c(0, 0))
    out <- opt$par
  } else {
    lambda.start <- (n - 1)/((alpha + 0.5*(1 - alpha))*sum(x^2))
    opt <- optimize(logllenet, c(0, 1e12), alpha=alpha, x=x)
    out <- c(alpha, opt$minimum)
  }
  return(out)
}

# function for estimation
grEBEN <- function(x, y, m, unpenalized=NULL, partitions=NULL, monotone=NULL, intercept=TRUE, 
                   alpha=0.05, lambda=NULL, psel=NULL, nfolds=NULL, recv=TRUE, posterior=FALSE, 
                   ELBO=TRUE, eps=0.001, maxiter=500, trace=TRUE) {
  
  # save argument list
  args <- list(x, y, m, unpenalized, partitions, monotone, intercept, alpha, lambda, psel, 
               nfolds, recv, posterior, ELBO, eps, maxiter, trace)
  names(args) <- formalArgs(grEBEN)
  
  # data characteristics (number of samples, variables etc)
  n <- nrow(x)
  xr <- x
  x <- cbind(unpenalized, xr)
  r <- ncol(xr)
  u <- ifelse(is.null(ncol(unpenalized)), 0, ncol(unpenalized))
  if(is.null(unpenalized)) {
    xu <- matrix(1, nrow=2)
  } else if(intercept){
    xu <- cbind(1, unpenalized)
  } else {
    xu <- unpenalized
  }
  p <- r + u
  ymat <- cbind(m - y, y)
  
  # if no penalty parameter lambda is given we estimate it by cross-validation
  cvll <- nzero <- numeric(0)
  if(is.null(lambda)) {
    if(trace) {cat("\r", "Estimating global lambda by cross-validation", sep="")}
    srt <- proc.time()[3]
    
    nfolds <- ifelse(is.null(nfolds), n, nfolds)
    rest <- n %% nfolds
    foldsize <- c(rep(n %/% nfolds + as.numeric(rest!=0), times=rest),
                  rep(n %/% nfolds, times=nfolds - rest))
    foldid <- sample(rep(1:nfolds, times=foldsize))
    
    cv.fit <- cv.glmnet(x, y, family="binomial", alpha=alpha, standardize=FALSE,
                        intercept=intercept, penalty.factor=c(rep(0, u), rep(1, r)),
                        dfmax=ifelse(is.null(psel), p + 1, psel + u), foldid=foldid,
                        grouped=FALSE)
    
    lambda <- cv.fit$lambda.min
    
    cvll <- c(cvll, cv.fit$cvm[which(cv.fit$lambda==cv.fit$lambda.min)])
    nzero <- c(nzero, cv.fit$nzero[which(cv.fit$lambda==cv.fit$lambda.min)])

    cv.time <- proc.time()[3] - srt
   
    if(trace) {cat("\n", "Global lambda estimated at ", round(lambda, 2), " in ", 
                   round(cv.time, 2), " seconds", sep="")}
  }
  lambda1 <- lambda*alpha*n*2
  lambda2 <- lambda*(1 - alpha)*n
  
  # assigning fixed (throughout algorithm) variables
  if(!is.list(partitions) & is.null(names(partitions))) {
    partitions <- list(partition1=partitions)
  } else if(is.null(names(partitions))) {
    names(partitions) <- paste("partition", 1:length(partitions), sep="")
  }
  partnames <- names(partitions)
  nparts <- length(partitions)
  if(is.null(monotone)) {
    monotone <- list(monotone=rep(FALSE, nparts), decreasing=rep(FALSE, nparts))
  }
  if(is.null(names(monotone))) {
    names(monotone) <- c("monotone", "decreasing")
  }
  sizes <- lapply(partitions, function(part) {rle(sort(part))$lengths})
  G <- lapply(partitions, function(part) {length(unique(part))})
  
  # auxiliary objects used in optimisation of penalty parameter
  intsec <- do.call("paste", c(partitions, sep=" "))
  uintsec <- unique(intsec)
  intsizes <- sapply(uintsec, function(int) {sum(intsec==int)})
  partsmat <- unique(do.call("cbind", partitions))
  partsind <- rep(1:nparts, times=unlist(G))
  
  # starting values for lambdag, lagrange multiplier s
  lambdagnew <- lambdag <- lambdagold <- lambdagseq <- lapply(G, function(gpart) {rep(1, gpart)})
  lambdamultvecold <- rep(1, r)
  s <- 0
  
  # fixed parameters (dont change with VB iterations)
  kappa <- y - m/2
  phi <- 0.25*lambda1^2/lambda2
  
  # starting values for the model parameters
  fit.start <- glmnet(x=x, y=ymat, family="binomial", alpha=0, lambda=lambda, standardize=FALSE, 
                      intercept=intercept, penalty.factor=c(rep(0, u), lambdamultvecold))
  startparam <- est_param3(xr, xu, kappa, m, n, p, as.numeric(predict(fit.start, newx=x, type="response")), 
                           phi, rep(phi, r), lambda2, lambdamultvecold, lambdamultvecold, intercept, 
                           !is.null(unpenalized), FALSE, FALSE, TRUE)
  dsigmaold <- as.numeric(startparam$dsigma)
  muold <- as.numeric(startparam$mu)
  ciold <- as.numeric(startparam$ci)
  chiold <- as.numeric(startparam$chi)
  
  # calculating the sums needed in optimisation routine
  sum1 <- sapply(uintsec, function(int) {
    ind <- which(intsec==int) + intercept;
    sum((dsigmaold[ind + u] + muold[ind + u]^2)*(1 + sqrt(phi/chiold[ind - intercept])))})
  
  # keep track of iterations and possibly evidence lower bound
  if(ELBO) {
    ELBOseq <- numeric(0)
    ELBOold <- Inf
  }
  niter2seq <- numeric(0)
  
  # outer loop of algorithm:
  opt.conv <- logical(0)
  vb.conv <- logical(0)
  conv1 <- FALSE
  iter1 <- 0
  if(trace) {cat("\n", "Estimating penalty multipliers by empirical Bayes", "\n", sep="")}
  srt <- proc.time()[3]
  while(!conv1 & (iter1 < maxiter)) {
    
    iter1 <- iter1 + 1
    
    # estimating new lambdag
    opt <- optim(par=c(s, log(unlist(lambdag, use.names=FALSE))), fn=fopt_groups, lambda2=lambda2,
                 nparts=nparts, partsind=partsind, partsmat=partsmat, sizes=unlist(sizes), G=G,
                 sum1=sum1, method="BFGS", control=list(maxit=1000))
    opt.conv <- c(opt.conv, opt$convergence==0)
    lambdagnew <- split(exp(opt$par[-1]), factor(rep(partnames, unlist(G)), levels=partnames))
    lambdagnew <- sapply(1:nparts, function(part) {
        if(monotone$monotone[part]) {
          return(as.numeric(pava(lambdagnew[[part]], sizes[[part]], monotone$decreasing[part])))
        } else {
          return(lambdagnew[[part]])
        }}, simplify=FALSE)
    if(any(monotone$monotone)) {
      lambdagnew <- sapply(1:nparts, function(part) {exp(log(lambdagnew[[part]]) - sum(sapply(1:nparts, function(part) {
        sum(sizes[[part]]*log(lambdagnew[[part]]))}))/(p*nparts))}, simplify=FALSE)
    }
    s <- opt$par[1]
    lambdagseq <- sapply(1:nparts, function(part) {cbind(lambdagseq[[part]], lambdagnew[[part]])}, simplify=FALSE)
    names(lambdagseq) <- partnames
    lambdamultvec <- apply(sapply(1:nparts, function(part) {lambdagnew[[part]][partitions[[part]]]}), 1, prod)
    
    # inner loop of algorithm:
    conv2 <- FALSE
    iter2 <- 0
    while(!conv2 & (iter2 < maxiter)) {
      iter2 <- iter2 + 1
      
      # estimating new model parameters
      newparam <- est_param3(xr, xu, kappa, m, n, p, ciold, phi, chiold, lambda2, lambdamultvec, lambdamultvecold,
                             intercept, !is.null(unpenalized), FALSE, ELBO, FALSE)
      
      dsigma <- as.numeric(newparam$dsigma)
      mu <- as.numeric(newparam$mu)
      ci <- as.numeric(newparam$ci)
      chi <- as.numeric(newparam$chi)
      if(ELBO) {
        ELBOnew <- newparam$elbo
        ELBOseq <- c(ELBOseq, ELBOnew)
      }
      
      # checking convergence of inner loop
      conv2 <- (max(c(abs((mu - muold)/ifelse(muold==0, muold + 0.00001, muold)), 
                      abs((dsigma - dsigmaold)/ifelse(dsigmaold==0, dsigmaold + 0.00001, dsigmaold)))) < eps)
      muold <- mu
      dsigmaold <- dsigma
      ciold <- ci
      chiold <- chi
      
    }
    vb.conv <- c(vb.conv, conv2)
    
    niter2seq <- c(niter2seq, iter2)
    
    # sum is needed in optimisation routine
    sum1 <- sapply(uintsec, function(int) {
      ind <- which(intsec==int) + intercept;
      sum((dsigmaold[ind + u] + muold[ind + u]^2)*(1 + sqrt(phi/chiold[ind - intercept])))})
    
    # checking convergence of outer loop:
    conv1 <- max(abs((unlist(lambdagnew) - unlist(lambdag))/unlist(lambdag))) < eps
    
    # updating lambdag for new iteration
    lambdagold <- lambdag
    lambdag <- lambdagnew
    lambdamultvecold <- lambdamultvec
    
    # printing progress
    if(trace) {
      cat("\r", "Penalty multipliers estimated at ", paste(partnames, lapply(lambdag, function(part) {
        paste(round(part, 2), collapse=", ")}), sep=": ", collapse=" and "), "      ", sep="")
    }
    
  }
  
  # if a posterior is desired, estimated it here
  if(posterior) {
    newparam <- est_param3(xr, xu, kappa, m, n, p, ci, phi, chi, lambda2, lambdamultvec, lambdamultvecold,
                           intercept, !is.null(unpenalized), TRUE, FALSE, FALSE)
    dsigma <- newparam$sigma
    mu <- as.numeric(newparam$mu)
    ci <- as.numeric(newparam$ci)
    chi <- as.numeric(newparam$chi)
  }
  
  # printing estimation time and final estimates
  eb.time <- proc.time()[3] - srt
  if(trace) {cat("\r", "Penalty multipliers estimated at ", paste(partnames, lapply(lambdag, function(part) {
    paste(round(part, 2), collapse=", ")}), sep=": ", collapse=" and "), " in ", 
    round(eb.time, 2), " seconds ", sep="")
  }
  
  # variable selection
  if(!is.null(psel)) {
    
    # setting the ratio of smallest to largest lambda and length of lambda sequence
    lambda.min.ratio <- ifelse(n < p, 0.01, 0.0001)
    nlambda <- 100
    
    # loop to find the number of variables
    lambdaseq <- NULL
    found <- inbet <- FALSE
    count <- 0
    while(!found) {
      suppressWarnings(fit.df <- glmnet(x, y, family="binomial", alpha=alpha, 
                                        nlambda=nlambda, lambda=lambdaseq, standardize=FALSE, 
                                        intercept=intercept, penalty.factor=c(rep(0, u), lambdamultvec),
                                        dfmax=ifelse(is.null(lambdaseq), psel + u, p + 1)))
      if(any(fit.df$df==psel)) {
        if(tail(fit.df$df, n=1L)==psel & !inbet) {
          
          lambdamax <- tail(fit.df$lambda[fit.df$df==psel], n=1L)
          lambdamin <- lambdamax*lambda.min.ratio
        } else {
          found <- TRUE
        }
      } else if(all(fit.df$df < psel)) {
        lambdamax <- tail(fit.df$lambda, n=1L)
        lambdamin <- lambdamax*lambda.min.ratio
      } else if(all(fit.df$df > psel)) {
        lambdamax <- max(abs(t(x) %*% (y - mean(y)*(1 - mean(y)))))/(n*alpha)
        lambdamin <- max(fit.df$lambda)
        inbet <- TRUE
        count <- count + 1
        if(count==5) {
          found <- TRUE
        }
      } else if(length(unique(lambdaseq))==1) {
        found <- TRUE
      } else {
        ind <- c(rle(fit.df$df < psel)$lengths[1], rle(fit.df$df < psel)$lengths[1] + 1)
        lambdamax <- max(fit.df$lambda[ind])
        lambdamin <- min(fit.df$lambda[ind])
        inbet <- TRUE
      }
      if(!found) {
        lambdaseq <- exp(seq(log(lambdamax), log(lambdamin), length.out=nlambda))
      }    
    }
    ind <- ifelse(any(fit.df$df==psel), tail(which(fit.df$df==psel), n=1L),
                  which.max(fit.df$df[which.min(abs(fit.df$df - psel))]))
    beta.sel <- as.numeric(coef(fit.df)[, ind])
  }
  
  
  # run penalized regression for prediction purposes
  fit.final <- glmnet(x=t(t(x)/c(rep(1, u), sqrt(lambdamultvec))), y=ymat, family="binomial", 
                      alpha=alpha, lambda=lambda, standardize=FALSE, intercept=intercept, 
                      penalty.factor=c(rep(0, u), lambdamultvec))
  beta <- as.numeric(coef(fit.final))
  fit.final.nogroups <- glmnet(x, y, family="binomial", alpha=alpha, lambda=lambda, 
                               standardize=FALSE, penalty.factor=c(rep(0, u), rep(1, r)))
  beta.nogroups <- as.numeric(coef(fit.final.nogroups))
  
  conv <- list(lambda.conv=conv1, opt.conv=opt.conv, vb.conv=vb.conv)
  if(ELBO & !is.null(psel)) {
    out <- list(mu=mu, sigma=dsigma, ci=ci, chi=chi, alpha=alpha, lambda=lambda, 
                lambdag=lambdagseq, ELBO=ELBOseq, nouteriter=iter1, ninneriter=niter2seq, 
                conv=conv, beta=beta, beta.nogroups=beta.nogroups, beta.sel=beta.sel, args=args)
  } else if(ELBO & is.null(psel)) {
    out <- list(mu=mu, sigma=dsigma, ci=ci, chi=chi, alpha=alpha, lambda=lambda, 
                lambdag=lambdagseq, ELBO=ELBOseq, nouteriter=iter1, ninneriter=niter2seq, 
                conv=conv, beta=beta, beta.nogroups=beta.nogroups, args=args)
  } else if(!ELBO & !is.null(psel)) {
    out <- list(mu=mu, sigma=dsigma, ci=ci, chi=chi, alpha=alpha, lambda=lambda, 
                lambdag=lambdagseq, nouteriter=iter1, ninneriter=niter2seq, 
                conv=conv, beta=beta, beta.nogroups=beta.nogroups, beta.sel=beta.sel, args=args)
  } else {
    out <- list(mu=mu, sigma=dsigma, ci=ci, chi=chi, alpha=alpha, lambda=lambda, 
                lambdag=lambdagseq, nouteriter=iter1, ninneriter=niter2seq, 
                conv=conv, beta=beta, beta.nogroups=beta.nogroups, args=args)
  }
  return(out)
  
}

# function to cross-validate the penalty parameters
cv.pen <- function(xr, y, unpenalized=NULL, intercept, psel=NULL, nfolds=NULL) {
  p <- ncol(xr)
  n <- nrow(xr)
  nfolds <- ifelse(is.null(nfolds), 10, nfolds)
  rest <- n %% nfolds
  foldsize <- c(rep(n %/% nfolds + as.numeric(rest!=0), times=rest),
                rep(n %/% nfolds, times=nfolds - rest))
  foldid <- sample(rep(1:nfolds, times=foldsize))
  r <- ncol(xr)
  u <- ifelse(is.null(ncol(unpenalized)), 0, ncol(unpenalized))
  x <- cbind(unpenalized, xr)
  seq.alpha <- seq(0.01, 0.99, length.out=50)
  seq.lam <- seq.df <- seq.cvll <- numeric(length(seq.alpha))
  for(a in 1:length(seq.alpha)) {
    cv.fit <- cv.glmnet(x, y, family="binomial", alpha=seq.alpha[a], standardize=FALSE,
                        intercept=intercept, penalty.factor=c(rep(0, u), rep(1, r)),
                        dfmax=ifelse(is.null(psel), p + 1, psel + u), foldid=foldid,
                        grouped=FALSE)
    ind <- which(cv.fit$lambda==cv.fit$lambda.min)
    seq.lam[a] <- cv.fit$lambda.min
    seq.df[a] <- cv.fit$nzero[ind]
    seq.cvll[a] <- cv.fit$cvm[ind]
  }
  lambda1bayes <- seq.lam*seq.alpha*n*2
  lambda2bayes <- seq.lam*(1 - seq.alpha)*n
  lambda1pen <- seq.lam*seq.alpha*n
  lambda2pen <- seq.lam*(1 - seq.alpha)*0.5*n
  
  out <- list(lambda1bayes=lambda1bayes, lambda2bayes=lambda2bayes, lambda1pen=lambda1pen,
              lambda2pen=lambda2pen, alpha=seq.alpha, lambda=seq.lam, cvll=seq.cvll)
  return(out)
  
}

grEBEN2 <- function(x, y, m, unpenalized=NULL, partitions=NULL, monotone=NULL, intercept=TRUE, 
                   alpha=0.05, lambda=NULL, psel=NULL, nfolds=NULL, recv=TRUE, posterior=FALSE, 
                   ELBO=TRUE, eps=0.001, maxiter=500, trace=TRUE) {
  
  # save argument list
  args <- list(x, y, m, unpenalized, partitions, monotone, intercept, alpha, lambda, psel, 
               nfolds, recv, posterior, ELBO, eps, maxiter, trace)
  names(args) <- formalArgs(grEBEN)
  
  # data characteristics (number of samples, variables etc)
  n <- nrow(x)
  xr <- x
  x <- cbind(unpenalized, xr)
  r <- ncol(xr)
  u <- ifelse(is.null(ncol(unpenalized)), 0, ncol(unpenalized))
  if(is.null(unpenalized)) {
    xu <- matrix(1, nrow=2)
  } else if(intercept){
    xu <- cbind(1, unpenalized)
  } else {
    xu <- unpenalized
  }
  p <- r + u
  ymat <- cbind(m - y, y)
  
  # if no penalty parameter lambda is given we estimate it by cross-validation
  nfolds <- ifelse(is.null(nfolds), n, nfolds)
  rest <- n %% nfolds
  foldsize <- c(rep(n %/% nfolds + as.numeric(rest!=0), times=rest),
                rep(n %/% nfolds, times=nfolds - rest))
  foldid <- sample(rep(1:nfolds, times=foldsize))
  cvll <- nzero <- numeric(0)
  if(is.null(lambda)) {
    if(trace) {cat("\r", "Estimating global lambda by cross-validation", sep="")}
    srt <- proc.time()[3]
    
    cv.fit <- cv.glmnet(x, y, family="binomial", alpha=alpha, standardize=FALSE,
                        intercept=intercept, penalty.factor=c(rep(0, u), rep(1, r)),
                        dfmax=ifelse(is.null(psel), p + 1, psel + u), foldid=foldid,
                        grouped=FALSE)
    
    lambda <- cv.fit$lambda.min
    
    cvll <- c(cvll, cv.fit$cvm[which(cv.fit$lambda==cv.fit$lambda.min)])
    nzero <- c(nzero, cv.fit$nzero[which(cv.fit$lambda==cv.fit$lambda.min)])
    
    cv.time <- proc.time()[3] - srt
    
    if(trace) {cat("\n", "Global lambda estimated at ", round(lambda, 2), " in ", 
                   round(cv.time, 2), " seconds", sep="")}
  }
  lambda1 <- lambda*alpha*n*2
  lambda2 <- lambda*(1 - alpha)*n
  
  # assigning fixed (throughout algorithm) variables
  if(!is.list(partitions) & is.null(names(partitions))) {
    partitions <- list(partition1=partitions)
  } else if(is.null(names(partitions))) {
    names(partitions) <- paste("partition", 1:length(partitions), sep="")
  }
  partnames <- names(partitions)
  nparts <- length(partitions)
  if(is.null(monotone)) {
    monotone <- list(monotone=rep(FALSE, nparts), decreasing=rep(FALSE, nparts))
  }
  if(is.null(names(monotone))) {
    names(monotone) <- c("monotone", "decreasing")
  }
  sizes <- lapply(partitions, function(part) {rle(sort(part))$lengths})
  G <- lapply(partitions, function(part) {length(unique(part))})
  
  
  # starting values for lambdag, lagrange multiplier s
  lambdagnew <- lambdag <- lambdagold <- lambdagseq <- lapply(G, function(gpart) {rep(1, gpart)})
  lambdamultvecold <- rep(1, r)
  s <- 0
  
  # fixed parameters (dont change with VB iterations)
  kappa <- y - m/2
  phi <- 0.25*lambda1^2/lambda2
  
  # starting values for the model parameters
  fit.start <- glmnet(x=x, y=ymat, family="binomial", alpha=0, lambda=lambda, standardize=FALSE, 
                      intercept=intercept, penalty.factor=c(rep(0, u), rep(1, r)))
  startparam <- est_param3(xr, xu, kappa, m, n, p, as.numeric(predict(fit.start, newx=x, type="response")), 
                           phi, rep(phi, r), lambda2, lambdamultvecold, lambdamultvecold, intercept, 
                           !is.null(unpenalized), FALSE, FALSE, TRUE)
  dsigmaold <- as.numeric(startparam$dsigma)
  muold <- as.numeric(startparam$mu)
  ciold <- as.numeric(startparam$ci)
  chiold <- as.numeric(startparam$chi)
  
  # calculating the sums needed in optimisation routine
  intsec <- do.call("paste", c(partitions, sep=" "))
  uintsec <- unique(intsec)
  intsizes <- sapply(uintsec, function(int) {sum(intsec==int)})
  sum1 <- sapply(uintsec, function(int) {
    ind <- which(intsec==int) + intercept;
    sum((dsigmaold[ind + u] + muold[ind + u]^2)*(1 + sqrt(phi/chiold[ind - intercept])))})
  partsmat <- unique(do.call("cbind", partitions))
  partsind <- rep(1:nparts, times=unlist(G))
  
  # keep track of iterations and possibly evidence lower bound
  if(ELBO) {
    ELBOseq <- numeric(0)
    ELBOold <- Inf
  }
  niter2seq <- numeric(0)
  
  # outer loop of algorithm:
  opt.conv <- logical(0)
  vb.conv <- logical(0)
  conv1 <- FALSE
  iter1 <- 0
  if(trace) {cat("\n", "Estimating penalty multipliers by empirical Bayes", "\n", sep="")}
  srt <- proc.time()[3]
  lambdaseq <- lambda
  while(!conv1 & (iter1 < maxiter)) {
    
    iter1 <- iter1 + 1
    
    # estimating new lambdag
    opt <- optim(par=c(s, log(unlist(lambdag, use.names=FALSE))), fn=fopt_groups, lambda2=lambda2,
                 nparts=nparts, partsind=partsind, partsmat=partsmat, sizes=unlist(sizes), G=G,
                 sum1=sum1, method="BFGS", control=list(maxit=1000))
    opt.conv <- c(opt.conv, opt$convergence==0)
    lambdagnew <- split(exp(opt$par[-1]), factor(rep(partnames, unlist(G)), levels=partnames))
    lambdagnew <- sapply(1:nparts, function(part) {
      if(monotone$monotone[part]) {
        return(as.numeric(pava(lambdagnew[[part]], sizes[[part]], monotone$decreasing[part])))
      } else {
        return(lambdagnew[[part]])
      }}, simplify=FALSE)
    if(any(monotone$monotone)) {
      lambdagnew <- sapply(1:nparts, function(part) {exp(log(lambdagnew[[part]]) - sum(sapply(1:nparts, function(part) {
        sum(sizes[[part]]*log(lambdagnew[[part]]))}))/(p*nparts))}, simplify=FALSE)
    }
    s <- opt$par[1]
    lambdagseq <- sapply(1:nparts, function(part) {cbind(lambdagseq[[part]], lambdagnew[[part]])}, simplify=FALSE)
    names(lambdagseq) <- partnames
    lambdamultvec <- apply(sapply(1:nparts, function(part) {lambdagnew[[part]][partitions[[part]]]}), 1, prod)
    
    # re-cross-validate lambda
    fit.lambda <- cv.glmnet(x=x, y=ymat, family="binomial", alpha=alpha, 
                            standardize=FALSE, intercept=intercept, 
                            penalty.factor=c(rep(0, u), lambdamultvec), 
                            foldid=foldid, grouped=FALSE)
    lambda <- fit.lambda$lambda.min
    lambdaseq <- c(lambdaseq, lambda)
    lambda1 <- lambda*alpha*n*2
    lambda2 <- lambda*(1 - alpha)*n
    
    # inner loop of algorithm:
    conv2 <- FALSE
    iter2 <- 0
    while(!conv2 & (iter2 < maxiter)) {
      iter2 <- iter2 + 1
      
      # estimating new model parameters
      newparam <- est_param3(xr, xu, kappa, m, n, p, ciold, phi, chiold, lambda2, lambdamultvec, lambdamultvecold,
                             intercept, !is.null(unpenalized), FALSE, ELBO, FALSE)
      
      dsigma <- as.numeric(newparam$dsigma)
      mu <- as.numeric(newparam$mu)
      ci <- as.numeric(newparam$ci)
      chi <- as.numeric(newparam$chi)
      if(ELBO) {
        ELBOnew <- newparam$elbo
        ELBOseq <- c(ELBOseq, ELBOnew)
      }
      
      # checking convergence of inner loop
      conv2 <- (max(c(abs((mu - muold)/ifelse(muold==0, muold + 0.00001, muold)), 
                      abs((dsigma - dsigmaold)/ifelse(dsigmaold==0, dsigmaold + 0.00001, dsigmaold)))) < eps)
      muold <- mu
      dsigmaold <- dsigma
      ciold <- ci
      chiold <- chi
      
    }
    vb.conv <- c(vb.conv, conv2)
    
    niter2seq <- c(niter2seq, iter2)
    
    # sum is needed in optimisation routine
    sum1 <- sapply(uintsec, function(int) {
      ind <- which(intsec==int) + intercept;
      sum((dsigmaold[ind + u] + muold[ind + u]^2)*(1 + sqrt(phi/chiold[ind - intercept])))})
    
    # checking convergence of outer loop:
    conv1 <- max(abs((unlist(lambdagnew) - unlist(lambdag))/unlist(lambdag))) < eps
    
    # updating lambdag for new iteration
    lambdagold <- lambdag
    lambdag <- lambdagnew
    lambdamultvecold <- lambdamultvec
    
    # printing progress
    if(trace) {
      cat("\r", "Penalty multipliers estimated at ", paste(partnames, lapply(lambdag, function(part) {
        paste(round(part, 2), collapse=", ")}), sep=": ", collapse=" and "), "      ", sep="")
    }
    
  }
  
  # if a posterior is desired, estimated it here
  if(posterior) {
    newparam <- est_param3(xr, xu, kappa, m, n, p, ci, phi, chi, lambda2, lambdamultvec, lambdamultvecold,
                           intercept, !is.null(unpenalized), TRUE, FALSE, FALSE)
    dsigma <- newparam$sigma
    mu <- as.numeric(newparam$mu)
    ci <- as.numeric(newparam$ci)
    chi <- as.numeric(newparam$chi)
  }
  
  # printing estimation time and final estimates
  eb.time <- proc.time()[3] - srt
  if(trace) {cat("\r", "Penalty multipliers estimated at ", paste(partnames, lapply(lambdag, function(part) {
    paste(round(part, 2), collapse=", ")}), sep=": ", collapse=" and "), " in ", 
    round(eb.time, 2), " seconds ", sep="")
  }
  
  # variable selection
  if(!is.null(psel)) {
    
    # setting the ratio of smallest to largest lambda and length of lambda sequence
    lambda.min.ratio <- ifelse(n < p, 0.01, 0.0001)
    nlambda <- 100
    
    # loop to find the number of variables
    lambdaseq <- NULL
    found <- inbet <- FALSE
    count <- 0
    while(!found) {
      suppressWarnings(fit.df <- glmnet(x, y, family="binomial", alpha=alpha, 
                                        nlambda=nlambda, lambda=lambdaseq, standardize=FALSE, 
                                        intercept=intercept, penalty.factor=c(rep(0, u), lambdamultvec),
                                        dfmax=ifelse(is.null(lambdaseq), psel + u, p + 1)))
      if(any(fit.df$df==psel)) {
        if(tail(fit.df$df, n=1L)==psel & !inbet) {
          
          lambdamax <- tail(fit.df$lambda[fit.df$df==psel], n=1L)
          lambdamin <- lambdamax*lambda.min.ratio
        } else {
          found <- TRUE
        }
      } else if(all(fit.df$df < psel)) {
        lambdamax <- tail(fit.df$lambda, n=1L)
        lambdamin <- lambdamax*lambda.min.ratio
      } else if(all(fit.df$df > psel)) {
        lambdamax <- max(abs(t(x) %*% (y - mean(y)*(1 - mean(y)))))/(n*alpha)
        lambdamin <- max(fit.df$lambda)
        inbet <- TRUE
        count <- count + 1
        if(count==5) {
          found <- TRUE
        }
      } else if(length(unique(lambdaseq))==1) {
        found <- TRUE
      } else {
        ind <- c(rle(fit.df$df < psel)$lengths[1], rle(fit.df$df < psel)$lengths[1] + 1)
        lambdamax <- max(fit.df$lambda[ind])
        lambdamin <- min(fit.df$lambda[ind])
        inbet <- TRUE
      }
      if(!found) {
        lambdaseq <- exp(seq(log(lambdamax), log(lambdamin), length.out=nlambda))
      }    
    }
    ind <- ifelse(any(fit.df$df==psel), tail(which(fit.df$df==psel), n=1L),
                  which.max(fit.df$df[which.min(abs(fit.df$df - psel))]))
    beta.sel <- as.numeric(coef(fit.df)[, ind])
  }
  
  
  # run penalized regression for prediction purposes
  fit.final <- glmnet(x=t(t(x)/c(rep(1, u), sqrt(lambdamultvec))), y=ymat, family="binomial", 
                      alpha=alpha, lambda=lambda, standardize=FALSE, intercept=intercept, 
                      penalty.factor=c(rep(0, u), lambdamultvec))
  beta <- as.numeric(coef(fit.final))
  fit.final.nogroups <- glmnet(x, y, family="binomial", alpha=alpha, lambda=lambda, 
                               standardize=FALSE, penalty.factor=c(rep(0, u), rep(1, r)))
  beta.nogroups <- as.numeric(coef(fit.final.nogroups))
  
  conv <- list(lambda.conv=conv1, opt.conv=opt.conv, vb.conv=vb.conv)
  if(ELBO & !is.null(psel)) {
    out <- list(mu=mu, sigma=dsigma, ci=ci, chi=chi, alpha=alpha, lambda=lambdaseq, 
                lambdag=lambdagseq, ELBO=ELBOseq, nouteriter=iter1, ninneriter=niter2seq, 
                conv=conv, beta=beta, beta.nogroups=beta.nogroups, beta.sel=beta.sel, args=args)
  } else if(ELBO & is.null(psel)) {
    out <- list(mu=mu, sigma=dsigma, ci=ci, chi=chi, alpha=alpha, lambda=lambdaseq, 
                lambdag=lambdagseq, ELBO=ELBOseq, nouteriter=iter1, ninneriter=niter2seq, 
                conv=conv, beta=beta, beta.nogroups=beta.nogroups, args=args)
  } else if(!ELBO & !is.null(psel)) {
    out <- list(mu=mu, sigma=dsigma, ci=ci, chi=chi, alpha=alpha, lambda=lambdaseq, 
                lambdag=lambdagseq, nouteriter=iter1, ninneriter=niter2seq, 
                conv=conv, beta=beta, beta.nogroups=beta.nogroups, beta.sel=beta.sel, args=args)
  } else {
    out <- list(mu=mu, sigma=dsigma, ci=ci, chi=chi, alpha=alpha, lambda=lambdaseq, 
                lambdag=lambdagseq, nouteriter=iter1, ninneriter=niter2seq, 
                conv=conv, beta=beta, beta.nogroups=beta.nogroups, args=args)
  }
  return(out)
  
}


# x=xtrain; y=ytrain; m=rep(1, length(ytrain)); unpenalized=utrain2;
# partitions=part.greben; monotone=NULL; intercept=TRUE;
# alpha=0.95; lambda=NULL; psel=psel; nfolds=10; posterior=FALSE;
# ELBO=FALSE; eps=0.001; maxiter=500; trace=TRUE
grEBEN3 <- function(x, y, m, unpenalized=NULL, partitions=NULL, monotone=NULL, intercept=TRUE, 
                    alpha=0.05, lambda=NULL, psel=FALSE, nfolds=NULL, posterior=FALSE, 
                    ELBO=FALSE, eps=0.001, maxiter=500, trace=TRUE) {
  
  # save argument list
  args <- list(x, y, m, unpenalized, partitions, monotone, intercept, alpha, 
               lambda, psel, nfolds, posterior, ELBO, eps, maxiter, trace)
  names(args) <- formalArgs(grEBEN3)
  
  # data characteristics (number of samples, variables etc)
  n <- nrow(x)
  xr <- x
  x <- cbind(unpenalized, xr)
  r <- ncol(xr)
  u <- ifelse(is.null(ncol(unpenalized)), 0, ncol(unpenalized))
  if(is.null(unpenalized)) {
    xu <- matrix(1, nrow=2)
  } else if(intercept){
    xu <- cbind(1, unpenalized)
  } else {
    xu <- unpenalized
  }
  p <- r + u
  ymat <- cbind(m - y, y)
  
  # if no penalty parameter lambda is given we estimate it by cross-validation
  nfolds <- ifelse(is.null(nfolds), n, nfolds)
  rest <- n %% nfolds
  foldsize <- c(rep(n %/% nfolds + as.numeric(rest!=0), times=rest),
                rep(n %/% nfolds, times=nfolds - rest))
  foldid <- sample(rep(1:nfolds, times=foldsize))
  cvll <- nzero <- numeric(0)
  if(is.null(lambda)) {
    if(trace) {cat("\r", "Estimating global lambda by cross-validation", sep="")}
    srt <- proc.time()[3]
    
    cv.fit <- cv.glmnet(x, y, family="binomial", alpha=alpha, standardize=FALSE,
                        intercept=intercept, penalty.factor=c(rep(0, u), rep(1, r)),
                        foldid=foldid, grouped=FALSE)
    
    lambda <- cv.fit$lambda.min
    
    cvll <- c(cvll, cv.fit$cvm[which(cv.fit$lambda==cv.fit$lambda.min)])
    nzero <- c(nzero, cv.fit$nzero[which(cv.fit$lambda==cv.fit$lambda.min)])
    
    cv.time <- proc.time()[3] - srt
    
    if(trace) {cat("\n", "Global lambda estimated at ", round(lambda, 2), " in ", 
                   round(cv.time, 2), " seconds", sep="")}
  }
  lambda1 <- lambda*alpha*n*2
  lambda2 <- lambda*(1 - alpha)*n
  
  # assigning fixed (throughout algorithm) variables
  if(!is.list(partitions) & is.null(names(partitions))) {
    partitions <- list(partition1=partitions)
  } else if(is.null(names(partitions))) {
    names(partitions) <- paste("partition", 1:length(partitions), sep="")
  }
  partnames <- names(partitions)
  nparts <- length(partitions)
  if(is.null(monotone)) {
    monotone <- list(monotone=rep(FALSE, nparts), decreasing=rep(FALSE, nparts))
  }
  if(is.null(names(monotone))) {
    names(monotone) <- c("monotone", "decreasing")
  }
  sizes <- lapply(partitions, function(part) {rle(sort(part))$lengths})
  G <- lapply(partitions, function(part) {length(unique(part))})
  
  # objects used in optimisation of penalty parameters
  intsec <- do.call("paste", c(partitions, sep=" "))
  uintsec <- unique(intsec)
  intsizes <- sapply(uintsec, function(int) {sum(intsec==int)})
  partsmat <- unique(do.call("cbind", partitions))
  partsind <- rep(1:nparts, times=unlist(G))
  
  # starting values for lambdag, lagrange multiplier s
  lambdagnew <- lambdag <- lambdagold <- lambdagseq <- 
    lapply(G, function(gpart) {rep(1, gpart)})
  lambdamultvecold <- rep(1, r)
  s <- 0
  
  # fixed parameters (dont change with VB iterations)
  kappa <- y - m/2
  phi <- 0.25*lambda1^2/lambda2
  
  # starting values for the model parameters
  fit.start <- glmnet(x=x, y=ymat, family="binomial", alpha=0, lambda=lambda, 
                      standardize=FALSE, intercept=intercept, 
                      penalty.factor=c(rep(0, u), lambdamultvecold))
  startparam <- est_param3(xr, xu, kappa, m, n, p, as.numeric(predict(fit.start, newx=x, type="response")), 
                           phi, rep(phi, r), lambda2, lambdamultvecold, lambdamultvecold, intercept, 
                           !is.null(unpenalized), FALSE, FALSE, TRUE)
  dsigmaold <- as.numeric(startparam$dsigma)
  muold <- as.numeric(startparam$mu)
  ciold <- as.numeric(startparam$ci)
  chiold <- as.numeric(startparam$chi)
  
  # calculating the sums needed in optimisation routine
  sum1 <- sapply(uintsec, function(int) {
    ind <- which(intsec==int) + intercept;
    sum((dsigmaold[ind + u] + muold[ind + u]^2)*(1 + sqrt(phi/chiold[ind - intercept])))})
  
  # keep track of iterations and possibly evidence lower bound
  if(ELBO) {
    ELBOseq <- numeric(0)
    ELBOold <- Inf
  }
  niter2seq <- numeric(0)
  
  # outer loop of algorithm:
  opt.conv <- logical(0)
  vb.conv <- logical(0)
  conv1 <- FALSE
  iter1 <- 0
  if(trace) {cat("\n", "Estimating penalty multipliers by empirical Bayes", "\n", sep="")}
  srt <- proc.time()[3]
  while(!conv1 & (iter1 < maxiter)) {
    
    iter1 <- iter1 + 1
    
    # estimating new lambdag
    opt <- optim(par=c(s, log(unlist(lambdag, use.names=FALSE))), fn=fopt_groups, lambda2=lambda2,
                 nparts=nparts, partsind=partsind, partsmat=partsmat, sizes=unlist(sizes), G=G,
                 sum1=sum1, method="BFGS", control=list(maxit=1000))
    opt.conv <- c(opt.conv, opt$convergence==0)
    lambdagnew <- split(exp(opt$par[-1]), factor(rep(partnames, unlist(G)), levels=partnames))
    lambdagnew <- sapply(1:nparts, function(part) {
      if(monotone$monotone[part]) {
        return(as.numeric(pava(lambdagnew[[part]], sizes[[part]], monotone$decreasing[part])))
      } else {
        return(lambdagnew[[part]])
      }}, simplify=FALSE)
    if(any(monotone$monotone)) {
      lambdagnew <- sapply(1:nparts, function(part) {exp(log(lambdagnew[[part]]) - sum(sapply(1:nparts, function(part) {
        sum(sizes[[part]]*log(lambdagnew[[part]]))}))/(p*nparts))}, simplify=FALSE)
    }
    s <- opt$par[1]
    lambdagseq <- sapply(1:nparts, function(part) {cbind(lambdagseq[[part]], lambdagnew[[part]])}, simplify=FALSE)
    names(lambdagseq) <- partnames
    lambdamultvec <- apply(sapply(1:nparts, function(part) {lambdagnew[[part]][partitions[[part]]]}), 1, prod)
    
    # inner loop of algorithm:
    conv2 <- FALSE
    iter2 <- 0
    while(!conv2 & (iter2 < maxiter)) {
      iter2 <- iter2 + 1
      
      # estimating new model parameters
      newparam <- est_param3(xr, xu, kappa, m, n, p, ciold, phi, chiold, lambda2, lambdamultvec, lambdamultvecold,
                             intercept, !is.null(unpenalized), FALSE, ELBO, FALSE)
      
      dsigma <- as.numeric(newparam$dsigma)
      mu <- as.numeric(newparam$mu)
      ci <- as.numeric(newparam$ci)
      chi <- as.numeric(newparam$chi)
      if(ELBO) {
        ELBOnew <- newparam$elbo
        ELBOseq <- c(ELBOseq, ELBOnew)
      }
      
      # checking convergence of inner loop
      conv2 <- (max(c(abs((mu - muold)/ifelse(muold==0, muold + 0.00001, muold)), 
                      abs((dsigma - dsigmaold)/ifelse(dsigmaold==0, dsigmaold + 0.00001, dsigmaold)))) < eps)
      muold <- mu
      dsigmaold <- dsigma
      ciold <- ci
      chiold <- chi
      
    }
    vb.conv <- c(vb.conv, conv2)
    
    niter2seq <- c(niter2seq, iter2)
    
    # sum is needed in optimisation routine
    sum1 <- sapply(uintsec, function(int) {
      ind <- which(intsec==int) + intercept;
      sum((dsigmaold[ind + u] + muold[ind + u]^2)*(1 + sqrt(phi/chiold[ind - intercept])))})
    
    # checking convergence of outer loop:
    conv1 <- max(abs((unlist(lambdagnew) - unlist(lambdag))/unlist(lambdag))) < eps
    
    # updating lambdag for new iteration
    lambdagold <- lambdag
    lambdag <- lambdagnew
    lambdamultvecold <- lambdamultvec
    
    # printing progress
    if(trace) {
      cat("\r", "Penalty multipliers estimated at ", paste(partnames, lapply(lambdag, function(part) {
        paste(round(part, 2), collapse=", ")}), sep=": ", collapse=" and "), "      ", sep="")
    }
    
  }
  
  # if a posterior is desired, estimated it here
  if(posterior) {
    newparam <- est_param3(xr, xu, kappa, m, n, p, ci, phi, chi, lambda2, lambdamultvec, lambdamultvecold,
                           intercept, !is.null(unpenalized), TRUE, FALSE, FALSE)
    dsigma <- newparam$sigma
    mu <- as.numeric(newparam$mu)
    ci <- as.numeric(newparam$ci)
    chi <- as.numeric(newparam$chi)
  }
  
  # printing estimation time and final estimates
  eb.time <- proc.time()[3] - srt
  if(trace) {cat("\r", "Penalty multipliers estimated at ", paste(partnames, lapply(lambdag, function(part) {
    paste(round(part, 2), collapse=", ")}), sep=": ", collapse=" and "), " in ", 
    round(eb.time, 2), " seconds ", sep="")
  }
  
  # variable selection and run penalized regression for prediction purposes
  if(!is.null(psel)) {
    lambda.min.ratio <- ifelse(n < p, 0.01, 0.0001)
    nlambda <- 100
    # lambdamax1 <- max(abs(t(xr) %*% (y - mean(y)*(1 - mean(y))))/
    #                     sqrt(lambdamultvec))/(n*ifelse(alpha==0, 0.001, alpha))
    fit.lambda1 <- glmnet(x, y, family="binomial", alpha=alpha,
                         standardize=FALSE, intercept=intercept,
                         penalty.factor=c(rep(0, u), lambdamultvec))
    fit.lambda2 <- glmnet(x, y, family="binomial", alpha=alpha,
                         standardize=FALSE, intercept=intercept,
                         penalty.factor=c(rep(0, u), rep(1, p)))
    lambdamax1 <- max(fit.lambda1$lambda)
    lambdamax2 <- max(fit.lambda2$lambda)
    # lambdamin1 <- max(lambdamax1*lambda.min.ratio, 1e-06)
    lambdamin1 <- lambdamin2 <- 1e-10
    
    fsel <- function(lambda, maxselec, alpha, groupreg) {
      if(groupreg) {
        fit.sel <- glmnet(x, y, family="binomial", alpha=alpha, lambda=lambda,
                          standardize=FALSE, intercept=intercept,
                          penalty.factor=c(rep(0, u), lambdamultvec))
      } else {
        fit.sel <- glmnet(x, y, family="binomial", alpha=alpha, lambda=lambda,
                          standardize=FALSE, intercept=intercept,
                          penalty.factor=c(rep(0, u), rep(1, r)))
      }
      return(fit.sel$df - u - maxselec)
    }
  
    sel.out <- sapply(psel, function(parsel) {
      lambda.sel <- uniroot(fsel, interval=c(lambdamin1, lambdamax1), 
                            maxiter=50, maxselec=parsel, alpha=alpha, 
                            groupreg=TRUE)$root
      fit.sel <- glmnet(x, y, family="binomial", alpha=alpha,
                        lambda=lambda.sel, standardize=FALSE, 
                        intercept=intercept,
                        penalty.factor=c(rep(0, u), lambdamultvec))
      return(list(beta=as.numeric(coef(fit.sel)), lambda=lambda.sel))}, 
      simplify=FALSE)
    beta <- sapply(sel.out, function(l) {return(l[[1]])})
    lambdaseq <- sapply(sel.out, function(l) {return(l[[2]])})

    sel.nogroups.out <- sapply(psel, function(parsel) {
      lambda.sel <- uniroot(fsel, interval=c(lambdamin1, lambdamax1), 
                            maxiter=50, maxselec=parsel, alpha=alpha, 
                            groupreg=FALSE)$root
      fit.sel <- glmnet(x, y, family="binomial", alpha=alpha,
                        lambda=lambda.sel, standardize=FALSE, 
                        intercept=intercept,
                        penalty.factor=c(rep(0, u), rep(1, r)))
      return(list(beta=as.numeric(coef(fit.sel)), lambda=lambda.sel))}, 
      simplify=FALSE)
    beta.nogroups <- sapply(sel.nogroups.out, function(l) {return(l[[1]])})
    lambdaseq.nogroups <- sapply(sel.nogroups.out, function(l) {return(l[[2]])})
    
  } else {
    fit.final <- glmnet(x, y=ymat, family="binomial", alpha=alpha, 
                        standardize=FALSE, intercept=intercept, 
                        penalty.factor=c(rep(0, u), lambdamultvec))
    beta <- as.matrix(coef(fit.final))
    lambdaseq <- fit.final$lambda
    
    fit.final.nogroups <- glmnet(x, y, family="binomial", alpha=alpha,
                                 standardize=FALSE, intercept=intercept,
                                 penalty.factor=c(rep(0, u), rep(1, r)))
    beta.nogroups <- as.matrix(coef(fit.final.nogroups))
    lambdaseq.nogroups <- fit.final.nogroups$lambda
  }
  
  # output list
  conv <- list(lambda.conv=conv1, opt.conv=opt.conv, vb.conv=vb.conv)
  if(ELBO) {
    out <- list(mu=mu, sigma=dsigma, ci=ci, chi=chi, alpha=alpha, 
                lambda.min=lambda, lambda=lambdaseq, 
                lambda.nogroups=lambdaseq.nogroups, lambdag=lambdagseq, 
                ELBO=ELBOseq, nouteriter=iter1, ninneriter=niter2seq, conv=conv, 
                beta=beta, beta.nogroups=beta.nogroups, args=args)
  } else {
    out <- list(mu=mu, sigma=dsigma, ci=ci, chi=chi, alpha=alpha, 
                lambda.min=lambda, lambda=lambdaseq, 
                lambda.nogroups=lambdaseq.nogroups, lambdag=lambdagseq, 
                nouteriter=iter1, ninneriter=niter2seq, conv=conv, beta=beta, 
                beta.nogroups=beta.nogroups, args=args)
  } 
  return(out)
  
}

# function to cross-validate the penalty parameters
cv.pen <- function(xr, y, unpenalized=NULL, intercept, psel=NULL, nfolds=NULL) {
  p <- ncol(xr)
  n <- nrow(xr)
  nfolds <- ifelse(is.null(nfolds), 10, nfolds)
  rest <- n %% nfolds
  foldsize <- c(rep(n %/% nfolds + as.numeric(rest!=0), times=rest),
                rep(n %/% nfolds, times=nfolds - rest))
  foldid <- sample(rep(1:nfolds, times=foldsize))
  r <- ncol(xr)
  u <- ifelse(is.null(ncol(unpenalized)), 0, ncol(unpenalized))
  x <- cbind(unpenalized, xr)
  seq.alpha <- seq(0.01, 0.99, length.out=50)
  seq.lam <- seq.df <- seq.cvll <- numeric(length(seq.alpha))
  for(a in 1:length(seq.alpha)) {
    cv.fit <- cv.glmnet(x, y, family="binomial", alpha=seq.alpha[a], standardize=FALSE,
                        intercept=intercept, penalty.factor=c(rep(0, u), rep(1, r)),
                        dfmax=ifelse(is.null(psel), p + 1, psel + u), foldid=foldid,
                        grouped=FALSE)
    ind <- which(cv.fit$lambda==cv.fit$lambda.min)
    seq.lam[a] <- cv.fit$lambda.min
    seq.df[a] <- cv.fit$nzero[ind]
    seq.cvll[a] <- cv.fit$cvm[ind]
  }
  lambda1bayes <- seq.lam*seq.alpha*n*2
  lambda2bayes <- seq.lam*(1 - seq.alpha)*n
  lambda1pen <- seq.lam*seq.alpha*n
  lambda2pen <- seq.lam*(1 - seq.alpha)*0.5*n
  
  out <- list(lambda1bayes=lambda1bayes, lambda2bayes=lambda2bayes, lambda1pen=lambda1pen,
              lambda2pen=lambda2pen, alpha=seq.alpha, lambda=seq.lam, cvll=seq.cvll)
  return(out)
  
}



# optimisation function for penalty multipliers in groups
fopt_groups <- function(par, lambda2, nparts, partsind, partsmat, sizes, G, sum1) {
  
  # s <- rep(par[1:nparts], times=unlist(G))
  # loglambdag <- par[-c(1:nparts)]
  s <- par[1]
  loglambdag <- par[-1]
  
  loglambdasum <- rowSums(sapply(1:nparts, function(part) {
    loglambdag[partsind==part][partsmat[, part]]}))
  
  partsum <- sum((0.5*lambda2*unlist(sapply(1:nparts, function(part) {
    tapply(exp(loglambdasum)*sum1, partsmat[, part], sum)})) + (s - 0.5)*sizes)^2)
  constr <- sum(loglambdag*sizes)^2
  magn <- partsum + constr
  return(magn)
  
}

# variable selection after the VBEM
select.grEBEN <- function(object, x, y, m, partitions, lambda1=NULL, lambda2=NULL, nsel, intercept, 
                          posterior, eps, maxiter, trace=TRUE, alphastart) {
  
  nparts <- length(partitions)
  sizes <- lapply(partitions, function(part) {rle(sort(part))$lengths})
  fit.grVBEM <- grVBEM2(x, y, m, partitions, lambda1, lambda2, intercept, monotone, 
                        posterior, eps, maxiter, trace, alphastart)
  lambdagvec <- exp(colSums(log(sapply(1:nparts, function(part) {
    rep(fit.grVBEM$lambdag[[part]][, fit.grVBEM$nouteriter + 1], times=sizes[[part]])}))))
  
  alpha <- lambda1/(2*lambda2 + lambda1)
  
  # selecting variables
  p <- ncol(x)
  n <- nrow(x)
  lambda.min.ratio <- ifelse(n < p, 0.01, 0.0001)
  nlambda <- 100
  
  
  lambdaseq <- NULL
  found <- inbet <- FALSE
  count <- 0
  while(!found) {
    suppressWarnings(fit.df <- glmnet(x, y, family="binomial", alpha=alpha, 
                                      nlambda=nlambda, lambda=lambdaseq, standardize=FALSE, 
                                      intercept=intercept, 
                                      dfmax=ifelse(is.null(lambdaseq), nsel, p + 1)))
    if(any(fit.df$df==nsel)) {
      if(tail(fit.df$df, n=1L)==nsel & !inbet) {
        
        lambdamax <- tail(fit.df$lambda[fit.df$df==nsel], n=1L)
        lambdamin <- lambdamax*lambda.min.ratio
      } else {
        found <- TRUE
      }
    } else if(all(fit.df$df < nsel)) {
      lambdamax <- tail(fit.df$lambda, n=1L)
      lambdamin <- lambdamax*lambda.min.ratio
    } else if(all(fit.df$df > nsel)) {
      lambdamax <- max(abs(t(x) %*% (y - mean(y)*(1 - mean(y)))))/(n*alpha)
      lambdamin <- max(fit.df$lambda)
      inbet <- TRUE
      count <- count + 1
      if(count==5) {
        found <- TRUE
      }
    } else if(length(unique(lambdaseq))==1) {
      found <- TRUE
    } else {
      ind <- c(rle(fit.df$df < nsel)$lengths[1], rle(fit.df$df < nsel)$lengths[1] + 1)
      lambdamax <- max(fit.df$lambda[ind])
      lambdamin <- min(fit.df$lambda[ind])
      inbet <- TRUE
    }
    if(!found) {
      lambdaseq <- exp(seq(log(lambdamax), log(lambdamin), length.out=nlambda))
    }    
  }
  ind <- ifelse(any(fit.df$df==nsel), tail(which(fit.df$df==nsel), n=1L),
                which.max(fit.df$df[which.min(abs(fit.df$df - nsel))]))
  indsel <- which(as.numeric(fit.df$beta[, ind])!=0)
  return(indsel)
}

# function to make predictions from estimated model
predict.grEBEN <- function(object, newx, unpenalized=NULL, 
                           type=c("VB", "penalized", "nogroups", "selection")) {
  newx <- cbind(unpenalized, newx)
  if(object$args$intercept) {newx = cbind(1, newx)}
  if(type=="VB") {
    lpred <- newx %*% as.matrix(object$mu)
    ppred <- exp(lpred)/(1 + exp(lpred))
  } else if(type=="penalized") {
    lpred <- newx %*% as.matrix(object$beta)
    ppred <- exp(lpred)/(1 + exp(lpred))
  } else if(type=="nogroups") {
    lpred <- newx %*% as.matrix(object$beta.nogroups)
    ppred <- exp(lpred)/(1 + exp(lpred))
  } else {
    lpred <- newx %*% as.matrix(object$beta.sel)
    ppred <- exp(lpred)/(1 + exp(lpred))
  }
  return(as.numeric(ppred))
}

predict.grEBEN2 <- function(object, newx, which=NULL, unpenalized=NULL, 
                            type=c("VB", "penalized", "nogroups", "selection")) {
  newx <- cbind(unpenalized, newx)
  if(object$args$intercept) {newx = cbind(1, newx)}
  if(type=="VB") {
    lpred <- newx %*% as.matrix(object$mu)
    ppred <- exp(lpred)/(1 + exp(lpred))
  } else if(type=="penalized") {
    lpred <- newx %*% as.matrix(object$beta[, which])
    ppred <- exp(lpred)/(1 + exp(lpred))
  } else if(type=="nogroups") {
    lpred <- newx %*% as.matrix(object$beta.nogroups[, which])
    ppred <- exp(lpred)/(1 + exp(lpred))
  } else {
    lpred <- newx %*% as.matrix(object$beta.sel[, which])
    ppred <- exp(lpred)/(1 + exp(lpred))
  }
  return(as.numeric(ppred))
}

