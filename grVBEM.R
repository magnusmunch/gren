##############################  preamble  #############################
# grMCEM implemented in R and C++                                     #
# version: 01                                                         #
# author: Magnus M?nch                                                #
# created: 14-03-2017                                                 #
# last edited: 28-03-2017                                             #
#######################################################################

###############################  notes  ###############################
# 28-03-2017: Added multiple partitions                               #
# 15-03-2017: First implementation of group-regularized elastic net   #
#             with lambda1 penalty multiplier, the square root of     #
#             lambda2 penalty multiplier                              #
#######################################################################

# paths
path.code <- as.character(ifelse(Sys.info()[1]=="Darwin","/Users/magnusmunch/Documents/PhD/EBEN/code/" ,"~/EBEN/code/"))

### libraries
library(Rcpp)
library(glmnet)
library(GRridge)

### functions
# source function for variational Bayes
sourceCpp(paste(path.code, "ENVB2.cpp", sep=""))

# functions to simulate elastic net prior model parameters
qtrunc <- function(p, spec, a = -Inf, b = Inf, ...) {
  
  tt <- p
  G <- get(paste("p", spec, sep = ""), mode = "function")
  Gin <- get(paste("q", spec, sep = ""), mode = "function")
  tt <- Gin(G(a, ...) + p*(G(b, ...) - G(a, ...)), ...)
  return(tt)
  
}

rtrunc <- function(n, spec, a = -Inf, b = Inf, ...){
  
  x <- u <- runif(n, min = 0, max = 1)
  x <- qtrunc(u, spec, a = a, b = b,...)
  return(x)
  
}

renbeta <- function(p, lambda1, lambda2) {
  
  if(lambda1==0) {
    beta <- rnorm(p, 0, sqrt(1/lambda2))
  } else if(lambda2==0) {
    u <- runif(p) - 0.5
    beta <- -2*sign(u)*log(1 - 2*abs(u))/lambda1
  } else {
    shape <- 0.5
    scale <- 8*lambda2/lambda1^2
    tau <- rtrunc(n=p, spec="gamma", a=1, b=Inf, shape=shape, scale=scale)
    tau <- ifelse(tau < 1, 1, tau)
    sigma <- ifelse(is.infinite(tau), 1/lambda2, 
                    (tau - 1)/(lambda2*tau))
    beta <- rnorm(p, 0, sqrt(sigma))
  }
  return(beta)
  
}

# function for estimation
grEBEN <- function(x, y, m, unpenalized=NULL, intercept=TRUE, partitions=NULL, codata=NULL, 
                   lambda1=NULL, lambda2=NULL, maxsel=NULL, posterior=FALSE, ELBO=TRUE, eps=0.001, 
                   maxiter=500, trace=TRUE) {
  
  # save argument list
  args <- list(x=x, y=y, m=m, unpenalized=unpenalized, intercept=intercept, partitions=partitions, 
               codata=codata, lambda1=lambda1, lambda2=lambda2, maxsel=maxsel, posterior=posterior, 
               ELBO=ELBO, eps=eps, maxiter=maxiter, trace=trace)
  
  # if no penalty parameters are given we estimate them by cross-validation
  if(is.null(lambda1) | is.null(lambda2)) {
    if(trace) {cat("\r", "Estimating global lambda1 and lambda2 by cross-validation", sep="")}
    srt <- proc.time()[3]
    opt.glob <- cv.pen(x, y, intercept)
    cv.time <- proc.time()[3] - srt
    lambda1 <- opt.glob$lambda1bayes[which.min(opt.glob$cvll)]
    lambda2 <- opt.glob$lambda2bayes[which.min(opt.glob$cvll)]
    if(trace) {cat("\n", "Global lambda1 and lambda2 estimated at ", round(lambda1, 2), " and ", 
                   round(lambda2, 2), " in ", round(cv.time, 2), " seconds", sep="")}
  } 
  n <- nrow(x)
  alpha <- lambda1/(lambda1 + 2*lambda2)
  lambda <- (0.5*lambda1 + lambda2)/n
  
  # assigning fixed (throughout algorithm) variables
  if(!is.null(partitions)) {
    if(!is.list(partitions) & is.null(names(partitions))) {
      partitions <- list(partition1=partitions)
    } else if(is.null(names(partitions))) {
      names(partitions) <- paste("partition", 1:length(partitions), sep="")
    }
    partnames <- names(partitions)
    nparts <- length(partitions)
    sizes <- lapply(partitions, function(part) {rle(sort(part))$lengths})
    G <- lapply(partitions, function(part) {length(unique(part))})
  } else {
    if(!is.list(codata)) {codata <- list(codata1=codata)}
    partitions <- lapply(codata, function(part) {return(c(1:p))})
    if(is.null(names(codata))) {
      names(partitions) <- paste("codata", 1:length(partitions), sep="")
    } else {
      names(partitions) <- names(codata)
    }
    partnames <- names(partitions)
    nparts <- length(partitions)
    sizes <- lapply(partitions, function(part) {return(rep(1, p))})
    G <- lapply(partitions, function(part) {return(p)})
  }
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
  kappa <- y - m/2
  phi <- 0.25*lambda1^2/lambda2
  
  # starting values for lambdag, lagrange multiplier s
  lambdagnew <- lambdag <- lambdagold <- lambdagseq <- lapply(G, function(gpart) {rep(1, gpart)})
  lambdamultvecold <- rep(1, r)
  # if(!is.null(codata)) {a <- 1}
  # if(!is.null(codata)) {a <- rep(1, 3)}
  if(!is.null(codata)) {a <- rep(1, 2)}
  s <- 0
  
  # starting values for the model parameters
  fit.start <- glmnet(x=x, y=ymat, family="binomial", alpha=0, lambda=(1 - alpha)*lambda + alpha*lambda, 
                      standardize=FALSE, intercept=intercept, penalty.factor=c(rep(0, u), rep(1, r)))
  startparam <- est_param3(xr, xu, kappa, m, n, p, as.numeric(predict(fit.start, newx=x, type="response")), 
                           phi, rep(phi, r), lambda2, lambdamultvecold, lambdamultvecold, intercept, 
                           !is.null(unpenalized), FALSE, FALSE, TRUE)
  dsigmaold <- as.numeric(startparam$dsigma)
  muold <- as.numeric(startparam$mu)
  ciold <- as.numeric(startparam$ci)
  chiold <- as.numeric(startparam$chi)
  
  cursec <- 1
  # calculating the sums needed in optimisation routine
  intsec <- do.call("paste", c(partitions, sep=" "))
  uintsec <- unique(intsec)
  sum1 <- sapply(1:length(uintsec), function(cursec) {
    ind <- which(intsec==uintsec[[cursec]]) + intercept;
    sum((dsigmaold[ind + u] + muold[ind + u]^2)*(1 + sqrt(phi/chiold[ind - intercept])))})
  intsizes <- as.numeric(t(table(as.data.frame(matrix(unlist(partitions), ncol=nparts)))))
  partsmat <- unique(matrix(unlist(partitions), ncol=nparts))
  partsind <- rep(1:nparts, times=unlist(G))
  
  # keep track of iterations and possibly evidence lower bound
  if(ELBO) {
    ELBOseq <- numeric(0)
    ELBOold <- Inf
  }
  niter2seq <- numeric(0)
  
  # outer loop of algorithm:
  conv1 <- FALSE
  iter1 <- 0
  if(trace) {cat("\n", "Estimating penalty multipliers by empirical Bayes", "\n", sep="")}
  srt <- proc.time()[3]
  while(!conv1 & (iter1 < maxiter)) {
    
    iter1 <- iter1 + 1
    
    # estimating new lambdag
    if(is.null(codata)) {
      opt <- optim(par=c(s, log(unlist(lambdag, use.names=FALSE))), fn=fopt_groups, lambda2=lambda2, 
                   nparts=nparts, partsind=partsind, partsmat=partsmat, sizes=unlist(sizes), G=G, 
                   sum1=sum1, method="Nelder-Mead", control=list(maxit=5000))
      
      lambdagnew <- split(exp(opt$par[-1]), factor(rep(partnames, unlist(G)), levels=partnames))
    } else {
      opt <- optim(par=c(s, a), fn=fopt_cont, codata=unlist(codata, use.names=FALSE), lambda2=lambda2, 
                   nparts=nparts, partsind=partsind, partsmat=partsmat, sizes=unlist(sizes), G=G, 
                   sum1=sum1, method="Nelder-Mead", control=list(maxit=5000))
      # a <- opt$par[2]
      # a <- opt$par[-1]
      a <- c(exp(opt$par[2]), opt$par[3])
      # lambdagnew <- split(unlist(codata)^a, factor(rep(partnames, unlist(G)), levels=partnames))
      # lambdagnew <- split(a[1]*exp(-a[2]*exp(-a[3]*unlist(codata))), 
      #                     factor(rep(partnames, unlist(G)), levels=partnames))
      lambdagnew <- split(a[1]*exp(a[2]*unlist(codata)), 
                          factor(rep(partnames, unlist(G)), levels=partnames))
      
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
      ELBOnew <- newparam$elbo
      ELBOseq <- c(ELBOseq, ELBOnew)
      
      # checking convergence of inner loop
      conv2 <- (max(c(abs((mu - muold)/ifelse(muold==0, muold + 0.00001, muold)), 
                      abs((dsigma - dsigmaold)/ifelse(dsigmaold==0, dsigmaold + 0.00001, dsigmaold)))) < eps)
      ELBOold <- ELBOnew
      muold <- mu
      dsigmaold <- dsigma
      ciold <- ci
      chiold <- chi
      
    }
    
    niter2seq <- c(niter2seq, iter2)
    
    # sum is needed in optimisation routine
    sum1 <- sapply(1:length(uintsec), function(cursec) {
      ind <- which(intsec==uintsec[[cursec]]) + intercept;
      sum((dsigmaold[ind + u] + muold[ind + u]^2)*(1 + sqrt(phi/chi[ind - intercept])))})
    
    # checking convergence of outer loop:
    conv1 <- max(abs((unlist(lambdagnew) - unlist(lambdag))/unlist(lambdag))) < eps
    
    # updating lambdag for new iteration
    lambdagold <- lambdag
    lambdag <- lambdagnew
    
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
  if(!is.null(maxsel)) {
    
  }
  
  # run penalized regression for prediction purposes
  fit.final <- glmnet(x=t(t(x)/c(rep(1, u), sqrt(lambdamultvec))), y=ymat, family="binomial", alpha=0, 
                      lambda=(1 - alpha)*lambda + alpha*lambda, standardize=FALSE, intercept=intercept, 
                      penalty.factor=c(rep(0, u), rep(1, r)))
  beta <- as.numeric(coef(fit.final))/c(rep(1, u + intercept), sqrt(lambdamultvec))
  
  if(ELBO) {
    out <- list(mu=mu, sigma=dsigma, ci=ci, chi=chi, lambda1=lambda1, lambda2=lambda2, 
                lambdag=lambdagseq, ELBO=ELBOseq, nouteriter=iter1, ninneriter=niter2seq, 
                conv=conv1, pen.beta=beta, args=args)
  } else {
    out <- list(mu=mu, sigma=dsigma, ci=ci, chi=chi, lambda1=lambda1, lambda2=lambda2, 
                lambdag=lambdagseq, nouteriter=iter1, ninneriter=niter2seq, 
                conv=conv1, pen.beta=beta, args=args)
  }
  return(out)
  
}

# function to cross-validate the penalty parameters
cv.pen <- function(x, y, intercept) {
  n <- nrow(x)
  seq.alpha <- seq(0.01, 0.99, length.out=50)
  seq.lam <- seq.df <- seq.cvll <- numeric(length(seq.alpha))
  for(a in 1:length(seq.alpha)) {
    cv.fit <- cv.glmnet(x, y, family="binomial", alpha=seq.alpha[a], standardize=FALSE,
                        intercept=intercept)
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
  
  loglambdamult <- rowSums(sapply(1:nparts, function(part) {
    loglambdag[partsind==part][partsmat[, part]]}))
  
  partsum <- sum((0.5*lambda2*unlist(sapply(1:nparts, function(part) {
    tapply(exp(loglambdamult)*sum1, partsmat[, part], sum)})) + (s - 0.5)*sizes)^2)
  constrsum <- sum(sapply(1:nparts, function(part) {
    sum(loglambdag[partsind==part]*sizes[partsind==part])^2}))
  magn <- partsum + constrsum
  return(magn)
  
}

# for continuous co-data
# lambdaj = cj^alpha
# fopt_cont <- function(par, codata, lambda2, nparts, partsind, partsmat, sizes, G, sum1) {
#   
#   s <- par[1]
#   alpha <- par[2]
#   
#   sum2 <- sum(log(codata))
#   magn <- sqrt((0.5*lambda2*sum(codata^alpha*log(codata)*sum1) + (s - 0.5)*sum2)^2 + alpha^2*sum2^2)
#     
#   return(magn)
# }

# lambdaj = alpha1*exp(-alpha2*exp(-alpha3*cj))
# fopt_cont <- function(par, codata, lambda2, nparts, partsind, partsmat, sizes, G, sum1) {
#   
#   p <- length(codata)
#   s <- par[1]
#   a1 <- par[2]
#   a2 <- par[3]
#   a3 <- par[4]
#   
#   ej <- exp(-a3*codata)
#   ej2 <- exp(-a2*ej)
#   der1 <- (s*p - p/2)/a1 + 0.5*lambda2*sum(ej2*sum1)
#   der2 <- 0.5*sum(ej) - 0.5*lambda2*a1*sum(ej2*sum1*ej) - s*sum(ej)
#   der3 <- 0.5*lambda2*a1*sum(ej2*sum1*a2*ej*codata) - 0.5*a2*sum(ej*codata) + s*a2*sum(ej*codata)
#   ders <- p*log(a1) - a2*sum(ej)
#   
#   magn <- ifelse(a1<=0, Inf, sqrt(der1^2 + der2^2 + der3^2 + ders^2))
#   return(magn)
# }

# lambdaj = alpha1*exp(alpha2*cj)
fopt_cont <- function(par, codata, lambda2, nparts, partsind, partsmat, sizes, G, sum1) {
  
  p <- length(codata)
  s <- par[1]
  loga1 <- par[2]
  a2 <- par[3]
  
  sum2 <- sum(codata)
  der1 <- 0.5*lambda2*sum(exp(loga1 + a2*codata)*sum1) + (s - 0.5)*p
  der2 <- 0.5*lambda2*sum(exp(loga1 + a2*codata)*sum1*codata) + (s - 0.5)*sum2
  ders <- p*loga1 + a2*sum2
  magn <- sqrt(der1^2 + der2^2 + ders^2)
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
predict.grEBEN <- function(object, newx, unpenalized=NULL, type=c("VB", "penalized")) {
  newx <- cbind(unpenalized, newx)
  if(object$args$intercept) {newx = cbind(1, newx)}
  if(type=="VB") {
    lpred <- newx %*% as.matrix(object$mu)
    ppred <- exp(lpred)/(1 + exp(lpred))
  } else {
    lpred <- newx %*% as.matrix(object$pen.beta)
    ppred <- exp(lpred)/(1 + exp(lpred))
  }
  return(as.numeric(ppred))
}
