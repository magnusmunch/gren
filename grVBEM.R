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
library(penalized)
library(nloptr)

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
grEBEN <- function(x, y, m, unpenalized=NULL, intercept=TRUE, partitions, lambda1=NULL, lambda2=NULL, 
                   maxsel=NULL, posterior=FALSE, ELBO=TRUE, eps=0.001, maxiter=500, trace=TRUE) {
  
  # save argument list
  args <- list(x=x, y=y, m=m, unpenalized=unpenalized, intercept=intercept, partitions=partitions, 
               lambda1=lambda1, lambda2=lambda2, maxsel=maxsel, posterior=posterior, ELBO=ELBO, eps=eps,
               maxiter=maxiter, trace=trace)
  
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
  if(!is.list(partitions) & is.null(names(partitions))) {
    partitions <- list(partition1=partitions)
  } else if(is.null(names(partitions))) {
    names(partitions) <- paste("partition", 1:length(partitions), sep="")
  }
  partnames <- names(partitions)
  nparts <- length(partitions)
  sizes <- lapply(partitions, function(part) {rle(sort(part))$lengths})
  G <- lapply(partitions, function(part) {length(unique(part))})
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
  
  cursec <- 3
  # calculating the sums needed in optimisation routine
  intsec <- do.call("paste", c(partitions, sep=" "))
  uintsec <- unique(intsec)
  sum1 <- sapply(1:length(uintsec), function(cursec) {
    ind <- which(intsec==uintsec[[cursec]]) + intercept;
    sum((dsigmaold[ind + u] + muold[ind + u]^2)*(1 + sqrt(phi/chi[ind - intercept])))})
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
    opt <- optim(par=c(s, log(unlist(lambdag, use.names=FALSE))), fn=grmagn3, lambda2=lambda2, 
                 nparts=nparts, partsind=partsind, partsmat=partsmat, sizes=unlist(sizes), G=G, 
                 sum1=sum1, method="Nelder-Mead", control=list(maxit=5000))
    s <- opt$par[1]
    lambdagnew <- split(exp(opt$par[-1]), factor(rep(partnames, unlist(G)), levels=partnames))
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

# optimisation function for penalty multipliers
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



# set.seed(123)
# n <- 200
# p <- 150
# G <- 3
# partitions <- list(groups=rep(1:G, each=p/G))
# x <- matrix(rnorm(n*p), ncol=p, nrow=n)
# lambda1 <- 1
# lambda2 <- 1
# lambdag <- c(0.2, 3, 1/(0.2*3))
# m <- rep(1, n)
# b0 <- rnorm(p + 1)
# sigma0 <- diag(rchisq(p + 1, 1))
# beta <- c(renbeta(p/G, lambda1*sqrt(lambdag[1]), lambda2*lambdag[1]), 
#           renbeta(p/G, lambda1*sqrt(lambdag[2]), lambda2*lambdag[2]),
#           renbeta(p/G, lambda1*sqrt(lambdag[3]), lambda2*lambdag[3]))
# y <- rbinom(n, m, exp(x %*% beta)/(1 + exp(x %*% beta)))
# 
# ntest <- 10
# xtest <- matrix(rnorm(n*p), ncol=p, nrow=ntest)
# 
# fit.test1 <- grEBEN(x, y, m, unpenalized=NULL, intercept=TRUE, partitions, lambda1=NULL, lambda2=NULL, 
#                     posterior=TRUE, ELBO=TRUE, eps=0.001, maxiter=500, trace=TRUE)
# 
# predpen.test1 <- predict.grEBEN(fit.test1, xtest, unpenalized=NULL, type="penalized")
# predvb.test1 <- predict.grEBEN(fit.test1, xtest, unpenalized=NULL, type="VB")
# 
# plot(fit.test1$pen.beta, fit.test1$mu)
# plot(predpen.test1, predvb.test1)


# grmagn1 <- function(par, lambda2, sizes, sum1) {
#   s <- par[1]
#   loglambdag <- par[-1]
#   magn <- sum((0.5*lambda2*sum1*exp(loglambdag) - 0.5*sizes + s*sizes)^2) + sum(sizes*loglambdag)^2
#   return(magn)
# }
# 
# grmagn2 <- function(par, lambda2, nparts, sizes, G, uintsec, sum1) {
#   
#   s <- par[1:nparts]
#   loglambdag <- par[-c(1:nparts)]
#   partind <- rep(1:nparts, times=unlist(G))
#   
#   ulambdagprod <- sapply(1:length(uintsec), function(cursec) {exp(sum(sapply(1:nparts, function(part) {
#     loglambdag[partind==part][as.numeric(strsplit(uintsec, split=" ")[[cursec]][part])]})))})
#   partsum <- sum(sapply(1:nparts, function(part) {sum((sapply(1:G[[part]], function(g) {
#     0.5*lambda2*sum(ulambdagprod*sum1*sapply(1:length(uintsec), function(cursec) {
#       as.numeric(strsplit(uintsec, split=" ")[[cursec]][part])==g}))}) - 0.5*sizes[[part]] + s[part]*sizes[[part]])^2)}))
#   constrsum <- sum(sapply(1:nparts, function(part) {sum(sizes[[part]]*loglambdag[partind==part])^2}))
#   magn <- partsum + constrsum
#   return(magn)
#   
# }


# grmagn3 <- function(par, lambda2, nparts, partsind, partsmat, sizes, intsizes, sum1) {
#   
#   s <- par[1]
#   loglambdag <- par[-1]
#   
#   loglambdamult <- rowSums(sapply(1:nparts, function(part) {
#     loglambdag[partsind==part][partsmat[, part]]}))
#   
#   partsum <- sum((0.5*lambda2*unlist(sapply(1:nparts, function(part) {
#     tapply(exp(loglambdamult)*sum1, partsmat[, part], sum)})) + (s - 0.5)*sizes*exp(loglambdag))^2)
#   constrsum <- sum(sizes*loglambdag)^2
#   magn <- partsum + constrsum
#   return(magn)
#   
# }



# # monotone optimization functions
# eval_f_list <- function(x, lambda2, sizes, sum1, sum2) {
#   comterm <- 0.5*lambda2*sum1*exp(cumsum(x))
#   return(list(objective=as.numeric(sum(comterm) - 0.5*sum(sum2*x)),
#               gradient=as.numeric(comterm - 0.5*sum2)))
# }
# 
# eval_g_list <- function(x, lambda2, sizes, sum1, sum2) {
#   return(list(constraints=as.numeric(sum(sum2*x)),
#               jacobian=as.numeric(sum2)))
# }
# 
# grmagn4 <- function(par, lambda2, sizes, sum1, sum2) {
#   s <- par[1]
#   loglambdag <- par[-1]
#   magn <- sqrt(sum((0.5*lambda2*sum1*exp(cumsum(loglambdag)) - 
#                       0.5*sum1*loglambdag + s*sum2)^2) + sum(sum2*loglambdag)^2)
#   return(ifelse(any(loglambdag[-1] <= 0), Inf, as.numeric(magn)))
# }

# monotone optimization with inequality constraints for ordering
# eval_f_list <- function(loglambdag, lambda2, sizes, sum1) {
#   return(list(objective=as.numeric(0.5*lambda2*sum(sum1*exp(loglambdag)) - 0.5*sum(sizes*loglambdag)),
#               gradient=as.numeric(0.5*lambda2*sum1*exp(loglambdag) - 0.5*sizes)))
# }
# 
# eval_g_ineq_list <- function(loglambdag, lambda2, sizes, sum1) {
#   jacobian <- cbind(diag(length(sizes) - 1), 0)
#   jacobian[(row(jacobian) - col(jacobian))==-1] <- -1
#   return(list(constraints=as.numeric(rev(diff(rev(loglambdag)))),
#               jacobian=jacobian))
# }
# 
# eval_g_eq_list <- function(loglambdag, lambda2, sizes, sum1) {
#   return(list(constraints=as.numeric(sum(sizes*loglambdag)),
#               jacobian=as.numeric(sizes)))
# }
# 
# grmagn4 <- function(par, lambda2, sizes, sum1) {
#   s <- par[1]
#   loglambdag <- par[-1]
#   magn <- sum((0.5*lambda2*sum1*exp(loglambdag) + s*sizes - 0.5*sizes)^2) + 
#     sum(sum1*loglambdag)^2
#   return(ifelse(any(diff(rev(loglambdag)) > 0), Inf, as.numeric(magn)))
# }

# # functions for Box-Cox transformed optimisation
# eval_f_bc <- function(x, cj, lambda2, p, sum1) {
#   alpha1 <- x[1]
#   alpha2 <- x[2]
#   
#   # cja2 <- cj + alpha2
#   # obj <- 0.5*lambda2*sum(sum1*((cja2)^alpha1 - 1)/alpha1) - 
#   #   0.5*sum(log(((cja2)^alpha1 - 1)/alpha1))
#   #   
#   # grad <- c(0.5*lambda2*sum((cja2)^alpha1*log(cja2)/alpha1 + 1/alpha1 - sum1*((cja2)^alpha1/alpha1^2 )) +
#   #             0.5*sum((cja2)^alpha1/(alpha1*((cja2)^alpha1 - 1)) - 1/(alpha1*((cja2)^alpha1 - 1)) - 
#   #                       (cja2)^alpha1*log(cja2)/((cja2)^alpha1 - 1)),
#   #           0.5*lambda2*sum(sum1*cja2^(alpha1 - 1)) -
#   #             0.5*sum(alpha1*(cja2)^(alpha1 - 1)/(cja2^(alpha1 - 1) - 1)))
#   obj <- 0.5*p*log(alpha1) + 0.5*lambda2*alpha2/alpha1*sum(sum1*cj^alpha1) + 
#     0.5*lambda2*alpha2/alpha1*sum(sum1) - 0.5*p*log(alpha2) - 0.5*sum(log(cj^alpha1 - 1))
#   grad <- c(0.5*p/alpha1 + 0.5*lambda2*alpha2/alpha1*sum(sum1*cj^alpha1*log(cj)) -
#               0.5*cj^alpha1*log(cj)/(cj^alpha1 - 1) - 0.5*lambda2*alpha2/alpha1*sum(sum1) - 
#               0.5*lambda2*alpha2/alpha1^2*sum(sum1*cj^alpha1))
#   return(list(objective=as.numeric(obj), gradient=as.numeric((grad))))
# }
# 
# eval_eq_bc <- function(x, cj, lambda2, p, sum1) {
#   alpha1 <- x[1]
#   alpha2 <- x[2]
#   # return(as.numeric(sum(log((cj + alpha2)^alpha1 - 1)) - p*log(alpha1)))
#   return(as.numeric(p*log(alpha2) - p*log(alpha1) + sum(log(cj^alpha1 - 1))))
# }
# 
# eval_jac_eq_bc <- function(x, cj, lambda2, p, sum1) {
#   alpha1 <- x[1]
#   alpha2 <- x[2]
#   
#   # cja2 <- cj + alpha2
#   # jac <- c(sum(cja2^alpha1*log(cja2)/(cja2^alpha1 - 1)) - p/alpha1, 
#   #          sum(alpha1*cja2^(alpha1 - 1)/(cja2^alpha1 - 1)))
#   jac <- c(sum(cj^alpha1*log(cj)/(cj^alpha1 - 1)) - p/alpha1, p/alpha2)
#   return(as.numeric(jac))
# }

# grVBEM <- function(x, y, m, partitions, lambda1=NULL, lambda2=NULL, intercept, monotone=FALSE, 
#                    eps, maxiter, trace=TRUE) {
#   
#   if(!is.list(partitions) & is.null(names(partitions))) {
#     partitions <- list(partition1=partitions)
#   } else if(is.null(names(partitions))) {
#     names(partitions) <- paste("partition", 1:length(partitions), sep="")
#   }
#   
#   partnames <- names(partitions)
#   
#   # assigning fixed (throughout algorithm) variables
#   nparts <- length(partitions)
#   sizes <- lapply(partitions, function(part) {rle(sort(part))$lengths})
#   G <- lapply(partitions, function(part) {length(unique(part))})
#   p <- ncol(x)
#   n <- nrow(x)
#   kappa <- y - m/2
#   
#   # if no penalty parameters are given we estimate them by cross-validation
#   if(is.null(lambda1) | is.null(lambda2)) {
#     if(trace) {cat("\r", "Estimating global lambda1 and lambda2 by cross-validation", sep="")}
#     srt <- proc.time()[3]
#     opt.glob <- cv.pen(x, y, intercept)
#     cv.time <- proc.time()[3] - srt
#     lambda1 <- opt.glob$lambda1[which.min(opt.glob$cvll)]
#     lambda2 <- opt.glob$lambda2[which.min(opt.glob$cvll)]
#     if(trace) {cat("\n", "Global lambda1 and lambda2 estimated at ", round(lambda1, 2), " and ", 
#         round(lambda2, 2), " in ", round(cv.time, 2), " seconds", sep="")}
#   }
#   
#   # in the multiplier setting phi does not change
#   phi <- 0.25*lambda1^2/lambda2
#   
#   # starting values for lambdag, lagrange multiplier s and alphag in case of monotonicity
#   lambdagnew <- lambdag <- lambdagold <- lambdagseq <- lapply(G, function(gpart) {rep(1, gpart)})
#   if(monotone) {
#     sum2 <- lapply(sizes, function(sizes) {rev(cumsum(rev(sizes)))})
#     # alphag <- sapply(1:nparts, function(part) {
#     #   alphag <- seq(1, 1.01, length.out=G[[part]]); 
#     #   return(alphag*exp(-sum(sum2[[part]]*log(alphag))/sum(sum2[[part]])))})
#     alphag <- lapply(G, function(G) {return(c(1, rep(1.01, G - 1)))})
#   } else {
#     s <- 0
#   }
#   
#   # starting values for mu and sigma
#   fit.pen <- penalized(y, x, unpenalized=formula(ifelse(intercept, "~1", "~0")), 
#                        model="logistic", lambda1=0, 
#                        lambda2=2*(lambda1 + lambda2), trace=FALSE)
#   #muold <- coef(fit.pen, which="all")
#   if(intercept) {
#     xadj <- cbind(1, x)
#   } else {
#     xadj <- x
#   }
#   b0 <- coef(fit.pen, which="all")
#   pred0 <- as.numeric(exp(xadj %*% b0)/(1 + exp(xadj %*% b0)))
#   w <- sqrt(pred0*(1 - pred0))
#   xw <- xadj*w
#   svdxw <- svd(xw)
#   d <- svdxw$d
#   v <- svdxw$v
#   invmat <- d^2/(d^2 + 4*(lambda1 + lambda2))^2
#   sigmaold <- t(t(v)*invmat) %*% t(v)
#   
#   # rest of the starting values follow from that
#   muold <- as.numeric(sigmaold %*% (t(xadj) %*% as.matrix(kappa)))
#   ci <- as.numeric(sqrt(colSums(t(xadj) * (sigmaold %*% t(xadj))) + (colSums(t(xadj)*muold))^2))
#   chi <- as.numeric(0.5*(lambda1 + lambda2)*(diag(sigmaold) + muold^2))[(intercept + 1):(intercept + p)]
#   
#   # sum is needed in optimisation routine
#   intsec <- do.call("paste", c(partitions, sep=" "))
#   uintsec <- unique(intsec)
#   sum1 <- sapply(1:length(uintsec), function(cursec) {
#     ind <- which(intsec==uintsec[[cursec]]) + intercept;
#     sum((diag(sigmaold)[ind] + muold[ind]^2)*(1 + sqrt(phi/chi[ind - intercept])))})
#   intsizes <- as.numeric(t(table(as.data.frame(matrix(unlist(partitions), ncol=nparts)))))
#   partsmat <- unique(matrix(unlist(partitions), ncol=nparts))
#   partsind <- rep(1:nparts, times=unlist(G))
#   
#   # keeping track of things:
#   lambdamult <- exp(rowSums(sapply(1:nparts, function(part) {
#     log(unlist(lambdag))[partsind==part][partsmat[, part]]})))
#   # lambdagprod <- sapply(1:length(uintsec), function(cursec) {prod(sapply(1:nparts, function(part) {
#   #   lambdag[[part]][as.numeric(strsplit(uintsec[cursec], split=" ")[[1]][part])]}))})
#   lowermllseq <- 0.5*sum(sapply(1:nparts, function(part) {sum(sizes[[part]]*log(lambdag[[part]]))})) - 
#     0.5*lambda2*sum(sum1*lambdamult)
#   # lowermllseq <- sum(sizes*log(lambdagold)) - 0.5*lambda2*sum(lambdagold*sum1)
#   niter2seq <- numeric(0)
#   
#   # outer loop of algorithm:
#   conv1 <- FALSE
#   iter1 <- 0
#   if(trace) {cat("\n", "Estimating penalty multipliers by empirical Bayes", "\n", sep="")}
#   srt <- proc.time()[3]
#   
#   while(!conv1 & (iter1 < maxiter)) {
#     
#     iter1 <- iter1 + 1
#     
#     # estimating new lambdag
#     # if(length(partitions)==1) {
#     #   # local_opts <- list(algorithm="NLOPT_LD_MMA", xtol_rel= 1.0e-7)
#     #   # opts <- list(algorithm="NLOPT_LD_AUGLAG", xtol_rel=1.0e-7, maxeval=1000,
#     #   #              local_opts=local_opts)
#     #   # opt <- nloptr(x0=log(lambdag), eval_f=grlowermll, eval_g_eq=grconstr,
#     #   #               opts=opts, lambda2=lambda2, sizes=as.numeric(sizes), sum1=sum1)
#     #   opt <- optim(par=c(s, log(unlist(lambdag, use.names=FALSE))), fn=grmagn1, lambda2=lambda2, 
#     #                sizes=sizes[[1]], sum1=sum1, method="BFGS", control=list(maxit=5000))
#     #   s <- opt$par[1]
#     #   lambdagnew <- list(exp(opt$par[-1]))
#     #   lambdagseq <- sapply(1:nparts, function(part) {cbind(lambdagseq[[part]], lambdagnew[[part]])}, simplify=FALSE)
#     #   names(lambdagseq) <- partnames  
#     # } else 
#     if(monotone) {
#       # opt <- optim(par=c(s, log(unlist(alphag, use.names=FALSE))), fn=grmagn_mono1, 
#       #              lambda2=lambda2, sizes=sizes[[1]], sum1=sum1[[1]], sum2=sum2[[1]], method="L-BFGS-B", 
#       #              lower=c(-Inf, -Inf, rep(0, unlist(G) - 1)),
#       #              control=list(maxit=5000))
#       unlist(alphag, use.names=FALSE)
#       opt <- nloptr(log(seq(1, 1.1, length.out=6)), eval_f=eval_f1, 
#                     eval_grad_f=eval_grad_f1, eval_g_eq=eval_g1, eval_jac_g_eq=eval_jac_g1, 
#                     lb=c(-Inf, rep(0, 5)), lambda2=lambda2, sizes=sizes[[1]], sum1=sum1[[1]], 
#                     sum2=sum2[[1]], 
#                     opts=list(algorithm="NLOPT_LD_AUGLAG", xtol_rel=sqrt(.Machine$double.eps), 
#                               print_level=2, maxeval=1000, 
#                               local_opts=list(algorithm="NLOPT_LD_MMA", 
#                                               xtol_rel=sqrt(.Machine$double.eps))))
#       opt$solution
#       alpha <- list(exp(opt$solution))
#       (lambdagnew <- list(cumprod(unlist(alpha, use.names=FALSE))))
#       lambdagseq <- sapply(1:nparts, function(part) {
#         cbind(lambdagseq[[part]], lambdagnew[[part]])}, simplify=FALSE)
#       names(lambdagseq) <- partnames
#     } else {
#       # opt <- optim(par=c(s, log(unlist(lambdag, use.names=FALSE))), fn=grmagn2, lambda2=lambda2, nparts=nparts, 
#       #              sizes=sizes, G=G, uintsec=uintsec, sum1=sum1, method="Nelder-Mead", control=list(maxit=5000))
#       # s <- opt$par[1:nparts]
#       # lambdagnew <- split(exp(opt$par[-c(1:nparts)]), factor(rep(partnames, unlist(G)), 
#       #                                                        levels=partnames))
#       # lambdagseq <- sapply(1:nparts, function(part) {cbind(lambdagseq[[part]], lambdagnew[[part]])}, simplify=FALSE)
#       # names(lambdagseq) <- partnames
#       
#       opt <- optim(par=c(s, log(unlist(lambdag, use.names=FALSE))), fn=grmagn3, lambda2=lambda2, nparts=nparts, 
#                    partsind=partsind, partsmat=partsmat, sizes=unlist(sizes), G=G, sum1=sum1, 
#                    method="Nelder-Mead", control=list(maxit=5000))
#       s <- opt$par[1]
#       lambdagnew <- split(exp(opt$par[-1]), factor(rep(partnames, unlist(G)), levels=partnames))
#       # s <- opt$par[1:nparts]
#       # lambdagnew <- split(exp(opt$par[-c(1:nparts)]), factor(rep(partnames, unlist(G)), levels=partnames))
#       lambdagseq <- sapply(1:nparts, function(part) {cbind(lambdagseq[[part]], lambdagnew[[part]])}, simplify=FALSE)
#       names(lambdagseq) <- partnames
#     }
#     
#     # inner loop of algorithm:
#     conv2 <- 0
#     iter2 <- 0
#     while(!conv2 & (iter2 < maxiter)) {
#       iter2 <- iter2 + 1
#       
#       # estimating new model parameters
#       lambdamultvec <- apply(sapply(1:nparts, function(part) {lambdagnew[[part]][partitions[[part]]]}), 1, prod)
#       newparam <- est_param(x, kappa, m, n, p, ci, rep(phi, p), chi, lambdamultvec*lambda2, intercept)
#       
#       sigma <- newparam$sigma
#       mu <- newparam$mu
#       ci <- newparam$ci
#       chi <- newparam$chi
#       
#       # checking convergence of inner loop
#       conv2 <- max(c(abs((mu - muold)/muold), abs(diag((sigma - sigmaold)/sigmaold)))) < eps
#       
#       muold <- mu
#       sigmaold <- sigma
#     }
#     
#     niter2seq <- c(niter2seq, iter2)
#     
#     # sum is needed in optimisation routine
#     sum1 <- sapply(1:length(uintsec), function(cursec) {
#       ind <- which(intsec==uintsec[[cursec]]) + intercept;
#       sum((diag(sigmaold)[ind] + muold[ind]^2)*(1 + sqrt(phi/chi[ind - intercept])))})
#     
#     # keeping track of lower bound on marginal log likelihood
#     lambdamult <- exp(rowSums(sapply(1:nparts, function(part) {
#       log(unlist(lambdagnew))[partsind==part][partsmat[, part]]})))
#     lowermllseq <- c(lowermllseq, 0.5*sum(sapply(1:nparts, function(part) {sum(sizes[[part]]*log(lambdag[[part]]))})) - 
#                        0.5*lambda2*sum(sum1*lambdamult))
#     # lowermllseq <- c(lowermllseq, sum(sizes*log(lambdag)) - 0.5*lambda2*sum(lambdag*sum1))
#     
#     # checking convergence of outer loop:
#     conv1 <- max(abs((unlist(lambdagnew) - unlist(lambdag))/unlist(lambdag))) < eps
#     
#     # updating lambdag for new iteration
#     lambdagold <- lambdag
#     lambdag <- lambdagnew
#     
#     # printing progress
#     if(trace) {cat("\r", "Penalty multipliers estimated at ", paste(partnames, lapply(lambdag, function(part) {
#       paste(round(part, 2), collapse=", ")}), sep=": ", collapse=" and "), "      ", sep="")}
#     # cat("\r", "Penalty multipliers estimated at ", paste(round(lambdag[-G], 2), collapse=", "), 
#     #     " and ", round(lambdag[G], 2), "      ", sep="")
#     
#   }
#   eb.time <- proc.time()[3] - srt
#   if(trace) {cat("\r", "Penalty multipliers estimated at ", paste(partnames, lapply(lambdag, function(part) {
#     paste(round(part, 2), collapse=", ")}), sep=": ", collapse=" and "), " in ", 
#     round(eb.time, 2), " seconds ", sep="")}
#   
#   out <- list(mu=mu, sigma=sigma, ci=ci, chi=chi, lambda1=lambda1, lambda2=lambda2, 
#               lambdag=lambdagseq, lowermll=lowermllseq, nouteriter=iter1, ninneriter=niter2seq, 
#               conv=conv1)
#   
#   return(out)
#   
# }

# grVBEM2 <- function(x, y, m, partitions, lambda1=NULL, lambda2=NULL, intercept, monotone, 
#                     iso=NULL, posterior, ELBO, eps, maxiter, trace=TRUE, alphastart=NULL) {
#   
#   if(!is.list(partitions) & is.null(names(partitions))) {
#     partitions <- list(partition1=partitions)
#   } else if(is.null(names(partitions))) {
#     names(partitions) <- paste("partition", 1:length(partitions), sep="")
#   }
#   
#   partnames <- names(partitions)
#   
#   # assigning fixed (throughout algorithm) variables
#   nparts <- length(partitions)
#   # if(!monotone) {
#   #   sizes <- lapply(partitions, function(part) {rle(sort(part))$lengths})
#   #   G <- lapply(partitions, function(part) {length(unique(part))})
#   # }
#   sizes <- lapply(partitions, function(part) {rle(sort(part))$lengths})
#   G <- lapply(partitions, function(part) {length(unique(part))})
#   p <- ncol(x)
#   n <- nrow(x)
#   kappa <- y - m/2
#   
#   # if no penalty parameters are given we estimate them by cross-validation
#   if(is.null(lambda1) | is.null(lambda2)) {
#     if(trace) {cat("\r", "Estimating global lambda1 and lambda2 by cross-validation", sep="")}
#     srt <- proc.time()[3]
#     opt.glob <- cv.pen(x, y, intercept)
#     cv.time <- proc.time()[3] - srt
#     lambda1 <- opt.glob$lambda1bayes[which.min(opt.glob$cvll)]
#     lambda2 <- opt.glob$lambda2bayes[which.min(opt.glob$cvll)]
#     if(trace) {cat("\n", "Global lambda1 and lambda2 estimated at ", round(lambda1, 2), " and ", 
#                    round(lambda2, 2), " in ", round(cv.time, 2), " seconds", sep="")}
#   }
#   
#   # in the multiplier setting phi does not change
#   phi <- 0.25*lambda1^2/lambda2
#   
#   # starting values for lambdag, lagrange multiplier s and possible approximate hessian B
#   # if(!monotone) {
#   #   lambdagnew <- lambdag <- lambdagold <- lambdagseq <- lapply(G, function(gpart) {rep(1, gpart)})
#   # } else {
#   #   lambdagnew <- lambdag <- lambdagold <- lambdagseq <- rep(list(numeric(p)), nparts)
#   # }
#   if(!monotone) {
#     lambdagnew <- lambdag <- lambdagold <- lambdagseq <- lapply(G, function(gpart) {rep(1, gpart)})
#   } else {
#     lambdagnew <- lambdag <- lambdagold <- lambdagseq <- sapply(1:nparts, function(part) {
#       logseq <- log(seq(0.9, 1.1, length.out=G[[part]])); return(exp(logseq - sum(sizes[[part]]*logseq)/p))},
#       simplify=FALSE)
#   }
#   # s <- rep(0, times=nparts)
#   if(monotone) {
#     # sum2 <- lapply(sizes, function(sizes) {rev(cumsum(rev(sizes)))})
#     # alphag <- sapply(1:nparts, function(part) {
#     #   alphag <- seq(1, 1.01, length.out=G[[part]]); 
#     #   return(alphag*exp(-sum(sum2[[part]]*log(alphag))/sum(sum2[[part]])))})
#     # alphag <- lapply(G, function(G) {return(c(1, seq(1.01, 3, length.out=G - 1)))})
#     alpha <- list(alphastart)
#   } else {
#     s <- 0
#   }
#   s <- 0
#   
#   # starting values for mu and sigma
#   fit.pen <- penalized(y, x, unpenalized=formula(ifelse(intercept, "~1", "~0")), 
#                        model="logistic", lambda1=0, 
#                        lambda2=2*(lambda1 + lambda2), trace=FALSE)
#   #muold <- coef(fit.pen, which="all")
#   if(intercept) {
#     xadj <- cbind(1, x)
#   } else {
#     xadj <- x
#   }
#   b0 <- coef(fit.pen, which="all")
#   pred0 <- as.numeric(exp(xadj %*% b0)/(1 + exp(xadj %*% b0)))
#   w <- sqrt(pred0*(1 - pred0))
#   xw <- xadj*w
#   svdxw <- svd(xw)
#   d <- svdxw$d
#   v <- svdxw$v
#   invmat <- d^2/(d^2 + 4*(lambda1 + lambda2))^2
#   sigmaold <- t(t(v)*invmat) %*% t(v)
#   
#   # rest of the starting values follow from that
#   muold <- as.numeric(sigmaold %*% (t(xadj) %*% as.matrix(kappa)))
#   ci <- as.numeric(sqrt(colSums(t(xadj) * (sigmaold %*% t(xadj))) + (colSums(t(xadj)*muold))^2))
#   chi <- as.numeric(0.5*(lambda1 + lambda2)*(diag(sigmaold) + muold^2))[(intercept + 1):(intercept + p)]
#   dsigmaold <- diag(sigmaold)
#   
#   # sum is needed in optimisation routine
#   intsec <- do.call("paste", c(partitions, sep=" "))
#   uintsec <- unique(intsec)
#   sum1 <- sapply(1:length(uintsec), function(cursec) {
#     ind <- which(intsec==uintsec[[cursec]]) + intercept;
#     sum((dsigmaold[ind] + muold[ind]^2)*(1 + sqrt(phi/chi[ind - intercept])))})
#   intsizes <- as.numeric(t(table(as.data.frame(matrix(unlist(partitions), ncol=nparts)))))
#   partsmat <- unique(matrix(unlist(partitions), ncol=nparts))
#   partsind <- rep(1:nparts, times=unlist(G))
#   # if(!monotone) {
#   #   intsec <- do.call("paste", c(partitions, sep=" "))
#   #   uintsec <- unique(intsec)
#   #   sum1 <- sapply(1:length(uintsec), function(cursec) {
#   #     ind <- which(intsec==uintsec[[cursec]]) + intercept;
#   #     sum((dsigmaold[ind] + muold[ind]^2)*(1 + sqrt(phi/chi[ind - intercept])))})
#   #   intsizes <- as.numeric(t(table(as.data.frame(matrix(unlist(partitions), ncol=nparts)))))
#   #   partsmat <- unique(matrix(unlist(partitions), ncol=nparts))
#   #   partsind <- rep(1:nparts, times=unlist(G))
#   # } else {
#   #  sum1 <- (dsigmaold[(intercept + 1):(intercept + p)] + 
#   #             muold[(intercept + 1):(intercept + p)]^2)*(1 + sqrt(phi/chi))
#   # }
#   
#   # remove intercept for iterations
#   if(intercept) {
#     muold <- muold[-1]
#     dsigmaold <- dsigmaold[-1]
#   }
#   
#   # keeping track of things:
#   # lambdamult <- exp(rowSums(sapply(1:nparts, function(part) {
#   #   log(unlist(lambdag))[partsind==part][partsmat[, part]]})))
#   # lambdagprod <- sapply(1:length(uintsec), function(cursec) {prod(sapply(1:nparts, function(part) {
#   #   lambdag[[part]][as.numeric(strsplit(uintsec[cursec], split=" ")[[1]][part])]}))})
#   # lowermllseq <- 0.5*sum(sapply(1:nparts, function(part) {sum(sizes[[part]]*log(lambdag[[part]]))})) - 
#   #   0.5*lambda2*sum(sum1*lambdamult)
#   # lowermllseq <- sum(sizes*log(lambdagold)) - 0.5*lambda2*sum(lambdagold*sum1)
#   # niter2seq <- numeric(0)
#   if(ELBO) {
#     ELBOseq <- numeric(0)
#     ELBOold <- Inf
#   }
#   niter2seq <- numeric(0)
#   
#   # outer loop of algorithm:
#   conv1 <- FALSE
#   iter1 <- 0
#   if(trace) {cat("\n", "Estimating penalty multipliers by empirical Bayes", "\n", sep="")}
#   srt <- proc.time()[3]
#   while(!conv1 & (iter1 < maxiter)) {
#     
#     iter1 <- iter1 + 1
#     
#     # estimating new lambdag
#     # if(length(partitions)==1) {
#     #   # local_opts <- list(algorithm="NLOPT_LD_MMA", xtol_rel= 1.0e-7)
#     #   # opts <- list(algorithm="NLOPT_LD_AUGLAG", xtol_rel=1.0e-7, maxeval=1000,
#     #   #              local_opts=local_opts)
#     #   # opt <- nloptr(x0=log(lambdag), eval_f=grlowermll, eval_g_eq=grconstr,
#     #   #               opts=opts, lambda2=lambda2, sizes=as.numeric(sizes), sum1=sum1)
#     #   opt <- optim(par=c(s, log(unlist(lambdag, use.names=FALSE))), fn=grmagn1, lambda2=lambda2, 
#     #                sizes=sizes[[1]], sum1=sum1, method="BFGS", control=list(maxit=5000))
#     #   s <- opt$par[1]
#     #   lambdagnew <- list(exp(opt$par[-1]))
#     #   lambdagseq <- sapply(1:nparts, function(part) {cbind(lambdagseq[[part]], lambdagnew[[part]])}, simplify=FALSE)
#     #   names(lambdagseq) <- partnames } else
#     if(monotone) {
#       # opt <- nloptr(log(unlist(alphag)), eval_f=eval_f_list,
#       #               eval_g_eq=eval_g1, eval_jac_g_eq=eval_jac_g1, 
#       #               lb=c(-Inf, rep(0, G[[1]] - 1)), lambda2=lambda2, sizes=sizes[[1]], 
#       #               sum1=sum1[[1]], sum2=sum2[[1]], 
#       #               opts=list(algorithm="NLOPT_LD_AUGLAG", xtol_abs=sqrt(.Machine$double.eps), 
#       #                         print_level=0, maxeval=1000, 
#       #                         local_opts=list(algorithm="NLOPT_LD_MMA", 
#       #                                         xtol_abs=sqrt(.Machine$double.eps))))
#       # alpha <- list(exp(opt$solution))
#       # lambdagnew <- list(cumprod(unlist(alpha, use.names=FALSE)))
#       # lambdagseq <- sapply(1:nparts, function(part) {
#       #   cbind(lambdagseq[[part]], lambdagnew[[part]])}, simplify=FALSE)
#       # names(lambdagseq) <- partnames
#       
#       # opt <- nloptr(unlist(alpha), eval_f=eval_f_bc,
#       #               eval_g_eq=eval_eq_bc, eval_jac_g_eq=eval_jac_eq_bc, 
#       #               lb=c(-Inf, max(codata)), cj=codata, lambda2=lambda2, p=p, sum1=sum1,
#       #               opts=list(algorithm="NLOPT_LD_AUGLAG", xtol_abs=sqrt(.Machine$double.eps), 
#       #                         print_level=0, maxeval=1000, 
#       #                         local_opts=list(algorithm="NLOPT_LD_MMA", 
#       #                                         xtol_abs=sqrt(.Machine$double.eps))))
#       # alpha <- list(opt$solution)
#       # if(unlist(alpha)[1]!=0) {
#       #   lambdagnew <- list(((codata + unlist(alpha)[2])^unlist(alpha)[1] - 1)/unlist(alpha)[1])
#       # } else {
#       #   lambdagnew <- list(log(codata + unlist(alpha[2])))
#       # }
#       
#       opt <- nloptr(log(unlist(lambdag, use.names=FALSE)), eval_f=eval_f_list,
#                     eval_g_eq=eval_g_eq_list, eval_g_ineq=eval_g_ineq_list,
#                     lambda2=lambda2, sum1=sum1, sizes=sizes[[1]],
#                     opts=list(algorithm="NLOPT_LD_AUGLAG", xtol_abs=sqrt(.Machine$double.eps),
#                               print_level=0, maxeval=5000,
#                               local_opts=list(algorithm="NLOPT_LD_MMA",
#                                               xtol_abs=sqrt(.Machine$double.eps))))
#       lambdagnew <- split(exp(opt$solution), factor(rep(partnames, unlist(G)), levels=partnames))
#       # opt <- optim(par=c(s, log(unlist(lambdag, use.names=FALSE))), fn=grmagn4, lambda2=lambda2, 
#       #              sizes=unlist(sizes), sum1=sum1, method="Nelder-Mead", 
#       #              control=list(maxit=5000, reltol=sqrt(.Machine$double.eps)))
#       # s <- opt$par[1]
#       # lambdagnew <- split(exp(opt$par[-1]), factor(rep(partnames, unlist(G)), levels=partnames))
#       lambdagseq <- sapply(1:nparts, function(part) {
#         cbind(lambdagseq[[part]], lambdagnew[[part]])}, simplify=FALSE)
#       names(lambdagseq) <- partnames
#       
#     } else {
#       # opt <- optim(par=c(s, log(unlist(lambdag, use.names=FALSE))), fn=grmagn2, lambda2=lambda2, nparts=nparts, 
#       #              sizes=sizes, G=G, uintsec=uintsec, sum1=sum1, method="Nelder-Mead", control=list(maxit=5000))
#       # s <- opt$par[1:nparts]
#       # lambdagnew <- split(exp(opt$par[-c(1:nparts)]), factor(rep(partnames, unlist(G)), 
#       #                                                        levels=partnames))
#       # lambdagseq <- sapply(1:nparts, function(part) {cbind(lambdagseq[[part]], lambdagnew[[part]])}, simplify=FALSE)
#       # names(lambdagseq) <- partnames
#       
#       opt <- optim(par=c(s, log(unlist(lambdag, use.names=FALSE))), fn=grmagn3, lambda2=lambda2, 
#                    nparts=nparts, 
#                    partsind=partsind, partsmat=partsmat, sizes=unlist(sizes), G=G, sum1=sum1, 
#                    method="Nelder-Mead", control=list(maxit=5000))
#       s <- opt$par[1]
#       lambdagnew <- split(exp(opt$par[-1]), factor(rep(partnames, unlist(G)), levels=partnames))
#       if(!is.null(iso)) { 
#         lambdagnew <- sapply(1:nparts, function(part) {
#           if(iso[[part]]) {
#             fit.iso <- biviso(cbind(lambdagnew[[part]], c(1:G[[part]]), sizes[[part]]));
#             curlambdag <- as.numeric(exp(log(fit.iso[, 1]) - sum(sizes[[part]]*log(fit.iso[, 1]))/
#                                            sum(sizes[[part]])))
#           } else {
#             curlambdag <- lambdagnew[[part]]
#           }
#           return(curlambdag)}, simplify=FALSE)
#       }
#       # s <- opt$par[1:nparts]
#       # lambdagnew <- split(exp(opt$par[-c(1:nparts)]), factor(rep(partnames, unlist(G)), levels=partnames))
#       lambdagseq <- sapply(1:nparts, function(part) {cbind(lambdagseq[[part]], lambdagnew[[part]])}, simplify=FALSE)
#       names(lambdagseq) <- partnames
#     }
#     lambdamult <- exp(rowSums(sapply(1:nparts, function(part) {
#       log(unlist(lambdagnew, use.names=FALSE))[partsind==part][partsmat[, part]]})))
#     
#     # inner loop of algorithm:
#     conv2 <- FALSE
#     iter2 <- 0
#     while(!conv2 & (iter2 < maxiter)) {
#       iter2 <- iter2 + 1
#       
#       ciold  <- ci
#       chiold <- chi
#       
#       # estimating new model parameters
#       lambdamultvec <- apply(sapply(1:nparts, function(part) {lambdagnew[[part]][partitions[[part]]]}), 1, prod)
#       
#       if(ELBO) {
#         newparam <- est_param(x, kappa, m, n, p, ci, phi, chi, lambdamultvec*lambda2, intercept)
#         
#         sigma <- newparam$sigma
#         dsigma <- diag(sigma)
#         mu <- as.numeric(newparam$mu)
#         ci <- as.numeric(newparam$ci)
#         chi <- as.numeric(newparam$chi)
#         
#         # calculation of ELBO
#         tempmat <- t(t(x)/(lambdamultvec + lambdamultvec*sqrt(phi/chiold))) %*% t(x)/lambda2
#         diag(tempmat) <- diag(tempmat) + 2*ciold/(m*tanh(ciold/2))
#         ldetsigma <- -determinant(tempmat, logarithm=TRUE)$modulus + sum(log(ciold)) -
#           sum(log(tanh(0.5*ciold))) - sum(log(lambdamultvec)) - 
#           sum(log(1 + 0.5*lambda1/sqrt(chiold*lambda2)))
#         ELBO <- sum((y - 0.5*m)*as.numeric(xadj %*% mu)) + sum(m*(0.5*ci - log(exp(ci) + 1))) +
#           0.25*sum(m*(ci - (rowSums((xadj %*% sigma)*xadj) + as.numeric(xadj %*% mu)^2)/ci)*tanh(0.5*ci)) +
#           0.5*ldetsigma + 0.5*sum(log(lambdamultvec)) - 0.25*lambda1*sum(sqrt(chi))/sqrt(lambda2) -
#           lambda2^(3/2)*sum(lambdamultvec*(0.5*lambda1/sqrt(lambda2) + 1/sqrt(chi))*
#                               (dsigma[(1 + intercept):(p + intercept)] + 
#                                  mu[(1 + intercept):(p + intercept)]^2))/lambda1
#         ELBOseq <- c(ELBOseq, ELBO)
#         
#         # checking convergence of inner loop
#         conv2 <- (max(c(abs((mu - muold)/ifelse(muold==0, muold + 0.00001, muold)), 
#                         abs((dsigma - dsigmaold)/ifelse(dsigmaold==0, dsigmaold + 0.00001, dsigmaold)))) < eps) |
#           abs(ELBO - ELBOold) < max(eps/1000, sqrt(.Machine$double.eps))
#         
#         ELBOold <- ELBO
#         muold <- mu
#         dsigmaold <- dsigma
#         
#       } else {
#         newparam <- est_param2(x, kappa, m, n, p, ci, phi, chi, lambdamultvec*lambda2, intercept=FALSE)
#         
#         dsigma <- as.numeric(newparam$dsigma)
#         mu <- as.numeric(newparam$mu)
#         ci <- as.numeric(newparam$ci)
#         chi <- as.numeric(newparam$chi)
#       
#         # checking convergence of inner loop
#         conv2 <- max(c(abs((mu - muold)/ifelse(muold==0, muold + 0.00001, muold)), 
#                        abs((dsigma - dsigmaold)/ifelse(dsigmaold==0, dsigmaold + 0.00001, dsigmaold)))) < eps
#       
#         muold <- mu
#         dsigmaold <- dsigma
#       
#       }
#     }
# 
#     niter2seq <- c(niter2seq, iter2)
#     
#     # sum is needed in optimisation routine
#     sum1 <- sapply(1:length(uintsec), function(cursec) {
#       ind <- which(intsec==uintsec[[cursec]]);
#       sum((dsigmaold[ind] + muold[ind]^2)*(1 + sqrt(phi/chi[ind])))})
#     
#     # keeping track of lower bound
#     
#     # lowermllseq <- c(lowermllseq, 0.5*sum(sapply(1:nparts, function(part) {sum(sizes[[part]]*log(lambdag[[part]]))})) -
#     #                    0.5*lambda2*sum(sum1*lambdamult))
#     # lowermllseq <- c(lowermllseq, sum(sizes*log(lambdag)) - 0.5*lambda2*sum(lambdag*sum1))
#     
#     
#     # checking convergence of outer loop:
#     conv1 <- max(abs((unlist(lambdagnew) - unlist(lambdag))/unlist(lambdag))) < eps
#     
#     # updating lambdag for new iteration
#     lambdagold <- lambdag
#     lambdag <- lambdagnew
#     
#     # printing progress
#     # if(trace & !monotone) {
#     #   cat("\r", "Penalty multipliers estimated at ", paste(partnames, lapply(lambdag, function(part) {
#     #   paste(round(part, 2), collapse=", ")}), sep=": ", collapse=" and "), "      ", sep="")
#     # } else if(trace) {
#     #   cat("\r", "Penalty model parameters estimated at ", 
#     #       paste(partnames, lapply(alpha, function(part) {
#     #         paste(round(part, 2), collapse=", ")}), sep=": ", collapse=" and "), "      ", sep="")
#     # }
#     if(trace) {
#       cat("\r", "Penalty multipliers estimated at ", paste(partnames, lapply(lambdag, function(part) {
#         paste(round(part, 2), collapse=", ")}), sep=": ", collapse=" and "), "      ", sep="")
#     }
#     
#   }
#   
#   if(posterior & !ELBO) {
#     newparam <- est_param(x, kappa, m, n, p, ci, rep(phi, p), chi, lambdamultvec*lambda2, intercept)
#     dsigma <- newparam$sigma
#     mu <- newparam$mu
#     ci <- newparam$ci
#     chi <- newparam$chi
#   }
#   eb.time <- proc.time()[3] - srt
#   if(trace) {cat("\r", "Penalty multipliers estimated at ", paste(partnames, lapply(lambdag, function(part) {
#     paste(round(part, 2), collapse=", ")}), sep=": ", collapse=" and "), " in ", 
#     round(eb.time, 2), " seconds ", sep="")}
#   
#   if(ELBO) {
#     out <- list(mu=mu, sigma=sigma, ci=ci, chi=chi, lambda1=lambda1, lambda2=lambda2, 
#                 lambdag=lambdagseq, ELBO=ELBOseq, nouteriter=iter1, ninneriter=niter2seq, 
#                 conv=conv1)
#   } else {
#     out <- list(mu=mu, sigma=dsigma, ci=ci, chi=chi, lambda1=lambda1, lambda2=lambda2, 
#                 lambdag=lambdagseq, nouteriter=iter1, ninneriter=niter2seq, 
#                 conv=conv1)
#   }
#   return(out)
#   
# }




# # monotonoicity by transforming the lambdas to alphas
# grVBEM2 <- function(x, y, m, partitions, lambda1=NULL, lambda2=NULL, intercept, monotone, 
#                     posterior, eps, maxiter, trace=TRUE, alphastart) {
#   
#   if(!is.list(partitions) & is.null(names(partitions))) {
#     partitions <- list(partition1=partitions)
#   } else if(is.null(names(partitions))) {
#     names(partitions) <- paste("partition", 1:length(partitions), sep="")
#   }
#   
#   partnames <- names(partitions)
#   
#   # assigning fixed (throughout algorithm) variables
#   nparts <- length(partitions)
#   sizes <- lapply(partitions, function(part) {rle(sort(part))$lengths})
#   G <- lapply(partitions, function(part) {length(unique(part))})
#   p <- ncol(x)
#   n <- nrow(x)
#   kappa <- y - m/2
#   
#   # if no penalty parameters are given we estimate them by cross-validation
#   if(is.null(lambda1) | is.null(lambda2)) {
#     if(trace) {cat("\r", "Estimating global lambda1 and lambda2 by cross-validation", sep="")}
#     srt <- proc.time()[3]
#     opt.glob <- cv.pen(x, y, intercept)
#     cv.time <- proc.time()[3] - srt
#     lambda1 <- opt.glob$lambda1[which.min(opt.glob$cvll)]
#     lambda2 <- opt.glob$lambda2[which.min(opt.glob$cvll)]
#     if(trace) {cat("\n", "Global lambda1 and lambda2 estimated at ", round(lambda1, 2), " and ", 
#                    round(lambda2, 2), " in ", round(cv.time, 2), " seconds", sep="")}
#   }
#   
#   # in the multiplier parametrization phi does not change
#   phi <- 0.25*lambda1^2/lambda2
#   
#   # starting values for lambdag, lagrange multiplier s and possible approximate hessian B
#   s <- 0
#   if(!monotone) {
#     lambdagnew <- lambdag <- lambdagold <- lambdagseq <- lapply(G, function(gpart) {
#       rep(1, gpart)})
#   } else {
#     # alphag <- sapply(1:nparts, function(part) {
#     #   alphag <- seq(1, 1.01, length.out=G[[part]]); 
#     #   return(alphag*exp(-sum(sum2[[part]]*log(alphag))/sum(sum2[[part]])))})
#     # alphag <- lapply(G, function(G) {return(c(1, seq(1.01, 3, length.out=G - 1)))})
#     alphag <- structure(list(alphastart), names=partnames)
#     lambdagnew <- lambdag <- lambdagold <- lambdagseq <- 
#       structure(list(exp(cumsum(log(unlist(alphag, use.names=FALSE))))), names=partnames)
#     sum2 <- lapply(sizes, function(sizes) {rev(cumsum(rev(sizes)))})
#   }
#   
#   # starting values for mu and sigma
#   fit.pen <- glmnet(x, y, family="binomial", intercept=intercept, alpha=0, standardize=FALSE,
#                     lambda=(0.5*lambda1 + lambda2)/n)
#   if(intercept) {
#     xadj <- cbind(1, x)
#   } else {
#     xadj <- x
#   }
#   b0 <- coef(fit.pen)
#   pred0 <- as.numeric(exp(xadj %*% b0)/(1 + exp(xadj %*% b0)))
#   w <- sqrt(pred0*(1 - pred0))
#   xw <- xadj*w
#   svdxw <- svd(xw)
#   d <- svdxw$d
#   v <- svdxw$v
#   invmat <- d^2/(d^2 + 4*(lambda1 + lambda2))^2
#   sigmaold <- t(t(v)*invmat) %*% t(v)
#   
#   # rest of the starting values follow from that
#   muold <- as.numeric(sigmaold %*% (t(xadj) %*% as.matrix(kappa)))
#   ci <- as.numeric(sqrt(colSums(t(xadj) * (sigmaold %*% t(xadj))) + (colSums(t(xadj)*muold))^2))
#   chi <- as.numeric(0.5*(lambda1 + lambda2)*(diag(sigmaold) + muold^2))[(intercept + 1):(intercept + p)]
#   dsigmaold <- diag(sigmaold)
#   
#   # sum is needed in optimisation routine
#   intsec <- do.call("paste", c(partitions, sep=" "))
#   uintsec <- unique(intsec)
#   sum1 <- sapply(1:length(uintsec), function(cursec) {
#     ind <- which(intsec==uintsec[[cursec]]) + intercept;
#     sum((dsigmaold[ind] + muold[ind]^2)*(1 + sqrt(phi/chi[ind - intercept])))})
#   intsizes <- as.numeric(t(table(as.data.frame(matrix(unlist(partitions), ncol=nparts)))))
#   partsmat <- unique(matrix(unlist(partitions), ncol=nparts))
#   partsind <- rep(1:nparts, times=unlist(G))
#   
#   # remove intercept for iterations
#   if(intercept) {
#     muold <- muold[-1]
#     dsigmaold <- dsigmaold[-1]
#   }
#   
#   # keeping track of things:
#   lambdamult <- exp(rowSums(sapply(1:nparts, function(part) {
#     log(unlist(lambdag, use.names=FALSE))[partsind==part][partsmat[, part]]})))
#   # lambdagprod <- sapply(1:length(uintsec), function(cursec) {prod(sapply(1:nparts, function(part) {
#   #   lambdag[[part]][as.numeric(strsplit(uintsec[cursec], split=" ")[[1]][part])]}))})
#   lowermllseq <- 0.5*sum(sapply(1:nparts, function(part) {sum(sizes[[part]]*log(lambdag[[part]]))})) - 
#     0.5*lambda2*sum(sum1*lambdamult)
#   # lowermllseq <- sum(sizes*log(lambdagold)) - 0.5*lambda2*sum(lambdagold*sum1)
#   niter2seq <- numeric(0)
#   
#   # outer loop of algorithm:
#   conv1 <- FALSE
#   iter1 <- 0
#   if(trace) {cat("\n", "Estimating penalty multipliers by empirical Bayes", "\n", sep="")}
#   srt <- proc.time()[3]
#   while(!conv1 & (iter1 < maxiter)) {
#     
#     iter1 <- iter1 + 1
#     
#     # estimating new lambdag
#     # if(length(partitions)==1) {
#     #   # local_opts <- list(algorithm="NLOPT_LD_MMA", xtol_rel= 1.0e-7)
#     #   # opts <- list(algorithm="NLOPT_LD_AUGLAG", xtol_rel=1.0e-7, maxeval=1000,
#     #   #              local_opts=local_opts)
#     #   # opt <- nloptr(x0=log(lambdag), eval_f=grlowermll, eval_g_eq=grconstr,
#     #   #               opts=opts, lambda2=lambda2, sizes=as.numeric(sizes), sum1=sum1)
#     #   opt <- optim(par=c(s, log(unlist(lambdag, use.names=FALSE))), fn=grmagn1, lambda2=lambda2, 
#     #                sizes=sizes[[1]], sum1=sum1, method="BFGS", control=list(maxit=5000))
#     #   s <- opt$par[1]
#     #   lambdagnew <- list(exp(opt$par[-1]))
#     #   lambdagseq <- sapply(1:nparts, function(part) {cbind(lambdagseq[[part]], lambdagnew[[part]])}, simplify=FALSE)
#     #   names(lambdagseq) <- partnames } else
#     if(monotone) {
#       # opt <- optim(par=c(s, log(unlist(alphag, use.names=FALSE))), fn=grmagn4, lambda2=lambda2, 
#       #              sizes=sizes[[1]], sum1=sum1[[1]], sum2=sum2[[1]],
#       #              method="L-BFGS-B", control=list(maxit=5000), 
#       #              lower=c(rep(-Inf, 2), rep(0, G[[1]] - 1)))
#       opt <- optim(par=c(s, log(unlist(alphag, use.names=FALSE))), fn=grmagn4, lambda2=lambda2, 
#                    sizes=sizes[[1]], sum1=sum1[[1]], sum2=sum2[[1]],
#                    method="Nelder-Mead", control=list(maxit=5000))
#       s <- opt$par[1]
#       alphag <- structure(list(exp(opt$par[-1])), names=partnames)
#       # opt <- nloptr(c(s, log(unlist(alphag))), eval_f=eval_f_list, eval_g_eq=eval_g_list,
#       #               lb=c(-Inf, rep(0, G[[1]] - 1)), lambda2=lambda2, sizes=sizes[[1]],
#       #               sum1=sum1[[1]], sum2=sum2[[1]],
#       #               opts=list(algorithm="NLOPT_LD_AUGLAG", xtol_abs=sqrt(.Machine$double.eps),
#       #                         print_level=0, maxeval=1000,
#       #                         local_opts=list(algorithm="NLOPT_LD_MMA",
#       #                                         xtol_abs=sqrt(.Machine$double.eps))))
#       # opt <- nloptr(log(unlist(alphag)), eval_f=eval_f_list, eval_g_eq=eval_g_list,
#       #               lb=c(-Inf, rep(0, G[[1]] - 1)), lambda2=lambda2, sizes=sizes[[1]],
#       #               sum1=sum1[[1]], sum2=sum2[[1]],
#       #               opts=list(algorithm="NLOPT_LD_AUGLAG", xtol_abs=sqrt(.Machine$double.eps),
#       #                         print_level=0, maxeval=1000,
#       #                         local_opts=list(algorithm="NLOPT_LD_MMA",
#       #                                         xtol_abs=sqrt(.Machine$double.eps))))
#       # alphag <- structure(list(exp(opt$solution)), names=partnames)
#       lambdagnew <- structure(list(cumprod(unlist(alphag, use.names=FALSE))), names=partnames)
#       lambdagseq <- structure(sapply(1:nparts, function(part) {
#         cbind(lambdagseq[[part]], lambdagnew[[part]])}, simplify=FALSE), names=partnames)
#     } else {
#       # opt <- optim(par=c(s, log(unlist(lambdag, use.names=FALSE))), fn=grmagn2, lambda2=lambda2, nparts=nparts, 
#       #              sizes=sizes, G=G, uintsec=uintsec, sum1=sum1, method="Nelder-Mead", control=list(maxit=5000))
#       # s <- opt$par[1:nparts]
#       # lambdagnew <- split(exp(opt$par[-c(1:nparts)]), factor(rep(partnames, unlist(G)), 
#       #                                                        levels=partnames))
#       # lambdagseq <- sapply(1:nparts, function(part) {cbind(lambdagseq[[part]], lambdagnew[[part]])}, simplify=FALSE)
#       # names(lambdagseq) <- partnames
#       
#       opt <- optim(par=c(s, log(unlist(lambdag, use.names=FALSE))), fn=grmagn3, lambda2=lambda2, nparts=nparts, 
#                    partsind=partsind, partsmat=partsmat, sizes=unlist(sizes), G=G, sum1=sum1, 
#                    method="Nelder-Mead", control=list(maxit=5000))
#       s <- opt$par[1]
#       lambdagnew <- split(exp(opt$par[-1]), factor(rep(partnames, unlist(G)), levels=partnames))
#       # s <- opt$par[1:nparts]
#       # lambdagnew <- split(exp(opt$par[-c(1:nparts)]), factor(rep(partnames, unlist(G)), levels=partnames))
#       lambdagseq <- sapply(1:nparts, function(part) {cbind(lambdagseq[[part]], lambdagnew[[part]])}, simplify=FALSE)
#       names(lambdagseq) <- partnames
#     }
#     
#     # inner loop of algorithm:
#     conv2 <- 0
#     iter2 <- 0
#     while(!conv2 & (iter2 < maxiter)) {
#       iter2 <- iter2 + 1
#       
#       # estimating new model parameters
#       lambdamultvec <- apply(sapply(1:nparts, function(part) {lambdagnew[[part]][partitions[[part]]]}), 1, prod)
#       newparam <- est_param2(x, kappa, m, n, p, ci, phi, chi, lambdamultvec*lambda2, intercept)
#       
#       dsigma <- as.numeric(newparam$dsigma)
#       mu <- as.numeric(newparam$mu)
#       ci <- as.numeric(newparam$ci)
#       chi <- as.numeric(newparam$chi)
#       
#       # checking convergence of inner loop
#       conv2 <- max(c(abs((mu - muold)/muold), abs((dsigma - dsigmaold)/dsigmaold))) < eps
#       
#       muold <- mu
#       dsigmaold <- dsigma
#     }
#     
#     niter2seq <- c(niter2seq, iter2)
#     
#     # sum is needed in optimisation routine
#     sum1 <- sapply(1:length(uintsec), function(cursec) {
#       ind <- which(intsec==uintsec[[cursec]]);
#       sum((dsigmaold[ind] + muold[ind]^2)*(1 + sqrt(phi/chi[ind])))})
#     
#     # keeping track of lower bound on marginal log likelihood
#     lambdamult <- exp(rowSums(sapply(1:nparts, function(part) {
#       log(unlist(lambdagnew))[partsind==part][partsmat[, part]]})))
#     lowermllseq <- c(lowermllseq, 0.5*sum(sapply(1:nparts, function(part) {sum(sizes[[part]]*log(lambdag[[part]]))})) - 
#                        0.5*lambda2*sum(sum1*lambdamult))
# 
#     # checking convergence of outer loop:
#     conv1 <- max(abs((unlist(lambdagnew) - unlist(lambdag))/unlist(lambdag))) < eps
#     
#     # updating lambdag for new iteration
#     lambdagold <- lambdag
#     lambdag <- lambdagnew
#     
#     # printing progress
#     cat("\r", "Penalty multipliers estimated at ", paste(partnames, lapply(lambdag, function(part) {
#       paste(round(part, 2), collapse=", ")}), sep=": ", collapse=" and "), "      ", sep="")
#     
#   }
#   
#   if(posterior) {
#     newparam <- est_param(x, kappa, m, n, p, ci, rep(phi, p), chi, lambdamultvec*lambda2, intercept)
#     dsigma <- newparam$sigma
#     mu <- newparam$mu
#     ci <- newparam$ci
#     chi <- newparam$chi
#   }
#   eb.time <- proc.time()[3] - srt
#   if(trace) {cat("\r", "Penalty multipliers estimated at ", paste(partnames, lapply(lambdag, function(part) {
#     paste(round(part, 2), collapse=", ")}), sep=": ", collapse=" and "), " in ", 
#     round(eb.time, 2), " seconds ", sep="")}
#   
#   out <- list(mu=mu, sigma=dsigma, ci=ci, chi=chi, lambda1=lambda1, lambda2=lambda2, 
#               lambdag=lambdagseq, lowermll=lowermllseq, nouteriter=iter1, ninneriter=niter2seq, 
#               conv=conv1)
#   
#   return(out)
#   
# }

# # Variational Bayes estimation of posterior
# VBen <- function(x, y, m, lambda1=NULL, lambda2=NULL, intercept, eps, maxiter, trace=TRUE) {
#   
#   # data characteristics
#   p <- ncol(x)
#   n <- nrow(x)
#   kappa <- y - m/2
#   
#   # if no penalty parameters are given we estimate them by cross-validation
#   cv.time <- NULL
#   if(is.null(lambda1) | is.null(lambda2)) {
#     if(trace) {cat("\r", "Estimating lambda1 and lambda2 by cross-validation", sep="")}
#     srt <- proc.time()[3]
#     opt.glob <- cv.pen(x, y, intercept)
#     cv.time <- proc.time()[3] - srt
#     lambda1 <- opt.glob$lambda1[which.min(opt.glob$cvll)]
#     lambda2 <- opt.glob$lambda2[which.min(opt.glob$cvll)]
#     if(trace) {cat("\n", "Lambda1 and lambda2 estimated at ", round(lambda1, 2), " and ", 
#                    round(lambda2, 2), " in ", round(cv.time, 2), " seconds", sep="")}
#   }
#   
#   # fixed parameters
#   phi <- 0.25*lambda1^2/lambda2
#   
#   # starting values for mu and sigma
#   fit.pen <- penalized(y, x, unpenalized=formula(ifelse(intercept, "~1", "~0")), 
#                        model="logistic", lambda1=0, 
#                        lambda2=2*(lambda1 + lambda2), trace=FALSE)
#   if(intercept) {
#     xadj <- cbind(1, x)
#   } else {
#     xadj <- x
#   }
#   b0 <- coef(fit.pen, which="all")
#   pred0 <- as.numeric(exp(xadj %*% b0)/(1 + exp(xadj %*% b0)))
#   w <- sqrt(pred0*(1 - pred0))
#   xw <- xadj*w
#   svdxw <- svd(xw)
#   d <- svdxw$d
#   v <- svdxw$v
#   invmat <- d^2/(d^2 + 4*(lambda1 + lambda2))^2
#   sigmaold <- t(t(v)*invmat) %*% t(v)
#   
#   # starting values
#   muold <- as.numeric(sigmaold %*% (t(xadj) %*% as.matrix(kappa)))
#   ci <- as.numeric(sqrt(colSums(t(xadj) * (sigmaold %*% t(xadj))) + (colSums(t(xadj)*muold))^2))
#   chi <- as.numeric(0.5*(lambda1 + lambda2)*(diag(sigmaold) + muold^2))[(intercept + 1):(intercept + p)]
#   
#   conv <- FALSE
#   iter <- 0
#   if(trace) {cat("\n", "Estimating model parameters by variational Bayes", "\n", sep="")}
#   srt <- proc.time()[3]
#   while(!conv & (iter < maxiter)) {
#     
#     iter <- iter + 1
#     
#     # estimating new model parameters
#     newparam <- est_param(x, kappa, m, n, p, ci, rep(phi, p), chi, rep(lambda2, p), intercept)
#     sigma <- newparam$sigma
#     mu <- newparam$mu
#     ci <- newparam$ci
#     chi <- newparam$chi
#       
#     # checking convergence
#     conv <- max(c(abs((mu - muold)/muold), abs(diag((sigma - sigmaold)/sigmaold)))) < eps
#       
#     muold <- mu
#     sigmaold <- sigma
#     
#     if(trace) {cat("\r", "Iteration ", iter, "  ", sep="")}
#   }
#   vb.time <- proc.time()[3] - srt
#   
#   # output
#   out <- list(mu=mu, sigma=sigma, ci=ci, chi=chi, lambda1=lambda1, lambda2=lambda2, 
#               niter=iter, conv=conv, cv.time=cv.time, vb.time=vb.time)
#   
#   return(out)
#    
# }
# 
# sourceCpp(paste(path.code, "ENgibbs.cpp", sep=""))
# MCMCen <- function(x, y, m, lambda1=NULL, lambda2=NULL, intercept, K=5000, trace=TRUE) {
#   
#   # data characteristics
#   p <- ncol(x)
#   n <- nrow(x)
#   
#   # if no penalty parameters are given we estimate them by cross-validation
#   cv.time <- NULL
#   if(is.null(lambda1) | is.null(lambda2)) {
#     if(trace) {cat("\r", "Estimating lambda1 and lambda2 by cross-validation", sep="")}
#     srt <- proc.time()[3]
#     opt.glob <- cv.pen(x, y, intercept)
#     cv.time <- proc.time()[3] - srt
#     lambda1 <- opt.glob$lambda1[which.min(opt.glob$cvll)]
#     lambda2 <- opt.glob$lambda2[which.min(opt.glob$cvll)]
#     if(trace) {cat("\n", "Lambda1 and lambda2 estimated at ", round(lambda1, 2), " and ", 
#                    round(lambda2, 2), " in ", round(cv.time, 2), " seconds", sep="")}
#   }
#   
#   # starting values for beta
#   fit.pen <- penalized(y, x, unpenalized=formula(ifelse(intercept, "~1", "~0")), 
#                        model="logistic", lambda1=0, 
#                        lambda2=2*(lambda1 + lambda2), trace=FALSE)
#   if(intercept) {
#     xadj <- cbind(1, x)
#   } else {
#     xadj <- x
#   }
#   b0 <- coef(fit.pen, which="all")
#   
#   # sampling from posterior
#   if(trace) {cat("\n", "Sampling from posterior", sep="")}
#   srt <- proc.time()[3]
#   gsamp <- gibbsC(x, y, m, n, p, lambda1, rep(lambda2, p), b0, intercept, K)
#   mcmc.time <- proc.time()[3] - srt
#   
#   # output
#   out <- list(beta=gsamp$beta, tau=gsamp$tau, omega=gsamp$omega, lambda1=lambda1, 
#               lambda2=lambda2, cv.time=cv.time, mcmc.time=mcmc.time)
#   
#   return(out)
#   
# }
# 
# 
#
# set.seed(9001)
# n <- 200
# p <- 200
# G <- 2
# x <- matrix(rnorm(n*p), ncol=p, nrow=n)
# lambda1 <- 1
# lambda2 <- 1
# lambdag <- c(0.2, 5)
# m <- rep(1, n)
# b0 <- rnorm(p + 1)
# sigma0 <- diag(rchisq(p + 1, 1))
# beta <- c(renbeta(p/G, sqrt(lambdag[1])*lambda1, lambdag[1]*lambda2), 
#           renbeta(p/G, sqrt(lambdag[2])*lambda1, lambdag[2]*lambda2))
# y <- rbinom(n, m, exp(x %*% beta)/(1 + exp(x %*% beta)))
# 
# test1.vb <- VBen(x, y, m, lambda1=1, lambda2=1, intercept=FALSE, eps=0.001, 
#                  maxiter=200, trace=TRUE)
# test1.mcmc <- MCMCen(x, y, m, lambda1=1, lambda2=1, intercept=FALSE, K=10000, trace=TRUE)
# 
# VBmu <- test1.vb$mu
# VBmcmc <- apply(test1.mcmc$beta, 1, mean)
# plot(VBmu, VBmcmc)
# 
# VBsigma <- diag(test1.vb$sigma)
# MCMCsigma <- apply(test1.mcmc$beta, 1, var)
# plot(VBsigma, MCMCsigma, col=rep(1:G, each=p/G))
# 
# VBmult <- VBsigma/exp(mean(log(VBsigma)))
# MCMCmult <- MCMCsigma/exp(mean(log(MCMCsigma)))
# plot(VBmult, MCMCmult)
# 
# VBgmsigmag <- sapply(1:G, function(g) {exp(mean(log(VBsigma[c(((g - 1)*p/G + 1):(g*(p/G)))])))})
# MCMCgmsigmag <- sapply(1:G, function(g) {exp(mean(log(MCMCsigma[c(((g - 1)*p/G + 1):(g*(p/G)))])))})
# 
# VBamsigmag <- sapply(1:G, function(g) {mean(VBsigma[c(((g - 1)*p/G + 1):(g*(p/G)))])})
# MCMCamsigmag <- sapply(1:G, function(g) {mean(MCMCsigma[c(((g - 1)*p/G + 1):(g*(p/G)))])})
# 
# (1/VBgmsigmag)/exp(mean(log(1/VBgmsigmag))) 
# (1/MCMCgmsigmag)/exp(mean(log(1/MCMCgmsigmag))) 
# 
# 1/VBamsigmag
# 1/MCMCamsigmag
# 
# dim(test1.mcmc$beta)
# cov(test1.mcmc$beta[1:100, ])
# cov(test1.mcmc$beta[1:100, ])
# 
# boxplot(cbind(VBsigma[1:100], VBsigma[101:200]))
# boxplot(cbind(MCMCsigma[1:100], MCMCsigma[101:200]))
# boxplot(cbind(beta[1:100], beta[101:200]))
# 
# exp(sum(log(VBmsigma))/3)/VBmsigma
# exp(sum(log(MCMCmsigma))/3)/MCMCmsigma
# 
# cbind(apply(test1.mcmc$beta, 1, var), diag(test1.vb$sigma))
# cbind(VBmult, MCMCmult)






