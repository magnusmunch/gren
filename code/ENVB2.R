##############################  preamble  #############################
# code belonging to ENVB2.pdf                                         #
# version: 02                                                         #
# author: Magnus Münch                                                #
# created: 10-10-2016                                                 #
# last edited: 31-10-2016                                             #
#######################################################################



###############################  notes  ###############################
# 31-10-2016: still have to update to include global l2 cv            #
# 27-10-2016: adjusted marginal likelihoods and their gradients       #
# 17-10-2016: Used woodbury identity for inverse, doesnt work yet,    #
#             I'm guessing it's due to numerical issues               #
# 14-10-2016: SVD is not accurate                                     #
# 10-10-2016: need to still adjust the gradients in lambda1 and       #
#             lambda2 seperately                                      #
# 10-10-2016: Version for multiple lambdas                            #
#######################################################################



### paths
path.data <- "C:/Users/Magnus/Documents/phd/data/"

### libraries
library(glmnet)
library(penalized)
library(mvtnorm)
library(GRridge)
library(pROC)

### functions
# these three functions are used in sampling from the elastic net prior
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

# below three functions for marginal likelihood calculations
marg.ll <- function(lambda, e.beta, v.beta, e.psi.inv, e.psi, G, sizes,
                    modmat) {
  # this functions calculates the marginal likelihood 
  # (the part proportional to lambda)
  # lambda is a vector with first element lambda1 and second lambda2
  # the e's are expectations from previous iteration
  
  lambda1 <- lambda[1:G]
  lambda2 <- lambda[(G + 1):(2*G)]
  lambda1vec <- rep(lambda1, times=sizes)
  lambda2vec <- rep(lambda2, times=sizes)
  part1 <- sum(sizes*log(lambda1))
  part2 <- 0.5*sum(lambda2vec*(v.beta + e.beta^2)*(1 + e.psi.inv))
  part3 <- sum(lambda1vec^2*(e.psi + 1)/(8*lambda2vec))
  part4 <- sum(sizes*pnorm(-lambda1/(sqrt(4*lambda2)), log.p=TRUE))
  ll <- part1 - part2 - part3 - part4
  return(ll)
  
}

marg.ll2 <- function(lambda2, lambda1, e.beta, v.beta, e.psi.inv, e.psi, sizes,
                     modmat) {
  # marginal likelihood in lambda2
  
  lambda1vec <- rep(lambda1, times=sizes)
  lambda2vec <- rep(lambda2, times=sizes)
  part1 <- 0.5*sum(lambda2vec*(v.beta + e.beta^2)*(1 + e.psi.inv))
  part2 <- sum(lambda1vec^2*(e.psi + 1)/(8*lambda2vec))
  part3 <- sum(sizes*pnorm(-lambda1/(sqrt(4*lambda2)), log.p=TRUE))
  ll <- -part1 - part2 -part3
  return(ll)
  
}

marg.ll1 <- function(lambda1, lambda2, e.psi, sizes, modmat) {
  # marginal likelihood in lambda1
  
  lambda1vec <- rep(lambda1, times=sizes)
  lambda2vec <- rep(lambda2, times=sizes)
  part1 <- sum(sizes*log(lambda1))
  part2 <- sum(lambda1vec^2*(e.psi + 1)/(8*lambda2vec))
  part3 <- sum(sizes*pnorm(-lambda1/(sqrt(4*lambda2)), log.p=TRUE))
  ll <- part1 - part2 - part3
  
  return(ll)
  
}

gr.marg.ll <- function(lambda, e.beta, v.beta, e.psi.inv, e.psi, G, sizes,
                       modmat) {
  # gradient of marginal likelihood
  lambda1 <- lambda[1:G]
  lambda2 <- lambda[(G + 1):(2*G)]
  part1.1 <- sizes/lambda1
  part1.2 <- lambda1*(t(as.matrix(e.psi + 1)) %*% modmat)/(4*lambda2)
  part1.3 <- sizes*dnorm(lambda1/sqrt(4*lambda2))/
    (sqrt(4*lambda2)*pnorm(-lambda1/sqrt(4*lambda2)))
  comp1 <- part1.1 - part1.2 + part1.3
  part2.1 <- lambda1^2*(t(as.matrix(e.psi + 1)) %*% modmat)/(8*lambda2^2)
  part2.2 <- 0.5*(t(as.matrix((v.beta + e.beta^2)*(1 + e.psi.inv))) %*% modmat)
  part2.3 <- sizes*lambda1*dnorm(lambda1/sqrt(4*lambda2))/
    (4*lambda2^(1.5)*pnorm(-lambda1/sqrt(4*lambda2)))
  comp2 <- part2.1 - part2.2 - part2.3
  return(c(comp1, comp2))
  
}

gr.marg.ll2 <- function(lambda2, lambda1, e.beta, v.beta, e.psi.inv, e.psi,
                        sizes, modmat) {
  # gradient of marginal likelihood in lambda2
  part1 <- lambda1^2*(t(as.matrix(e.psi + 1)) %*% modmat)/(8*lambda2^2)
  part2 <- 0.5*t(as.matrix((v.beta + e.beta^2)*(e.psi.inv + 1))) %*% modmat
  part3 <- sizes*lambda1*dnorm(lambda1/sqrt(4*lambda2))/
    (4*lambda2^(1.5)*pnorm(-lambda1/sqrt(4*lambda2)))
  gr <- part1 - part2 - part3
  return(gr)
  
}

gr.marg.ll1 <- function(lambda1, lambda2, e.psi, sizes, modmat) {
  # gradient of marginal likelihood in lambda1
  part1 <- sizes/lambda1
  part2 <- lambda1*(t(matrix(e.psi + 1)) %*% modmat)/(4*lambda2)
  part3 <- sizes*dnorm(lambda1/sqrt(4*lambda2))/
    (sqrt(4*lambda2)*pnorm(-lambda1/sqrt(4*lambda2)))
  gr <- part1 - part2 + part3
  return(gr)
  
}

# the fitting function
envb2 <- function(x, y, groups, lambda1=NULL, lambda2=NULL, mustart=NULL, 
                  sigmastart=NULL, model=c("binomial", "multinomial"), opt.prior=TRUE, 
                  fix.global.lambda2=FALSE, inv=c(NULL, "svd", "woodbury"), maxiter=1000, 
                  epsilon=1e-07, trace=TRUE) {
  
  #input: data x and y, y is a matrix with columns of category counts;
  #       groups is a numeric vector of length p, with values from 
  #       1 to number of groups; penalty parameters, lambda1 and lambda2;
  #       starting values, mustart and sigmastart, in binomial these are 
  #       a vector and matrix respectively, in multinomial a matrix and
  #       list of matrices, respectively; model; opt.prior determines 
  #       whether prior is optimised, if this is TRUE both optimised,
  #       if vector then lambda's corresponding to 1 or TRUE in vector
  #       optimised; fix.global.lambda2=FALSE implies that lambda2 is 
  #       optimised using MML, otherwise global lambda2 is estimated by CV,
  #       and multipliers are estimated by MML; inv is the matrix 
  #       inversion method used, NULL means no special one; maximum number 
  #       of iterations; tolerance for convergence check; if trace is TRUE, 
  #       intermediate interation number and lambda values are printed
  
  # checking the inputs and issueing warnings and/or errors
  errors <- c("x is not numeric", "groups is not numeric", 
              "length(groups) should be equal to ncol(x)",
              "lambda1 is not NULL or numeric",
              "lambda2 is not NULL or numeric", 
              "mustart and/or sigmastart are not NULL or numeric",
              "model is not binomial or multinomial",
              "opt.prior should be 0, 1, c(0, 1) or c(1, 0), or the logical equivalents thereof",
              "inv should be NULL, svd or woodbury",
              "maxiter is not an integer", "epsilon is not numeric",
              "trace is not a logical")
  checks <- c(is.numeric(x), is.numeric(groups), length(groups)==ncol(x),
              is.null(lambda1) | is.numeric(lambda1), 
              is.null(lambda2) | is.numeric(lambda2), 
              is.null(mustart) | is.numeric(mustart) | is.null(sigmastart) |
                is.numeric(sigmastart), model %in% c("binomial", "multinomial"),
              opt.prior==0 | opt.prior==1 | opt.prior==c(1, 0) | opt.prior==c(0, 1),
              is.null(inv) | inv %in% c("svd", "woodbury"), 
              round(maxiter, 0)==maxiter, is.numeric(epsilon), 
              is.logical(trace))
  if(sum(!checks)==1) {
    mes <- errors[!checks]
    stop(mes)
  } else if(sum(!checks)==2) {
    mes <- paste(errors[min(which(!checks))], ", in addition: ", 
                 errors[max(which(!checks))], sep="")
    stop(mes)
  } else if(sum(!checks) > 2) {
    seler <- errors[which(!checks)[which(!checks)!=min(which(!checks))]]
    mes2 <- paste("(", 2:sum(!checks), ") ", seler, collapse=", ", 
                  sep="")
    mes <- paste(errors[min(which(!checks))], ", in addition: ",
                 mes2, sep="")
    stop(mes)
  }
  
  p <- ncol(x)
  n <- nrow(x)
  sizes <- rle(groups)$lengths
  G <- length(unique(groups))
  if(sum(opt.prior) > 0) {
    # calculate the model matrix for the groups (needed in optim)
    modmat <- matrix(0, ncol=G, nrow=p)
    modmat <- sapply(1:G, function(g) {as.numeric(groups==g)})
  }
  
  # check whether y is a numeric or matrix with one column,
  # in that case interpreted as binary logistic and column of target 0 added
  if(ifelse(!is.numeric(y), dim(y)[2] < 2, is.null(dim(y)))) {
    y <- cbind(1 - y, y)
  }
  C <- ncol(y) # number of categories (2 in binomial)
  
  # these are all fixed quantities, only calculate once
  b <- m <- rowSums(y)
  kappa <- y[, 2:C] - m/2 # this is a matrix now (for multinomial case)
  lambda <- 0.5
  
  # these quantities change only if opt.prior=TRUE 
  # we also initialize a sequence in which we store the hyperparameters
  # over the iterations
  # we first check whether starting values are provided otherwise use cv
  if(is.null(lambda1) | is.null(lambda2)) {
    
    if(is.null(lambda1) & is.null(lambda2)) {
      
      fit <- cv.glmnet(x, y, family=model, alpha=0, intercept=FALSE,
                       standardize=FALSE)
      lambda1old <- lambda1seq <- lambda1 <- rep_len(fit$lambda.min/2)
      lambda2old <- lambda2seq <- lambda2 <- rep_len(fit$lambda.min/2)
      mu <- muold <- as.numeric(coef(fit, s="lambda.min"))
      
    } else {
      
      fit <- glmnet(x, y, family=model, alpha=0, lambda=(!is.null(lambda2))*lambda2*2.5 + 
                      (!is.null(lambda1))*lambda1*1.5, intercept=FALSE)
      lambda1old <- lambda1seq <- lambda1 <- rep_len(fit$lambda/2)
      lambda2old <- lambda2seq <- lambda2 <- rep_len(fit$lambda/2)
      
    }  
    
    if(model=="binomial") {
      # if the model is binomial these are the starting values
      phat <- exp(x %*% mu)/(1 + exp(x %*% mu))
      w <- sqrt(phat*(1 - phat))
      Xw <- x*w
      # use the first lambda, since theyre all the same
      invmat <- solve(crossprod(Xw) + 2*(lambda1[1] + lambda2[1]/2)*diag(p))
      sigma <- sigmaold <- invmat %*% crossprod(Xw) %*% invmat
    }   
  } else {
    
    # if provided but length < G, we recycle
    lambda1old <- lambda1seq <- lambda1 <- rep_len(lambda1, length.out=G)
    lambda2old <- lambda2seq <- lambda2 <- rep_len(lambda2, length.out=G)
    
  }
  lambda1vec <- rep(lambda1old, times=sizes)
  lambda2vec <- rep(lambda2old, times=sizes)
  phi <- lambda1vec^2/(4*lambda2vec) # this is vector of length p
  
  # check whether starting values provided, otherwise estimate with glmnet
  # if we estimated lambda1 and/or lambda2, then starting values are 
  # estimated earlier
  if(model=="binomial") {
    if(is.null(mustart)) {
      if((length(unique(lambda1)) > 1) | (length(unique(lambda2)) > 1)) {
        # if one of the two is a vector use grridge
        part <- CreatePartition(as.factor(groups))
        fit <- grridge(t(x), as.factor(y), list(part=part), ~0)
        mu <- muold <- as.numeric(fit$betas)
      } else {
        # if both are length 1 (or all the same), use glmnet (faster)
        fit <- glmnet(x, y, family=model, lambda=mean(lambda1vec) + mean(lambda2vec)/2, 
                      alpha=0, intercept=FALSE, standardize=FALSE)
        mu <- muold <- as.numeric(fit$beta)
      }
    } else {
      mu <- muold <- mustart
    } 
    if(is.null(sigmastart)) {
      xmu <- x %*% mu
      phat <- as.numeric(exp(xmu)/(1 + exp(xmu)))
      W <- diag(sqrt(phat*(1 - phat)))
      Xw <- W %*% x
      svdxw <- svd(Xw)
      U <- svdxw$u
      V <- svdxw$v
      d <- svdxw$d
      td <- min(n, p)
      invmat <- 1/(d^2 + 2*(mean(lambda1vec) + mean(lambda2vec)/2))
      part1 <- invmat^2*d^2
      sigma <- t(t(V)*part1) %*% t(V)
    } else {
      sigma <- sigmaold <- sigmastart
    }
    
    ci <- ciold <- as.numeric(sqrt(colSums(t(x) * (sigma %*% t(x))) + (colSums(t(x)*mu))^2))
    chi <- chiold <- as.numeric(lambda2vec*(diag(sigma) + mu^2))
    
  } else {
    # THIS DOES NOT WORK CURRENTLY
    if(is.null(mustart)) {
      fit <- glmnet(x, y, family=model, lambda=mean(lambda1vec) + lambda2/2, 
                    alpha=0, intercept=FALSE)
      mu <- muold <- do.call(cbind, lapply(fit$beta, function(b) {
        as.matrix(b)}))[, -C] # we keep the last category fixed to 0
    } else {
      mu <- muold <- mustart
    }
    if(is.null(sigmastart)) {
      sigma <- sigmaold <- vector(mode="list", length=C - 1)
      for(c in 1:(C - 1)) {
        # store the sigma matrices in a list
        W <- diag(sqrt(phat[, c]*(1 - phat[, c])))
        Xw <- W %*% x
        invmat <- solve(crossprod(Xw) + 2*(lambda1 + lambda2/2)*diag(p))
        sigma[[c]] <- sigmaold[[c]] <- invmat %*% crossprod(Xw) %*% invmat
      }
    } else {
      sigma <- sigmaold <- sigmastart
    }
    
    # these are matrices now
    ci <- ciold <- sqrt(do.call(cbind, lapply(sigma, function(s) {
      diag(x %*% s %*% t(x))})) + (x %*% mu)^2)
    chi <- chiold <- lambda2*(do.call(cbind, lapply(sigma, diag)) + mu^2)
    
    R <- matrix(NA, nrow=n, ncol=C - 1) # initialize this for later
    
  }
  
  conv <- FALSE
  niter <- 0
  
  while(!conv & (niter < maxiter)) {
    
    if(model=="binomial") {
      
      # these are updating equations for the binomial model 
      if(!is.null(inv)) {
        # faster inversion methods
        h <- (1 + sqrt(phi/chiold))*lambda2vec
        om <- (0.5*m/ciold)*tanh(ciold/2)
        omsq <- sqrt(om)
        hinvsq <- 1/sqrt(h)
        if(inv=="svd") {
          # using the svd of Om^(1/2)*X*H^(-1/2)
          svdx <- svd(t(x * omsq) * hinvsq)
          U <- svdx$v
          V <- svdx$u
          d <- svdx$d
          hinvsqV <- V*hinvsq
          sigma <- hinvsqV %*% (t(hinvsqV)/(d^2 + 1))
        } else {
          # using the Woodbury identity 
          V <- x*omsq
          U <- t(V)
          hinv <- 1/h
          vhinv <- t(t(V)*hinv)
          sigma <- diag(hinv) - 
            (hinv*U) %*% solve(diag(n) + vhinv %*% U) %*% vhinv
        }
      } else {
        Om <- diag((0.5*m/ciold)*tanh(ciold/2))
        Z <- diag(sqrt(phi/chiold))
        sigma <- solve(t(x) %*% Om %*% x + lambda2vec*diag(p) + lambda2vec*Z)
      }
      mu <- as.numeric(sigma %*% (t(x) %*% kappa))
      ci <- as.numeric(sqrt(colSums(t(x) * (sigma %*% t(x))) + (colSums(t(x)*mu))^2))
      chi <- as.numeric(lambda2vec*(diag(sigma) + mu^2))
      
      if(sum(opt.prior) > 0) {
        # if we want to calculate MML estimates for hyperparameters,
        # this step is included in every iteration
        
        # recalculate phi, since it depends on lambda1 and lambda2
        phi <- lambda1vec^2/(4*lambda2vec)
        
        # fixed parameters needed in function optimisation
        e.beta <- mu
        v.beta <- diag(sigma)
        e.psi.inv <- sqrt(phi/chi)
        e.psi <- 1/phi + sqrt(chi/phi)
        
        # optimisation routine to find MML
        # if opt.prior=TRUE then optimise both, if opt.prior is a vector, 
        # then we optimise only lambda1 if opt.prior=c(1, 0) or 
        # c(TRUE, FALSE), only lambda2 optimisation if 
        # opt.prior=c(0, 1) or c(FALSE, TRUE)
        if((sum(opt.prior)/length(opt.prior))==1) {
          # optimise both (not recommended)
          lambdaold <- c(lambda1old, lambda2old)
          opt.rout <- optim(par=lambdaold, fn=marg.ll, gr=gr.marg.ll,
                            method="L-BFGS-B", lower=rep(0.001, 2*G), 
                            upper=rep(Inf, 2*G), 
                            control=list(fnscale=-1), e.beta=e.beta, 
                            v.beta=v.beta, e.psi.inv=e.psi.inv, e.psi=e.psi, G=G,
                            sizes=sizes, modmat=modmat)
          lambda1 <- opt.rout$par[1:G]
          lambda2 <- opt.rout$par[(G + 1):(2*G)]
        } else if(opt.prior[2]) {
          #optimise lambda2 (recommended)
          opt.rout <- optim(par=lambda2old, fn=marg.ll2, gr=gr.marg.ll2,
                            method="L-BFGS-B", lower=rep(0.001, G), 
                            upper=rep(Inf, G),
                            control=list(fnscale=-1), lambda1=lambda1old,
                            e.beta=e.beta, v.beta=v.beta, e.psi.inv=e.psi.inv, 
                            e.psi=e.psi, sizes=sizes, modmat=modmat)
          lambda2 <- opt.rout$par
        } else {
          # optimse lambda1
          opt.rout <- optim(par=lambda1old, fn=marg.ll1, gr=gr.marg.ll1,
                            method="L-BFGS-B", lower=rep(0.001, G), 
                            upper=rep(Inf, G),
                            control=list(fnscale=-1), lambda2=lambda2old,
                            e.psi=e.psi, sizes=sizes, modmat=modmat)
          lambda1 <- opt.rout$par
        }
        
      }
      
      # check the convergence of all parameters (including hyperparameters,
      # which are fixed if we set opt.prior=FALSE)
      conv <- max(abs(c(sigma - sigmaold, mu - muold))) < epsilon
      niter <- niter + 1
      
      # if we want the iteration trace, print current iteration
      if(trace) {
        cat("\r", "Iteration: ", niter, ", ", 
            paste("lambda1=[", paste(round(lambda1, 2), collapse=", "), 
                  sep=""), "], ", 
            paste("lambda2=[", paste(round(lambda2, 2), collapse=", "), 
                  sep=""), "]", sep="")
      }
      
      # update old parameters to new ones
      sigmaold <- sigma
      muold <- mu
      ciold <- ci
      chiold <- chi
      
      # update old hyperparameters to new ones
      lambda1old <- lambda1
      lambda2old <- lambda2
      lambda1vec <- rep(lambda1old, times=sizes)
      lambda2vec <- rep(lambda2old, times=sizes)
      
      # update sequence of hyperparameters over the iterations
      lambda1seq <- cbind(lambda1seq, lambda1)
      lambda2seq <- cbind(lambda2seq, lambda2)
      
    } else {
      
      # THIS IS NOT FUNCTIONAL MOMENTARILY
      # and here in the multinomial setting
      for(c in 1:(C - 1)) {
        R[, c] <- log(sum(1)) ## left of here with multinomial
        Om <- diag((0.5*m/ciold[, c])*tanh(ciold[, c]/2))
        Z <- diag(sqrt(phi/chiold[, c]))
        sigma[[c]] <- solve(t(x) %*% Om %*% x + lambda2*diag(p) + lambda2*Z)
        mu[, c] <- sigma[[c]] %*% (t(x) %*% kappa[, c] + t(x) %*% Om)
      }
      
      ci <- ciold <- sqrt(do.call(cbind, lapply(sigma, function(s) {
        diag(x %*% s %*% t(x))})) + (x %*% mu)^2)
      chi <- chiold <- lambda2*(do.call(cbind, lapply(sigma, diag)) + mu^2)
      
    }
  }
  
  # create output list with components:
  # niter: number of iterations of algorithm
  # conv: did the algorithm (parameters) converge or not?
  # sigma: the final covariance matrix of the beta
  # mu: the final vector of expectations of the beta
  # c: the final vector of c parameters for the omega posterior part
  # chi: final vector of chi parameters for the GIG posterior part
  # lambda1 and lambda2: final hyperparameter estimates (they are the same
  #                      as input if opt.prior=FALSE)
  out <- list(niter=niter, conv=conv, sigma=sigma, mu=mu, c=ci, chi=chi,
              lambda1=lambda1seq, lambda2=lambda2seq)
  return(out)
  
}

# test 1
set.seed(123)
n <- 50
p <- 200
G <- 2
groups <- rep(1:G, each=p/G)
lambda1 <- rep(0.5, G)
lambda2 <- rep(1, G)
x <- matrix(rnorm(n*p), ncol=p, nrow=n)
beta <- c(sapply(1:G, function(g) {
  renbeta(p/G, lambda1[g], lambda2[g])}))
prob <- exp(x %*% beta)/(1 + exp(x %*% beta))
y <- as.numeric(runif(n) < prob)
y <- cbind(1 - y, y)

test1 <- envb2(x, y, groups, lambda1=lambda1, lambda2=1, 
               mustart=NULL, sigmastart=NULL, 
               model="binomial", opt.prior=c(0, 1), inv=NULL, 
               maxiter=1000, epsilon=1e-07, trace=TRUE)
cred.int <- cbind(test1$mu - qnorm(0.975)*sqrt(diag(test1$sigma)),
                  test1$mu + qnorm(0.975)*sqrt(diag(test1$sigma)))
test2 <- optL2(y[, 2], x, unpenalized=~0, model="logistic", 
               lambda1=mean(lambda1), standardize=FALSE)
test3 <- cv.glmnet(x, y[, 2], alpha=0.5/(2*(test2$lambda + 0.25)),
                   family="binomial", standardize=FALSE, intercept=FALSE,
                   type.measure="deviance")

# comparing outputs
plot(test1$mu, test2$fullfit@penalized)
plot(test1$mu, test3$glmnet.fit$beta
     [, which(test3$lambda==test3$lambda.min)])
plot(test2$fullfit@penalized, test3$glmnet.fit$beta
     [, which(test3$lambda==test3$lambda.min)])
plot(beta, test1$mu)
plot(beta, test2$fullfit@penalized)
plot(beta, test3$glmnet.fit$beta
     [, which(test3$lambda==test3$lambda.min)])

var(beta - test1$mu)
var(beta - test2$fullfit@penalized)
var(beta - test3$glmnet.fit$beta
    [, which(test3$lambda==test3$lambda.min)])
round(cbind(beta=beta, penalized=test2$fullfit@penalized, 
            glmnet=test3$glmnet.fit$beta
            [, which(test3$lambda==test3$lambda.min)], ENVB=test1$mu),
      digits=5)



# test 2
set.seed(1)
n <- 50
p <- 500
G <- 2
groups <- rep(1:G, each=p/G)
lambda1 <- rep(0.5, G)
lambda2 <- c(1, 10)
x <- matrix(rnorm(n*p), ncol=p, nrow=n)
beta <- c(sapply(1:G, function(g) {
  renbeta(p/G, lambda1[g], lambda2[g])}))
prob <- exp(x %*% beta)/(1 + exp(x %*% beta))
y <- as.numeric(runif(n) < prob)
y <- cbind(1 - y, y)

test1 <- envb2(x, y, groups, lambda1=lambda1, lambda2=mean(lambda2), 
               mustart=NULL, sigmastart=NULL, 
               model="binomial", opt.prior=TRUE, inv=NULL, 
               maxiter=1000, epsilon=0.001, trace=TRUE)
test2 <- optL2(y[, 2], x, unpenalized=~0, model="logistic",
               lambda1=mean(lambda1), standardize=FALSE)
test3 <- cv.glmnet(x, y[, 2], alpha=0.5/(2*(test2$lambda + 0.25)),
                   family="binomial", standardize=FALSE, intercept=FALSE,
                   type.measure="deviance")

# comparing outputs
plot(test1$mu, test2$fullfit@penalized)
plot(test1$mu, test3$glmnet.fit$beta
     [, which(test3$lambda==test3$lambda.min)])
plot(test2$fullfit@penalized, test3$glmnet.fit$beta
     [, which(test3$lambda==test3$lambda.min)])
plot(beta, test1$mu)
plot(beta, test2$fullfit@penalized)
plot(beta, test3$glmnet.fit$beta
     [, which(test3$lambda==test3$lambda.min)])

var(beta - test1$mu)
var(beta - test2$fullfit@penalized)
var(beta - test3$glmnet.fit$beta
    [, which(test3$lambda==test3$lambda.min)])
round(cbind(penalized=test2$fullfit@penalized, 
            glmnet=test3$glmnet.fit$beta
            [, which(test3$lambda==test3$lambda.min)], ENVB=test1$mu),
      digits=5)
ci <- cbind(qnorm(0.40, mean=test1$mu, sd=sqrt(diag(test1$sigma))),
            qnorm(0.40, mean=test1$mu, sd=sqrt(diag(test1$sigma)),
                  lower.tail=FALSE))
sel <- which(ci[, 1]*ci[, 2] > 0)
bsel <- rep(0, p)
bsel[sel] <- test1$mu[sel]

# comparing selections
var(beta - test2$fullfit@penalized)
var(beta - test3$glmnet.fit$beta
    [, which(test3$lambda==test3$lambda.min)])
var(beta - bsel)

# testing on new data
xtest <- matrix(rnorm(n*p), ncol=p, nrow=n)
probtest <- exp(xtest %*% beta)/(1 + exp(xtest %*% beta))
ytest <- as.numeric(runif(n) < probtest)
ytest <- cbind(1 - ytest, ytest)

estENVB <- optL2(y[, -1], x[, sel], unpenalized=~0, lambda1=0)
estpenalized <- optL2(y[, -1], x[, test2$fullfit@penalized!=0],
                      unpenalized=~0, lambda1=0)
estglmnet <- penalized(y[, -1], x[, test3$glmnet.fit$beta[, which(test3$lambda==test3$lambda.min)]
                                  !=0], unpenalized=~0, lambda1=0,
                       lambda2=1e-18, model="logistic")

xbENVB <- as.matrix(xtest[, sel]) %*% 
  as.matrix(estENVB$fullfit@penalized)
xbpenalized <- as.matrix(xtest[, test2$fullfit@penalized!=0]) %*% 
  as.matrix(estpenalized$fullfit@penalized)
xbglmnet <- as.matrix(xtest[, test3$glmnet.fit$beta[, which(test3$lambda==test3$lambda.min)]
                            !=0]) %*% 
  as.matrix(estglmnet@penalized)
predENVB <- exp(xbENVB)/(1 + exp(xbENVB))
predpenalized <- exp(xbpenalized)/(1 + exp(xbpenalized))
predglmnet <- exp(xbglmnet)/(1 + exp(xbglmnet))

(auc1 <- auc(as.numeric(predENVB > 0.5), ytest[, -1]))
(auc2 <- auc(as.numeric(predpenalized > 0.5), ytest[, -1]))
(auc3 <- auc(as.numeric(predglmnet > 0.5), ytest[, -1]))
length(sel)
sum(test2$fullfit@penalized!=0)
sum(test3$glmnet.fit$beta[, which(test3$lambda==test3$lambda.min)]!=0)

# test 3
set.seed(123)
n <- 50
p <- 200
G <- 2
groups <- rep(1:G, each=p/G)
lambda1 <- rep(0.5, G)
lambda2 <- c(100, 10)
rho <- 0.8
sigma <- matrix(0, ncol=p, nrow=p)
sigma[1:(p/G), 1:(p/G)] <- sigma[(p/G + 1):p, (p/G + 1):p] <- rho
diag(sigma) <- 1
x <- rmvnorm(n, sigma=sigma)
beta <- c(sapply(1:G, function(g) {
  renbeta(p/G, lambda1[g], lambda2[g])}))
prob <- exp(x %*% beta)/(1 + exp(x %*% beta))
y <- as.numeric(runif(n) < prob)
y <- cbind(1 - y, y)

test1 <- envb2(x, y, groups, lambda1=lambda1, lambda2=mean(lambda2), 
               mustart=NULL, sigmastart=NULL, 
               model="binomial", opt.prior=c(0, 1), inv=NULL, 
               maxiter=1000, epsilon=0.001, trace=TRUE)
test2 <- optL2(y[, 2], x, unpenalized=~0, model="logistic",
               lambda1=mean(lambda1), standardize=FALSE)
test3 <- cv.glmnet(x, y[, 2], alpha=0.5/(2*(test2$lambda + 0.25)),
                   family="binomial", standardize=FALSE, intercept=FALSE,
                   type.measure="deviance")

# comparing outputs
plot(test1$mu, test2$fullfit@penalized)
plot(test1$mu, test3$glmnet.fit$beta
     [, which(test3$lambda==test3$lambda.min)])
plot(test2$fullfit@penalized, test3$glmnet.fit$beta
     [, which(test3$lambda==test3$lambda.min)])
plot(beta, test1$mu)
plot(beta, test2$fullfit@penalized)
plot(beta, test3$glmnet.fit$beta
     [, which(test3$lambda==test3$lambda.min)])



# test 4
set.seed(123)
n <- 50
p <- 200
G <- 2
groups <- rep(1:G, each=p/G)
lambda1 <- rep(0.5, G)
lambda2 <- c(100, 10)
rho <- 0.8
sigma <- matrix(0, ncol=p, nrow=p)
sigma[1:(p/G), 1:(p/G)] <- sigma[(p/G + 1):p, (p/G + 1):p] <- rho
diag(sigma) <- 1
x <- rmvnorm(n, sigma=sigma)
beta <- c(sapply(1:G, function(g) {
  renbeta(p/G, lambda1[g], lambda2[g])}))
prob <- exp(x %*% beta)/(1 + exp(x %*% beta))
y <- as.numeric(runif(n) < prob)
y <- cbind(1 - y, y)

test1.1 <- envb2(x, y, groups, lambda1=lambda1, lambda2=mean(lambda2), 
                 mustart=NULL, sigmastart=NULL, 
                 model="binomial", opt.prior=TRUE, inv=NULL, 
                 maxiter=1000, epsilon=0.001, trace=TRUE)
test1 <- penalized(y[, 2], x, unpenalized=~0, model="logistic",
                   lambda1=rep(test1.1$lambda1[, test1.1$niter + 1], 
                               each=p/G), 
                   lambda2=rep(test1.1$lambda2[, test1.1$niter + 1], 
                               each=p/G))
test2 <- optL2(y[, 2], x, unpenalized=~0, model="logistic",
               lambda1=mean(lambda1), standardize=FALSE)
var(beta - test1@penalized)
var(beta - test2$fullfit@penalized)


test3 <- optL2(y[, 2], x[, which(test1@penalized!=0)], unpenalized=~0, 
               model="logistic", lambda1=0)
test4 <- optL2(y[, 2], x[, which(test2$fullfit@penalized!=0)], 
               unpenalized=~0, model="logistic", lambda1=0)

xtest <- rmvnorm(n, sigma=sigma)
probtest <- exp(xtest %*% beta)/(1 + exp(xtest %*% beta))
ytest <- as.numeric(runif(n) < probtest)
ytest <- cbind(1 - ytest, ytest)

betaest1 <- betaest2 <- rep(0, p)
betaest1[which(test1@penalized!=0)] <- test3$fullfit@penalized
betaest2[which(test2$fullfit@penalized!=0)] <- test4$fullfit@penalized
var(beta - betaest1)
var(beta - betaest2)
pred1 <- exp(xtest %*% betaest1)/(1 + exp(xtest %*% betaest1))
pred2 <- exp(xtest %*% betaest2)/(1 + exp(xtest %*% betaest2))
auc(as.numeric(pred1 > 0.5), ytest[, -1])
auc(as.numeric(pred2 > 0.5), ytest[, -1])




# data 1
data(dataVerlaat)
x <- apply(t(as.matrix(datcenVerlaat)), 2, function(x) {
  (x - mean(x))/sd(x)})
y <- respVerlaat
cutoffs <- quantile(pvalFarkas, probs=c(0.2, 0.4, 0.6, 0.8))
groups <- 1 + rowSums(sapply(1:4, function(g) {pvalFarkas >= cutoffs[g]}))
x <- x[, order(groups)]
groups <- sort(groups)

test1.1 <- envb2(x, y, groups, lambda1=1, lambda2=1, 
                 mustart=NULL, sigmastart=NULL, 
                 model="binomial", opt.prior=TRUE, inv="svd", 
                 maxiter=1000, epsilon=0.001, trace=TRUE)
test1 <- penalized(y[, 2], x, unpenalized=~0, model="logistic",
                   lambda1=rep(test1.1$lambda1[, test1.1$niter + 1], 
                               each=p/G), 
                   lambda2=rep(test1.1$lambda2[, test1.1$niter + 1], 
                               each=p/G))
test2 <- optL2(y[, 2], x, unpenalized=~0, model="logistic",
               lambda1=mean(lambda1), standardize=FALSE)
var(beta - test1@penalized)
var(beta - test2$fullfit@penalized)


test3 <- optL2(y[, 2], x[, which(test1@penalized!=0)], unpenalized=~0, 
               model="logistic", lambda1=0)
test4 <- optL2(y[, 2], x[, which(test2$fullfit@penalized!=0)], 
               unpenalized=~0, model="logistic", lambda1=0)

xtest <- rmvnorm(n, sigma=sigma)
probtest <- exp(xtest %*% beta)/(1 + exp(xtest %*% beta))
ytest <- as.numeric(runif(n) < probtest)
ytest <- cbind(1 - ytest, ytest)

betaest1 <- betaest2 <- rep(0, p)
betaest1[which(test1@penalized!=0)] <- test3$fullfit@penalized
betaest2[which(test2$fullfit@penalized!=0)] <- test4$fullfit@penalized
var(beta - betaest1)
var(beta - betaest2)
pred1 <- exp(xtest %*% betaest1)/(1 + exp(xtest %*% betaest1))
pred2 <- exp(xtest %*% betaest2)/(1 + exp(xtest %*% betaest2))
auc(as.numeric(pred1 > 0.5), ytest[, -1])
auc(as.numeric(pred2 > 0.5), ytest[, -1])












p <- 100
n <- 30
lambda <- 3
W <- matrix(rexp(n^2), ncol=n, nrow=n)
X <- matrix(rnorm(n*p), nrow=n, ncol=p)
mu <- rnorm(p)
Xw <- W %*% X
invmat1 <- solve(t(Xw) %*% Xw + lambda*diag(p))
sigma <- invmat1 %*% t(Xw) %*% Xw %*% invmat1

svdxw <- svd(Xw)
U <- svdxw$u
V <- svdxw$v
d <- svdxw$d
invmat2 <- 1/(d^2 + lambda)
part1 <- (invmat2^2)*(d^2)
sigma2 <- t(t(V)*part1) %*% t(V)
sigma3 <- V %*% diag(1/(d^2 + lambda)) %*% diag(d^2) %*% 
  diag(1/(d^2 + lambda)) %*% t(V)

sigma1[1:5, 1:5]
sigma2[1:5, 1:5]
sigma3[1:5, 1:5]



library(HyperbolicDist)
library(SuppDists)

p <- 0.5
lambda1 <- 0.5
lambda2 <- 0.7
beta <- 0.4
a <- lambda1^2/(4*lambda2)
b <- lambda2*beta^2
lambda <- b
mu <- 2*abs(beta)*lambda2/lambda1
x1 <- rgig(100000, c(p, a, b))
x2 <- 1/rinvGauss(100000, nu=mu, lambda=lambda)

hist(x1, breaks=50, freq=FALSE)
hist(x2, breaks=50, add=TRUE, freq=FALSE)