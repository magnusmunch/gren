##############################  preamble  #############################
# code belonging to ENVB2.pdf                                         #
# version: 01                                                         #
# author: Magnus Münch                                                #
# created: 27-09-2016                                                 #
# last edited: 27-10-2016                                             #
#######################################################################



###############################  notes  ###############################
# 27-10-2016: added test case                                         #
# 10-10-2016: Version for one lambda1 and/or one lambda2              #
# 06-10-2016: works for one lambda1 and one lambda2 now. Still needs  #
#             needs extension to multiple lambdas and extensive tests #
# 04-10-2016: still have to figure out multinomial, so not functional #
#             right now                                               #
#######################################################################



### paths

### libraries
library(glmnet)
library(penalized)
library(mvtnorm)

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
marg.ll <- function(lambda, e.beta, v.beta, e.psi.inv, e.psi, p) {
  # this functions calculates the marginal likelihood 
  # (the part proportional to lambda)
  # lambda is a vector with first element lambda1 and second lambda2
  # the e's are expectations from previous iteration
  
  lambda1 <- lambda[1]
  lambda2 <- lambda[2]
  ll <- - 0.5*lambda2*sum((v.beta + e.beta^2)*(1 + e.psi.inv)) + 
    p*log(lambda1) - p*pnorm(-lambda1/sqrt(4*lambda2), log.p=TRUE) -
    lambda1^2/(8*lambda2)*sum(e.psi + 1)
  
  return(ll)
  
}

marg.ll2 <- function(lambda2, e.beta, v.beta, e.psi.inv, e.psi, p, lambda1) {
  # marginal likelihood in lambda2

  ll <- -0.5*lambda2*sum((v.beta + e.beta^2)*(1 + e.psi.inv)) - 
    p*pnorm(-lambda1/sqrt(4*lambda2), log.p=TRUE) -
    lambda1^2/(8*lambda2)*sum(e.psi + 1)
  
  return(ll)
  
}

marg.ll1 <- function(lambda1, e.psi, p, lambda2) {
  # marginal likelihood in lambda1
  
  ll <- p*log(lambda1) -p*pnorm(-lambda1/sqrt(4*lambda2), log.p=TRUE) -
    lambda1^2/(8*lambda2)*sum(e.psi + 1)
  
  return(ll)
  
}

gr.marg.ll <- function(lambda, e.beta, v.beta, e.psi.inv, e.psi, p) {
  # gradient of marginal likelihood
  lambda1 <- lambda[1]
  lambda2 <- lambda[2]
  comp1 <- p/lambda1 - lambda1/(4*lambda2)*sum(e.psi + 1) - 
    p/sqrt(4*lambda2)*dnorm(lambda1/sqrt(4*lambda2))/
    pnorm(-lambda1/sqrt(4*lambda2))
  comp2 <- lambda1^2/(8*lambda2^2)*sum(e.psi + 1) - 
    0.5*sum((v.beta + e.beta^2)*(e.psi.inv + 1)) - 
    p*lambda1/(4*lambda2^(1.5))*dnorm(lambda1/sqrt(4*lambda2))/
    pnorm(-lambda1/sqrt(4*lambda2))
  return(c(comp1, comp2))
  
}

gr.marg.ll2 <- function(lambda2, e.beta, v.beta, e.psi.inv, e.psi, p, 
                        lambda1) {
  # gradient of marginal likelihood in lambda2
  gr <- lambda1^2/(8*lambda2^2)*sum(e.psi + 1) - 
    0.5*sum((v.beta + e.beta^2)*(e.psi.inv + 1)) - 
    p*lambda1/(4*lambda2^(1.5))*dnorm(lambda1/sqrt(4*lambda2))/
    pnorm(-lambda1/sqrt(4*lambda2))
  return(gr)
  
}

gr.marg.ll1 <- function(lambda1, e.psi, p, 
                        lambda2) {
  # gradient of marginal likelihood in lambda1
  gr <- p/lambda1 - lambda1/(4*lambda2)*sum(e.psi + 1) - 
    p/sqrt(4*lambda2)*dnorm(lambda1/sqrt(4*lambda2))/
    pnorm(-lambda1/sqrt(4*lambda2))
  return(gr)
  
}

# the fitting function
envb2.1 <- function(x, y, lambda1, lambda2, mustart=NULL, sigmastart=NULL,
                    model=c("binomial", "multinomial"), opt.prior=TRUE, 
                    maxiter, epsilon) {
  
  #input: data x and y, y is a matrix with columns of category counts; 
  #       penalty parameters, lambda1 and lambda2; starting values, 
  #       mustart and sigmastart, in binomial these are 
  #       a vector and matrix respectively, in multinomial a matrix and
  #       list of matrices, respectively; model; opt.prior determines 
  #       whether prior is optimised, if this is TRUE both optimised,
  #       if vector then lambda's corresponding to 1 or TRUE in vector
  #       optimised; maximum number of 
  #       iterations; tolerance for convergence check
  
  p <- ncol(x)
  n <- nrow(x)
  
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
  lambda1old <- lambda1seq <- lambda1
  lambda2old <- lambda2seq <- lambda2
  phi <- lambda1^2/(4*lambda2)
  
  # check whether starting values provided, otherwise estimate with glmnet
  if(model=="binomial") {
    if(is.null(mustart)) {
      fit <- glmnet(x, y, family=model, lambda=lambda1 + lambda2/2, 
                    alpha=0, intercept=FALSE)
      mu <- muold <- as.numeric(fit$beta)
    } else {
      mu <- muold <- mustart
    } 
    if(is.null(sigmastart)) {
      phat <- exp(x %*% mu)/(1 + exp(x %*% mu))
      W <- diag(sqrt(phat*(1 - phat)))
      Xw <- W %*% x
      invmat <- solve(crossprod(Xw) + 2*(lambda1 + lambda2/2)*diag(p))
      sigma <- sigmaold <- invmat %*% crossprod(Xw) %*% invmat
    } else {
      sigma <- sigmaold <- sigmastart
    }
    
    ci <- ciold <- as.numeric(sqrt(diag(x %*% sigma %*% t(x)) + (x %*% mu)^2))
    chi <- chiold <- as.numeric(lambda2*(diag(sigma) + mu^2))
    
  } else {
    if(is.null(mustart)) {
      fit <- glmnet(x, y, family=model, lambda=lambda1 + lambda2/2, 
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
      Om <- diag((0.5*m/ciold)*tanh(ciold/2))
      Z <- diag(sqrt(phi/chiold))
      sigma <- solve(t(x) %*% Om %*% x + lambda2old*diag(p) + lambda2old*Z)
      mu <- sigma %*% t(x) %*% kappa
      ci <- as.numeric(sqrt(diag(x %*% sigma %*% t(x)) + (x %*% mu)^2))
      chi <- as.numeric(lambda2old*(diag(sigma) + mu^2))
    
      if(sum(opt.prior) > 0) {
        # if we want to calculate MML estimates for hyperparameters,
        # this step is included in every iteration
        
        # recalculate phi, since it depends on lambda1 and lambda2
        phi <- lambda1old^2/(4*lambda2old)
        
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
          lambdavec <- c(lambda1old, lambda2old)
          opt.rout <- optim(par=lambdavec, fn=marg.ll, gr=gr.marg.ll,
                            method="L-BFGS-B", lower=c(0.001, 0.001), 
                            upper=c(Inf, Inf), control=list(fnscale=-1), 
                            e.beta=e.beta, v.beta=v.beta, 
                            e.psi.inv=e.psi.inv, e.psi=e.psi, p=p)
          lambda1 <- opt.rout$par[1]
          lambda2 <- opt.rout$par[2]
        } else if(opt.prior[2]) {
          #optimise lambda2
          opt.rout <- optim(par=lambda2old, fn=marg.ll2, gr=gr.marg.ll2,
                            method="Brent", lower=0.001, upper=10000,
                            control=list(fnscale=-1), e.beta=e.beta, 
                            v.beta=v.beta, e.psi.inv=e.psi.inv, e.psi=e.psi,
                            p=p, lambda1=lambda1old)
          lambda2 <- opt.rout$par
        } else {
          # optimse lambda1
          opt.rout <- optim(par=lambda1old, fn=marg.ll1, gr=gr.marg.ll1,
                            method="Brent", lower=0.001, upper=10000,
                            control=list(fnscale=-1), e.psi=e.psi,
                            p=p, lambda2=lambda2old)
          lambda1 <- opt.rout$par
        }
        
      }
      
      # check the convergence of all parameters (including hyperparameters,
      # which are fixed if we set opt.prior=FALSE)
      conv <- max(c(abs(sigma - sigmaold), abs(mu - muold), abs(ci - ciold), 
                    abs(chi - chiold), abs(lambda1 - lambda1old),
                    abs(lambda2 - lambda2old))) < epsilon
      niter <- niter + 1
    
      # update old parameters to new ones
      sigmaold <- sigma
      muold <- mu
      ciold <- ci
      chiold <- chi
      
      # update old hyperparameters to new ones
      lambda1old <- lambda1
      lambda2old <- lambda2
      
      # update sequence of hyperparameters over the iterations
      lambda1seq <- c(lambda1seq, lambda1)
      lambda2seq <- c(lambda2seq, lambda2)
      
    } else {
      
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
n <- 200
p <- 20
sigma <- matrix(0, ncol=p, nrow=p)
sigma[1:(p/2), 1:(p/2)] <- sigma[(p/2 + 1):p, (p/2 + 1):p] <- 0.9
diag(sigma) <- 1
x <- matrix(rmvnorm(n, sigma=sigma), ncol=p, nrow=n)
beta <- rnorm(p)
prob <- exp(x %*% beta)/(1 + exp(x %*% beta))
y <- as.numeric(runif(n) < prob)

test1 <- envb2.1(x, y, lambda1=1, lambda2=0.05, mustart=rnorm(p), 
               sigmastart=diag(p), model="binomial", opt.prior=c(0, 1),
               maxiter=1000, epsilon=0.001)
test2 <- optL2(y, x, unpenalized=~0, model="logistic", 
               lambda1=test1$lambda1[test1$niter + 1], standardize=FALSE)
plot(test1$mu, test2$fullfit@penalized)
windows()
par(mfrow=c(2, 1))
plot(beta, test1$mu)
plot(beta, test2$fullfit@penalized)
par(mfrow=c(1, 1))

var(beta - test1$mu)
var(beta - test2$fullfit@penalized)
cbind(penalized=test2$fullfit@penalized, ENVB=test1$mu)

image(x=1:20, y=1:20, z=test1$sigma)

data.frame(method=c("ENVB", "penalized"),
           lambda1=c(test1$lambda1[test1$niter + 1], 
                     test2$fullfit@lambda1),
           lambda2=c(test1$lambda2[test1$niter + 1],
                     test2$fullfit@lambda2))



# test 2
set.seed(123)
n <- 200
p <- 20
lambda1 <- 1
lambda2 <- 1
sigma <- matrix(0, ncol=p, nrow=p)
sigma[1:(p/2), 1:(p/2)] <- sigma[(p/2 + 1):p, (p/2 + 1):p] <- 0.9
diag(sigma) <- 1
x <- matrix(rmvnorm(n, sigma=sigma), ncol=p, nrow=n)
beta <- renbeta(p, lambda1, lambda2)
prob <- exp(x %*% beta)/(1 + exp(x %*% beta))
y <- as.numeric(runif(n) < prob)

test1 <- envb2.1(x, y, lambda1=1, lambda2=0.05, mustart=rnorm(p), 
                 sigmastart=diag(p), model="binomial", opt.prior=c(0, 1),
                 maxiter=1000, epsilon=0.001)

test2 <- optL2(y, x, unpenalized=~0, model="logistic", 
               lambda1=test1$lambda1[test1$niter + 1], standardize=FALSE)
plot(test1$mu, test2$fullfit@penalized)
windows()
par(mfrow=c(2, 1))
plot(beta, test1$mu)
plot(beta, test2$fullfit@penalized)
par(mfrow=c(1, 1))

var(beta - test1$mu)
var(beta - test2$fullfit@penalized)
cbind(penalized=test2$fullfit@penalized, ENVB=test1$mu)

image(x=1:20, y=1:20, z=test1$sigma)

data.frame(method=c("ENVB", "penalized"),
           lambda1=c(test1$lambda1[test1$niter + 1], 
                     test2$fullfit@lambda1),
           lambda2=c(test1$lambda2[test1$niter + 1],
                     test2$fullfit@lambda2))




# test 2
set.seed(123)
n <- 50
p <- 100
lambda1 <- 1
lambda2 <- 1
sigma <- matrix(0, ncol=p, nrow=p)
sigma[1:(p/2), 1:(p/2)] <- sigma[(p/2 + 1):p, (p/2 + 1):p] <- 0.9
diag(sigma) <- 1
x <- matrix(rmvnorm(n, sigma=sigma), ncol=p, nrow=n)
beta <- renbeta(p, lambda1, lambda2)
prob <- exp(x %*% beta)/(1 + exp(x %*% beta))
y <- as.numeric(runif(n) < prob)

test1 <- envb2.1(x, y, lambda1=1, lambda2=0.05, mustart=rnorm(p), 
                 sigmastart=diag(p), model="binomial", opt.prior=c(1, 1),
                 maxiter=1000, epsilon=0.001)
