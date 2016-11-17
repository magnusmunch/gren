##############################  preamble  #############################
# code belonging to abstract_SMPGD_2017_V01.pdf                       #
# version: 01                                                         #
# author: Magnus Münch                                                #
# created: 15-11-2016                                                 #
# last edited: 15-11-2016                                             #
#######################################################################

###############################  notes  ###############################
#                                                                     #
#######################################################################

### paths
path.results <- "C:/Users/Magnus/Documents/phd/ENVB/abstract_SMPGD_2017/results/"

### libraries
library(glmnet)
library(penalized)
library(mvtnorm)
library(GRridge)
library(pROC)

### functions
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

marg.ll1.2g <- function(lambda, e.beta, v.beta, e.psi.inv, e.psi, G, sizes, modmat) {
  # marginal likelihood in lambda1 and lambda2g
  
  lambda1 <- lambda[1]
  lambda2 <- lambda[2:(G + 1)]
  lambda2vec <- rep(lambda2, times=sizes)
  part1 <- p*log(lambda1)
  part2 <- 0.5*sum(lambda2vec*(v.beta + e.beta^2)*(1 + e.psi.inv))
  part3 <- lambda1^2*sum((e.psi + 1)/(8*lambda2vec))
  part4 <- sum(sizes*pnorm(-lambda1/(sqrt(4*lambda2)), log.p=TRUE))
  ll <- part1 - part2 - part3 - part4
  
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

gr.marg.ll1.2g <- function(lambda, e.beta, v.beta, e.psi.inv, e.psi, G, sizes, modmat) {
  # gradient of marginal likelihood in lambda1 and lambda2g
  lambda1 <- lambda[1]
  lambda2 <- lambda[2:(G + 1)]
  part1.1 <- p/lambda1
  part1.2 <- lambda1*sum((t(as.matrix(e.psi + 1)) %*% modmat)/lambda2)/4
  part1.3 <- sum(sizes*dnorm(lambda1/sqrt(4*lambda2))/
    (sqrt(4*lambda2)*pnorm(-lambda1/sqrt(4*lambda2))))
  comp1 <- part1.1 - part1.2 + part1.3
  part2.1 <- lambda1^2*(t(as.matrix(e.psi + 1)) %*% modmat)/(8*lambda2^2)
  part2.2 <- 0.5*as.numeric((t(as.matrix((v.beta + e.beta^2)*(1 + e.psi.inv))) %*% modmat))
  part2.3 <- lambda1*sizes*dnorm(lambda1/sqrt(4*lambda2))/
    (4*lambda2^(1.5)*pnorm(-lambda1/sqrt(4*lambda2)))
  comp2 <- part2.1 - part2.2 - part2.3
  return(c(comp1, comp2))
  
}

# the fitting function
envb2 <- function(x, y, groups, lambda1=NULL, lambda2=NULL, mustart=NULL, 
                  sigmastart=NULL, model=c("binomial", "multinomial"), opt.prior=FALSE, 
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
  #       or one lambda1 and multiple lambda2's with c(2, 1)
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
              "opt.prior should be 0, 1, c(0, 1), c(1, 0) and the logical equivalents thereof, or c(2, 1)",
              "inv should be NULL, svd or woodbury",
              "maxiter is not an integer", "epsilon is not numeric",
              "trace is not a logical")
  checks <- c(is.numeric(x), is.numeric(groups), length(groups)==ncol(x),
              is.null(lambda1) | is.numeric(lambda1), 
              is.null(lambda2) | is.numeric(lambda2), 
              is.null(mustart) | is.numeric(mustart) | is.null(sigmastart) |
                is.numeric(sigmastart), model %in% c("binomial", "multinomial"),
              opt.prior==0 | opt.prior==1 | opt.prior==c(1, 0) | opt.prior==c(0, 1) | opt.prior==c(2, 1),
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
      
      fit <- optL2(y[, 2], x, unpenalized=~0, lambda1=0, model="logistic")
      lambda1old <- lambda1seq <- lambda1 <- rep_len(fit$lambda/2, G)
      lambda2old <- lambda2seq <- lambda2 <- rep_len(fit$lambda/2, G)
      mustart <- as.numeric(fit$fullfit@penalized)
      
    } else {
      
      lambda1old <- lambda1seq <- lambda1 <- rep_len(lambda1*(!is.null(lambda1)) + lambda2*(!is.null(lambda2)), G)
      lambda2old <- lambda2seq <- lambda2 <- rep_len(lambda1*(!is.null(lambda1)) + lambda2*(!is.null(lambda2)), G)
      
    }  
    
  } else {
    
    # if provided but length < G, we recycle
    lambda1old <- lambda1seq <- lambda1 <- rep_len(lambda1, length.out=G)
    lambda2old <- lambda2seq <- lambda2 <- rep_len(lambda2, length.out=G)
    
  }
  lambda1vec <- rep(lambda1old, times=sizes)
  lambda2vec <- rep(lambda2old, times=sizes)
  phi <- lambda1vec^2/(4*lambda2vec) # this is vector of length p
  
  # check whether starting values provided, otherwise estimate with penalized
  # if we estimated lambda1 and/or lambda2, then starting values are 
  # estimated earlier
  if(model=="binomial") {
    if(is.null(mustart)) {

      fit <- penalized(y[, 2], x, model="logistic", unpenalized=~0, lambda1=0, 
                       lambda2=lambda1vec + lambda2vec/2, standardize=FALSE)
      mu <- muold <- as.numeric(fit@penalized)

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
      invmat <- 1/(d^2 + 2*mean(lambda1vec) + mean(lambda2vec))
      part1 <- invmat^2*d^2
      sigma <- sigmaold <- t(t(V)*part1) %*% t(V)
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
    
    niter <- niter + 1
    
    # if we want the iteration trace, print current iteration
    if(trace) {
      cat("\r", "Iteration: ", niter, ", ", 
          paste("lambda1=[", paste(round(lambda1, 2), collapse=", "), 
                sep=""), "], ", 
          paste("lambda2=[", paste(round(lambda2, 2), collapse=", "), 
                sep=""), "]", sep="")
    }
    
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
      
      # check the convergence of the model parameters
      conv <- max(abs(c(sigma - sigmaold, mu - muold))) < epsilon
      
      # update old parameters to new ones
      sigmaold <- sigma
      muold <- mu
      ciold <- ci
      chiold <- chi
      
      if((sum(opt.prior) > 0) & !conv) {
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
        } else if(opt.prior[1]==2) {
          # optimise one lambda1 and multiple lambda2
          lambdaold <- c(lambda1old[1], lambda2old)
          opt.rout <- optim(par=lambdaold, fn=marg.ll1.2g, gr=gr.marg.ll1.2g,
                            method="L-BFGS-B", lower=rep(0.001, G + 1), 
                            upper=rep(Inf, G + 1), control=list(fnscale=-1), 
                            e.beta=e.beta, v.beta=v.beta, e.psi.inv=e.psi.inv, 
                            e.psi=e.psi, G=G, sizes=sizes, modmat=modmat)
          lambda1 <- opt.rout$par[1]
          lambda2 <- opt.rout$par[2:(G + 1)]
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
        
        # update old hyperparameters to new ones
        lambda1old <- lambda1
        lambda2old <- lambda2
        lambda1vec <- ifelse(opt.prior[1]==2, rep(lambda1old, times=p), 
                             rep(lambda1old, sizes))
        lambda2vec <- rep(lambda2old, times=sizes)
        
        # update sequence of hyperparameters over the iterations
        lambda1seq <- cbind(lambda1seq, lambda1)
        lambda2seq <- cbind(lambda2seq, lambda2)
        
      } else {
        conv <- conv1
      }
      
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

# function to sparsify coefficient vector in simulations
sparsify <- function(vec, frac){
  if(frac==0){
    return(vec)
  } else {
    N0 <- round(frac*length(vec))
    vecnew <- vec
    vecnew[1:N0] <- 0
    return(vecnew)
  }
}

### simulations
set.seed(123)
n <- 100           # Nb of observations    
ntest <- 1000          ## Nb of test set observations     
p <- 200            # Nb of variables per group
G <- 10            # Nb of groups
meanBeta <- 0.01   # Beta variances per group are VarBeta*(1:G); use this for CorX=0.5
CorX <- 0.5       # correlation within variable block
Nblock <- 10*G    # number of correlation blocks
settings <- c(ntrain=n, p=p, G=G, meanBeta=meanBeta, CorX=CorX, Nblock=Nblock)

nrep <- 1  #number of repeats per simulation setting
facvec <- c(1.3, 1.6, 2)   #tunes how much weaker each next group is. The '2' means that the second group is twice as weak as the first, etc
fractvec <- c(0, 0.7, 0.9) #tunes the sparsity per group. E.g. 0.9 means that 9/10 betas in a group are set to 0

aucmat <- c()
briermat <- c()
msemat <- c()
for(fac in facvec){
  for(fract in fractvec){
    for(reptit in 1:nrep){
      print(paste("fac=", fac))
      print(paste("fract=", fract))
      print(paste("repeat=", reptit))
      reps <- rev(sapply(0:(G - 1), function(i) {fac^(-i)}))
      meansB <- rep(reps, each=p)*meanBeta/mean(rep(reps, each=p))
      Beta <- meansB
      Beta <- rev(sparsify(Beta, frac=fract))
      
      ### FITTING THE MODELS
      pblock <- G*p/Nblock
      grs <- rep(1:G,each=p)
      P <- G*p #Complete number of variables
      X <- Reduce(cbind, lapply(1:Nblock, function(z) {
        matrix(rep(rnorm(n, sd=sqrt(CorX/(1 - CorX))), times=pblock), n, pblock)})) + matrix(rnorm(n*G*p), n, G*p)
      X <- t((t(X) - apply(t(X), 1, mean))/apply(t(X), 1, sd))
      
      lpred <- X %*% Beta 
      logisticintercept <- 0
      prob <- 1/(1 + exp(-(lpred + logisticintercept)))
      Y <- rbinom(length(prob), 1, prob)
      
      # ENVB
      vbSim <- envb2(X, Y, grs, lambda1=1, lambda2=1, mustart=NULL, sigmastart=NULL, 
                     model="binomial", opt.prior=c(2, 1), inv="woodbury", maxiter=1000, epsilon=1e-06, 
                     trace=TRUE)
      
      # GRridge
      groups <- CreatePartition(grs, grsize=p, uniform=T, decreasing=F)
      partsim <- list(grouping=groups)
      grSim <- grridge(t(X), Y, unpenal=~0, partsim, savepredobj="all", innfold=10, method="stable")
      
      # calculating mse
      grMse <- var((grSim$betas - Beta)^2)
      vbMse <- var((vbSim$mu - Beta)^2)
      mses <- c(fac, fract, repit, grMse, vbMse)
      
      ### TESTING THE MODELS
      # making the test data
      Xtest <- Reduce(cbind, lapply(1:Nblock, function(z) {
        matrix(rep(rnorm(ntest, sd=sqrt(CorX/(1 - CorX))), times=pblock), ntest, pblock)})) + 
        matrix(rnorm(ntest*G*p), ntest, G*p)
      lpredtest <- Xtest %*% Beta 
      logisticintercept <- 0
      
      probtest <- 1/(1 + exp(-(lpredtest + logisticintercept)))
      Ytest <- rbinom(length(probtest), 1, probtest)
      
      # making predictions
      vbPred <- 1/(1 + exp(-(Xtest %*% vbSim$mu)))
      grPred <- predict.grridge(grSim, t(Xtest))
      
      # calculating brier residuals and auc
      grBrier <- mean((grPred[, 2] - probtest)^2)
      vbBrier <- mean((vbPred - probtest)^2)
      briers <- c(fac, fract, reptit, grBrier, vbBrier)
      
      cutoffs <- rev(seq(0, 1, by=0.005))
      grRoc <- GRridge::roc(probs=as.numeric(grPred[, 2]), true=Ytest, cutoffs) #ridge, sel
      vbRoc <- GRridge::roc(probs=as.numeric(vbPred), true=Ytest, cutoffs)
      aucs <- c(fac, fract, reptit, GRridge::auc(grRoc), GRridge::auc(vbRoc))
      
      msemat <- rbind(msemat, mses)
      print("mses")
      print(msemat)
      briermat <- rbind(briermat, briers)
      print("briers")
      print(briermat)
      aucmat <- rbind(aucmat, aucs)
      print("aucs")
      
    }
  }   
}

test1 <- matrix(c(1:5), ncol=5, nrow=5)
t(t(test1)*c(1:5))













