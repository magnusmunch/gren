### gives back the number of selected features in grpreg for specific lambda
froot <- function(lambda, psel, inp) {
  # fit the model with specific lambda
  fit <- do.call("grpreg", c(inp, lambda=lambda))
  
  # calculate the number of selected features
  out <- psel - sum(fit$beta[-1, ]!=0)
  return(out)
}

## grpreg that selects specific number of features
sel.grpreg <- function(X, y, group=1:ncol(X), 
                       penalty=c("grLasso", "grMCP", "grSCAD", "gel", "cMCP"), 
                       family=c("gaussian", "binomial", "poisson"),
                       lambda.min={if (nrow(X) > ncol(X)) 1e-4 else .05},
                       log.lambda = TRUE, alpha=1, eps=1e-4, max.iter=10000, 
                       dfmax=p, gmax=length(unique(group)), 
                       gamma=ifelse(penalty == "grSCAD", 4, 3), tau = 1/3, 
                       group.multiplier=rep(1, length(unique(group))), 
                       warn=TRUE, returnX = FALSE, psel) {
  
  # store the arguments
  # args <- as.list(sys.call())
  inp <- lapply(as.list(match.call(definition=sys.function(),
                                   call=sys.call()))[-1], eval)
  psel <- inp$psel
  inp <- inp[-length(inp)]

  # find range of lambdas using built-in grpreg functions
  bilevel <- strtrim(penalty, 2) != "gr"
  yy <- grpreg:::newY(y, family)
  XG <- grpreg:::newXG(X, group, group.multiplier, attr(yy, "m"), bilevel)
  init.lambda <- grpreg:::setupLambda(XG$X, yy, XG$g, family, penalty, alpha, 
                                      lambda.min, log.lambda, 2, XG$m)
  
  # check whether smallest lambda gives largest psel, otherwise, this psel
  # not possible
  init.fit <- do.call("grpreg", c(inp, lambda=init.lambda, nlambda=2))
  pselfit <- psel[psel <= max(colSums(as.matrix(init.fit$beta[-1, ])!=0))]
  
  # if no model has fewer selected features than any psel, we return initial
  # model
  if(length(pselfit) > 0) {
    # find the lambdas that give the specified psels
    sel.lambda <- numeric(length(pselfit))
    for(i in 1:length(pselfit)) {
      fit.root <- suppressWarnings(uniroot(froot, interval=init.lambda, 
                                           psel=pselfit[i], inp=inp))
      sel.lambda[i] <- fit.root$root
    }
  
    # final fit with the found lambdas
    inp[["lambda"]] <- sel.lambda
    inp[["nlambda"]] <- max(length(sel.lambda), 2)
  
    fit <- do.call("grpreg", inp)
  } else {
    fit <- init.fit
  }
  return(fit)
  
}

n <- 100
p <- 20
x <- matrix(rnorm(n*p), ncol=p, nrow=n)
beta <- -9.5:9.5
y <- rbinom(n, 1, as.numeric(1/(1 + exp(-x %*% beta))))
m <- rep(1, n)
alpha <- 0.5
lambda <- 1
lambdamult <- rep(1, p)
intercept <- FALSE
control <- list(nsamples=1000)

fit.mcmc <- mcmc.gren(x, y, m=rep(1, nrow(x)), 
                      alpha=0.5, lambda=1, lambdamult=rep(1, p), 
                      intercept=FALSE, control=list(nsamples=1000))
plot(c(beta), apply(fit.mcmc$samples$beta, 1, mean))


# mcmc version of VB steps
# devtools::install_version("BayesLogit", "0.2-0")
mcmc.gren <- function(x, y, m=rep(1, nrow(x)), alpha=0.5, lambda=NULL, 
                      lambdamult, intercept=TRUE, control=list(nsamples=1000)) {
  
  n <- nrow(x)
  p <- ncol(x)
  
  lambda1 <- alpha*lambda
  lambda2 <- 0.5*(1 - alpha)*lambda
  
  # create matrices for samples
  seq.beta <- matrix(NA, ncol=control$nsamples, nrow=p + intercept)
  seq.psi <- matrix(NA, ncol=control$nsamples, nrow=p)
  seq.omega <- matrix(NA, ncol=control$nsamples, nrow=n)
  
  # starting values
  beta <- as.numeric(coef(glmnet::glmnet(x, y, "binomial", alpha=0, 
                                         lambda=lambda/n, intercept=intercept)))
  if(!intercept) {beta <- beta[-(p + 1)]}
  psi <- rep(1, p)
  Sigma <- matrix(NA, ncol=p + intercept, nrow=p + intercept)
  
  # creating mcmc samples
  for(k in 1:control$nsamples) {
    h <- lambda2*lambdamult + lambda2*lambdamult/psi
    hinvxtr <- t(x)/h
    if(intercept) {
      omega <- BayesLogit::rpg(n, m, abs(as.numeric(cbind(1, x) %*% beta)))
      
      # calculating Sigma with block inversion
      omtrx <- t(matrix(omega)) %*% x
      Dinv <- diag(1/h) - hinvxtr %*% solve(diag(1/omega) + x %*% hinvxtr) %*% 
        t(hinvxtr)
      omtrxDinv <- omtrx %*% Dinv
      val <- 1/as.numeric(sum(omega) - omtrxDinv %*% t(omtrx))
      Sigma[1, 1] <- val
      Sigma[1, 2:(p + 1)] <- - val*omtrxDinv
      Sigma[2:(p + 1), 1] <- t(Sigma[1, 2:(p + 1)])
      Sigma[2:(p + 1), 2:(p + 1)] <- Dinv - Sigma[2:(p + 1), 1] %*% omtrxDinv
      beta <- as.numeric(rmvnorm(1, Sigma %*% rbind(1, t(x)) %*% (y - m/2), 
                                 Sigma))
      psi <- 1/statmod::rinvgauss(p, lambda/(2*lambda*sqrt(lambdamult)*
                                               abs(beta[-1])),
                                  rep(lambda1^2/(4*lambda2), p))
    } else {
      omega <- BayesLogit::rpg(n, m, abs(as.numeric(x %*% beta)))
      Sigma <- diag(1/h) - hinvxtr %*% solve(diag(1/omega) + x %*% hinvxtr) %*% 
        t(hinvxtr)
      beta <- as.numeric(rmvnorm(1, Sigma %*% t(x) %*% (y - m/2), Sigma))
      psi <- 1/statmod::rinvgauss(p, lambda/(2*lambda*sqrt(lambdamult)*abs(beta)),
                                  rep(lambda1^2/(4*lambda2), p))
    }
    
    seq.omega[, k] <- omega
    seq.beta[, k] <- beta
    seq.psi[, k] <- psi
  }
  
  out <- list(samples=list(omega=seq.omega, beta=seq.beta, psi=seq.psi))
  return(out)
  
}

                