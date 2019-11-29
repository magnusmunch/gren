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

<<<<<<< HEAD
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


=======
>>>>>>> 31db1da7fc892df40e05ce4287d82f38e942befc
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

<<<<<<< HEAD
=======
vb.gren <- function(x, y, m=rep(1, nrow(x)), alpha=0.5, lambda=NULL, 
                    lambdamult, intercept=TRUE, 
                    control=list(maxit=100, epsilon=0.001)) {
  
  n <- nrow(x)
  p <- ncol(x)
  
  lambda1 <- alpha*lambda
  lambda2 <- 0.5*(1 - alpha)*lambda
  
  # generate starting values
  phi <- 0.25*lambda1^2/lambda2
  fit.start <- glmnet::glmnet(x, y, "binomial", alpha=0, 
                              lambda=lambda/n, intercept=intercept)
  pred.start <- as.numeric(predict(fit.start, newx=x, type="response"))
  startparam <- gren:::est_param(x, matrix(1, nrow=2), y - 0.5*m, m, n, p, 
                                 pred.start, phi, rep(phi, p), lambda2, 
                                 lambdamult, lambdamult, intercept, FALSE, 
                                 FALSE, FALSE, TRUE)
  dsigmaold <- as.numeric(startparam$dsigma)
  muold <- as.numeric(startparam$mu)
  ciold <- as.numeric(startparam$ci)
  chiold <- as.numeric(startparam$chi)
  
  conv <- FALSE
  iter <- 0
  while(!conv & (iter < control$maxit)) {
    iter <- iter + 1
    
    # estimating new model parameters
    newparam <- gren:::est_param(x, matrix(1, nrow=2), y - 0.5*m, m, n, p, 
                                 ciold, phi, chiold, lambda2, lambdamult, 
                                 lambdamult, intercept, FALSE, FALSE, 
                                 FALSE, FALSE)
    
    dsigma <- as.numeric(newparam$dsigma)
    mu <- as.numeric(newparam$mu)
    ci <- as.numeric(newparam$ci)
    chi <- as.numeric(newparam$chi)
    
    # checking convergence of inner loop
    conv <- max(c(abs((mu - muold)/ifelse(muold==0, muold + 0.00001, muold)), 
                  abs((dsigma - dsigmaold)/
                        ifelse(dsigmaold==0, dsigmaold + 0.00001, 
                               dsigmaold)))) < control$epsilon
    
    # updating vb parameters
    muold <- mu
    dsigmaold <- dsigma
    ciold <- ci
    chiold <- chi
  }
  
  # obtaining the full posterior
  newparam <- gren:::est_param(x, matrix(1, nrow=2), y - 0.5*m, m, n, p, 
                               ciold, phi, chiold, lambda2, lambdamult, 
                               lambdamult, intercept, FALSE, TRUE, 
                               FALSE, FALSE)
  Sigma <- newparam$sigma
  mu <- as.numeric(newparam$mu)
  ci <- as.numeric(newparam$ci)
  chi <- as.numeric(newparam$chi)
  
  # calculating the linear reponse variances for beta
  A <- 1/(lambda2*lambdamult*(1 + 0.5*lambda1/sqrt(lambda2*chi)))
  if(intercept) {
    B <- mu[-1]^2*(1 + 0.5*lambda1/sqrt(lambda2*chi))/
      (lambda2*lambdamult*(1 - 0.25*lambda1*chi^(-3/2)/sqrt(lambda2))^2)
    C <- -0.25/chi^2 + lambda1*chi^(-3/2)/(16*sqrt(lambda2)) + 0.5*lambda2*
      lambdamult*(1 + 3*lambda1*chi^(-5/2)/(8*sqrt(lambda2)))*
      (mu[-1]^2 + dsigma[-1])
    D <- 0.25*lambda2^2*lambdamult^2*
      (1 - 0.25*lambda1*chi^(-3/2)/sqrt(lambda2))/diag(solve(Sigma))[-1]^2
  } else {
    B <- mu^2*(1 + 0.5*lambda1/sqrt(lambda2*chi))/
      (lambda2*lambdamult*(1 - 0.25*lambda1*chi^(-3/2)/sqrt(lambda2))^2)
    C <- -0.25/chi^2 + lambda1*chi^(-3/2)/(16*sqrt(lambda2)) + 0.5*lambda2*
      lambdamult*(1 + 3*lambda1*chi^(-5/2)/(8*sqrt(lambda2)))*(mu^2 + dsigma)
    D <- 0.25*lambda2^2*lambdamult^2*
      (1 - 0.25*lambda1*chi^(-3/2)/sqrt(lambda2))/diag(solve(Sigma))^2
  }
  LRdsigma <- A + A/(B*C - B*D - 1)
  
  out <- list(posterior=list(sigma=Sigma, mu=mu, ci=ci, chi=chi, 
                             LRdsigma=LRdsigma))
  return(out)
  
}






>>>>>>> 31db1da7fc892df40e05ce4287d82f38e942befc
                