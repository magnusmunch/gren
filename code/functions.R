### function to find a specific number of features
froot <- function(lambda, psel, x, y, part, alpha, type) {
  if(type=="sgl") {
    fit <- SGL(list(x=x, y=y), part, type="logit", alpha=alpha, standardize=FALSE,
               lambdas=lambda, nlam=1, min.frac=0.001)
    out <- psel - sum(fit$beta!=0)
  } else if(type=="cmcp") {
    fit <- grpreg(x, y, part, "cMCP", "binomial", lambda=lambda, alpha=alpha)
    out <- psel - sum(fit$beta[-1, ]!=0)
  } else if(type=="gel") {
    fit <- grpreg(x, y, part, "gel", "binomial", lambda=lambda, alpha=alpha)
    out <- psel - sum(fit$beta[-1, ]!=0)
  }
  return(out)
}

sel.grpreg <- function(X, y, group=1:ncol(X), 
                       penalty=c("grLasso", "grMCP", "grSCAD", "gel", "cMCP"), 
                       family=c("gaussian", "binomial", "poisson"),
                       nlambda=100, lambda, 
                       lambda.min={if (nrow(X) > ncol(X)) 1e-4 else .05},
                       log.lambda = TRUE, alpha=1, eps=1e-4, max.iter=10000, 
                       dfmax=p, gmax=length(unique(group)), 
                       gamma=ifelse(penalty == "grSCAD", 4, 3), tau = 1/3, 
                       group.multiplier, warn=TRUE, returnX = FALSE) {
  
  
  
}
sapply(psel, function(csel) {
  fit.root <- uniroot(froot, range(test.init$lambda), psel=csel, x, y, part3,
                         alpha=0.05, type="gelasso", maxiter=100);
  fit <- grpreg(x, y, part3, "gel", "binomial", lambda=test.root$root,
                       alpha=0.05);
    return(sum(test.fit$beta[-1, ]!=0))})
  
  
  
  
  