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

                