cv.gren <- function(x, y, m=rep(1, nrow(x)), partitions, keep.pred=TRUE, 
                    fix.lambda=FALSE, nfolds.out=nrow(x), foldid.out=NULL,
                    type.measure=c("auc", "deviance", "class.error"), ...) {
  
  argum <- formals(gren)
  
  # change some input for convenience
  if(is.data.frame(x)) {
    names.xr <- colnames(x)
    x <- as.matrix(x)
  } else {
    names.xr <- NULL
  }
  if(is.data.frame(unpenalized)) {
    names.xu <- colnames(unpenalized)
    if(intercept) {names.xu <- c("Intercept", names.xu)}
    unpenalized <- as.matrix(unpenalized)
  } else {
    names.xu <- NULL
  }
  if(is.vector(partitions) & is.atomic(partitions)) {
    partitions <- list(partition1=partitions)
  }
  
  if(!is.numeric(x) | !is.numeric(y) | !is.numeric(m)) {
    stop("only numerical input data is supported at the moment")
  } else if(!(is.null(ncol(y)) | ncol(y)==2)) {
    stop("y is either a vector or matrix with two columns")
  } else if(!is.list(partitions) & !is.numeric(partitions)) {
    stop("partitions should be either a list of partitions or one partition as 
         a numeric vector")
  } else if(!is.logical(keep.pred)) {
    stop("keep.pred should be either TRUE or FALSE")
  } else if(!is.logical(fix.lambda)) {
    stop("fix.lambda should be either TRUE or FALSE")
  } else if(!is.logical(keep.pred)) {
    stop("keep.pred should be either TRUE or FALSE")
  } else if(!is.character(type.measure)) {
    stop("type.measure should be one or more of: auc, deviance, class.error")
  } else if(!(all(type.measure %in% c("auc", "deviance", "class.error")))) {
    stop("type.measure should be one or more of: auc, deviance, class.error")
  } else if(any(m!=1) & any(type.measure %in% c("auc", "class.error"))) {
    stop("auc or misclassification error are currently only possible with 
         binary classification problems")
  }
  
  # extra check for data
  if(nrow(x)!=nrow(ymat)) {
    stop("number of observations in y not equal to the number of 
         observations in x")
  } 
  
  # determine the folds
  if(is.null(foldid.out)) {
    rest <- n %% nfolds.out
    foldsize <- c(rep(n %/% nfolds.out + as.numeric(rest!=0), times=rest),
                  rep(n %/% nfolds.out, times=nfolds.out - rest))
    foldid.out <- sample(rep(1:nfolds.out, times=foldsize))
  } else if(length(foldid.out)!=n | !all(foldid.out %in% c(1:n))) {
    stop("foldid.out must be a vector of length n, containing only whole 
         numbers from 1 to nrow(x)")
  }
  nfolds.out <- length(unique(foldid.out))
  
  # if fix.lambda=TRUE, we estimate the global lambda only once
  if(fix.lambda) {
    if(trace) {cat("\r", "Estimating fixed global lambda by cross-validation", 
                   sep="")}
    srt <- proc.time()[3]
    r <- ncol(x)
    u <- ifelse(is.null(unpenalized), 0, ncol(unpenalized) - intercept)
    cv.fit <- cv.glmnet(x, ymat, family="binomial", alpha=alpha, 
                        standardize=FALSE, intercept=intercept, 
                        foldid=foldid.out, grouped=FALSE,
                        penalty.factor=c(rep(0, u, rep(1, r))))
    lambda <- cv.fit$lambda.min
    cv.time <- proc.time()[3] - srt
    
    if(trace) {cat("\n", "Global lambda estimated at ", round(lambda, 2), 
                   " in ", round(cv.time, 2), " seconds", sep="")
    }
  }
  
  for(k in 1:nfolds.out) {
    ytrain <- ymat[foldid.out!=k, ]
    xtrain <- x[foldid.out!=k, ]
    xtest <- as.matrix(x[foldid.out==k, ], ncol=r)
    utrain <- unpenalized
  }
  
  
}