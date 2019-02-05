# function to make predictions with the gren object
predict.gren <- function(object, newx, unpenalized, s=NULL, 
                         type=c("VB", "freq")) {
  if(is.null(object$freq) & type=="freq") {
    stop("object should contain a frequentist model fit")
  } else if(!is.null(s) & !is.numeric(s)) {
    stop("s is either NULL or a numeric")
  } else if(!(type %in% c("VB", "freq"))) {
    stop("type is either VB or freq")
  }
  
  x <- cbind(unpenalized, newx)
  if(type=="VB") {
    if(object$args$intercept) {
      x <- cbind(1, x)
    } 
    prob <- as.numeric(1/(1 + exp(-x %*% as.matrix(object$vb.post$mu))))
  } else {
    prob <- predict(object$freq, x, s=s, type="response")
  }
  return(prob)
}