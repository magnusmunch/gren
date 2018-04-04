# retrieve coefficients from gren estimated model
coef.gren <- function(object, s=NULL, type="freq") {
  if(!is.null(s) & !is.numeric(s)) {
    stop("s is either NULL or a numeric")
  } else if(is.null(s)) {
    s <- object$lambda
  } else if(!any(type %in% c("VB", "freq"))) {
    stop("type is either VB or freq")
  }
  
  if(type=="VB") {
    coefs <- list(mu=object$vb.post$mu, sigma=object$vb.post$sigma)
  } else {
    coefs <- coef(object$freq, s=s)
  }
  return(coefs)
}