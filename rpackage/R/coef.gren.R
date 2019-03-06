# retrieve coefficients from gren estimated model
coef.gren <- function(object, s=NULL, type=c("groupreg", "regular")) {
  if(class(object)!="gren") {
    stop("object should be a gren fit") 
  } else if (is.null(s) & !is.numeric(s)) {
    stop("s is either NULL or a numeric")
  } else (!any(type %in% c("groupreg", "regular"))) {
    stop("type is groupreg or regular")
  }
  
  if(type=="groupreg") {
    coefs <- coef(object$freq.model$groupreg, s=s)
  } else {
    coefs <- coef(object$freq.model$regular, s=s)
  } 
  return(coefs)
}