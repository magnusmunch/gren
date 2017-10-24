# library(devtools)
# devtools::check("/Users/magnusmunch/Downloads/myglmnet")
# install.packages("/Users/magnusmunch/Downloads/glmnet_2.0-13.tar.gz", repos=NULL, type="source")

library(glmnet)
library(penalized)

set.seed(123)
nreps <- 100
cormat <- matrix(NA, ncol=2, nrow=nreps)
for(r in 1:nreps) {
  n <- 100
  p <- 200
  x <- matrix(rnorm(n*p), ncol=p, nrow=n)
  beta <- rnorm(p, 0, 0.1)
  y <- rbinom(n, 1, as.numeric(exp(x %*% beta)/(1 + exp(x %*% beta))))
  
  lambda1 <- 0.5
  lambda2 <- 0.5
  lambdag <- exp(seq(-2, 2, length.out=p))
  fit.penal1 <- penalized(y, x, ~1, lambda1=sqrt(lambdag)*lambda1, lambda2=lambdag*lambda2, model="logistic")
  fit.penal2 <- penalized(y, x, ~1, lambda1=lambdag*lambda1, lambda2=lambdag*lambda2, model="logistic")
  fit.glmnet <- glmnet(x, y, family="binomial", standardize=FALSE, thresh=1e-10,
                       alpha=lambda1/(2*lambda2 + lambda1), lambda=(2*lambda2 + lambda1)/n,
                       penalty.factor=lambdag)
  
  fit.lm1 <- lm(as.numeric(coef(fit.glmnet)) ~ as.numeric(coef(fit.penal1, which="all")))
  fit.lm2 <- lm(as.numeric(coef(fit.glmnet)) ~ as.numeric(coef(fit.penal2, which="all")))
  
  cormat[r, ] <- c(cor(as.numeric(coef(fit.glmnet)), as.numeric(coef(fit.penal1, which="all"))),
                   cor(as.numeric(coef(fit.glmnet)), as.numeric(coef(fit.penal2, which="all"))))
  
}

apply(cormat, 2, median)
apply(cormat, 2, sd)
