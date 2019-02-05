##############################  preamble  #############################
# code to test the Rcpp functions                                     #
# version: 01                                                         #
# author: Magnus Münch                                                #
# created: 25-11-2016                                                 #
# last edited: 05-12-2016                                             #
#######################################################################

###############################  notes  ###############################
#                                                                     #
#######################################################################

### paths
path.code <- "C:/Users/Magnus/Documents/phd/ENVB/code/"

### libraries
library(Rcpp)
library(penalized)
library(mvtnorm)
library(GeneralizedHyperbolic)

### functions
# c++ version
sourceCpp(paste(path.code, "ENVB2.cpp", sep=""))
sourceCpp(paste(path.code, "ENgibbs.cpp", sep=""))

# R version
est_param2 <- function(x, kappa, m, n, p, ciold, phi, chiold, lambda2, intercept) {
  
  h <- (1 + sqrt(phi/chiold))*lambda2
  om <- (0.5*m/ciold)*tanh(ciold/2)
  omsq <- sqrt(om)
  hinv <- 1/h
  sigma <- matrix(0, ncol=p + 1, nrow=p + 1)
  xaug <- cbind(1, x)
  
  # using the Woodbury identity 
  if(intercept) {
    invtrom <- 1/sum(om)
    Omadj <- diag(om) - invtrom*as.matrix(om) %*% t(as.matrix(om))
    xhinv <- t(t(x)*hinv)
    Ainv <- diag(hinv) - t(xhinv) %*% Omadj %*% solve(diag(n) + x %*% t(xhinv) %*% Omadj) %*% xhinv
    ainvxom <- Ainv %*% t(x) %*% as.matrix(om)
    sigma[1, 1] <- invtrom + invtrom^2*t(as.matrix(om)) %*% x %*% ainvxom
    sigma[2:(p + 1), 1] <- sigma[1, 2:(p + 1)] <- -invtrom*ainvxom
    sigma[2:(p + 1), 2:(p + 1)] <- Ainv
    mu <- as.numeric(sigma %*% (t(xaug) %*% kappa))
    ci <- as.numeric(sqrt(colSums(t(xaug) * (sigma %*% t(xaug))) + (colSums(t(xaug)*mu))^2))
    chi <- as.numeric(lambda2*(diag(sigma[-1, -1]) + mu[-1]^2))
  } else {
    V <- x*omsq
    vhinv <- t(t(V)*hinv)
    sigma <- diag(hinv) - t(vhinv) %*% solve(diag(n) + vhinv %*% t(V)) %*% vhinv
    mu <- as.numeric(sigma %*% (t(x) %*% kappa))
    ci <- as.numeric(sqrt(colSums(t(x) * (sigma %*% t(x))) + (colSums(t(x)*mu))^2))
    chi <- as.numeric(lambda2*(diag(sigma) + mu^2))
  }
  
  return(list(sigma=sigma, mu=mu, ci=ci, chi=chi))
  
}

### testing
# testing whether it works 
p <- 1000
n <- 100
m <- rep(1, n)
x <- matrix(rnorm(n*p), ncol=p, nrow=n)
y <- rbinom(n, 1, 0.5)
kappa <- y - m/2
ciold <- rchisq(n, 1)
phi <- chiold <- rchisq(p, 1)
lambda2 <- seq(0.1, 2, length.out=p)
test1 <- est_param(x, kappa, m, n, p, ciold, phi, chiold, lambda2, FALSE)
test2 <- est_param(x, kappa, m, n, p, ciold, phi, chiold, lambda2, TRUE)
test3 <- est_param2(x, kappa, m, n, p, ciold, phi, chiold, lambda2, FALSE)
test4 <- est_param2(x, kappa, m, n, p, ciold, phi, chiold, lambda2, TRUE)
sum(round(test1$sigma, 8)!=round(test3$sigma, 8))
sum(round(test2$sigma, 8)!=round(test4$sigma, 8))

# testing timing versus R estimation
n <- 100
(pseq <- floor(exp(seq(6, 8.5, 0.1))))
out <- numeric(3)
for(p in pseq) {
  m <- rep(1, n)
  x <- matrix(rnorm(n*p), ncol=p, nrow=n)
  y <- rbinom(n, 1, 0.5)
  kappa <- y - m/2
  ciold <- rchisq(n, 1)
  phi <- chiold <- rchisq(p, 1)
  lambda2 <- seq(0.1, 2, length.out=p)
  
  time1 <- system.time(
    test1 <- est_param(x, kappa, m, n, p, ciold, phi, chiold, lambda2, TRUE)
  )
  time2 <- system.time(
    test2 <- est_param2(x, kappa, m, n, p, ciold, phi, chiold, lambda2, TRUE)
  )  
  out <- cbind(out, c(time1[3], time2[3], 
                      all(round(test1$sigma, 7)==round(test2$sigma, 7))))
  
}
res <- cbind(pseq, t(out)[-1, ])
plot(res[, 1], res[, 2], type="l", ylim=range(res[, 2:3]))
lines(res[, 1], res[, 3], col=2)

### testing Gibbs sampler with intercept
set.seed(123)
n <- 50
p <- 100
x <- matrix(rnorm(n*p), ncol=p, nrow=n)
lambda1 <- 1
lambda2 <- 1
beta <- seq(-0.1, 0.1, length.out=p)
m <- rep(1, n)
y <- rbinom(n, m, exp(x %*% beta)/(1 + exp(x %*% beta)))

sourceCpp(paste(path.code, "ENgibbs.cpp", sep=""))
sam <- gibbsC(x, y, m, n, p, lambda1, rep(lambda2, p), b0=c(0, beta), intercept=TRUE, 1)
sam <- replicate(1000, tryCatch({
  gibbsC(x, y, m, n, p, lambda1, rep(lambda2, p), b0=c(0, beta), intercept=TRUE, 1)$beta},
  error={function(war) {return(rep(NA, p + 1))}}))
n.na <- apply(sam, 1, function(r) {sum(is.na(r))})
beta.est <- apply(sam, 1, function(b) {d <- density(b, na.rm=TRUE); return(d$x[which.max(d$y)])})
test2 <- glm(y ~ x, family="binomial")
plot(c(0, beta), test2$coefficients)
plot(c(0, beta), beta.est)
plot(test2$coefficients, beta.est)

sourceCpp(paste(path.code, "ENgibbs.cpp", sep=""))
test1 <- rbetaC(x, y - 0.5, n, p, rchisq(p, 1), rchisq(n, 1), rep(lambda2, p), TRUE)
chol(test1)






