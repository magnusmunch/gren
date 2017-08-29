##############################  preamble  #############################
# simulations to check proportionality of multiplier estimates        #
# version: 01                                                         #
# author: Magnus M?nch                                                #
# created: 15-03-2017                                                 #
# last edited: 15-03-2017                                             #
#######################################################################

###############################  notes  ###############################
#                                                                     #
#######################################################################

### paths
path.code <- as.character(ifelse(Sys.info()[1]=="Darwin","/Users/magnusmunch/Documents/PhD/EBEN/code/" ,"~/EBEN/code/"))
path.res <- as.character(ifelse(Sys.info()[1]=="Darwin","/Users/magnusmunch/Documents/PhD/EBEN/results/" ,"~/EBEN/results/"))
path.data <- as.character(ifelse(Sys.info()[1]=="Darwin","/Users/magnusmunch/Documents/PhD/EBEN/data/hepB/" ,"~/EBEN/data/hepB/"))
path.graph <- "/Users/magnusmunch/Documents/PhD/EBEN/graphs/"

### libraries
library(Rcpp)
library(glmnet)
library(penalized)
library(GRridge)

### function to estimate variational Bayes parameters
VBen <- function(x, y, b0, lambda1, lambda2, lambdagvec, eps) {
  
  n <- nrow(x)
  p <- ncol(x)
  phi <- 0.25*lambda1^2/lambda2
  kappa <- y - 0.5
  
  # starting values
  pred0 <- as.numeric(exp(x %*% b0)/(1 + exp(x %*% b0)))
  w <- sqrt(pred0*(1 - pred0))
  xw <- x*w
  svdxw <- svd(xw)
  d <- svdxw$d
  v <- svdxw$v
  invmat <- d^2/(d^2 + 4*(lambda1 + lambda2))^2
  sigmaold <- t(t(v)*invmat) %*% t(v)
  muold <- as.numeric(sigmaold %*% (t(x) %*% as.matrix(kappa)))
  ci <- as.numeric(sqrt(colSums(t(x) * (sigmaold %*% t(x))) + (colSums(t(x)*muold))^2))
  chi <- as.numeric(0.5*(lambda1 + lambda2)*(diag(sigmaold) + muold^2))
  
  conv <- FALSE
  while(!conv) {
    
    # estimating new model parameters
    newparam <- est_param(x, kappa, m, n, p, ci, rep(phi, p), chi, lambdagvec*lambda2, 
                          intercept=FALSE)
    
    sigma <- newparam$sigma
    mu <- newparam$mu
    ci <- newparam$ci
    chi <- newparam$chi
    
    # checking convergence of inner loop
    conv <- max(c(abs((mu - muold)/muold), abs(diag((sigma - sigmaold)/sigmaold)))) < eps
    
    muold <- mu
    sigmaold <- sigma
  }
  
  out <- list(mu=mu, sigma=sigma, chi=chi, ci=ci, phi=phi)
  return(out)
}

### source grENVB functions
source(paste(path.code, "grVBEM.R", sep=""))

### compile MCMC sampler
sourceCpp(paste(path.code, "ENgibbs.cpp", sep=""))

### the simulation
p <- 500
n <- 100
G <- 2
groups <- rep(1:G, each=p/G)
beta <- c(rep(0.5, p/G), rep(1, p/G))
nreps <- 100
m <- rep(1, n)

mat.sum2 <- matrix(NA, nrow=nreps, ncol=2*G)
mat.ratio2 <- matrix(NA, nrow=nreps, ncol=2*G)
lambda1 <- 1
lambda2 <- 1
for(r in 1:nreps) {

  set.seed(r + 100)
  print(paste("Iteration: ", r, sep=""))
  x <- matrix(rnorm(n*p), ncol=p, nrow=n)
  y <- rbinom(n, 1, as.numeric(exp(x %*% as.matrix(beta))/(1 + exp(x %*% as.matrix(beta)))))

  fit.pen <- penalized(y, x, unpenalized=~0, model="logistic", lambda1=0, 
                       lambda2=2*(lambda1 + lambda2), trace=FALSE)
  b0 <- coef(fit.pen, which="all")
  
  lambdagvec <- c(rep(0.5, p/G), rep(2, p/G))
  fit.grVBen <- VBen(x, y, b0, lambda1, lambda2, lambdagvec, eps=0.001)
  fit.grGibbs <- gibbsC(x, y, m, n, p, lambda1, lambda2*lambdagvec, b0, intercept=FALSE, K=10000)
  
  # hist(fit.grGibbs$beta[1, -c(1:500)], freq=FALSE)
  # curve(dnorm(x, mean=fit.grVBen$mu[1], sd=fit.grVBen$sigma[1, 1]), add=TRUE)
  # 
  # plot(apply(fit.grGibbs$beta[, -c(1:500)], 1, var), diag(fit.grVBen$sigma))
  # lm(apply(fit.grGibbs$beta[, -c(1:500)], 1, var) ~ diag(fit.grVBen$sigma))
  # plot(apply(fit.grGibbs$beta[, -c(1:500)], 1, mean), as.numeric(fit.grVBen$mu))
  # lm(apply(fit.grGibbs$beta[, -c(1:500)], 1, mean) ~ as.numeric(fit.grVBen$mu))
  
  var.gibbs <- sapply(1:p, function(j) {var(fit.grGibbs$beta[j, -c(1:500)])})
  var.vb <- diag(fit.grVBen$sigma)
  
  sum.gibbs <- sapply(1:p, function(j) {
    mean(fit.grGibbs$beta[j, -c(1:500)]^2*fit.grGibbs$tau[j, -c(1:500)]/
           (fit.grGibbs$tau[j, -c(1:500)] - 1))})
  sum.vb <- (diag(fit.grVBen$sigma) + as.numeric(fit.grVBen$mu))*(1 + sqrt(as.numeric(fit.grVBen$phi/fit.grVBen$chi)))
  
  mat.sum2[r, ] <- c(sum(sum.gibbs[1:(p/G)]), sum(sum.gibbs[(p/G + 1):p]),
                     sum(sum.vb[1:(p/G)]), sum(sum.vb[(p/G + 1):p]))
  mat.ratio2[r, ] <- c(sum(var.gibbs[1:(p/G)]), sum(var.gibbs[(p/G + 1):p]),
                       sum(var.vb[1:(p/G)]), sum(var.vb[(p/G + 1):p]))
}

save(mat.ratio2, mat.sum2, file=paste(path.res, "grVBEM_sim4_res2.Rdata", sep=""))

cbind(mat.sum2[, 1]/mat.sum2[, 2], mat.sum2[, 3]/mat.sum2[, 4])
t.test(mat.ratio1[, 1]/mat.ratio1[, 2], mat.ratio1[, 3]/mat.ratio1[, 4])

# colnames(mat.mult1) <- c(paste(rep("grEBEN", 5), 1:5, sep=""),
#                         paste(rep("GRridge", 5), 1:5, sep=""))
# save(mat.mult1, file=paste(path.res, "grVBEM_sim3_res1.Rdata", sep=""))

# load(paste(path.res, "grVBEM_sim3_res1.Rdata", sep=""))


p <- 5
beta <- seq(0, 5)




mat <- matrix(0, ncol=3, nrow=3)
mat[1, 1] <- 3
mat[2, 1] <- mat[1, 2] <- mat[2, 2] <- 2
mat[3, 1] <- mat[3, 2] <- mat[3, 3] <- mat[2, 3] <- mat[1, 3] <- 1
solve(mat)

x <- matrix(rnorm(3), ncol=1)
t(x) %*% mat %*% x


