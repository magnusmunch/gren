# rfiles <- c("coef.gren", "denet", "fopt_groups", "gren", "cv.gren", 
#             "predict.gren", "qtau", "renet", "rtau")
# cfiles <- c("est_param")
# library(Rcpp)
# Rcpp.package.skeleton(name="gren", author="Magnus M. Münch", 
#                       maintainer="Magnus M. Münch", email="m.munch@vumc.nl", 
#                       example_code=FALSE, module=FALSE,
#                       path="/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/rpackages",
#                       code_files=paste("/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/code/", 
#                                        rfiles, ".R", sep=""),
#                       cpp_files=paste("/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/code/", 
#                                       cfiles, ".cpp", sep=""))

# R CMD build "/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/rpackages/gren

# local
install.packages("/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/rpackages/gren_0.0.0.9000.tar.gz",
                 repos=NULL)
library(gren)
vignette("gren")

# # github
# library(devtools)
# install_github("magnusmunch/gren")
# library(gren)
# vignette("gren")

# test package
## Create data
p <- 1000
n <- 100
set.seed(2018)
x <- matrix(rnorm(n*p), ncol=p, nrow=n)
beta <- c(rnorm(p/2, 0, 0.1), rnorm(p/2, 0, 1))
m <- rep(1, n)
y <- rbinom(n, m, as.numeric(1/(1 + exp(-x %*% as.matrix(beta)))))
partitions <- list(groups=rep(c(1, 2), each=p/2))
fit.gren <- gren(x, y, m, partitions=partitions, lambda=NULL)
pred.gren <- predict(fit.gren, x)
coef.gren <- coef(fit.gren)
str(fit.gren)
summary(fit.gren$freq)
str(fit.gren$freq.model)
plot(fit.gren$freq.model$lambda, fit.gren$freq.model$dev.ratio, type="l")

# fit.cv.gren <- cv.gren(x, y, m, partitions=partitions, fix.lambda=TRUE, intercept=TRUE)
xl <- matrix(rnorm(1000), ncol=100, nrow=10)
fit.lm <- lm(y ~ x)
plot(fit.lm)




# test glmnet and myglmnet
# install.packages("glmnet")
install.packages("/Users/magnusmunch/Documents/OneDrive/PhD/EBEN/rpackages/myglmnet_2.0-13.tar.gz",
                 repos=NULL)
library(myglmnet)
library(glmnet)

p <- 1000
n <- 100
set.seed(2018)
x <- scale(matrix(rnorm(n*p), ncol=p, nrow=n))
beta <- c(rnorm(p/2, 0, 0.1), rnorm(p/2, 0, 1))
m <- rep(1, n)
y <- rbinom(n, m, as.numeric(1/(1 + exp(-x %*% as.matrix(beta)))))

alpha <- 0.5
pen.fac <- exp(seq(-1, 1, length.out=p))*p/sum(exp(seq(-1, 1, length.out=p)))

fit.glmnet1 <- glmnet::glmnet(x, y, family="binomial", alpha=alpha, 
                              standardize=FALSE, intercept=TRUE, 
                              penalty.factor=pen.fac)
fit.glmnet2 <- glmnet::glmnet(x, y, family="binomial", alpha=alpha, 
                              standardize=FALSE, intercept=TRUE)
fit.myglmnet1 <- myglmnet::glmnet(x, y, family="binomial", alpha=alpha, 
                                  standardize=FALSE, intercept=TRUE,
                                  penalty.factor=pen.fac)

lam.glmnet1.max <- max(abs(t(x) %*% (y - mean(y)*(1 - mean(y))))/
                         (n*alpha*pen.fac.rescale))
lam.glmnet2.max <- max(as.numeric(abs(t(x) %*% (y - mean(y)*(1 - mean(y)))))/
                         (n*alpha))
lam.myglmnet1.max <- max(abs(t(x) %*% (y - mean(y)*(1 - mean(y))))/
                          (n*alpha*sqrt(pen.fac.rescale)))

c(max(fit.glmnet1$lambda), lam.glmnet1.max)
c(max(fit.glmnet2$lambda), lam.glmnet2.max)
c(max(fit.myglmnet1$lambda), lam.myglmnet1.max)


fit.glmnet3 <- glmnet::glmnet(x, y, family="binomial", alpha=1, 
                              standardize=FALSE, intercept=TRUE,
                              penalty.factor=pen.fac)
fit.myglmnet2 <- myglmnet::glmnet(x, y, family="binomial", alpha=1, 
                                  standardize=FALSE, intercept=TRUE,
                                  penalty.factor=pen.fac)

identical(fit.glmnet3[-12], fit.myglmnet2[-12])
identical(myglmnet::glmnet, glmnet::glmnet)

