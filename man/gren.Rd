\name{gren}
\alias{gren}
\title{
Group-regularized logistic elastic net regression
}
\description{
Function that estimates a group-regularized elastic net model.
}
\usage{
gren(x, y, m, unpenalized=NULL, partitions=NULL, alpha=0.5, lambda=NULL, 
     intercept=TRUE, monotone=NULL, psel=TRUE, posterior=FALSE, nfolds=nrow(x), 
     foldid=NULL, trace=TRUE, 
     init=list(lambdag=NULL, mu=NULL, sigma=NULL, chi=NULL, 
               ci=NULL),
     control=list(epsilon=0.001, maxit=500, maxit.opt=1000, maxit.vb=100))
}

\arguments{
\item{x}{
feature data as either \code{matrix} or \code{data.frame}.}
\item{y}{
response as either a \code{numeric} with binomial/binary successes of length \code{nrow(x)} or a \code{matrix} of \code{nrow(x)} rows and two columns, where the first column contains the binomial/binary failures and the second column the binomial/binary successes.}
\item{m}{
\code{numeric} of length \code{nrow(x)} that contains the number of Bernoulli trials.}
\item{unpenalized}{
Optional \code{matrix} or \code{data.frame} of unpenalized covariates of \code{nrow(x)} rows.}
\item{partitions}{
\code{list} that contains the (possibly multiple) partitions of the data. Every \code{list} object corresponds to one partition, where every partition is a \code{numeric} of length \code{ncol(x)} containg the group ids of the features.}
\item{alpha}{
proportion of L1 penalty as a \code{numeric} of length 1.}
\item{lambda}{
global penalty parameter. The default \code{NULL} will result in estimation by cross-validation.}
\item{intercept}{
\code{logical} to indicate whether an intercept should be included.}
\item{monotone}{
\code{list} of two \code{logical} vectors of length \code{length(partitions)}. The first one \code{monotone} indicates whether the corresponding partition's penalty parameters should be monotonically estimates, the second vector \code{decreasing} indicates whether the penalty parameters are decreasing with group number.}
\item{psel}{
either a \code{numeric} vector that indicates the number of features to select or a \code{logical}. If \code{TRUE} feature selection is done by letting \code{\link{glmnet}} determine the penalty parameter sequence.}
\item{posterior}{
if \code{TRUE}, the full variational Bayes posterior is returned.}
\item{nfolds}{
\code{numeric} of length 1 with the number of folds used in the cross-validation of the global \code{lambda}. The default is \code{nrow(x)}.}
\item{foldid}{
optional \code{numeric} vector of length \code{nrow(x)} with the fold assignments of the observations.}
\item{trace}{
if \code{TRUE}, progress of the algorithm is printed.}
\item{init}{
optional \code{list} containing the starting values of the iterative algorithm.}
\item{control}{
a \code{list} of algorithm control parameters. See Details for more information.}
}

\details{
This is the main function of the package that estimates a group-regularized elastic net regression. The elastic net penalty's proportion of L1-norm penalisation is determined by \code{alpha}. \code{alpha} close to 0 implies more ridge-like penalty, while \code{alpha} close to 1 implies lasso-like penalty. The algorithm is a two-step procedure: first, a global lambda penalty is estimates by cross-validation. Next, the groupwise lambda multipliers are estimates by an EM algorithm. The EM algorithm consists of: i) an expectation step in which the expected marginal likelihood of the penalty multipliers is iteratively approximated by a variational Bayes EM algorithm and ii) a maximisation step in which the approximate expected marginal likelihood is maximised with respect to the penalty multipliers. After convergence of the algorithm an (optional) frequentist elastic net model is fit using the estimated penalty multipliers by setting \code{psel=TRUE} or by setting \code{psel} to a \code{numeric} vector. \\ 

The user may speed up the procedure by specifying initial values for the EM algorithm in \code{init}. \code{init} is a \code{list} that contains initial values for \deqn{\lambda_g} in a \code{list} of length \code{length(partitions)}, initial values for the \deqn{\mu_j}, \deqn{\chi_j}, and \deqn{c_i} in \code{numeric} vectors of lengths \code{ncol(x)}, \code{ncol(x)}, and \code{nrow(x)}, respectively. The initial values for the \deqn{\Sigma_{ij}} are given in a \code{matrix} with \code{ncol(x)} rows and columns. \\

\code{control}
}

\value{
Function returns a list of output.
}

\references{
TBA
}

\author{
Magnus M. Münch <m.munch@vumc.nl>
}

\seealso{
\code{\link{predict.gren}},\code{\link{coef.gren}},\code{\link{cv.gren}}
}

\examples{
## Create data
p <- 1000
n <- 100
set.seed(2018)
x <- matrix(rnorm(n*p), ncol=p, nrow=n)
beta <- c(rnorm(p/2, 0, 0.1), rnorm(p/2, 0, 1))
m <- rep(1, n)
y <- rbinom(n, m, as.numeric(1/(1 + exp(-x \%*\% as.matrix(beta)))))
partitions <- list(groups=rep(c(1, 2), each=p/2))

## estimate model
fit.gren <- gren(x, y, m, partitions=partitions)
}