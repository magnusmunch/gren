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
     control=list(epsilon=0.001, maxit=500, maxit.opt=1000, maxit.vb=100))
}

\arguments{
\item{x}{
data as either \code{matrix} or \code{data.frame}.}
}

\details{
This is the main function of the package.
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