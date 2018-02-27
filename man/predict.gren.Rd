\name{predict.gren}

\alias{predict.gren}
\alias{coef.gren}

\title{
Predictions from gren model.}

\description{
Create predictions from new data using a fitted gren model.}

\usage{
predict.gren(object, newx, unpenalized=NULL, s=NULL, type=c("VB", "freq"))}

\arguments{
  \item{object}{
A fitted \code{gren} model.}}

\details{
This is the predict function of the package.}

\value{
Function returns a \code{numeric} with predicted probabilities.}

\references{
TBA}

\author{
Magnus M. Münch <m.munch@vumc.nl>}

\seealso{
\code{\link{gren}},\code{\link{coef.gren}}}

\examples{
## Create data
p <- 1000
n <- 100
set.seed(2018)
x <- matrix(rnorm(n*p), ncol=p, nrow=n)
beta <- c(rnorm(p/2, 0, 0.1), rnorm(p/2, 0, 1))
m <- rep(1, n)
y <- rbinom(n, m, as.numeric(1/(1 + exp(-x \%*\% as.matrix(beta))))
partitions <- list(groups=rep(c(1, 2), each=p/2))

## estimate model
fit.gren <- gren(x, y, m, partitions=partitions)

## create new data
xnew <- matrix(rnorm(n*p), ncol=p, nrow=n)

## create predictions
preds <- predict.gren(fit.gren, xnew, type="VB")}