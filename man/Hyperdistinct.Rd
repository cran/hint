\name{Hyperdistinct}
\alias{Hyperdistinct}
\alias{dhydist}
\alias{phydist}
\alias{qhydist}
\alias{rhydist}

\title{
Drawing Distinct Categories from a Single Urn
}
\description{
Density, distribution function, quantile function and random generation for the distribution of distinct categories drawn from a single urn in which there are duplicates in q of the categories.
}
\usage{
dhydist(n, a, q, range = NULL, log = FALSE)
phydist(n, a, q, vals, upper.tail = TRUE, log.p = FALSE)
qhydist(p, n, a, q, upper.tail = TRUE, log.p = FALSE)
rhydist(num = 5, n, a, q)
}

\arguments{
  \item{n}{
An integer specifying the number of categories in the urn.
}
  \item{a}{
An integer specifying the number of balls drawn from the urn.
}
  \item{q}{
An integer specifying the number of categories in the urn which have duplicate members.
}
  \item{p}{
A probability between 0 and 1.
}
  \item{num}{
An integer specifying the number of random numbers to generate. Defaults to 5.
}
  \item{range, vals}{
A vector of integers specifying the intersection sizes for which probabilities (\code{dhydist}) or cumulative probabilites (\code{phydist}) should be computed (can be a single number). If \code{range} is NULL (default) then probabilities will be returned over the entire range of possible values.
}
  \item{log, log.p}{
Logical. If TRUE, probabilities p are given as log(p). Defaults to FALSE.
}
  \item{upper.tail}{
Logical. If TRUE, probabilities are P(X >= c), else P(X <= c). Defaults to TRUE.
}

}
\details{
The distribution of the number of distinct categories drawn when sampling without replacement from a single urn containing duplicates in \eqn{q}{q} of its \eqn{n}{n} categories is given by
\deqn{P(X=c) = {q \choose a-c} \sum^{q}_{j=0} {q-a+c \choose j}{n-a+c-j \choose 2c-a-j}/{n+q \choose a}}{P(X = c) = choose(q,a-c)*sum_j choose(q-a+c,j)*choose(n-a+c-j,2c-a-j)/choose(n+q,a)}
When all of the \eqn{n}{n} categories contain duplicates, this can be expressed in a closed form:
\deqn{P(X=c) = {n \choose c}{c \choose a-c}2^{2c-a} /{2n \choose a}}{P(X = c) = choose(n,c)*choose(c,a-c)*2^(2c-a) / choose(2n,a)}


}
\value{
\code{dhydist}, \code{phydist}, and \code{qhydist} return a data frame with two columns: c, the number of distinct categories drawn, and p, the associated p-values. \code{rhydist} returns an integer vector of random samples based on the distribution of distinct categories when sampling from a single urn containing \eqn{q}{q} duplicates in \eqn{n}{n} categories.
}
\references{
Kalinka, A.T. (2013). The probability of drawing intersections: extending the hypergeometric distribution. \href{http://arxiv.org/abs/1305.0717}{arXiv.1305.0717}
}
\author{
Alex T. Kalinka \email{alex.t.kalinka@gmail.com}
}

\seealso{
\code{\link{Hyperintersection}}, \code{\link{plotDistr}}.
}
\examples{
## Generate the distribution of distinct categories drawn from a single urn.
dd <- dhydist(20, 10, 12)
## Restrict the range of intersections.
dd <- dhydist(20, 10, 12, range = 5:10)
## Generate cumulative probabilities.
pp <- phydist(29, 15, 8, vals = 5)
pp <- phydist(29, 15, 8, vals = 2, upper.tail = FALSE)
## Extract quantiles:
qq <- qhydist(0.15, 23, 12, 10)
## Generate random samples based on this distribution.
rr <- rhydist(num = 10, 18, 9, 12)
}

\keyword{distribution}
\keyword{univar}


