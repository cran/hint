\name{Binomialintersection}
\alias{Binomialintersection}
\alias{dbint}
\alias{pbint}
\alias{qbint}
\alias{rbint}

\title{
The Binomial Intersection Distribution
}
\description{
Density, distribution function, quantile function and random generation for the binomial intersection distribution.
}
\usage{
dbint(n, A, range = NULL, log = FALSE)
pbint(n, A, vals, upper.tail = TRUE, log.p = FALSE)
qbint(p, n, A, upper.tail = TRUE, log.p = FALSE)
rbint(num = 5, n, A)
}

\arguments{
  \item{n}{
An integer specifying the number of categories in the urns.
}
  \item{A}{
A vector of integers specifying the numbers of balls drawn from each urn. The length of the vector equals the number of urns.
}
  \item{p}{
A probability between 0 and 1.
}
  \item{num}{
An integer specifying the number of random numbers to generate. Defaults to 5.
}
  \item{range, vals}{
A vector of integers specifying the intersection sizes for which probabilities (\code{dhint}) or cumulative probabilites (\code{phint}) should be computed (can be a single number). If \code{range} is NULL (default) then probabilities will be returned over the entire range of possible values.
}
  \item{log, log.p}{
Logical. If TRUE, probabilities p are given as log(p). Defaults to FALSE.
}
  \item{upper.tail}{
Logical. If TRUE, probabilities are P(X >= v), else P(X <= v). Defaults to TRUE.
}

}
\details{
The binomial intersection distribution is given by
\deqn{  P(X = v|N) = {b \choose v} \left(\prod_{i=1}^{N-1} p_{i}\right)^{v} \left(1 - \prod_{i=1}^{N-1} p_{i}\right)^{b-v} }{ P(X = v|N) = choose(b,v)*((prod_{i=1}^{N-1} p_i)^v) * (1 - prod_{i=1}^{N-1} p_i)^(b-v) }
where b gives the sample size which is smallest. This is an approximation for the hypergeometric intersection distribution when \eqn{n}{n} is large and \eqn{b}{b} is small relative to the samples taken from the \eqn{N-1}{N-1} other urns.
}
\value{
\code{dbint}, \code{pbint}, and \code{qbint} return a data frame with two columns: v, the intersection size, and p, the associated p-values. \code{rbint} returns an integer vector of random samples based on the binomial intersection distribution.
}
\references{
Kalinka, A.T. (2013). The probability of drawing intersections: extending the hypergeometric distribution. \href{http://arxiv.org/abs/1305.0717}{arXiv.1305.0717}
}
\author{
Alex T. Kalinka \email{alex.t.kalinka@gmail.com}
}

\seealso{
\code{\link{Hyperintersection}}, \code{\link{hint.test}}, \code{\link{hint.dist.test}}, \code{\link{plotDistr}}, \code{\link{Hyperdistinct}}.
}
\examples{
## Generate the distribution of intersections sizes:
dd <- dbint(20, c(10, 12, 11, 14))
## Restrict the range of intersections.
dd <- dbint(20, c(10, 12), range = 0:5)
## Generate cumulative probabilities.
pp <- pbint(29, c(15, 8), vals = 5)
pp <- pbint(29, c(15, 8), vals = 2, upper.tail = FALSE)
## Extract quantiles:
qq <- qbint(0.15, 23, c(12, 10))
## Generate random samples from Binomial intersection distributions.
rr <- rbint(num = 10, 18, c(9, 14))
}

\keyword{distribution}
\keyword{univar}


