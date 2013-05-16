\name{Hyperintersection}
\alias{Hyperintersection}
\alias{dhint}
\alias{phint}
\alias{qhint}
\alias{rhint}

\title{
The Hypergeometric Intersection Family of Distributions
}
\description{
Density, distribution function, quantile function and random generation for the hypergeometric intersection distribution.
}
\usage{
dhint(n, a, b, q = 0, range = NULL, log = FALSE, verbose = TRUE)
phint(n, a, b, q = 0, vals, upper.tail = TRUE, log.p = FALSE)
qhint(p, n, a, b, q = 0, upper.tail = TRUE, log.p = FALSE)
rhint(num = 5, n, a, b, q = 0)
}

\arguments{
  \item{n}{
An integer specifying the number of categories in the urns.
}
  \item{a}{
An integer specifying the number of balls drawn from the first urn.
}
  \item{b}{
An integer specifying the number of balls drawn from the second urn.
}
  \item{q}{
An integer specifying the number of categories in the second urn which have duplicate members. If \code{q} is 0 (default) then the symmetrical, singleton case is computed, otherwise the asymmetrical, duplicates case is computed (see Details).
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
  \item{verbose}{
Logical. If TRUE, progress of calculation in the asymmetric, duplicates case is printed to the screen.
}
  \item{upper.tail}{
Logical. If TRUE, probabilities are P(X >= v), else P(X <= v). Defaults to TRUE.
}

}
\details{
The hypergeometric intersection distributions describe the distribution of intersection sizes when sampling without replacement from two separate urns in which reside balls belonging to the same n object categories. In the simplest case when there is exactly one ball in each category in each urn (symmetrical, singleton case), then the distribution is hypergeometric:
\deqn{P(X=v)=\frac{{a \choose v}{n-a \choose b-v}}{{n \choose b}}}{P(X=v) = (choose(a,v)*choose(n-a,b-v))/choose(n,b)}
If, however, we allow duplicates in \eqn{q \leq n}{q <= n} of the categories in the second urn, then the distribution of intersection sizes is described by the following variant of the hypergeometric:
\deqn{P(X=v) = \sum_{m=0}^{\alpha} \sum_{l=0}^{\beta} \sum_{j=0}^{l} {n-q \choose v-l} {q \choose l} {q-l \choose m} {n-v-q+l \choose a-v-m}  {l \choose j} {n+q-a-m-j \choose b-v} / {n \choose a}{n+q \choose b}}{P(X=v) = sum_m sum_l sum_j choose(n-q,v-l)*choose(q,l)*choose(q-l,m)*choose(n-v-q+l,a-v-m)*choose(l,j)*choose(n+q-a-m-j,b-v)/ choose(n,a)*choose(n+q,b)}

}
\value{
\code{dhint}, \code{phint}, and \code{qhint} return a data frame with two columns: v, the intersection size, and p, the associated p-values. \code{rhint} returns an integer vector of random samples based on the hypergeometric intersection distribution.
}
\references{
Kalinka, A.T. (2013). The probability of drawing intersections: extending the hypergeometric distribution. \href{http://arxiv.org/abs/1305.0717}{arXiv.1305.0717}
}
\author{
Alex T. Kalinka \email{alex.t.kalinka@gmail.com}
}

\seealso{
\code{\link{hint.test}}, \code{\link{hint.dist.test}}, \code{\link{plotDistr}}, \code{\link{Hyperdistinct}}.
}
\examples{
## Generate the distribution of intersections sizes without duplicates:
dd <- dhint(20, 10, 12)
## Restrict the range of intersections.
dd <- dhint(20, 10, 12, range = 0:5)
## Allow duplicates in q of the categories in the second urn:
dd <- dhint(35, 15, 11, 22, verbose = FALSE)
## Generate cumulative probabilities.
pp <- phint(29, 15, 8, vals = 5)
pp <- phint(29, 15, 8, vals = 2, upper.tail = FALSE)
pp <- phint(29, 15, 8, 23, vals = 2)
## Extract quantiles:
qq <- qhint(0.15, 23, 12, 10)
qq <- qhint(0.15, 23, 12, 10, 18)
## Generate random samples from Hypergeometric intersection distributions.
rr <- rhint(num = 10, 18, 9, 14)
rr <- rhint(num = 10, 22, 11, 17, 12)
}

\keyword{distribution}
\keyword{univar}


