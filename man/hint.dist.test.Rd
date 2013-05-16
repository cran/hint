\name{hint.dist.test}
\alias{hint.dist.test}

\title{
To Test for Significant Distances between pairs of Intersection Sizes
}
\description{
Tests whether the absolute distance between two intersection sizes would be expected by chance, i.e. whether they fall into opposite tails of their respective Hypergeometric Intersection distributions.
}
\usage{
hint.dist.test(d, n1, a1, b1, n2, a2, b2, q1 = 0, q2 = 0, alternative = "greater")
}

\arguments{
  \item{d}{
A positive integer specifying the observed distance to be tested.
}
  \item{n1}{
An integer specifying the number of categories in the urns for the first distribution.
}
  \item{a1}{
An integer specifying the number of balls drawn from the first urn for the first distribution.
}
  \item{b1}{
An integer specifying the number of balls drawn from the second urn for the first distribution.
}
  \item{n2}{
An integer specifying the number of categories in the urns for the second distribution.
}
  \item{a2}{
An integer specifying the number of balls drawn from the first urn for the second distribution.
}
  \item{b2}{
An integer specifying the number of balls drawn from the second urn for the second distribution.
}
  \item{q1}{
An integer specifying the number of categories with duplicates in the second urn of the first distribution. If 0 then the symmetric, singleton case is computed, otherwise the asymmetric, duplicates case is computed (see \code{\link{Hyperintersection}}).
}
  \item{q2}{
An integer specifying the number of categories with duplicates in the second urn of the second distribution. If 0 then the symmetric, singleton case is computed, otherwise the asymmetric, duplicates case is computed (see \code{\link{Hyperintersection}}).
}
  \item{alternative}{
A characer string specifying the hypothesis to be tested. Can be one of \code{"greater"}, \code{"less"}, or \code{"two.sided"}.
}
}
\details{
The distribution of absolute distances between two hypergeometric intersection sizes is given by
\deqn{P(X=d) = \sum_{\{v_{1},v_{2}\}_{i} \in D_{d}}^{|D_{d}|} P(v_{1_i}|n_{1},a_{1},b_{1},...)\cdot P(v_{2_i}|n_{2},a_{2},b_{2},...) }{P(X=d) = sum_{v1,v2} P(v1)*P(v2) }
where \eqn{D_{d}}{D_d} is the set of pairs of intersection sizes, \eqn{\{v_{1},v_{2}\}}{(v_1,v_2)}, with absolute differences of size \eqn{d}{d}.
}
\value{
An object of class \code{hint.test}, which is a list containing the following components:
\item{parameters }{An integer vector giving the parameter values.}
\item{p.value }{A numerical value giving the p-value associated with the test.}
\item{alternative }{A character string naming the hypothesis that was tested.}
}
\references{
Kalinka, A.T. (2013). The probability of drawing intersections: extending the hypergeometric distribution. \href{http://arxiv.org/abs/1305.0717}{arXiv.1305.0717}
}
\author{
Alex T. Kalinka \email{alex.t.kalinka@gmail.com}
}


\seealso{
\code{\link{phint}}, \code{\link{hint.test}}, \code{\link{plotDistr}}, \code{\link{plot.hint.test}}.
}
\examples{
hint.dist.test(4, 20, 11, 9, 24, 14, 17, alternative = "greater")

}

\keyword{distribution}
\keyword{htest}


