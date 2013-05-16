\name{hint.test}
\alias{hint.test}

\title{
A Hypergeometric Intersection Test
}
\description{
Apply the hypergeometric intersection test to categorical data to test for enrichment or depletion of intersections between two samples.
}
\usage{
hint.test(cats, draw1, draw2, alternative = "greater")
}

\arguments{
  \item{cats}{
A data frame or matrix with 3 columns; the first gives the category identifier, and the second and third give the number of balls belonging to this category in the first and second urns respectively.
}
  \item{draw1}{
A vector of objects corresponding to the categories given in \code{cats} drawn from the first urn.
}
  \item{draw2}{
A vector of objects corresponding to the categories given in \code{cats} drawn from the second urn.
}
  \item{alternative}{
A characer string specifying the hypothesis to be tested. Can be one of \code{"greater"}, \code{"less"}, or \code{"two.sided"}.
}
}
\details{
The hypergeometric intersection distributions describe the distribution of intersection sizes when sampling without replacement from two separate urns in which reside balls belonging to the same n object categories (see \code{\link{Hyperintersection}}).
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
\code{\link{phint}}, \code{\link{hint.dist.test}}, \code{\link{plotDistr}}, \code{\link{plot.hint.test}}.
}
\examples{
## Build data input:
dd <- data.frame(letters[1:20], rep(1,20), rep(1,20))
hint.test(dd, letters[1:9], letters[4:14], alternative = "greater")

dd <- data.frame(letters[1:20], rep(1,20), c(rep(1,4),rep(2,16)))
hint.test(dd, letters[1:9], letters[4:14], alternative = "less")

}

\keyword{distribution}
\keyword{htest}

