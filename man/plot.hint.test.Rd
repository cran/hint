\name{plot.hint.test}
\alias{plot.hint.test}

\title{
Plot hint.test Objects.
}
\description{
This function visualises the results of a Hypergeometric Intersection test.
}
\usage{
\method{plot}{hint.test}(x, ...)
}

\arguments{
  \item{x}{
An object of class \code{hint.test}.
}
  \item{\dots}{
Additional arguments to be passed to \code{plot}.
}
}
\details{
Plots the relevant Hypergeometric Intersection distribution as a segment plot, and highlights the region where the observed statistic falls, i.e. the region from which the probability is computed (two.sided tests are visualised in one tail, the one with the smallest density). This can be especially useful for pedagogical purposes. 
}
\value{
Plots to the current device.
}

\author{
Alex T. Kalinka \email{alex.t.kalinka@gmail.com}
}

\seealso{
\code{\link{plotDistr}}, \code{\link{hint.test}}, \code{\link{hint.dist.test}} \code{\link{dhint}}.
}
\examples{
## Build hint.test object:
dd <- data.frame(letters[1:20], rep(1,10), rep(1,10))
tt <- hint.test(dd, letters[1:9], letters[7:15], alternative = "greater")

## Plot the distribution with the test region highlighted in red.
plot(tt, col = "blue")

}

\keyword{distribution}
\keyword{htest}

