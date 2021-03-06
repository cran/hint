\name{plotDistr}
\alias{plotDistr}

\title{
Plot a Distribution
}
\description{
Plot a distribution or visualise the result of a hypothesis test.
}
\usage{
plotDistr(distr, col = "black", test.col = "red", xlim = NULL, ylim = NULL, 
          xlab = "Intersection size (v)", ylab = "Probability", add = FALSE, ...)
}

\arguments{
  \item{distr}{
A data frame or matrix in which the first column gives random variable values, and the second gives probabilities. Can also be a vector (in which case random variables of \code{0:length(distr)} will be automatically assigned, or an object of class \code{hint.test}.
}
  \item{col}{
A character string naming the colour to use for the distribution. Defaults to \code{"black"}.
}
  \item{test.col}{
A character string naming the colour to use for the region in which the cumulative probability of the hypothesis test was derived (if it exists). Defaults to \code{"red"}.
}
  \item{xlim}{
A vector of two numbers giving the range for the x-axis. If NULL (default), then this is determined by the maximum and minimum values in \code{distr}.
}
  \item{ylim}{
A vector of two numbers giving the range for the y-axis. If NULL (default), then this is determined by the maximum and minimum values in \code{distr}.
}
  \item{xlab}{
A character string giving a label for the x-axis. Deafults to \code{"Intersection size (v)"}.
}
  \item{ylab}{
A character string giving a label for the y-axis. Deafults to \code{"Probability"}.
}
  \item{add}{
Logical. Whether the plot will be added to an existing plot or not. Defaults to FALSE.
}
  \item{\dots}{
Additional arguments to be passed to \code{plot}.
}
}
\details{
Visualising the results of a hypothesis test may often be of interest, but can be especially useful for pedagogical purposes.
}
\value{
Plots to the current device.
}

\author{
Alex T. Kalinka \email{alex.t.kalinka@gmail.com}
}


\seealso{
\code{\link{plot.hint.test}}, \code{\link{hint.test}}, \code{\link{hint.dist.test}}.
}
\examples{
dd <- dhint(25, c(12, 9))
plotDistr(dd)
}

\keyword{distribution}
\keyword{htest}

