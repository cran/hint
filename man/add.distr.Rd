\name{add.distr}
\alias{add.distr}

\title{
Add Distributions to a Plot
}
\description{
This function will add one or more distributions or hypothesis tests to an existing plot.
}
\usage{
add.distr(..., cols = "blue", test.cols = "red")
}

\arguments{
  \item{\dots}{
One or more distributions or objects of class \code{hint.test}.
}
  \item{cols}{
A character string vector naming the colours of the distributions. If \code{length(cols)} is less than the number of distributions, the colours will be recycled. Defaults to \code{"blue"}.
}
  \item{test.cols}{
A character string vector naming the colours to use for the regions in which the cumulative probability of the hypothesis test was derived (if it exists). If \code{length(test.cols)} is less than the number of distributions, the colours will be recycled. Defaults to \code{"red"}.
}

}

\value{
Plots to the current device.
}

\author{
Alex T. Kalinka \email{alex.t.kalinka@gmail.com}
}


\seealso{
\code{\link{plotDistr}}, \code{\link{plot.hint.test}}.
}
\examples{
dd <- dhint(25, c(12, 9))
plotDistr(dd)
dd1 <- dhint(12, c(4, 7))
dd2 <- dhint(8, c(4, 2))
add.distr(dd1, dd2, cols = c("blue", "darkgreen"))
}

