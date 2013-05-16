\name{overlay.sim}
\alias{overlay.sim}

\title{
Overlay Simulation Results on Distribution
}
\description{
This function adds simulated data to a distribution plot.
}
\usage{
overlay.sim(sim, breaks, col = "red", pch = 1, lwd = 1)
}

\arguments{
  \item{sim}{
A vector of integers (simulated data). See \code{\link{sim.hypint}} and \code{\link{sim.hydist}}.
}
  \item{breaks}{
A vector of integers specifying the x-axis range for which simulated data values should be plotted.
}
  \item{col}{
A character string naming the colour of the data points. Defaults to \code{"red"}.
}
  \item{pch}{
A number specifying the point type. Defaults to 1.
}
  \item{lwd}{
A number specifying the line width. Defaults to 1.
}
}

\value{
Simulated data will be added to an existing plot.
}
\references{
Kalinka, A.T. (2013). The probability of drawing intersections: extending the hypergeometric distribution. \href{http://arxiv.org/abs/1305.0717}{arXiv.1305.0717}
}
\author{
Alex T. Kalinka \email{alex.t.kalinka@gmail.com}
}


\seealso{
\code{\link{sim.hypint}}, \code{\link{sim.hydist}}.}
\examples{
## Generate exact distribution and plot:
dd <- dhint(108, c(54,76,89))
plot(dd[,1], dd[,2], col="blue")
## Simulate drawing from 3 urns (sims should be larger).
nn <- sim.hydist(108, c(54,76,89), sims = 10)
## Add data to existing distribution plot:
overlay.sim(nn, breaks = 0:54)
}

