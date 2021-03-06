\name{sim.hypint}
\alias{sim.hypint}

\title{
Simulating Intersection Sizes
}
\description{
This function simulates the drawing of balls from two urns to generate distributions of intersection sizes.
}
\usage{
sim.hypint(n, A, sims = 10000, Na = NULL)
}

\arguments{
  \item{n}{
An integer specifying the number of categories in the urns.
}
  \item{A}{
A vector of integers specifying the samples to be taken from each urn. The length of the vector equals the number of urns.
}
  \item{sims}{
An integer specifying the number of replicate draws for the simulation. Defaults to 10,000.
}
  \item{Na}{
A list of integer vectors specifying the number of balls in each category for each urn.
}
}
\details{
The simulation allows any number of balls in each of the \code{n} categories via the \code{Na} argument.
}
\value{
An integer vector of length \code{sims} giving intersection sizes for all replicates of the simulation.
}
\references{
Kalinka, A.T. (2013). The probability of drawing intersections: extending the hypergeometric distribution. \href{http://arxiv.org/abs/1305.0717}{arXiv.1305.0717}
}
\author{
Alex T. Kalinka \email{alex.t.kalinka@gmail.com}
}


\seealso{
\code{\link{Hyperintersection}}, \code{\link{hint.test}}.
}
\examples{
nn <- sim.hypint(108, c(54, 73), sims = 10)
}

