\name{sim.hydist}
\alias{sim.hydist}

\title{
Simulating the Drawing of Distinct Categories
}
\description{
This function simulates the drawing of balls from a single urn containing q duplicates in n categories to generate distributions of the number of distinct categories drawn.
}
\usage{
sim.hydist(n, a, sims = 10000, Na = rep(2, n))
}

\arguments{
  \item{n}{
An integer specifying the number of categories in the urn.
}
  \item{a}{
An integer specifying the number of balls drawn from the urn.
}
  \item{sims}{
An integer specifying the number of replicate draws for the simulation. Defaults to 10,000.
}
  \item{Na}{
A vector of integers specifying the number of balls in each category of the urn.
}
}
\details{
The simulation allows any number of balls in each of the \code{n} categories via the \code{Na} argument.
}
\value{
An integer vector of length \code{sims} giving the numbers of distinct categories for all replicates of the simulation.
}
\references{
Kalinka, A.T. (2013). The probability of drawing intersections: extending the hypergeometric distribution. \href{http://arxiv.org/abs/1305.0717}{arXiv.1305.0717}
}
\author{
Alex T. Kalinka \email{alex.t.kalinka@gmail.com}
}


\seealso{
\code{\link{Hyperdistinct}}.
}
\examples{
## All n categories contain duplicates:
nn <- sim.hydist(108, 74, sims = 10)
## Duplicates in 78 of the categories:
nn <- sim.hydist(108, 74, sims = 10, Na = c(rep(1,(108-78)),rep(2,78)))
}

