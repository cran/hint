\name{print.hint.test}
\alias{print.hint.test}

\title{
Print a Summary of a hint.test Object
}
\description{
Prints a summary of a \code{hint.test} object to the screen.
}
\usage{
\method{print}{hint.test}(x, ...)
}

\arguments{
  \item{x}{
An object of class \code{hint.test}.
}
  \item{\dots}{
Additional arguments to be passed to \code{print}.
}
}

\value{
Prints summary output to the screen.
}

\author{
Alex T. Kalinka \email{alex.t.kalinka@gmail.com}
}


\seealso{
\code{\link{hint.test}}, \code{\link{hint.dist.test}}.
}
\examples{
## Build test object.
dd <- data.frame(letters[1:20], rep(1,10), rep(1,10))
tt <- hint.test(dd, letters[1:9], letters[7:15], alternative = "greater")

## Print summary:
print(tt)

}

