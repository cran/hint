\name{hint-package}
\alias{hint-package}
\alias{hint}
\docType{package}
\title{
The \code{hint} package
}
\description{
Tools for hypothesis testing based on Hypergeometric Intersection distributions.
}
\details{
Hypergeometric Intersection distributions are a broad group of distributions that describe the probability of picking intersections when drawing independently from two (or more) urns containing variable numbers of balls belonging to the same \code{n} categories.
}


\author{
Alex T. Kalinka \email{alex.t.kalinka@gmail.com}
}
\references{
Kalinka, A.T. (2013). The probability of drawing intersections: extending the hypergeometric distribution. \href{http://arxiv.org/abs/1305.0717}{arXiv.1305.0717}
}

\keyword{ package }
\keyword{distribution}
\keyword{htest}
\keyword{univar}

\seealso{
\code{\link{Hyperintersection}}, \code{\link{hint.test}}, \code{\link{Binomialintersection}}, \code{\link{Hyperdistinct}}.
}
\examples{
## Generate a hypergeometric intersection distribution and plot it.
dd <- dhint(110, c(65, 71))
plotDistr(dd)

## Allow duplicates in the second urn and add the distributions to the existing plot.
dd <- dhint(110, c(65, 71), 82)
dd2 <- dhint(59, c(44, 28), 47)
add.distr(dd, dd2, cols = c("blue", "darkgreen"))

## Test whether the intersection of categories in data sets is significant and plot the result:
dd <- data.frame(letters[1:20], rep(1,20), rep(1,20))
tt <- hint.test(dd, letters[1:9], letters[4:15], alternative = "greater")
plot(tt)

## Allow duplicates in the second urn.
dd <- data.frame(letters[1:20], rep(1,20), c(rep(1,4),rep(2,16)))
tt <- hint.test(dd, letters[1:9], letters[9:14], alternative = "less")
plot(tt)

}


