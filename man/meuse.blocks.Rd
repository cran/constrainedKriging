\name{meuse.blocks}
\alias{meuse.blocks}
\title{Meuse block}
\description{
\code{meuse.blocks} is an object of the class
\code{SpatialPolygonsDataFrame} that contains the coordinates of 259
artificially defined by putting a grid with 150 m mesh width over the the
flood plain area of the river Meuse, near the village Stein.  The 259 x 2
data frame contains the covariate \code{dist}
and the attribute \code{M} for each block.}
\format{
The object contains the following slots:
\describe{
\item{polygons}{an object of the class \code{SpatialPolygons} that contains the coordinates
  of the 259 blocks}
\item{data}{a 259 x 2 data frame contains:
\describe{
\item{dist}{mean Euclidean distance of the blocks from the river, normalized to the interval [0;1]}
\item{M}{the (spatial) variance of the mean distance between points, uniformly distributed within the blocks, and the river}
}} }
}
\usage{
data(meuse.blocks)
}
\author{Christoph Hofer \email{christoph.hofer@alumni.ethz.ch}}
\seealso{\link{preCKrige} and \link{CKrige}}
\keyword{datasets}
\examples{
\dontrun{
data(meuse.blocks)
summary(meuse.blocks)
### show the shape of the 259 blocks
plot(meuse.blocks)
### plot maps of the covariate \code{dist} and attribute \code{M}
spplot(meuse.blocks)
}
}
