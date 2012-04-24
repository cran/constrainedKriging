\name{plot.preCKrigePolygons}
\alias{plot.preCKrigePolygons}
\title{Plotting polygon configurations and their approximated areas}
%
\description{Plotting method for objects of the class
\code{preCKrige.polygons}. The plot show the block neighbourhood
configuration for one polygon (block) of a \code{preCKrige.polygons} object
as well as its area approximation by the pixels.
}
%
\usage{
\S3method{plot}{preCKrigePolygons}(x, index, \ldots)
}
%
\arguments{
\item{x}{object of the class \code{preCKrigePolygons}.  In general
the output object of a \code{\link{preCKrige}} function call.}
\item{index}{numeric value, list index of the desired polygon (block) of the polygon list
  \code{x@polygons}}
\item{\ldots}{further plotting parameters}
}
%
\author{Christoph Hofer \email{christoph.hofer@alumni.ethz.ch}}

\seealso{ \code{\link{preCKrige}} }
%
\examples{
\dontrun{
### load data
data(meuse,meuse.blocks)

### plot blocks
plot(meuse.blocks)

### compute the approximated block variance of each block in
### meuse.blocks without the definition of neighbours blocks (default)
preCK_1  <- preCKrige(newdata = meuse.blocks,
    model = covmodel("exponential", 0.05, 0.15, scale = 192.5),
    pwidth = 75, pheight = 75)

### plot block approximation of block 59
plot(preCK_1, 59)


### define neighbours
if(require(spdep))
{
neighbours <- poly2nb(meuse.blocks)
class(neighbours)
### neighbours should be an object of the class "list"
class(neighbours) <- "list"
### compute the approximated block variance-covariance matrices of each block in
### meuse.blocks without the defined block neighbours
preCK_2 <- preCKrige(newdata = meuse.blocks, neighbours = neighbours,
    model = covmodel("exponential", 0.05, 0.15, scale = 192.5),
    pwidth = 75, pheight = 75)

### plot block approximation of block 59 and its
### block neighbours
plot(preCK_1, 59)
}
if(!require(spdep))
{
cat("Please, install the package spdep to excute this example.\n")
}

}
}
\keyword{dplot}