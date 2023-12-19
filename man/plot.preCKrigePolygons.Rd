% 2023-12-15 A. Papritz correction of spelling errors and minor revision of
% text
% cd R.user.home/constrainedKriging/constrainedKriging_02-7/pkg/man/
% R CMD Rdconv -t html -o bla.html plot.preCKrigePolygons.Rd ; open bla.html;
% R CMD Rd2pdf --force plot.preCKrigePolygons.Rd;
% hunspell -d en_GB-ise,en_GB-ize,en_GB-large -t -i mac  plot.preCKrigePolygons.Rd
\name{plot.preCKrigePolygons}
\alias{plot.preCKrigePolygons}
\title{Plotting a Polygon Neighbourhood Configuration}
%
\description{Plotting method for objects of the class
\dQuote{\code{preCKrige.polygons}}.  The plot shows the polygon
neighbourhood configuration for one polygon (block) in a
\code{preCKrige.polygons} object as well as its representation by the
pixels.}
%
\usage{
\S3method{plot}{preCKrigePolygons}(x, index, \ldots)
}
%
\arguments{

\item{x}{an object of the class \dQuote{\code{preCKrigePolygons}}.  In
general the output object of a \code{\link{preCKrige}} function call.}

\item{index}{a numeric scalar with the index of the desired polygon
(block) in the list of polygons \code{x@polygons}.}

\item{\ldots}{further plotting parameters.}
}
%

\value{No return value, called for side effects.}

\references{
  Hofer, C. and Papritz, A. (2011).  constrainedKriging: an R-package
  for customary, constrained and covariance-matching constrained point
  or block Kriging.  \emph{Computers & Geosciences}. \bold{37}, 1562--1569,
  \doi{10.1016/j.cageo.2011.02.009}.
}

\author{Christoph Hofer, \email{christoph.hofer@alumni.ethz.ch}}

\seealso{ \code{\link{preCKrige}} and \code{\link{preCKrigePolygons}}.}
%
\examples{
\donttest{
### load data
data(meuse, package = "sp")
data(meuse.blocks)

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
if(!requireNamespace("spdep", quietly = TRUE)){
  stop("install package spdep to run example")
}
neighbours <- spdep::poly2nb(meuse.blocks)
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
plot(preCK_2, 59)

}
}
\keyword{dplot}
