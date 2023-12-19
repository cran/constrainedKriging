% 2023-12-15 A. Papritz correction of spelling errors and minor revision of
% text
% cd R.user.home/constrainedKriging/constrainedKriging_02-7/pkg/man/
% R CMD Rdconv -t html -o bla.html meuse.blocks.Rd ; open bla.html;
% R CMD Rd2pdf --force meuse.blocks.Rd;
% hunspell -d en_GB-ise,en_GB-ize,en_GB-large -t -i mac  meuse.blocks.Rd
\name{meuse.blocks}
\alias{meuse.blocks}
\title{Meuse block}

\description{ \code{meuse.blocks} is an object of the class
\dQuote{\code{SpatialPolygonsDataFrame}} that contains the coordinates
of 259 artificially defined (mostly square) blocks obtained by intersecting
a grid with 150 m mesh width over the the flood plain area of the river
Meuse, near the village Stein.  The 259 x 2 data frame contains the
covariate
\code{dist} and the attribute \code{M} (spatial variance of
\code{sqrt(dist)}) for each block.}

\format{
The object contains the following slots:
\describe{

\item{data}{a 259 x 2 data frame contains:

\describe{

\item{dist}{mean Euclidean distance of the blocks from the
river, normalized to the interval [0;1].}

\item{M}{the (spatial) variance of \code{sqrt(dist)}
 within the blocks, see \cite{Hofer & Papritz (2011)}.}
}
}

\item{polygons}{an object of the class \code{SpatialPolygons} that
contains the coordinates of the 259 blocks, see
\code{\link[sp]{SpatialPolygons}}.}


\item{plotOrder}{see \code{\link[sp]{SpatialPolygons}}.}


\item{bbox}{see \code{\link[sp]{SpatialPolygons}}.}


\item{proj4string}{see \code{\link[sp]{SpatialPolygons}}.}

}
}

\usage{
data(meuse.blocks)
}
\author{Christoph Hofer, \email{christoph.hofer@alumni.ethz.ch}}

\references{
  Hofer, C. and Papritz, A. (2011).  constrainedKriging: an R-package
  for customary, constrained and covariance-matching constrained point
  or block Kriging.  \emph{Computers & Geosciences}. \bold{37}, 1562--1569,
  \doi{10.1016/j.cageo.2011.02.009}.
}

\seealso{\link{preCKrige} and \link{CKrige}}
\keyword{datasets}
\examples{
data(meuse.blocks)
summary(meuse.blocks)
### show the shape of the 259 blocks
plot(meuse.blocks)
### plot maps of the covariate dist and attribute M
spplot(meuse.blocks)
}
