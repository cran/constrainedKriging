% 2023-12-15 A. Papritz correction of spelling errors and major revision of
% text, extended examples
% cd R.user.home/constrainedKriging/constrainedKriging_02-7/pkg/man/
% R CMD Rdconv -t html -o bla.html internal.Rd ; open bla.html;
% R CMD Rd2pdf --force internal.Rd;
% hunspell -d en_GB-ise,en_GB-ize,en_GB-large -t -i mac  internal.Rd
\name{internal-functions}
\alias{f.point.block.cov}
\alias{f.point.block.cov}
\alias{K}

\title{Internal functions exported for use in dependent packages}

\description{The internal functions \code{K} and \code{f.point.block.cov}
are not meant to be called by users directly.  They are exported because
other R packages that import from \pkg{constrainedKriging} use them.}

\usage{
K(dist, model)

f.point.block.cov(pixconfig, locations, model)
}

\arguments{
\item{dist}{a numeric vector with distances.}

\item{model}{an object of class \dQuote{\code{covmodel}} that defines an
isotropic covariance model.}

\item{pixconfig}{a list of lists  with the information about the pixels used
for the covariance approximation of the polygon neighbourhood
configurations, see respective information in section \emph{Value} of
\code{\link{preCKrige}}.}

\item{locations}{a numeric matrix with two columns that contain the
coordinates of the observations.}
}

\value{
\code{K} computes a vector with autocorrelations.

\code{f.point.block.cov} computes a (\var{n} \eqn{times}{x} \var{m} )
matrix with the autocovariances between \var{n} observations and 
\var{m} prediction target blocks.
}

\author{Christoph Hofer, \email{christoph.hofer@alumni.ethz.ch}}

\keyword{internal}
