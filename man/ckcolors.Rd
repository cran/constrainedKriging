% 2023-12-15 A. Papritz correction of spelling errors and major revision of
% text, extended examples
% cd R.user.home/constrainedKriging/constrainedKriging_02-7/pkg/man/
% R CMD Rdconv -t html -o bla.html ckcolors.Rd ; open bla.html;
% R CMD Rd2pdf --force ckcolors.Rd;
% hunspell -d en_GB-ise,en_GB-ize,en_GB-large -t -i mac  ckcolors.Rd
\name{ck.colors}
\alias{ck.colors}

\title{Colour Palette}

\description{Create a vector of \var{n} contiguous colours based on
\code{\link[grDevices]{hsv}}.}
%
\usage{
ck.colors(n)
}

\arguments{
 \item{n}{a positive integer scalar with the numbers of colours in the 
 palette.}
}

\value{A character vector with \var{n} hex colour codes.}

\author{Christoph Hofer, \email{christoph.hofer@alumni.ethz.ch}}

\examples{
ck.colors(3L)
}

\keyword{methods}
