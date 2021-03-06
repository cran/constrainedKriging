\name{covmodel}
\alias{covmodel}
\alias{covmodellist}
\alias{print.covmodel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Create isotropic covariance model
}
\description{Function to generate isotropic covariance models, or add an
isotropic covariance
model to an existing isotropic model}
\usage{
covmodel(modelname, mev, nugget,variance, scale,
parameter, add.covmodel)

\S3method{print}{covmodel}(x, \ldots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{modelname}{character vector, name of the covariance model, e.g. "exponential",
"spherical", "gauss". A call of \code{covmodel()} without a function argument
 displays a table with all available models and their parameters. Check
the \code{CovarianceFct} in the \pkg{RandomFields} package for detailed
information about the covariance functions.
}
\item{mev}{numeric value, variance of the measurement error}
\item{nugget}{numeric value, variance of microstructure white noise process (range smaller than the data support)
}
  \item{variance}{numeric value, partial sill of the variogram model
}
  \item{scale}{numeric value, scale parameter of the variogram model
}
  \item{parameter}{numeric vector of covariance parameters, missing for some model
like \code{nugget}, \code{spherical} or \code{gauss} or
}
\item{add.covmodel}{ object of the class \code{covmodel} that is
added to the covariance model defined by \code{modelname} (see
examples)}
\item{x}{ a covariance model generated by \code{covmodel}}
\item{\ldots}{further printing arguments}
}

\note{The names and parametrisation of the covariance model originate
from the \code{CovarianceFct} in the \pkg{RandomFields} package.  The
values of the arguments \code{mev}, \code{nugget}, \code{variance} and \code{scale}
are by default = 0.

Please, be aware that you only can generate spatial isotropic covariance
models, Time-Space models or so called (hypermodels) are not implemented.
}

\value{
an object of the class \code{covmodel} that define a covariance model.
}
\author{
Christoph Hofer \email{christoph.hofer@alumni.ethz.ch}
}

\examples{
\dontrun{
# table with all available covariance models and their
# parameters
covmodel()

# exponential model without a measurement error and without a nugget,
# partial sill = 10, scale  parameter = 15
covmodel(modelname = "exponential", variance = 10, scale = 15)

# exponential model with a measurement error ( mev = 0.5) and a
# nugget (nugget = 2.1), exponential partial  sill (variance = 10)
# and scale parameter = 15
covmodel(modelname = "exponential", mev  = 0.5, nugget = 2.1,
variance = 10, scale = 15)
}
}
%
\keyword{models}

