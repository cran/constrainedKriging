\name{CKrige}
\docType{methods}
\alias{CKrige}
\alias{CKrige.points}
\alias{CKrige.polygons}
\alias{CKrige-methods}
\alias{CKrige,formula,data.frame,formula,preCKrigePoints-method}
\alias{CKrige,formula,data.frame,formula,preCKrigePolygons-method}
\alias{print.CKrige.exout.polygons}
\alias{summary.CKrige.exout.polygons}
\alias{print.CKrige.exout.points}
\alias{summary.CKrige.exout.points}
\title{Spatial interpolation (2D) by constrained, covariance-matching 
constrained and universal, global Point or Block Kriging
}
\description{Function for constrained, covariance-matching 
constrained and universal (external drift kriging) point or block (of any shape) kriging in a global 
neighbourhood and for isotropic covariance models.}
%
\usage{
CKrige( formula, data, locations, object, \ldots)

\S4method{CKrige}{formula,data.frame,formula,preCKrigePolygons}(formula, data, locations, object, method = 2, ex.out = F)

\S4method{CKrige}{formula,data.frame,formula,preCKrigePoints}(formula, data, locations, object, method = 2, ex.out = F)
}
%
\arguments{
  \item{formula}{formula of the linear regression model in the form \code{response ~ terms of covariates}, for ordinary kriging use the 
	formula \code{response ~ 1}}
   \item{data}{a data frame with the values of the covariates, the names of the
   covariates used in the formula object must match the column
   names of \code{data}.}
  \item{locations}{a \code{formula} object that describe the coordinates 
	of the data locations (e.g \code{ ~ x+y })}
  \item{object}{either an object of the class \dQuote{\code{preCKrigePolygons}} for block kriging or of the class
	\dQuote{\code{preCKrigePoints}} for point kriging. In general the output object of the \code{\link{preCKrige}} function.}
\item{\ldots}{two further arguments to control the spatial 
interpolation method and the output}
  \item{method}{numeric value to choose the kriging method \code{1} 
	universal (external drift), \code{method = 2} constrained, 
	\code{method = 3} covariance-matching constrained kriging. By default, \code{method = 2}}
  \item{ex.out}{logical value, if \code{ex.out} is set \code{TRUE} \code{CKrige} returns an extended output with additional information,
  see details for more informations,  by default  \code{ex.out = F} }
}
\details{The \code{CKrige} function depends always on a \code{\link{preCKrige}}
output object that contains the parameter of the isotropic covariance model as well
as the covariates of the prediction targets.}
%
\value{By default, \code{CKrige} returns an object of
the class \code{SpatialPointsDataFrame} or \cr  \code{SpatialPolygonsDataFrame}
depending whether the input object for the \code{object} argument is of the
class  \dQuote{\code{preCKrigePoints}} or \dQuote{\code{preCKrigePolygons}}.\cr


The data frame of the returned object contains the following columns independent of the selected kriging method:
\item{prediction}{numeric vector with the kriging prediction of the chosen method}
\item{prediction.se}{numeric vector with the root mean square error (kriging standard error)}

The data frame contains 3 additional columns with constrained kriging
parameters, if the argument \code{method = 2} of the \code{CKrige} function:
\item{sqrt.P}{numeric vector with sqrt( Var[ target point or block ] -
Var[ fitted values ] )}
\item{sqrt.Q}{numeric vector with sqrt( Var[ universal kriging predictor ] -
Var[ fitted values ] )}
\item{K}{numeric vector with \code{sqrt.P / sqrt.Q}}

The data frame contains 3 additional columns with covariance-matching
constrained kriging parameters, if the argument
\code{method = 3} of the \code{CKrige} function:
\item{P1.11}{numeric vector, first element of the matrix P1 = ( Cov[target
point or block] - Cov[fitted values] )^(1/2)}
\item{Q1.11}{numeric vector, first element of the matrix Q1 = (
Cov[universal kriging predictor] - Cov[fitted values] )^(1/2)}
\item{K.11}{numeric vector, first element of the matrix K = O1^{-1}P1[1,1]}\cr\cr

The \code{CKrige} function returns a list with the following components if
the argument \code{ex.out = T} and the argument \code{method} is either \code{1} or \code{2}:
\item{object}{either an object of the class \code{SpatialPolygonsDataFrame}
or \cr \code{SpatialPointsDataFrame} as described above}
\item{krig.method}{numeric scalar, number of the chosen kriging method 1, 2
or 3.}
\item{parameter}{list with 2 components.  First component \code{beta.coef}
is the vector with the Generalized Least Square coefficients of the linear
regression and the second component \code{cov.beta} contains the covariance
matrix of the Generalized Least Square coefficients.}
\item{sk.weights}{if argument \code{method = 1} or \code{method = 2} \code{sk.weights} is a matrix with the simple kriging weights.  The ith column
contains the simple kriging weights of the ith prediction target object. If the argument \code{method = 3} the list component \code{sk.weights} is a list.
Each list component contains the matrix with the  simple kriging weights of the prediction target and its defined neighbours.}
\item{inv.Sigma}{matrix, inverse covariance matrix of the data}
\item{residuals}{numeric vector with the Generalized Least Square residuals
of the linear regression.}\cr

The list of the extended output contains the additional component \code{CMCK.par} if the argument \code{method = 3}. The
\code{CMCK.par} component is a list of lists with CMCK parameters, in particular  \code{P1} list of the P1 matrices, 
\code{Q1} list of the Q1 matrices and  \code{K} list of the K matrices.} 
%
\references{
See main help page of the \pkg{\link{constrainedKriging}} package.}
\author{Christoph Hofer, \email{christoph.hofer@alumni.ethz.ch}
}
%
\note{\code{print} and \code{summary} methods for the different \code{CKrige} output objects are available.}
\seealso{ \code{\link{preCKrige}} }
\examples{
# load data
data(meuse,meuse.blocks)

# approximation of block variance
# pixel area = 75m x 75m
# exponential covariance function with measurement error = 0, nugget = 0.05,
# part. sill =  0.15 and range parameter = 192.5
preCK <- preCKrige(newdata = meuse.blocks,
         model = covmodel("exponential", 0, 0.05, 0.15, scale = 192.5),
	 pwidth = 75, pheight = 75)
    
# block prediction by constrained kriging
ck <- CKrige(formula = log(zinc) ~ sqrt(dist),
    data = meuse, locations = ~ x + y, object = preCK, ex.out = TRUE)
    
# backtransformation to the original scale for the CK prediction
beta <- ck$parameter$beta.coef
c1 <-  0.2
c2  <- beta[2]^2 * meuse.blocks@data$M
c3  <- unlist( preCK@covmat )
ck$object@data$Zn <- exp(ck$object@data$prediction + 0.5*(c1 + c2 - c3))

# factor: the factor multiplied with the corresponding block mean
# prediction gives the upper bound of the 97.5% confidence interval
# of the prediction
ck$object@data$factor <- exp(ck$object@data$prediction
                         + 1.96*ck$object@data$prediction.se)/ck$object@data$Zn

# plots with spplot, generic function in the sp package

# plot of the constrained block kriging prediction
# function ck.colors(n) create a vector of n colors
ticks <- seq(0, 1850, by = 185) 
spplot(ck$object, zcol ="Zn" , col.regions = ck.colors(10), at = ticks,
    colorkey= list(labels = list(at = ticks, labels = ticks)))
    
# plot of the factor to get the upper bound of the 97.5% confidence interval
spplot(ck$object, zcol = "factor", col.regions = ck.colors(10), at = seq(1,3, by = 0.2))
}

\keyword{methods}