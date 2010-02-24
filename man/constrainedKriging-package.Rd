\name{constrainedKriging-package}
\alias{constrainedKriging-package}
\alias{constrainedKriging}
\docType{package}
\title{constrainedKriging a package for nonlinear spatial predictions with local change of support}
\description{The constrainedKriging package provides functions for
tow-dimensional spatial interpolation by
constrained, covariance-matching constrained and universal (external drift)
kriging for points or block of any shape for data with a nonstationary mean
function and an isotropic weakly stationary variogram.  The linear spatial
interpolation methods, constrained and covariance-matching constrained
kriging, provide approximately unbiased prediction for nonlinear target
values under change of support.

In principle, the package provides two user functions, preCKrige
and CKrige to calculate spatial prediction in two
steps:

1. Call of \code{\link{preCKrige}} to calculate the variance-covariance
   matrices for defined sets of points or polygons (blocks).
   
2. Call of \code{\link{CKrige}} by using the output of
   preCKrige to calculate the spatial interpolation by one of the three kriging methods.
}
%
\details{
The constrained kriging predictor proposed by Cressie (1993) and the
covariance-matching constrained kriging predictor proposed by Aldworth and
Cressie (2003) are linear in the data like the universal kriging predictor.
However, the constrained kriging predictor satisfies in addition to the
unbiasedness constraint of universal kriging a second constraint that
matches the variances of the predictions to the variances of the prediction
target (either a point value or a block mean).
%
The covariance-matching constrained kriging predictor matches for a set of
blocks both the variances and covariances of predictions and prediction
targets ( points or block means) and is the extended version of the
constrained kriging predictor.
%
Like constrained kriging, covariance-matching constrained kriging is less
biased than universal kriging for nonlinear predictions and exactly
unbiased if the spatial variable is Gaussian.

The formulas of the three kriging predictors are given below:
%%% formula kriging predictors
\itemize{
\item{universal kriging predictor for 1 block: \deqn{\widehat{Y}_{\mathrm{UK}}(B) =
\mathbf{x}(B)^{\prime}\widehat{\beta}_{\mathrm{GLS}} +
\mathbf{ c }^{\prime} \boldsymbol{\Sigma} ( \mathbf{ Z } -
\mathbf{X}\widehat{\beta}_{\mathrm{GLS}} ) ,}{Y_UK(B) =
x(B)'beta.hat_(GLS) + c'Sigma( Z - Xbeta.hat_(GLS) ),}}
\item{constrained kriging predictor for 1 block: \deqn{\widehat{Y}_{\mathrm{CK}}(B) =
\mathbf{x}(B)^{\prime}\widehat{\beta}_{\mathrm{GLS}} + K
\mathbf{ c }^{\prime} \boldsymbol{\Sigma} ( \mathbf{ Z } -
\mathbf{X}\widehat{\beta}_{\mathrm{GLS}}) ,}{Y_CK(B) = x(B)'beta.hat_(GLS)
+ K_(CK)c'Sigma( Z -
Xbeta.hat_(GLS) ),}}
\item{covariance-matching constrained kriging predictor for m blocks: \deqn{\widehat{
\mathbf{ Y } }_{\mathrm{CMCK}} =
\mathbf{X}_{m}\widehat{\beta}_{\mathrm{GLS}} + \mathbf{ K }^{\prime}
\mathbf{C}^{\prime}\boldsymbol{\Sigma} (\mathbf{ Z } -
\mathbf{X}\widehat{\beta}_{\mathrm{GLS}}),}{Y_CMCK =
X_m'beta.hat_(GLS) +
K_(CMCK)'C'Sigma( Z - X'beta.hat_(GLS) ),}}} where \eqn{\mathbf{Z} =
c(Z(s_{1}), \ldots,Z(s_{n}))^{^\prime}}{Z = (Z(s_1), ..., Z(s_n) )'} is the
vector with the data; \eqn{\mathbf{s} = (x,y)^{\prime}}{s = (x,y)'}
indicates a location in the survey domain and
\eqn{Y(B)}{Y(B)} is the block mean value of the block area \eqn{B}{B};
\eqn{\mathbf{X} = (\mathbf{x}_{1}(B), \ldots,
\mathbf{x}_{n}(B))^{\prime}}{X =
(x(B_1), ..., x(B_m))'} is the design matrix of the data and
\eqn{\mathbf{X}_{m}}{X_m} is the design matrix of the target blocks;
\eqn{\widehat{\boldsymbol{\beta}}_{\mathrm{GLS}}}{beta.hat_(GLS)} is the
vector with the generalised least square estimate of the linear regression
coefficients;
\eqn{\mathbf{C} = (\mathbf{c}_{1},
\ldots, \mathbf{c}_{m})}{C = (c_1, ..., c_m)} is a (n,m)-matrix that
contains the covariances between the \eqn{m}{m} prediction targets (point
or blocks) and the
\eqn{n}{n} data points;
\eqn{\boldsymbol{\Sigma}}{Sigma} is the covariance matrix of the data;
the scalar \deqn{K = ( \mathrm{Var}[ Y(B)] - \mathrm{Var}[
\mathbf{x}(B)\widehat{\boldsymbol{\beta}}_{\mathrm{GLS}}] )^{\frac{1}{2}} /
( \mathrm{Var}[ \widehat{Y}_{UK}(B)  ]- \mathbf{ Var }[
\mathbf{x}(B)\widehat{\boldsymbol{\beta}}_{\mathrm{GLS}}] )^{\frac{1}{2}} =
(P/Q)^{\frac{1}{2}}}{ K_(CK) = ( Var[ Y(B) ] - Var[ x(B)'beta.hat_(GLS) ] )^0.5
/ ( Var[ Y(B)_UK ] - Var[ x(B)'beta.hat_(GLS) ] )^0.5 = P^0.5 / Q^0.5 }
% K CMCK
and the (m, m)-matrix \deqn{\mathbf{K} = \mathbf{Q}_{1}^{-1}\mathbf{P_{1}},}{K_(CMCK) =
Q1^(-1)P1,} where the (m,m)-matrix \deqn{\mathbf{P}_{1} = ( \mathrm{Cov}[\mathbf{Y},
\mathbf{Y}^{\prime}] -
\mathrm{Cov}[\mathbf{X}_{m}\widehat{\boldsymbol{\beta}}_{\mathrm{GLS}},
(\mathbf{X}_{m}\widehat{\boldsymbol{\beta}}_{\mathrm{GLS}})^{\prime}] )^{\frac{1}{2}}
}{P1 = Cov[ Y, Y' ] - Cov[ X_mbeta.hat_(GLS), X_mbeta.hat_(GLS)' ]}
and the (m, m)-matrix
\deqn{\mathbf{Q}_{1} = (\mathrm{Cov}[\mathbf{Y}_{\mathrm{UK}},
\mathbf{Y}_{\mathrm{UK}}^{\prime}] -
\mathrm{Cov}[\mathbf{X}_{m}\widehat{\boldsymbol{\beta}}_{\mathrm{GLS}},
(\mathbf{X}_{m}\widehat{\boldsymbol{\beta}}_{\mathrm{GLS}})^{\prime}])^{\frac{1}{2}}}{Q1
= Cov[ Y_(UK), Y_(UK)' ] - Cov[ X_mbeta.hat_(GLS), (X_mbeta.hat_(GLS))' ]}


The mean square prediction error (MSEP) of the three predictors are:
\itemize{
\item{MSPE of the universal kriging predictor for 1 block:}
\deqn{ MSPE[ \widehat{Y}_{UK}(B) ] = \mathrm{ Cov }[ (
\widehat{Y}_{UK}(B) - Y(B))( \widehat{Y}_{UK}(B) - Y(B) )^{\prime} ]}{MSPE[ Y(B)_(UK) ] = Var[ ( Y(B)_(UK) - Y(B) ), (
Y(B)_(UK) - Y(B) )' ],}
\item{MSPE of the constrained kriging predictor for 1 block:}
\deqn{ MSPE[ \widehat{\mathbf{Y}}_{CK} ] = MSPE[ \widehat{\mathbf{Y}}_{UK}
] + (P^{\frac{1}{2}} + Q^{\frac{1}{2}})^{2},}{MSPE[ Y(B)_(CK) ] = MSPE[
Y(B)_(UK) ] +( P^0.5 + Q^0.5 )^2,}
\item{MSPE of the covariance-matching constrained kriging predictor for m blocks:}
\deqn{ MSPE[ \widehat{\mathbf{Y}}_{CMCK} ] = MSPE[
\widehat{\mathbf{Y}}_{UK} ] +
(\mathbf{P}_{1}-\mathbf{Q}_{1})(\mathbf{P}_{1}-\mathbf{Q}_{1}).}{MSPE[
Y_(CMCK) ] = MSPE[ Y_(UK) ] +( P1 + Q1 ) (P1 + Q1 ).}
}
}% end details

\author{
Christoph Hofer \email{christoph.hofer@allumni.ethz.ch}
}
\references{%
Aldworth, J. and Cressie, N. (2003). Prediction of nonlinear spatial 
functionals. \emph{Journal of Statistical Planning and Inference}, \bold{112}, 3--41

Cressie, N. (1993). Aggregation in geostatistical problems. In 
A. Soares, editor, \emph{Geostatistics Troia 92}, \bold{1}, pages 25--36, 
Dordrecht. Kluwer Academic Publishers.

Hofer, C. and Papritz, A. (2010). Predicting threshold exceedance by local block means
in soil pollution surveys. \emph{Mathematical Geosciences}.  \bold{42}, 631-656, doi: 10.1007/s11004-010-9287-4

Hofer, C. and Papritz, A. (in preparation). constrainedKriging: an R-package for customary, constrained and
covariance-matching constrained point or block kriging. \emph{Computers & Geosciences}.
}
%
%
\keyword{ package }