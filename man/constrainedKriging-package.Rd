% 2023-12-09 A. Papritz correction of spelling errors and minor revision of
% text
% cd R.user.home/constrainedKriging/constrainedKriging_02-7/pkg/man/
% R CMD Rdconv -t html -o bla.html constrainedKriging-package.Rd ; open bla.html;
% R CMD Rd2pdf --force constrainedKriging-package.Rd;
% hunspell -d en_GB-ise,en_GB-ize,en_GB-large -t -i mac  constrainedKriging-package.Rd
\name{constrainedKriging-package}
\alias{constrainedKriging-package}
\alias{constrainedKriging}
\docType{package}

\title{Nonlinear Spatial Kriging Predictions under Change of Support}

\description{The package \pkg{constrainedKriging} provides functions for
spatial interpolation by \emph{constrained},\cr
\emph{covariance-matching constrained} and \emph{universal}
(\emph{external drift}) \emph{Kriging} for points or blocks of any shape in
a two-dimensional domain from data with a non-stationary mean function and
an isotropic weakly stationary variogram.  The linear spatial interpolation
methods constrained and covariance-matching constrained Kriging provide
approximately unbiased predictions for non-linearly transformed target
values under change of support.

The package provides two main user functions, \code{\link{preCKrige}} and
\code{\link{CKrige}}, to calculate spatial predictions in two steps:

   \enumerate{ \item \code{\link{preCKrige}} computes the
   variance-covariance matrices for sets of prediction points or
   prediction blocks (polygons).

   \item \code{\link{CKrige}} computes from the output of
   \code{\link{preCKrige}} spatial predictions by one of three Kriging
   methods mentioned above.
}}

%
\details{ The constrained Kriging predictor introduced by
\cite{Cressie (1993)} and the covariance-matching constrained Kriging
predictor proposed by \cite{Aldworth and Cressie (2003)} are linear in the
data like the universal Kriging predictor.  However, the constrained
Kriging predictor satisfies in addition to the unbiasedness constraint of
universal Kriging a second constraint that matches the variances of the
predictions to the variances of the prediction targets (either point values
or block means).  The covariance-matching constrained Kriging predictor
matches for a set of points or blocks both the variances and covariances of
predictions and prediction targets and is the extended version of the
constrained Kriging predictor.  Like constrained Kriging,
covariance-matching constrained Kriging is less biased than universal
Kriging for predictions of non-linearly transformed functionals of a spatial
variable and exactly unbiased if the variable is Gaussian.  We summarize
the formulae of the three Kriging methods below for predicting block means
from point observations, but analogous equations could be given for
problems that do not involve change of support.

For a set of \eqn{m} blocks, \eqn{B_1, \ldots, B_m}, the covariance-matching
constrained Kriging predictor is given by
\deqn{\widehat{\mathbf{Y}}_{\mathrm{CMCK}} =
  \mathbf{X}_{m}\widehat{\boldsymbol{\beta}}_{\mathrm{GLS}} +
  \mathbf{K}^{\prime}
  \mathbf{C}^{\prime}\boldsymbol{\Sigma}^{-1} (\mathbf{Z} -
  \mathbf{X}\widehat{\boldsymbol{\beta}}_{\mathrm{GLS}}),
}{
  Yhat_CMCK = X_m beta.hat_GLS + K' C' Sigma^{-1} (Z - X beta.hat_GLS),
}
where
\eqn{
  \mathbf{Y} = (Y(B_1), \ldots, Y(B_m))^\prime
}{Y = (Y(B_1), \ldots, Y(B_m))'} is the set of block means to be predicted,
with \eqn{Y(B_i)} the mean value of
\eqn{Y} averaged over the block \eqn{B_i};
\eqn{
  \mathbf{Z} = (Z(\mathbf{s}_{1}), \ldots, Z(\mathbf{s}_{n}))^{^\prime}
}{Z = (Z(s_1), \ldots, Z(s_n))'} is the vector with data
\eqn{
  Z(\mathbf{s}_{i}) = Y(\mathbf{s}_{i}) + \epsilon_i
}{Z(s_i) = Y(s_i) + \epsilon_i} where the response
\eqn{Y(\mathbf{s}_{i})}{Y(s_i)} is possibly contaminated by measurement error
\eqn{\epsilon_i};
\eqn{\mathbf{s} = (x,y)^{\prime}}{s = (x,y)'} indicates a location in
the survey domain;
\eqn{
  \mathbf{X} = (\mathbf{x}(\mathbf{s}_{1}), \ldots,
  \mathbf{x}(\mathbf{s}_{n}))^{\prime}
}{X = (x(s_1), \ldots, x(s_n))'}
is the design matrix of the data and
\eqn{
  \mathbf{X}_m = (\mathbf{x}(B_{1}), \ldots, \mathbf{x}(B_{m}))^{\prime}
}{X_m = (x(B_1), \ldots, x(B_m))'}  the design matrix of the target blocks;
\eqn{\widehat{\boldsymbol{\boldsymbol{\beta}}}_{\mathrm{GLS}}}{beta.hat_GLS}
is the vector with the generalised least square estimate of the
linear regression coefficients;
\eqn{
  \mathbf{C} = (\mathbf{c}_{1}, \ldots, \mathbf{c}_{m})
}{C = (c_1, \ldots, c_m)} is a \eqn{(n\times m)}-matrix that
contains the covariances between the \eqn{n} data points and the \eqn{m}
prediction targets;
\eqn{\boldsymbol{\Sigma}}{\Sigma} is the \eqn{(n\times n)}-covariance
matrix of the data;
and \eqn{\mathbf{K}}{K} is the \eqn{(m \times m)}-matrix
\deqn{
  \mathbf{K} = \mathbf{Q}_{1}^{-1}\mathbf{P_{1}},
}{K = Q_1^{-1} P_1,}
where \eqn{\mathbf{P}_{1}}{P_1} is the \eqn{(m \times m)}-matrix
\deqn{
  \mathbf{P}_{1} = ( \mathrm{Cov}[\mathbf{Y}, \mathbf{Y}^{\prime}] -
    \mathrm{Cov}[\mathbf{X}_{m}\widehat{\boldsymbol{\boldsymbol{\beta}}}_{\mathrm{GLS}},
   (\mathbf{X}_{m}\widehat{\boldsymbol{\boldsymbol{\beta}}}_{\mathrm{GLS}})^{\prime}] )^{\frac{1}{2}}
}{P1 = Cov[Y, Y'] - Cov[X_m beta.hat_GLS, X_m beta.hat_GLS' ]}
and \eqn{\mathbf{Q}_{1}}{Q_1} the \eqn{(m \times m)}-matrix
\deqn{
  \mathbf{Q}_{1} = (\mathrm{Cov}[\widehat{\mathbf{Y}}_{\mathrm{UK}},
  \widehat{\mathbf{Y}}_{\mathrm{UK}}^{\prime}] -
  \mathrm{Cov}[\mathbf{X}_{m}\widehat{\boldsymbol{\boldsymbol{\beta}}}_{\mathrm{GLS}},
  (\mathbf{X}_{m}\widehat{\boldsymbol{\boldsymbol{\beta}}}_{\mathrm{GLS}})^{\prime}])^{\frac{1}{2}}
}{Q1 = Cov[Y_(UK), Y_(UK)'] - Cov[X_m beta.hat_GLS, (X_m beta.hat_GLS)']}
with
\deqn{\widehat{\mathbf{Y}}_{\mathrm{UK}} =
  \mathbf{X}_{m}\widehat{\boldsymbol{\beta}}_{\mathrm{GLS}} +
  \mathbf{C}^{\prime}\boldsymbol{\Sigma}^{-1} (\mathbf{Z} -
  \mathbf{X}\widehat{\boldsymbol{\beta}}_{\mathrm{GLS}}),
}{
  Yhat_UK = X_m beta.hat_GLS +  C' Sigma^{-1} (Z - X beta.hat_GLS),
}
denoting the universal Kriging predictor of \eqn{\mathbf{Y}}.

For \eqn{m = 1} \eqn{\widehat{\mathbf{Y}}_{\mathrm{CMCK}}}{Yhat_CMCK} reduces
to the constrained Kriging predictor
\deqn{\widehat{Y}_{\mathrm{CK}}(B_1) =
  \mathbf{x}(B_1)^{\prime}\widehat{\boldsymbol{\beta}}_{\mathrm{GLS}} + K
  \mathbf{c_1}^{\prime} \boldsymbol{\Sigma}^{-1} ( \mathbf{Z} -
  \mathbf{X}\widehat{\boldsymbol{\beta}}_{\mathrm{GLS}}) ,
}{
  Yhat_CK(B_1) = x(B_1)' beta.hat_GLS + K c_1' Sigma^{-1} (Z - X beta.hat_GLS),
}
with the scalar
\deqn{
  K = (\mathrm{Var}[Y(B_1)] -
  \mathrm{Var}[\mathbf{x}(B_1)^{\prime}\widehat{\boldsymbol{\boldsymbol{\beta}}}_{\mathrm{GLS}}] )^{\frac{1}{2}} /
  (\mathrm{Var}[\widehat{Y}_{\mathrm{UK}}(B_1)] - \mathrm{Var}[
  \mathbf{x}(B_1)^{\prime}\widehat{\boldsymbol{\boldsymbol{\beta}}}_{\mathrm{GLS}}] )^{\frac{1}{2}} =
  (P/Q)^{\frac{1}{2}}.
}{K = (
  Var[Y(B_1)] - Var[x(B_1)' beta.hat_GLS])^0.5
  / (Var[Y(B_1)_UK] - Var[x(B_1)' beta.hat_GLS])^0.5 = P^0.5 / Q^0.5.
}

The mean square prediction error (MSEP) of
\eqn{\widehat{\mathbf{Y}}_{\mathrm{CMCK}}}{Yhat_CMCK} is given by
\deqn{
  \mathrm{MSPE}[\widehat{\mathbf{Y}}_{\mathrm{CMCK}}] = \mathrm{MSPE}[
  \widehat{\mathbf{Y}}_{\mathrm{UK}} ] +
  (\mathbf{P}_{1}-\mathbf{Q}_{1})(\mathbf{P}_{1}-\mathbf{Q}_{1}).
}{MSPE[Yhat_CMCK] = MSPE[ Yhat_(UK)] +(P1 + Q1) (P1 + Q1 )}
and of \eqn{\widehat{Y}_{\mathrm{CK}}(B_1)}{Yhat_CK(B_1)} by
\deqn{
  \mathrm{MSPE}[\widehat{Y}_{\mathrm{CK}}(B_{1})] = \mathrm{MSPE}[
  \widehat{Y}_{\mathrm{UK}}(B_{1})] + (P^{\frac{1}{2}} - Q^{\frac{1}{2}})^{2},
}{MSPE[Yhat_CK(B_1)] = MSPE[Yhat_UK(B_1)] +( P^0.5 - Q^0.5 )^2,}
where the MSEP of universal Kriging is given by
\deqn{
  \mathrm{MSPE}[\widehat{\mathbf{Y}}_\mathrm{UK}] =
  \mathrm{Cov}[\mathbf{Y}, \mathbf{Y}^{\prime}] -
  \mathbf{C}^{\prime}\boldsymbol{\Sigma}^{-1}\mathbf{C} +
  (\mathbf{X}_{m}^{\prime} - \mathbf{X}^{\prime}\boldsymbol{\Sigma}^{-1}\mathbf{C})^{\prime}
  (\mathbf{X}^{\prime}\boldsymbol{\Sigma}^{-1}\mathbf{X})^{-1}
  (\mathbf{X}_{m}^{\prime} - \mathbf{X}^{\prime}\boldsymbol{\Sigma}^{-1}\mathbf{C}).
}{Cov[Y, Y'] - C' \Sigma^{-1} C +
  (X_m - X' \Sigma^{-1} C)' (X' \Sigma^{-1} C) (X_m - X' \Sigma^{-1} C).
}
}

\author{
Christoph Hofer, \email{christoph.hofer@alumni.ethz.ch}
}

\references{%
Aldworth, J. and Cressie, N. (2003). Prediction of non-linear spatial
functionals. \emph{Journal of Statistical Planning and Inference},
\bold{112}, 3--41, \doi{10.1016/S0378-3758(02)00321-X}.

Cressie, N. (1993). Aggregation in geostatistical problems. In
A. Soares, editor, \emph{Geostatistics Troia 92}, \bold{1}, pages 25--36,
Dordrecht. Kluwer Academic Publishers, \doi{10.1007/978-94-011-1739-5_3}.

Hofer, C. and Papritz, A. (2010).  Predicting threshold exceedance by
local block means in soil pollution surveys.  \emph{Mathematical
Geosciences}.  \bold{42}, 631--656, \doi{10.1007/s11004-010-9287-4}

Hofer, C. and Papritz, A. (2011).  constrainedKriging: an R-package
for customary, constrained and covariance-matching constrained point
or block Kriging.  \emph{Computers & Geosciences}. \bold{37}, 1562--1569,
\doi{10.1016/j.cageo.2011.02.009}

}
%
%
\keyword{ package }
