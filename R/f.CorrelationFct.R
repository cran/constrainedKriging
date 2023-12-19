##   ##############################################################################
### f.CorrelationFct

f.CorrelationFct <- function( x, model, scale, param = NULL ){

  ## function computes correlation of some stationary autocorrelation
  ## models available in the RandomFields package version 2.  the function
  ## is used as a substitute for CovarianceFct{RandomFields} the
  ## parametrization of the models follow the help page of
  ## CovarianceFct{RandomFields}.

  ## arguments:
  ## x          	    lag distances
  ## model            name of autocorrelation model
  ## scale            range parameter
  ## param            extra parameters of autocorrelation model

  ## function returns vector with correlations

  ## 2023-11-20 Andreas Papritz
  ## 2023-12-10 AP correction of error for nugget model occurring with scale = 0

### - scale lag distances

  if(scale > 0.){
    x <- x / scale
  }

### - compute autocorrelations

  result <- switch(
    model[1],

#### --- bessel
    bessel = {

      A <- unname( param[1] )

      res <- rep( 1., length( x ) )
      sel <- x > 0.
      xx <- x[sel]

      res[sel] <- ( 2.^A * gamma(1.+A) * xx^(-A) * besselJ( xx, A ) )
      res

    },

#### --- cauchy
    cauchy = {

      A <- unname( param[1] )
      ( 1. + x^2 )^(-A)

    },

#### --- cauchytbm (mit \gamma = 3)
    cauchytbm = {

      A  <- unname( param[1] )
      B  <- unname( param[2] )

      ( 1. + ( 1. - B/3.) * x^A ) * ( 1. + x^A )^(-B / (A - 1.) )

    },

#### --- circular (compact support)
    circular = {

      res <- rep( 0., length( x ) )
      sel <- x < 1.
      xx <- x[sel]

      res[sel] <- 1. -2./pi * ( xx * sqrt( 1. - xx^2 ) + asin( xx ) )
      res

    },

#### --- constant
    constant = {

      rep( 1., length( x ) )

    },

#### --- cubic (compact support)
    cubic = {

      res <- rep( 0., length( x ) )
      sel <- x < 1.
      xx <- x[sel]

      res[sel] <- 1. - 7*xx^2 + 8.75*xx^3 - 3.5*xx^5 + 0.75*xx^7
      res

    },

#### --- dampedcosine
    dampedcosine = {

      A <- unname( param[1] )
      exp( -A * x ) * cos( x )

    },

#### --- exponential
    exponential = {

      exp( -x )

    },

#### --- gauss
    gauss = {

       exp( -x^2 )

    },

#### --- spherical (compact support)
    spherical = {

      res <- rep( 0., length( x ) )
      sel <- x < 1.
      xx <- x[sel]

      res[sel] <- 1. -1.5 * xx + 0.5 * xx^3
      res

    },

#### --- gencauchy
    gencauchy = {

      A <- unname( param[1] )
      B <- unname( param[2] )

      ( 1. + x^A )^(-B/A)

    },

#### --- gengneiting (compact support)
    gengneiting = {

      A <- unname( param[1] ) # n
      B <- unname( param[2] ) # alpha

      res <- rep( 0., length( x ) )
      sel <- x < 1.
      xx <- x[sel]

      res[sel] <- if( identical( as.integer(A), 1L ) ){

        BB <- B + 1.
        ( 1. + BB * xx ) * ( 1. - xx)^BB

      } else if( identical( as.integer(A), 2L ) ){

        BB <- B + 2.
        ( 1. + BB * xx + (BB^2 - 1.) * xx^2 / 3. ) * (1. - xx)^BB

      } else if( identical( as.integer(A), 3L ) ){

        BB <- B + 3
        (
          1 + BB * xx + (2. * BB^2 - 3.) * xx^2 / 5. +
          (BB^2 - 4.) * BB * xx^3 / 15.
        ) * (1. - xx)^BB

      } else {
        stop( "gengneiting model undefined for 'n' !%in% 1:3" )
      }
      res

    },

#### --- gneiting (compact support)
    gneiting = {

      res <- rep( 0., length( x ) )
      sel <- x < 1. / 0.301187465825
      xx <- 0.301187465825 * x[sel]

      res[sel] <- (1. + 8.* xx + 25. * xx^2 + 32. * xx^3 ) * (1. - xx )^8
      res

    },

#### --- hyperbolic
    hyperbolic = {

      A <- unname( param[1] ) # nu
      B <- unname( param[2] ) # lambda
      CC <- unname( param[3] ) # delta

      xx <- CC^2 + x^2
      xx^(0.5 * B) * besselK(A * sqrt(xx), B) / ( CC^B * besselK(A * CC, B) )

    },

#### --- lgd1 (cf Gneiting, T. / Schlather, M.
  ## Stochastic models that separate fractal dimension and the Hurst effect,
  ## 2004, SIAM Review , Vol. 46, No. 2, p. 269-282
    lgd1 = {

      A <- unname( param[1] )
      B <- unname( param[2] )

      res <- rep( NA_real_, length( x ) )
      sel <- x <= 1.

      res[sel]  <- 1. - B / (A + B) * x[sel]^A
      res[!sel] <- A / (A + B) * x[!sel]^(-B)
      res

    },

#### --- nugget
    nugget = {

      res <- rep( 0., length( x ) )
      sel <- x <= 0.
      res[sel] <- 1.
      res

    },

#### --- penta (compact support)
    penta = {

      res <- rep( 0., length( x ) )
      sel <- x < 1.
      xx <- x[sel]

      res[sel] <- 1. - 22./3. * xx^2 + 33. * xx^4 - 38.5 * xx^5 +
        16.5 * xx^7 - 5.5 * xx^9 + 5./6. * xx^11
      res

    },

#### --- power (compact support)
    power = {

      A <- unname( param[1] )

      res <- rep( 0., length( x ) )
      sel <- x < 1.
      xx <- x[sel]

      res[sel] <- (1 - xx)^A
      res

    },

#### --- wave
    wave = {

      res <- rep( 1., length( x ) )
      sel <- x > 0.
      xx <- x[sel]

      res[sel] <- sin( xx ) / xx
      res

    },

#### --- qexponential
    qexponential = {

      A <- unname( param[1] )

      ( 2. * exp(-x) - A * exp(-2.*x) ) / (2. - A)

    },

#### --- matern
    matern = {

      A <- unname( param[1] )

      res <- rep( 1., length( x ) )
      sel <- x > 0.
      xx <- sqrt(2. * A) * x[sel]

      res[sel] <- 2.^(1. - A) / gamma(A) * xx^A * besselK( xx, A )
      res

    },

#### --- whittle
    whittle = {

      A <- unname( param[1] )

      res <- rep( 1., length( x ) )
      sel <- x > 0.
      xx <- x[sel]

      res[sel] <- 2.^(1 - A) / gamma(A) * xx^A * besselK( xx, A )
      res

    },

#### --- stable
    stable = {

      A <- unname( param[1] )

      exp( -x^A )

    },

    stop( "correlation model ", model, " undefined" )

  )

  return(result)

}
