#####################################################################
#                                                                   #
#   R Funktion um Punkt Rechteck Kovarianz der Pixel zu berechnen   #
#   ch: 23-02-2010                                                  #
#                                                                   #
#########################################################################
####
#### point block covariance in C
####
f.point.rec.covmat<- function(
  pxcoord  = 0,
  pycoord = 0,
  crxcoord = 0,
  crycoord = 0,
  rxwidth = 0,
  rywidth = 0,
  epsabs =.Machine$double.eps,
  epsrel = epsabs,
  subdivisions = 10000,
  param = list()
)
{
  # Name des Packages in dem sich die Funktion befindet
  thispackage  <- "constrainedKriging"
  #n.mod =AnzahlderKovrianzmodelle,nugggt,exponetiell
  n.mod <- length( param )
  stopifnot(n.mod > 0)

  # Nugget Parameter auf Null setzen
  param <- lapply( param, function(x){

      if( is.null(x$parameter) ){ x$parameter <- 0 }
      return(x)

    }
  )

  # Anzahl Parameter pro Kovarianzmodell
  n.mod.param <- unlist( lapply( param, function(x){ length(x$parameter) } ) )


  # model name wird codiert:
  # nugget =1, exponetial =2, usw.

  param <- lapply(
    param,
    function(x)
    {
      if( x$model== "nugget" ){ x$model <- 0; x$parameter <- 0 }

      if( x$model== "exponential" ){ x$model <- 1; x$parameter = c(1) }

      if( x$model == "stable" ){ x$model <- c(1) }

      if( x$model== "gauss" ){ x$model <- 1 ; x$parameter = c(2)  }

      if( x$model== "spherical" ){ x$model <- 2}

      if( x$model== "bessel" ){ x$model <- 4 }

      if( x$model== "cauchy" ){ x$model <- 5 }

      if( x$model== "cauchytbm" ){ x$model <- 6 }

      if( x$model== "circular" ){ x$model <- 7 }

      if( x$model== "constant" ){ x$model <- 8 }

      if( x$model== "cubic" ){ x$model <- 9 }

      if( x$model== "dampedcosine" ){ x$model <- 10 }

      if( x$model== "gencauchy" ){ x$model <- 11 }

      if( x$model== "gengneiting" && x$parameter == 1 ){ x$model <- 12 }

      if( x$model== "gengneiting" && x$parameter == 2 ){ x$model <- 13 }

      if( x$model== "gengneiting" && x$parameter == 3 ){ x$model <- 14 }

      if( x$model== "gneiting" ){ x$model <- 15 }

      if( x$model== "hyperbolic" ){ x$model <- 16 }

      if( x$model== "penta" ){ x$model <- 17 }

      if( x$model== "lgd1" ){ x$model <- 18 }

      if( x$model== "power" ){ x$model <- 19 }

      if( x$model== "wave" ){ x$model <- 20 }

      if( x$model== "qexponential" ){ x$model <- 21 }

      if( x$model == "matern" && x$parameter == 0.5 ){ x$model <- c(1) ; x$parameter = c(1) }
      if( x$model == "matern" && x$parameter != 0.5 ){ x$model <- c(3) }

      if( x$model == "whittle"){ x$model <- 22 }



      return(x)
    }# end of function
  ) # end of lapply

  #Die Liste param wird in ein Vektor umgewandelt mit den folgenden Elementen
  # Modell Nummer, Sill, Range, Parameter, Modell Nummer, Sill, Range, Parameter, usw...
  # z.B c( 0, 7, 0, 0, 1, 23,23, 4, 0.5)
  #  0 = Nugget, variance = 7, scale = 0, parameter = 0, 1 = Exponential, variance = 23, scale = 23.4, parameter = 0.5 (kappa)
  #

  covpar <- vector()
  j = 1
  for(i in 1:length(param)){
    covpar[j] <- param[[i]][names(param[[i]]) == "model"][[1]]
    j= j+1
    covpar[j] <- param[[i]][names(param[[i]]) == "variance"][[1]]
    j = j+1
    covpar[j] <- param[[i]][names(param[[i]]) == "scale"][[1]]
    j = j+1
    covpar[j:( j + (n.mod.param[i]-1) )] <- param[[i]][names(param[[i]]) == "parameter"][[1]]
    j = j + n.mod.param[i]
  }

  # number of support.points
  if( is.null( dim(pxcoord) ) ){
    npoints <- length(pxcoord)
  } else {
    npoints <- dim(pxcoord)[1]
  }
  # number of pixel
  if( is.null( dim(crxcoord) ) ){
    npixel <- length(crxcoord)
  } else {
    npixel <- dim(crxcoord)[1]
  }

  limit  <-  subdivisions
  lenw  <-  4*limit

  t.result <- .C(
    "PointRectCov",
    as.double(pxcoord), as.double(pycoord), as.double(crxcoord), as.double(crycoord),
    as.double(rxwidth), as.double(rywidth), as.double(npoints), as.double(npixel),
    as.double(epsabs), as.double(epsrel), result = double( npoints * npixel), abserr = as.double(1),
    as.integer(limit), as.integer(lenw + 1), as.integer(lenw), as.integer(lenw +1),
    as.integer(limit+1), integer(limit+1), double(lenw+1),  as.double(covpar),
    as.integer(n.mod), as.integer(n.mod.param),
    DUP = TRUE, NAOK=FALSE, PACKAGE = thispackage
  )
  return( t.result)
} # end of function

