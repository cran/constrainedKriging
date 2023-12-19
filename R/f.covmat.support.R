#######################################
#                                    #
#   Berechnung der Kovarianzmatrix   #
#   ch: 22-02-2010                   #
#                                    #
######################################

# 2023-11-20 A. Papritz subsitution of function
#                       CovarianceFct{RandomsFields} by f.CorrelationFct
# 2023-12-10 AP correction of error for nugget model occurring with scale = 0

f.covmat.support <- function(model, locations)
### purpose: calculate the the n X n  support covarianzmatrix
### arguments:
###  model= list with the covarianz parameter,
###                       1. element = coords of the center of the lower left pixel
###  locations  = Koordinaten der St?tzpunkte (1D Vektor oder n X
###  			, oder n X 3 Matrix
###
### output: n X n Kovarianzmatrix der Distanzen
{

  t.distmat <- f.row.dist( locations, locations )

  t.part.covmat.list <- lapply(
    model,
    function(x, t.distmat){

      if(x$model == "mev"){
        x$model = "nugget"
      }

#       # using function CovarianceFct{RandomsFields}
#
#       t.part.covmat <- x$variance * CovarianceFct(
#         x = c( t.distmat ),
#         model = x$model,
#         param = c(mean = NA,
#           variance = 1,
#           nugget = 0,
#           scale = x$scale,
#           parameter = x$parameter)
#       )

      # using substitute function f.CorrelationFct

      t.part.covmat <- x$variance * f.CorrelationFct(
        x = c( t.distmat ),
        scale = x$scale,
        model = x$model,
        param = x$parameter
      )

#       cat("\n\n", x$model, "  parameter", x$parameter,"\n")
#       print(range(t.part.covmat.new - t.part.covmat))


      dim( t.part.covmat ) <- dim( t.distmat )
      return( t.part.covmat )
    }
    , t.distmat
  )


  t.covmat.support <- matrix(
    rowSums( matrix(
        unlist(t.part.covmat.list),
        ncol = length(t.part.covmat.list) )),
    ncol = dim(t.part.covmat.list[[1]])[2]
  )

  rm(t.distmat)
  return(t.covmat.support)

}
