## computeV call the function K if the class attribut in computeV is
## set as "special" for more ?computeV of the  R-package spatialCovariance
K <- function( dist, model )
{

# 2023-11-20 A. Papritz subsitution of function
#                       CovarianceFct{RandomsFields} by f.CorrelationFct
# 2023-12-10 AP correction of error for nugget model occurring with scale = 0

  #   # using function CovarianceFct{RandomsFields}
  #
  #   t.cov <- CovarianceFct( x = dist, model = model$model,
  #     param = c(NA, variance = 1,
  #       nugget = 0, scale = model$scale, model$parameter )
  #   )

  # using substitute function f.CorrelationFct

  t.cov <- f.CorrelationFct(
    x = c(dist),
    scale = model$scale,
    model = model$model,
    param = model$parameter
  )

  return( t.cov )
}
# # #
# # #
# # #
f.pp.cov <- function(t.dist, model)
### purpose: calculate the  spatial covariance for a
###          given distance vector with the cov model in model
###
### arguments:
###            t.dist = numeric vector with distances
###            model = list with covariance parameter
###                        one element of the list =
###                        $model = cov model
###                        $scale = range parameter
###                        $parameter = vector with the cov model parameter
###                        $variance = sill matrix of the cov model
###
### author: ch.hofer
### date: 20.2.2006
{

# 2023-11-20 A. Papritz subsitution of function
#                       CovarianceFct{RandomsFields} by f.CorrelationFct
# 2023-12-10 AP correction of error for nugget model occurring with scale = 0


  t.part.cov.list <- lapply(
    model,
    function(model, t.dist){

      if( model$model != "special" ){

        #         # using function CovarianceFct{RandomsFields}
        #
        #         t.part.cor <- CovarianceFct(
        #           x = t.dist,
        #           model = model$model,
        #           param = c(
        #             mean = NA,
        #             variance = 1,
        #             nugget = 0,
        #             scale.par = model$scale,
        #             cov.par = model$parameter
        #           )
        #         )

        # using substitute function f.CorrelationFct
        t.part.cor <- f.CorrelationFct(
          x = c(t.dist),
          model = model$model,
          scale = model$scale,
          param = model$parameter
        )


      } else {
        t.part.cor <- K(
          dist = t.dist,
          model = model
        )
      }
      return( t.part.cor * model$variance )
    },
    t.dist = t.dist
  )

  t.part.cov.mat <- matrix( unlist( t.part.cov.list ), ncol = length( model ) )
  rm( t.part.cov.list )
  return( rowSums( t.part.cov.mat ) )

} ## end function f.pp.cov
