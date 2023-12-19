f.pixconfig.sfc <- function(
  polygons,
  neighbours,
  pixgrid,
  n = 1
)
## purpose: this function build the pixel representation for all polygons
## plus optionally their neighbours
##
## arguments:
##   polygons     simple feature geometry list-column with the polygons
##   neighbours   list with each element is a vector with the polygons
##                indices of the neighbours
##   pixgrid      information about pixel grid for the largest bounding box
##                (output of the function f.pixelgrid.sfc)
##   n            number of repetitions of pixel representations
##
## author: Ch. Hofer
## date: 14.12.2006

## 2023-11-23 A. Papritz implementation for simple feature geometry
##                       list-column (sfc object); slight revision of code;
{

  #
  pixgridInfo = pixgrid$pixgridInfo

  # logical vector defining polygons with area smaller than pixel area
  sa.polygons <- st_area(polygons) < pixgridInfo$rowwidth * pixgridInfo$colwidth

  # centroids of polygons
  polygon.centroids <- t(sapply(
    st_centroid( polygons ),
    st_coordinates
  ))

  ### the left lower corner of the pixel grid is choosed randomly therefore
  ### one have to set a seed to get reproducable results
  # # # if( n == 1 ){
  # # #  print( "Pixelgrid is fix!!!\n" )
  # # # }
  # # # else
  # # # {
  # # # cat( "The mean of ", n, "random Pixelgrids are calculated  !!!\n" )
  # # # }

  # list of indices with target polygon plus optional neighbours
  polygons.config <- lapply(
    1L:length( polygons ),
    function(i, neighbours)
    {
      t.poly.config <- c(i, neighbours[[i]])
    },
    neighbours = neighbours
  )

  # vector of displacement of origin of pixel grid
  if( n > 1){
    # random displacement when several repetition of pixel representation
    delta.x.pix <-  runif( 1, -0.5*pixgridInfo$colwidth, pixgridInfo$colwidth )
    delta.y.pix <-  runif( 1, -0.5*pixgridInfo$rowwidth, pixgridInfo$rowwidth )
  } else {
    # fixed displacement when just one pixel representation is used
    delta.x.pix <-  0.5 * pixgridInfo$colwidth
    delta.y.pix <-  0.5 * pixgridInfo$rowwidth
  }

  # compute pixel representation of polygons plus their optional neighbours
  # loop over all polygons

  pixconfig <- lapply(
    polygons.config,
    function(
      ith.poly,
      polygons, sa.polygons, polygon.centroids,
      pixgridInfo, t.delta.x.pix, t.delta.y.pix
    ){

      pixel.config <- f.build.pixel.sfc(
        posindex = ith.poly,
        polygons = polygons[ ith.poly ],
        sa.polygons = sa.polygons[ ith.poly ],
        polygon.centroids = polygon.centroids[ ith.poly, , drop = FALSE],
        pixgridInfo  = pixgridInfo,
        delta.x.pix = delta.x.pix,
        delta.y.pix = delta.y.pix
      )

      return( pixel.config )
    },
    polygons,
    sa.polygons,
    polygon.centroids,
    pixgridInfo,
    delta.x.pix,
    delta.y.pix
  )

  return( pixconfig )
} ## end function

