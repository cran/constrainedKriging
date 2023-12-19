f.build.pixel.sfc  <- function(
  posindex,
  polygons,
  sa.polygons,
  polygon.centroids,
  pixgridInfo,
  delta.x.pix,
  delta.y.pix)

# function computes the pixel representation of a polygon plus optional
# neighbours (= neighbourhood), pixels with centers in small area polygons
# are assigned to the polygon that has the largest intersection with the pixel

# arguments
#   posindex    numeric vector with indices of polygon and its optional neighbours
#   polygons    simple feature geometry list-column with the polygon and its optional neighbours
#   sa.polygons logical vector defining which of  the polygon and its optional
#               neighbours have small area (=area < pixel area)
#   polygon.centroids numeric matrix with centroids of polygon and its optional neighbours
#   pixgridInfo definition of pixel grid (computed by precompute{spatialCovariance}
#   delta.x.pix displacement of x-coordinate of origin of pixel grid
#   delta.y.pix displacement of y-coordinate of origin of pixel grid

# author of original code based on gpc.poly: Ch. Hofer date: 14.12.2006
# 2023-11-24 A. Papritz implementation for simple feature geometry
#                       list-column (sfc) object and further revision of code
{

  # parameters of pixel grid
  t.rowwidth <- pixgridInfo$rowwidth
  t.colwidth <- pixgridInfo$colwidth
  t.mesh.n.row <- pixgridInfo$nrows
  t.mesh.n.col <- pixgridInfo$ncol


  # bounding box of target polygons plus optionally neighbours
  t.grid.bbox <- matrix( st_bbox(polygons), nrow = 2 )

  ## print(t.grid.bbox)
  ## da das Zentrum des Pixels bestimmt ob ein Polygon approximiert wird
  ## wird der Startpunkt des Gitters nach links unten verschoben. Andernfalls
  ## haette die untere linke Ecke nicht die gleiche Wahrscheinlichkeit ein Pixel
  ## punkt zu besitzten

  # shift origin of pixel grid to left and bottom by half of pixel width
  # and height
  t.grid.bbox[1,1] <- t.grid.bbox[1,1] - 0.5 * t.colwidth
  t.grid.bbox[2,1] <- t.grid.bbox[2,1] - 0.5 * t.rowwidth

  # generate matrix with coordinates of pixel centers (optionally randomly displaced)
  t.grid.coords <- pixgridInfo$indices.preLimits
  pixcenter <- cbind(
    t.grid.bbox[1,1] + delta.x.pix + 0.5 * t.colwidth + t.grid.coords$bx,
    t.grid.bbox[2,1] + delta.y.pix + 0.5 * t.rowwidth + t.grid.coords$ax)


  # rearrage pixel centers such that they are sorted columwise
  if( dim(pixcenter)[1] > 1){
    pixcenter <- pixcenter[ order( pixcenter[,1] ), ]
  }

  # define pix center if NULL
  if( is.null(pixcenter) ){
    pixcenter <- matrix(pixcenter, ncol = 2)
  }

  # convert to pixel centers to simple feature geometry list-column (sfc
  # object)
    pixcenter.sfc <- st_sfc(lapply(
    1:NROW(pixcenter),
    function(i, pixcenter) st_point(pixcenter[i, ]),
    pixcenter = pixcenter
  ))

  # intersect pixelcenters with polygons, result is a matrix with entries
  # equal to one if pixel center lies in a polygon and zero otherwise
  pix.in.poly <- matrix(as.numeric(
    st_intersects(pixcenter.sfc, polygons, sparse = FALSE)
  ), nrow = length(pixcenter.sfc) )

  # number of pixel centers per polygon
  no.pix.in.poly <- colSums(pix.in.poly)

  # set elements of sa.polygons to TRUE for polygons without a pixelcenter
  sa.polygons[ no.pix.in.poly == 0 & sa.polygons == FALSE ] <- TRUE

  # assign pixels with centers in small area polygons to other polygons

  if( any( sa.polygons & no.pix.in.poly > 0 ) ){

    # generate pixel polygons

    pixel.sfc <- st_sfc(
      apply(
        pixcenter, 1,
        function(x, cw, rw){

          # define corners of pixel
          pixcorner <- rbind(
            cbind(t.l.x <- x[1] - 0.5*cw, t.l.y <- x[2] - 0.5*rw), ### lower left
            cbind(t.r.x <- x[1] + 0.5*cw, t.l.y), ## lower right
            cbind(t.r.x, t.r.y <- x[2] + 0.5*rw), ### upper right
            cbind(t.l.x, t.r.y) ## upper left
          )

          # convert to sfg polygon object
          pixpolygon <- st_polygon(list(rbind(pixcorner, pixcorner[1, ])))

          return(pixpolygon)

        },
        cw = t.colwidth,
        rw = t.rowwidth
      )
    )


    # compute changes in pix.in.poly due to small polygons

    t.pix.in.poly <- lapply(
      which(sa.polygons & no.pix.in.poly > 0),
      function(ipoly, pix.in.poly, polygons, pixel.sfc){
        res <- lapply(
          which(pix.in.poly[, ipoly] > 0),
          function(ipix, ipoly, pix.in.poly, polygons, pixel.sfc){

            # intersect pixel with center in current small area polygon with all
            # polygons

            # compute area of intersections of pixel with polygons

            t.inter.area <- rep(0., length(polygons))
            sel <- st_intersects(pixel.sfc[ipix], polygons, sparse = FALSE)
            t.inter.area[sel] <- st_area(
              st_intersection( pixel.sfc[ipix], polygons)
            )

            # index of polygon with maximum intersection area

            iipoly <- which.max(t.inter.area)[1]

            # adjust pix.in.poly and no.pix.in.poly

            pix.in.poly[ipix, ipoly]   <- 0
            pix.in.poly[ipix, iipoly]  <- 1

            # return result

            c(ipix, pix.in.poly[ipix, ])

          }, ipoly = ipoly,
          pix.in.poly = pix.in.poly, polygons = polygons, pixel.sfc = pixel.sfc
        )
      },
      pix.in.poly = pix.in.poly, polygons = polygons, pixel.sfc = pixel.sfc
    )

    # convert to matrix

    t.pix.in.poly <- matrix(
      unlist(t.pix.in.poly), ncol = length(polygons) + 1L, byrow = TRUE
    )

    # update pix.in.poly by changes

    pix.in.poly[t.pix.in.poly[, 1], ] <- t.pix.in.poly[, -1]

  }

  # return result

  list(
    pixcenter = pixcenter,
    rowwidth = t.rowwidth,
    colwidth = t.colwidth,
    nrows = t.mesh.n.row,
    ncols = t.mesh.n.col,
    #pixel.poly.dim = dim(pix.in.poly),
    no.pix.in.poly = no.pix.in.poly,
    sa.polygons = sa.polygons,
    polygon.centroids = matrix(polygon.centroids, ncol = 2),
    posindex = posindex,
    #which.poly = which.poly,
    pix.in.poly = pix.in.poly
  )


} ## end function

