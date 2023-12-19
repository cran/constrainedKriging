f.bbox.information.sfc <- function( polygons, neighbours )
## purpose: find the largest bounding box of a polygon and its neighbours
## (= neighbourhood) and find the polygon with the smallest area
##
## arguments:
##   polygons     simple feature geometry list-column with the polygons
##   neighbours   list with each element is a vector with the polygons
##                indices of the neighbours
##
## wholes in polygons are not considered !!!!!!!!
##
##
## author: Ch.  Hofer date: 29.1.2007

## 2023-11-23 A. Papritz implementation for simple feature geometry
##                       list-column (sfc object); slight revision of code;
##                       function returns now width and height of largest
##                       bounding box (and no longer maxima of bounding box
##                       widths and heights)
##
{

# compute bounding box for all polygons plus their neighbours

  t.bbox.list <- lapply(
    1:length(polygons),
    function(i, polygons, neighbours){

       t.bbox.surrounding <- matrix(
        st_bbox(polygons[c(i, neighbours[[i]])]),
        nrow = 2, byrow = TRUE
      )

      t.dx.neighbourhood <- diff(t.bbox.surrounding[, 1])
      t.dy.neighbourhood <- diff(t.bbox.surrounding[, 2])

      # return width, height, and area of bounding box of ith
      # neighbourhood and area of ith polygon
      return(
        list(
          dx.neighbourhood = t.dx.neighbourhood,
          dy.neighbourhood = t.dy.neighbourhood,
          area.polygon = st_area(polygons[i]),
          area.neighbourhood = t.dx.neighbourhood * t.dy.neighbourhood
        )
      )
    },
    polygons = polygons, neighbours = neighbours
  )

  # extract indices of largest and smallest elements in t.bbox.list

  t.bbox.matrix <- matrix( unlist( t.bbox.list), ncol = 4, byrow = TRUE )
  t.which.max.bbox <- apply(t.bbox.matrix, 2, which.max)
  t.which.min.bbox <- apply(t.bbox.matrix, 2, which.min)

  # return indices of largest neighbourhood and of smallest polygon along
  # with smallest polygon area and width and height of largest bounding box
  return(	list(
      largest.box.index = t.which.max.bbox[4],
      min.area.index = t.which.min.bbox[3],
      min.area = t.bbox.matrix[t.which.min.bbox[3],3],
      t.xy.range = cbind(
        t.bbox.matrix[t.which.max.bbox[1] ,1],
        t.bbox.matrix[t.which.max.bbox[2] ,2]
      )
    )
  )
} #end function
