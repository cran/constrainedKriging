# # # Definition of the preCKrigePolygons class
# # # Christoph Hofer, 30-03-2010
# setClass(
#   Class = "preCKrigePolygons",
#   representation( covmat = "list", se.covmat = "list", pixconfig = "list",
#     pixcovmat = "matrix",  model = "list", data = "data.frame",
#     polygons = "list" ),
#   prototype = list(),
#   sealed = TRUE,
#   S3methods = FALSE
# )

## 2023-01-25 A. Papritz revised class definition
preCKrigePolygons <- setClass(
  "preCKrigePolygons",
  slots = c(
    covmat = "list",
    se.covmat = "list",
    pixconfig = "list",
    pixcovmat = "matrix",
    model = "list",
    data = "data.frame",
    polygons = "list"
  )
)
