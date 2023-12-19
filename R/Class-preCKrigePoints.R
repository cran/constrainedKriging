# # # Definition of the preCKrigePoints class
# # # Christoph Hofer, 30-03-
# setClass(
#   Class = "preCKrigePoints",
#   representation( covmat = "list", posindex = "list", model = "list",
#     data = "data.frame", coords = "matrix" ),
#   prototype = list(),
#   sealed = TRUE,
#   S3methods =FALSE
# )

## 2023-01-25 A. Papritz revised class definition
preCKrigePoints <- setClass(
  "preCKrigePoints",
  slots = c(
    covmat = "list",
    posindex = "list",
    model = "list",
    data = "data.frame",
    coords = "matrix"
  )
)
