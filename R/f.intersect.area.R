f.intersect.area <- function( polygons, gpc.pixel)

### purpose: function returns a list with the intersect area between
###          the pixels and the polygons
### arguments:
###            polygons = polygon list
###            gpc.pixel = gpc.poly object of the pixel
###
### author: ch.hofer
###  date: 27.2.2007
{
t.intersect.area.list <- lapply( gpc.pixel,
                                function( gpc.pixel, polygons){

                                    ### old return( area.poly( intersect( gpc.pixel, polygons ) ) )
                                    area <- try(
                                                area.poly(
                                                          try( intersect( gpc.pixel, polygons ), silent = T )
                                                          ),
                                                          silent = T
                                                )
                                    if( is( area ) == "try-error" ){ area <- 0}
                                    return( area )
                                },
                                polygons
			    )
return( as.vector( t.intersect.area.list, mode = "numeric"))
}
