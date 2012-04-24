.CONSTRAINEDKRIGING_CACHE <- new.env(FALSE, parent=globalenv())

.onAttach <- function(lib, pkg) {
    assign("gpclib", FALSE, envir=.CONSTRAINEDKRIGING_CACHE)
    packageStartupMessage(paste("\nType 'help(constrainedKriging)' for a brief introduction to constrainedKriging.\n\n",
	    "\n\tNote: constrainedKrigig depends on the package gpclib,\n",
	    "\twhich has a restricted license. The non-commercial use of the\n",
	    "\tgpclib software is free but the commercial use is prohibited.\n \n",
	      "\tType 'gpc.copyright()' for details.\n"
	      ))
}

.onUnload <- function(libpath) {
    rm(.CONSTRAINEDKRIGING_CACHE)
}
