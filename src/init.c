// 2023-11-30 A. Papritz registration of C function PointRectCov
// cf section 5.4 Manual Writing R Extensions, Version 4.3.2

#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include <Rinternals.h>


/* .C calls */
extern void PointRectCov(
                         double *, double *, double *, double *, 
                         double *, double *, double *, double *, 
                         double *, double *, double *, double *, 
                         int *, int *, int *, int *, 
                         int *, int *, double *, double *, 
                         int *, int *);

static R_NativePrimitiveArgType PointRectCov_type[] = {
         REALSXP, REALSXP, REALSXP, REALSXP, 
         REALSXP, REALSXP, REALSXP, REALSXP, 
         REALSXP, REALSXP, REALSXP, REALSXP, 
         INTSXP, INTSXP, INTSXP, INTSXP, 
         INTSXP, INTSXP, REALSXP, REALSXP, 
         INTSXP, INTSXP
};

static const R_CMethodDef CEntries[] = {
    {"PointRectCov", (DL_FUNC) &PointRectCov, 22, PointRectCov_type},
    {NULL, NULL, 0, NULL}
};

void R_init_constrainedKriging(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
