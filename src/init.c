#include "mean.circular.h"
#include "median.circular.h"
#include "medianHL.circular.h"
#include "weighted.mean.circular.h"
#include "minuspipluspi.h"
#include "rvonmises.h"
#include "distance.h"
#include "dwrappednormal.h"
#include "mle.wrappednormal.h"
#include <R_ext/RS.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

#define C_DEF(name, n)  {#name, (DL_FUNC) &name, n}

static const R_CMethodDef CEntries[]  = {
    C_DEF(MeanCircularRad, 3),
    C_DEF(MedianCircularRad, 5),
    C_DEF(MedianHLCircularRad, 5),
    C_DEF(MedianHLCircularPropRad, 5),
    C_DEF(sampleReplace, 4),
    C_DEF(sampleNoReplace, 5),
    C_DEF(MinusPiPlusPiRad, 2),
    C_DEF(WeightedMeanCircularRad, 4),
    C_DEF(rvm, 4),
    C_DEF(R_angularseparation, 5),
    C_DEF(R_chord, 5),
    C_DEF(R_geodesic, 5),
    C_DEF(R_correlation, 5),
    C_DEF(R_distance, 6),
    {NULL, NULL, 0}
};

static const R_FortranMethodDef FortEntries[] = {
    {"dwrpnorm", (DL_FUNC) &F77_NAME(dwrpnorm),  7},
    {"mlewrpno", (DL_FUNC) &F77_NAME(mlewrpno),  10},    
    {NULL, NULL, 0}
};

void attribute_visible R_init_circular(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, FortEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
