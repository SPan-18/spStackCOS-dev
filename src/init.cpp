#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include "spStackCOS.h"

static const R_CallMethodDef CallEntries[] = {
  {"idist",                   (DL_FUNC) &idist,                   6},
  {"R_cholRankOneUpdate",     (DL_FUNC) &R_cholRankOneUpdate,     6},
  {"R_cholRowDelUpdate",      (DL_FUNC) &R_cholRowDelUpdate,      4},
  {"R_cholRowBlockDelUpdate", (DL_FUNC) &R_cholRowBlockDelUpdate, 5},
  {"spLMexact",               (DL_FUNC) &spLMexact,               14},
  {"spLMexactLOO",            (DL_FUNC) &spLMexactLOO,            16},
  {"sptLMexactCOS",           (DL_FUNC) &sptLMexactCOS,           26},
  {"sptLMexactTimeAgg",       (DL_FUNC) &sptLMexactTimeAgg,       13},
  {"sptPointPredictTimeAgg",  (DL_FUNC) &sptPointPredictTimeAgg,  16}
};

void
#ifdef HAVE_VISIBILITY_ATTRIBUTE
__attribute__ ((visibility ("default")))
#endif
  R_init_sp(DllInfo *dll)
  {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
  }
