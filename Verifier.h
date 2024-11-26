#ifndef BLNS_VERIFIER_H
#define BLNS_VERIFIER_H

#include "params.h"
#include "Holder.h"
#include "ISIS.h"


int  V_Verify(const VP_STRUCT& VP, const CRS_Data2& crs, const mat_ZZ_p& B_f);

#endif