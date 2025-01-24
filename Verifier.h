#ifndef BLNS_VERIFIER_H
#define BLNS_VERIFIER_H

#include "params.h"
#include "Holder.h"
#include "ISIS.h"


int  V_Verify(const VP_t& VP, const string& inputStr, const CRS2_t& crs, const mat_zz_p& B_f);

#endif