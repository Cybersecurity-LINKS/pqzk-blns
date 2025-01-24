#ifndef BLNS_COM_H
#define BLNS_COM_H

#include <sstream>

#include "params.h"
#include "LHC.h"
#include "Hash.h"
#include "Squares.h"

typedef struct
{
    vec_zz_pX       t_A, t_y, t_g, w;
    Vec<vec_zz_pX>  com_1, com_2;
    vec_ZZ          z_3;
    vec_zz_pX       h;
    zz_pX           t, f0;
    vec_ZZX         z_1,  z_2;
    Vec<vec_ZZX>    op_1, op_2;
    int             valid;
} PROOF_Com;

#include "Issuer.h"


void       Preprocessing_Com(vec_ZZ& s1, const vec_ZZ& s, const ZZ B_goth2);
PROOF_Com  Prove_Com( const string inputStr, const CRS_Data& crs, const IPK_STRUCT& ipk, const mat_ZZ& P0, const vec_ZZ& u0, const ZZ B_goth2, const vec_ZZ& w0);
int        Verify_Com(const string inputStr, const CRS_Data& crs, const IPK_STRUCT& ipk, const mat_ZZ& P0, const vec_ZZ& u0, const ZZ B_goth2, const PROOF_Com& Pi);

#endif