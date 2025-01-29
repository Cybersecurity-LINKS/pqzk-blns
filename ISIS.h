#ifndef BLNS_ISIS_H
#define BLNS_ISIS_H

#include <sstream>

#include "params.h"
#include "LHC.h"
#include "Hash.h"
#include "Squares.h"
#include "Issuer.h"


typedef struct
{
    vec_zz_pX       t_A, t_y, t_g, w;    
    vec_ZZ          z_3;
    vec_zz_pX       h;
    zz_pX           t, f0;
    vec_ZZX         z_1, z_2;
    int             valid;
} PROOF_I_t;


void Preprocessing_ISIS(vec_ZZ& s1, vec_ZZ& r1, const vec_ZZ& s0, const ZZ& B_goth_s2, const vec_ZZ& r0, const ZZ& B_goth_r2);
void Prove_ISIS(PROOF_I_t& Pi, const string& inputStr, const CRS_t& crs, const IPK_t& ipk, const mat_ZZ& P0, const mat_ZZ& C0, const vec_zz_p& mex, const mat_ZZ& B0, const vec_ZZ& Bounds, const ZZ& aux, const Vec<vec_ZZ>& w0);
int  Verify_ISIS(const string& inputStr, const CRS_t& crs, const IPK_t& ipk, const mat_ZZ& P0, const mat_ZZ& C0, const vec_zz_p& mex, const mat_ZZ& B0, const vec_ZZ& Bounds, const ZZ& aux, const PROOF_I_t& Pi );

#endif