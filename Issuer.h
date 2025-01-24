#ifndef BLNS_ISSUER_H
#define BLNS_ISSUER_H

#include "params.h"
#include "Lattice.h"

typedef struct
{
    zz_pX       a1;
    vec_zz_pX   a2;
    vec_zz_pX   c0;
    vec_zz_pX   c1;
} IPK_t;

#include "Com.h"

 
void    I_KeyGen(IPK_t& ipk, mat_L& isk);
void    I_VerCred(vec_ZZ& s_0, vec_ZZX& w, ZZ& x, const string& inputStr, const CRS2_t& crs, const mat_zz_p& B_f, const IPK_t& ipk, const mat_L& isk, const Vec<string>& attrs_prime, const zz_pX& u, const PROOF_C_t& Pi);

#endif