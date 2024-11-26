#ifndef BLNS_ISSUER_H
#define BLNS_ISSUER_H

#include "params.h"
#include "Lattice.h"
#include "Com.h"


typedef struct
{
    ZZ_pX       a1;
    vec_ZZ_pX   a2;
    vec_ZZ_pX   c0;
    vec_ZZ_pX   c1;
} IPK_STRUCT;

 
void    I_KeyGen(IPK_STRUCT& ipk, mat_L& isk);
void    I_VerCred(vec_ZZ& s_0, vec_ZZX& w, ZZ& x, const CRS_Data2& crs, const mat_ZZ_p& B_f, const IPK_STRUCT& ipk, const mat_L& isk, const Vec<string>& attrs_prime, const ZZ_pX& u, const PROOF_Com& Pi);

#endif