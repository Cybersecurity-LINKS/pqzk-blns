#ifndef BLNS_HOLDER_H
#define BLNS_HOLDER_H

#include "params.h"
#include "Issuer.h"
#include "Hash.h"
#include "Com.h"
#include "ISIS.h"


typedef struct
{
    vec_ZZX     m;
    vec_ZZX     r;
} STATE_STRUCT;

typedef struct
{
    vec_ZZX     s;
    vec_ZZX     r;
    ZZ          x;
    int         valid;
} CRED_STRUCT;

typedef struct
{
    // vec_ZZ      cp;
    IPK_STRUCT  ipk;
    Vec<string> attrs_prime;
    vec_ZZ      idx;
    PROOF_ISIS  pi;
    int         valid;
} VP_STRUCT;


void    H_Init(CRS_Data2& crs, Vec<string>& attrs, const string inputStr);

void    H_VerCred1(zz_pX& u, PROOF_Com& Pi, STATE_STRUCT& state, const CRS_Data2& crs, const IPK_STRUCT& ipk, const Vec<string>& attrs);
void    H_VerCred2(CRED_STRUCT& cred, const IPK_STRUCT& ipk, const mat_zz_p& B_f, const vec_ZZ& s_0, const vec_ZZX& w, const ZZ& x, const STATE_STRUCT& state);
void    H_VerPres(VP_STRUCT& VP, const CRED_STRUCT& cred, const CRS_Data2& crs, const IPK_STRUCT& ipk, const mat_zz_p& B_f, const Vec<string>& attrs);

#endif