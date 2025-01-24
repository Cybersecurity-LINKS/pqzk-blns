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
} STATE_t;

typedef struct
{
    vec_ZZX     s;
    vec_ZZX     r;
    ZZ          x;
    int         valid;
} CRED_t;

typedef struct
{
    // vec_ZZ      cp;
    IPK_t  ipk;
    Vec<string> attrs_prime;
    vec_ZZ      idx;
    PROOF_I_t   pi;
    int         valid;
} VP_t;


void    H_Init(CRS2_t& crs, Vec<string>& attrs, const string& inputStr);

void    H_VerCred1(zz_pX& u, PROOF_C_t& Pi, STATE_t& state, const string& inputStr, const CRS2_t& crs, const IPK_t& ipk, const Vec<string>& attrs);
void    H_VerCred2(CRED_t& cred, const IPK_t& ipk, const mat_zz_p& B_f, const vec_ZZ& s_0, const vec_ZZX& w, const ZZ& x, const STATE_t& state);
void    H_VerPres(VP_t& VP, const CRED_t& cred, const string& inputStr, const CRS2_t& crs, const IPK_t& ipk, const mat_zz_p& B_f, const Vec<string>& attrs);

#endif