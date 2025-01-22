#ifndef BLNS_RANDOM_H
#define BLNS_RANDOM_H

#include "params.h"
#include "Utils.h"


void NTRU_TrapGen(zz_pX& a1, mat_L& B);

void preGSampler(vec_ZZ& v, const mat_L& B, const RR& sigma, const vec_ZZ& c);
void GSampler(vec_ZZ& s, vec_ZZX& w, const zz_pX h, const vec_zz_pX a, const mat_L& B, const RR sigma, const zz_pX d);

ZZ   ZSampler(const RR sigma, const RR c);
ZZX  polySampler(const RR sigma);

#endif
