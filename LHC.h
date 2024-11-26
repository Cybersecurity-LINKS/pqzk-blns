#ifndef BLNS_LHC_H
#define BLNS_LHC_H

#include "params.h"
#include "Lattice.h"
#include "Utils.h"


void LHC_Com(Vec<vec_ZZ_pX>& com, Vec<vec_ZZX>& st, const int index, const Mat<ZZ_pX>& A_i, const Mat<ZZ_pX>& B_i, const vec_ZZX s, const vec_ZZX y);

int Rej(const int index, const vec_ZZX z, const vec_ZZX v, const RR s, const RR M);
int Rej_v_ZZ(const vec_ZZ z, const vec_ZZ v, const RR s, const RR M);
int Rej_v_ZZX(const vec_ZZX z, const vec_ZZX v, const RR s, const RR M);

Vec<vec_ZZX> LHC_Open(const int index, const ZZX c, const Vec<vec_ZZX>& st);

int LHC_Verify(const int index, const Mat<ZZ_pX>& A_i, const Mat<ZZ_pX>& B_i, const Vec<vec_ZZ_pX>& com, const ZZX c, const vec_ZZX z, const Vec<vec_ZZX>& op);

#endif