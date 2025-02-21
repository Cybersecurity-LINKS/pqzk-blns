// Copyright 2025 Fondazione LINKS

// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at

//     http://www.apache.org/licenses/LICENSE-2.0

// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef BLNS_LHC_H
#define BLNS_LHC_H

#include "params.h"
#include "Lattice.h"
#include "Utils.h"


void LHC_Com(Vec<vec_zz_pX>& com, Vec<vec_ZZX>& st, const long& index, const Mat<zz_pX>& A_i, const Mat<zz_pX>& B_i, const vec_ZZX& s, const vec_ZZX& y) ;

long Rej(const long& index, const vec_ZZX& z, const vec_ZZX& v, const RR& s, const RR& M);
long Rej_v_ZZ(const vec_ZZ& z, const vec_ZZ& v, const RR& s, const RR& M);
long Rej_v_zzp(const vec_zz_p& z, const vec_zz_p& v, const long& q, const RR& s, const RR& M);
long Rej_v_ZZX(const vec_ZZX& z, const vec_ZZX& v, const RR& s, const RR& M);
long Rej_v_zzpX(const vec_zz_pX& z, const vec_zz_pX& v, const long& q, const RR& s, const RR& M);

void LHC_Open(Vec<vec_ZZX>& op, const long& index, const ZZX& c, const Vec<vec_ZZX>& st);

long LHC_Verify(const long& index, const Mat<zz_pX>& A_i, const Mat<zz_pX>& B_i, const Vec<vec_zz_pX>& com, const ZZX& c, const vec_ZZX& z, const Vec<vec_ZZX>& op);

#endif