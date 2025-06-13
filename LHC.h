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

typedef struct
{
    vec_zz_pX   t_1, t_2, w_1, w_2;
} LHC_COM_t;

typedef struct
{
    vec_zz_pX   e_1, e_2, e_3, f_1, f_2, f_3;
} LHC_ST_t;

typedef struct
{
    vec_zz_pX   z_1, z_2, z_3;
    long        valid = 0;
} LHC_OP_t;


void LHC_Com(LHC_COM_t& com, LHC_ST_t& st, const long& index, const Mat<zz_pX>& A_i, const Mat<zz_pX>& B_i, const vec_zz_pX& s, const vec_zz_pX& y, const long& m);

long Rej_v_zzp(const vec_zz_p& z, const vec_zz_p& v, const long& q, const RR& s, const RR& M);
long Rej_v_zzpX(const vec_zz_pX& z, const vec_zz_pX& v, const long& q, const RR& s, const RR& M);

void LHC_Open(LHC_OP_t& op, const long& index, const zz_pX& c, const LHC_ST_t& st, const long& m);

long LHC_Verify(const long& index, const Mat<zz_pX>& A_i, const Mat<zz_pX>& B_i, const LHC_COM_t& com, const zz_pX& c, const vec_zz_pX& z, const LHC_OP_t& op, const long& m);

#endif