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
    long            valid;
} PROOF_C_t;

#include "Issuer.h"


void Preprocessing_Com(vec_ZZ& s1, const vec_ZZ& s, const ZZ& B_goth2);
void Prove_Com(PROOF_C_t& Pi, const string& inputStr, const CRS_t& crs, const IPK_t& ipk, const mat_ZZ& P0, const vec_ZZ& u0, const ZZ& B_goth2, const vec_ZZ& w0);
long Verify_Com(const string& inputStr, const CRS_t& crs, const IPK_t& ipk, const mat_ZZ& P0, const vec_ZZ& u0, const ZZ& B_goth2, const PROOF_C_t& Pi);

#endif