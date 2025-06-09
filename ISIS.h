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

#ifndef BLNS_ISIS_H
#define BLNS_ISIS_H

#include "params.h"
#include "LHC.h"
#include "Hash.h"
#include "Squares.h"

typedef struct
{
    long            valid = 0;
    vec_zz_pX       t_A, t_y, t_g, w;    
    vec_zz_p        z_3;
    vec_zz_pX       h;
    zz_pX           t, f0;
    vec_zz_pX       z_1, z_2;
} PROOF_I_t;


void Preprocessing_ISIS(vec_zz_p& s, vec_zz_p& r, const vec_ZZ& s0, const vec_ZZ& r0, const ZZ& B_goth_s2, const ZZ& B_goth_r2);
void Prove_ISIS(uint8_t** Pi_ptr, const uint8_t* seed_crs, const CRS_t& crs, const uint8_t* ipk_bytes, const mat_zz_p& P, const mat_zz_p& C, const vec_zz_p& mex, const mat_zz_p& B_f, const vec_ZZ& Bounds, const long& aux, const Vec<vec_ZZ>& w0);
long Verify_ISIS(const uint8_t* seed_crs, const CRS_t& crs, const uint8_t* ipk_bytes, const mat_zz_p& P, const mat_zz_p& C, const vec_zz_p& mex, const mat_zz_p& B_f, const vec_ZZ& Bounds, const long& aux, uint8_t** Pi_ptr);

#endif