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

#include "params.h"
#include "LHC.h"
#include "Hash.h"
#include "Squares.h"

typedef struct
{
    long            valid = 0;
    vec_zz_pX       t_A, t_y, t_g, w;
    LHC_COM_t       com_1, com_2;
    vec_zz_p        z_3;
    vec_zz_pX       h;
    zz_pX           t, f0;
    vec_zz_pX       z_1,  z_2;
    LHC_OP_t        op_1, op_2;
} PROOF_C_t;

#include "Issuer.h"


void Preprocessing_Com(vec_zz_p& s, const vec_ZZ& s0, const ZZ& B_goth2);
void Prove_Com(uint8_t** Pi_ptr, const unsigned char* seed_crs, const CRS_t& crs, const uint8_t* seed_ipk, const mat_zz_p& P, const vec_zz_p& u0, const ZZ& B_goth2, const vec_ZZ& w0);
long Verify_Com(const unsigned char* seed_crs, const CRS_t& crs, const uint8_t* seed_ipk, const mat_zz_p& P, const vec_zz_p& u0, const ZZ& B_goth2, uint8_t** Pi_ptr);

#endif