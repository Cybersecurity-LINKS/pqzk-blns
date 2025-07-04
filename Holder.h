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
    long        valid = 0;
} CRED_t;

typedef struct
{
    // vec_ZZ   cp;
    const uint8_t* ipk_bytes;
    Vec<string> attrs_prime;
    vec_ZZ      idx;
    uint8_t*    Pi;
    long        valid = 0;
} VP_t;


void    H_Init(CRS2_t& crs, mat_zz_p& B_f, uint8_t* seed_crs, Vec<string>& attrs, const long &num_idx_hid);

void    H_VerCred1(RHO1_t& Rho1, STATE_t& state, const uint8_t* seed_crs, const CRS2_t& crs, const uint8_t* ipk_bytes, const Vec<string>& attrs, const vec_UL &idx_pub);
void    H_VerCred2(CRED_t& cred, const uint8_t* ipk_bytes, const mat_zz_p& B_f, uint8_t** Rho2_ptr, const STATE_t& state);
void    H_VerPres(VP_t& VP, const CRED_t& cred, const uint8_t* seed_crs, const CRS2_t& crs, const uint8_t* ipk_bytes, const mat_zz_p& B_f, Vec<string>& attrs, const vec_UL &idx_pub);

void    H_VerCred_Plain(CRED_t& cred, const uint8_t* ipk_bytes, const mat_zz_p& B_f, uint8_t** Rho_ptr, Vec<string>& attrs);

#endif