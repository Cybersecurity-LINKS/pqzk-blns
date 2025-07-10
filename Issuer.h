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

#ifndef BLNS_ISSUER_H
#define BLNS_ISSUER_H

#include "params.h"
#include "Lattice.h"
#include "Com.h"
#include "serialize.h"
#include "Utils.h"

typedef struct
{
    uint8_t*    u;
    uint8_t*    Pi;
} RHO1_t;


void    I_KeyGen(uint8_t** ipk_prt, ISK_t& isk);
void    CompleteIPK(IPK_t& ipk, const uint8_t* ipk_bytes);

void    I_VerCred(uint8_t** Rho2_ptr, const uint8_t* seed_crs, const CRS2_t& crs, const mat_zz_p& B_f, const uint8_t* ipk_bytes, const ISK_t& isk, const Vec<string>& attrs_prime, const vec_UL &idx_pub, RHO1_t& Rho1);

void    I_VerCred_Plain(uint8_t** Rho_ptr, const mat_zz_p& B_f, const uint8_t* ipk_bytes, const ISK_t& isk, Vec<string>& attrs);

#ifdef USE_REVOCATION
    void I_UpdateSign(uint8_t** Rho2_ptr, const mat_zz_p& B_f, const uint8_t* ipk_bytes, const ISK_t& isk, const uint8_t* u_bytes, const string& old_timestamp, const string& new_timestamp);
#endif

#endif