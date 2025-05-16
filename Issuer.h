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

typedef struct
{
    zz_pX       a1;
    vec_zz_pX   a2;
    vec_zz_pX   c0;
    vec_zz_pX   c1;
} IPK_t;

typedef struct
{
    uint8_t*    u;
    uint8_t*    Pi;
} RHO1_t;

#include "Com.h"

 
void    I_KeyGen(IPK_t& ipk, ZZX& f, ZZX& g, ZZX& F, ZZX& G);
void    I_VerCred(uint8_t** Rho2_ptr, const unsigned char* seed_crs, const CRS2_t& crs, const mat_zz_p& B_f, const IPK_t& ipk, const ZZX& f, const ZZX& g, const ZZX& F, const ZZX& G, const Vec<string>& attrs_prime, RHO1_t& Rho1);

#endif