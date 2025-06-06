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

#ifndef BLNS_RANDOM_H
#define BLNS_RANDOM_H

#include "params.h"
#include "Utils.h"

typedef struct
{
    zz_pX       a1;
    vec_zz_pX   a2;
    vec_zz_pX   c0;
    vec_zz_pX   c1;
    uint8_t     seed_ipk[SEED_LEN];
} IPK_t; // Issuer Public Key

typedef struct
{
    ZZX         f;
    ZZX         g;
    ZZX         F;
    ZZX         G;
} ISK_t; // Issuer Secret Key


#ifdef ENABLE_FALCON
void Falcon_keygen(zz_pX& a1, ISK_t& isk);
void Falcon_GSampler(vec_ZZ& s, vec_ZZX& w, const zz_pX& h, const vec_zz_pX& a, const ISK_t& isk, const zz_pX& d);
#endif

void NTRU_TrapGen(zz_pX& a1, ISK_t& isk);

void preGSampler(vec_ZZ& v, const mat_L& B, const vec_ZZ& c);
void GSampler(vec_ZZ& s, vec_ZZX& w, const zz_pX& h, const vec_zz_pX& a, const mat_L& B, const zz_pX& d);

void ZSampler(ZZ& x, const double& sigma, const double& c);
void polySampler(ZZX& s, const double& sigma);
void polySampler_hat(zz_pX& s, const double& sigma);

#endif
