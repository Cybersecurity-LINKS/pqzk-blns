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

#ifdef ENABLE_FALCON
void Falcon_keygen(zz_pX& a1, ZZX& f, ZZX& g, ZZX& F, ZZX& G);
#endif

void NTRU_TrapGen(zz_pX& a1, ZZX& f, ZZX& g, ZZX& F, ZZX& G);

void preGSampler(vec_ZZ& v, const mat_L& B, const double& sigma, const vec_ZZ& c);
void GSampler(vec_ZZ& s, vec_ZZX& w, const zz_pX& h, const vec_zz_pX& a, const mat_L& B, const double& sigma, const zz_pX& d);

void ZSampler(ZZ& x, const double& sigma, const double& c);
void polySampler(ZZX& s, const double& sigma);
void polySampler_hat(zz_pX& s, const double& sigma);

#endif
