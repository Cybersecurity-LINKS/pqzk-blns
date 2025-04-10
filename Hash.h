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

#ifndef BLNS_HASH_H
#define BLNS_HASH_H

#include <string.h>

#include "params.h"
#include "shake128.h"


typedef      Vec<Mat<zz_pX>>   CRS_t;
typedef  Vec<Vec<Mat<zz_pX>>>  CRS2_t;
typedef  shake128_state_t      HASH_STATE_t; 


HASH_STATE_t* Hash_Init(const uint8_t* v, const size_t len);
void Hash_Update(HASH_STATE_t *state, const uint8_t* v, const size_t len);
HASH_STATE_t* Hash_Copy(const HASH_STATE_t *state0);

void Hash_zz_pX(zz_pX& out_poly, HASH_STATE_t *state, const long& n_coeffs, const size_t& b_coeffs);
void Hash_v_zz_p(vec_zz_p& out_vec, HASH_STATE_t *state, const long& n_elems, const size_t& b_num);
void Hash_R_goth(vec_zz_p& out, HASH_STATE_t *state, const long& n_elems);
void Hash_ZZ_xi0(ZZ& out, HASH_STATE_t *state, const size_t& b_num);

void Hcrs(CRS2_t& crs, mat_zz_p& B_f, const unsigned char* seed_crs);

void HCom1(mat_zz_p& R_goth, const HASH_STATE_t *state0);
void HCom2(mat_zz_p& gamma, const HASH_STATE_t *state0);
void HCom3(vec_zz_pX& mu, const HASH_STATE_t *state0);
void HCom4(zz_pX& c, const HASH_STATE_t *state0);

void HISIS1(mat_zz_p& R_goth, const HASH_STATE_t *state0);
void HISIS2(mat_zz_p& gamma, const HASH_STATE_t *state0);
void HISIS3(vec_zz_pX& mu, const HASH_STATE_t *state0);
void HISIS4(zz_pX& c, const HASH_STATE_t *state0);

void HM(vec_ZZ& m_i, const string& a_i);

#endif