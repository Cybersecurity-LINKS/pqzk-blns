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
#include <openssl/evp.h>

#include "params.h"


typedef      Vec<Mat<zz_pX>>   CRS_t;
typedef  Vec<Vec<Mat<zz_pX>>>  CRS2_t;


long     CustomHash(const long  x, const size_t out_len);

EVP_MD_CTX* Hash_Init(const string& inputStr);

void Hash_zz_pX(zz_pX& out_poly, EVP_MD_CTX *mdctx, const long& n_coeffs, const size_t& b_coeffs);
void Hash_v_zz_p(vec_zz_p& out_vec, EVP_MD_CTX *mdctx, const long& n_elems, const size_t& b_num);
void Hash_R_goth(vec_L& out, EVP_MD_CTX *mdctx, const long& n_elems);
void Hash_ZZ_xi0(ZZ& out, EVP_MD_CTX *mdctx, const size_t& b_num);

void Hcrs(CRS2_t& crs, const string& inputStr);

void HCom1(mat_L& R_goth, const string& inputStr);
void HCom2(mat_zz_p& gamma, const string& inputStr);
void HCom3(vec_zz_pX& mu, const string& inputStr);
void HCom4(ZZX& c, const string& inputStr);

void HISIS1(mat_L& R_goth, const string& inputStr);
void HISIS2(mat_zz_p& gamma, const string& inputStr);
void HISIS3(vec_zz_pX& mu, const string& inputStr);
void HISIS4(ZZX& c, const string& inputStr);

void HM(vec_ZZ& m_i, const string& a_i);

#endif