#ifndef BLNS_HASH_H
#define BLNS_HASH_H

#include <string.h>
#include <openssl/evp.h>

#include "params.h"


typedef      Vec<Mat<ZZ_pX>>   CRS_Data;
typedef  Vec<Vec<Mat<ZZ_pX>>>  CRS_Data2;


long int            CustomHash(const long int x, const size_t out_len);

EVP_MD_CTX*         Hash_Init(const string inputStr);
ZZ_pX               Hash_ZZ_pX(EVP_MD_CTX *mdctx, const long n_coeffs, const size_t b_coeffs);
vec_ZZ_p            Hash_v_ZZp(EVP_MD_CTX *mdctx, const long n_elems, const size_t b_num);
vec_ZZ              Hash_bits(EVP_MD_CTX *mdctx, const long n_elems);
ZZ                  Hash_ZZ_xi0(EVP_MD_CTX *mdctx, const size_t b_num);

CRS_Data2           Hcrs(const string inputStr);

Vec<Mat<vec_ZZ>>    HCom1(const string inputStr);
mat_ZZ_p            HCom2(const string inputStr);
vec_ZZ_pX           HCom3(const string inputStr);
ZZX                 HCom4(const string inputStr);

Vec<Mat<vec_ZZ>>    HISIS1(const string inputStr);
mat_ZZ_p            HISIS2(const string inputStr);
vec_ZZ_pX           HISIS3(const string inputStr);
ZZX                 HISIS4(const string inputStr);

vec_ZZ              HM(const string a_i);

#endif