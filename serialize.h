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

#include <NTL/lzz_pX.h>

using namespace std;
using namespace NTL;


// compute size functions
size_t calc_ser_size(const long r, const long c);
size_t calc_ser_size_minbyte(const long r, const long c, int nbits);
size_t calc_ser_size_minbits(const long r, const long c, int nbits);
size_t calc_ser_size_vec_poly(const long n, const long d);

// serialize functions
void serialize_mat_zz_p(uint8_t* v, const size_t s, const long r, const long c, const mat_zz_p& m);
void serialize_minbyte_mat_zz_p(uint8_t* v, const size_t s, const long r, const long c, const int nbits, const mat_zz_p& m);
void serialize_minbits_mat_zz_p(uint8_t* v, const size_t s, const long r, const long c, const int nbits, const mat_zz_p& m);
void serialize_vec_poly_zz_pX(uint8_t* v, const size_t s, const long n, const long d, const vec_zz_pX& p);

// deserialize functions
void deserialize_mat_zz_p(mat_zz_p& m, const long r, const long c, const uint8_t* v, const size_t s);
void deserialize_minbyte_mat_zz_p(mat_zz_p& m, const long r, const long c, const int nbits, const uint8_t* v, const size_t s);
void deserialize_minbits_mat_zz_p(mat_zz_p& m, const long r, const long c, const int nbits, const uint8_t* v, const size_t s);
void deserialize_vec_poly_zz_pX(vec_zz_pX& p, const long n, const long d, const uint8_t* v, const size_t s);

