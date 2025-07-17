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

#include "params.h"


// compute size functions
size_t calc_ser_size(const long r, const long c);
size_t calc_ser_size_minbyte(const long r, const long c, int nbits);
size_t calc_ser_size_minbits(const long r, const long c, int nbits);
size_t calc_ser_size_vec_poly(const long n, const long d);
size_t calc_ser_size_vec_ZZX(const long n, const long d);
size_t calc_ser_size_poly(long d);
size_t calc_ser_size_vec_zz_p(long l);
size_t calc_ser_size_vec_UL(long l);
size_t calc_ser_size_vec_ZZ(long l);
size_t calc_ser_size_ZZ(void);
size_t calc_ser_size_big_ZZ(long nbits);
size_t calc_ser_size_vec_poly_minbyte(const long n, const long d, const int nbits);
size_t calc_ser_size_poly_minbyte(const long d, const int nbits);
size_t calc_ser_size_vec_zz_p_minbyte(const long l, const int nbits);
size_t calc_ser_size_vec_ZZ_minbyte(const long l, const int nbits);
size_t calc_ser_size_vec_poly_minbits(const long n, const long d, const int nbits);
size_t calc_ser_size_poly_minbits(const long d, const int nbits);
size_t calc_ser_size_vec_zz_p_minbits(const long l, const int nbits);



// serialize functions
void serialize_mat_zz_p(uint8_t* v, const size_t s, const long r, const long c, const mat_zz_p& m);
void serialize_minbyte_mat_zz_p(uint8_t* v, const size_t s, const long r, const long c, const int nbits, const mat_zz_p& m);
void serialize_minbits_mat_zz_p(uint8_t* v, const size_t s, const long r, const long c, const int nbits, const mat_zz_p& m);
void serialize_vec_poly_zz_pX(uint8_t* v, const size_t s, const long n, const long d, const vec_zz_pX& p);
void serialize_vec_ZZX(uint8_t* v, const size_t s, const long n, const long d, const vec_ZZX& p);
void serialize_poly_zz_pX(uint8_t* v, const size_t s, const long d, const zz_pX& p);
void serialize_vec_zz_p(uint8_t* v, const size_t s, const long d, const vec_zz_p& p);
void serialize_vec_UL(uint8_t* v, const size_t s, const long d, const vec_UL& p);
void serialize_vec_ZZ(uint8_t* v, const size_t s, const long d, const vec_ZZ& p);
void serialize_ZZ(uint8_t* v, const size_t s, const ZZ& p);
void serialize_big_ZZ(uint8_t* v, const size_t s, const ZZ& p);
void serialize_minbyte_vec_poly_zz_pX(uint8_t* v, const size_t s, const long n, const long d, const int nbits, const vec_zz_pX& p);
void serialize_minbyte_poly_zz_pX(uint8_t* v, const size_t s, const long d, const int nbits, const zz_pX& p);
void serialize_minbyte_vec_zz_p(uint8_t* v, const size_t s, const long l, const int nbits, const vec_zz_p& p);
void serialize_minbyte_vec_ZZ(uint8_t* v, const size_t s, const long l, const int nbits, const vec_ZZ& p);
void serialize_minbits_vec_poly_zz_pX(uint8_t* v, const size_t s, const long n, const long d, const int nbits, const vec_zz_pX& p);
void serialize_minbits_poly_zz_pX(uint8_t* v, const size_t s, const long d, const int nbits, const zz_pX& p);
void serialize_minbits_vec_zz_p(uint8_t* v, const size_t s, const long l, const int nbits, const vec_zz_p& p);


// deserialize functions
void deserialize_mat_zz_p(mat_zz_p& m, const long r, const long c, const uint8_t* v, const size_t s);
void deserialize_minbyte_mat_zz_p(mat_zz_p& m, const long r, const long c, const int nbits, const uint8_t* v, const size_t s);
void deserialize_minbits_mat_zz_p(mat_zz_p& m, const long r, const long c, const int nbits, const uint8_t* v, const size_t s);
void deserialize_vec_poly_zz_pX(vec_zz_pX& p, const long n, const long d, const uint8_t* v, const size_t s);
void deserialize_vec_ZZX(vec_ZZX& p, const long n, const long d, const uint8_t* v, const size_t s);
void deserialize_poly_zz_pX(zz_pX& p, const long d, const uint8_t* v, const size_t s);
void deserialize_vec_zz_p(vec_zz_p& p, const long d, const uint8_t* v, const size_t s);
void deserialize_vec_UL(vec_UL& p, const long d, const uint8_t* v, const size_t s);
void deserialize_vec_ZZ(vec_ZZ& p, const long d, const uint8_t* v, const size_t s);
void deserialize_ZZ(ZZ& p, const uint8_t* v, const size_t s);
void deserialize_big_ZZ(ZZ& p, const uint8_t* v, const size_t s);
void deserialize_minbyte_vec_poly_zz_pX(vec_zz_pX& p, const long n, const long d, const int nbits, const uint8_t* v, const size_t s);
void deserialize_minbyte_poly_zz_pX(zz_pX& p, const long d, const int nbits, const uint8_t* v, const size_t s);
void deserialize_minbyte_vec_zz_p(vec_zz_p& p, const long l, const int nbits, const uint8_t* v, const size_t s);
void deserialize_minbyte_vec_ZZ(vec_ZZ& p, const long l, const int nbits, const uint8_t* v, const size_t s);
void deserialize_minbits_vec_poly_zz_pX(vec_zz_pX& p, const long n, const long d, const int nbits, const uint8_t* v, const size_t s);
void deserialize_minbits_poly_zz_pX(zz_pX& p, const long d, const int nbits, const uint8_t* v, const size_t s);
void deserialize_minbits_vec_zz_p(vec_zz_p& p, const long l, const int nbits, const uint8_t* v, const size_t s);
