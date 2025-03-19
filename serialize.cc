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

#include "serialize.h"

#include <cmath>

#include <iostream>
#include <sstream> 
#include <fstream>


size_t calc_ser_size(const long r, const long c) {
    return ((r*c)*sizeof(long));
}

size_t calc_ser_size_minbyte(const long r, const long c, int nbits) {
    int b = static_cast<int>(std::ceil(nbits / 8.0));
    return (r*c*b);
}

size_t calc_ser_size_minbits(const long r, const long c, int nbits) {
    return static_cast<int>(std::ceil((r*c*nbits) / 8.0));;
}



void serialize_mat_zz_p(uint8_t* v, const size_t s, const long r, const long c, const mat_zz_p& m) {
    long i, j;
    
    long* data_ptr = reinterpret_cast<long*>(v);
    for (i = 0; i < r; i++) {
        for (j = 0; j < c; j++) {
            data_ptr[i * c + j] = conv<long>(m[i][j]);
        }
    }
}


void deserialize_mat_zz_p(mat_zz_p& m, const long r, const long c, const uint8_t* v, const size_t s) {
    long i, j;
    
    const long* data_ptr = reinterpret_cast<const long*>(v);
    for (long i = 0; i < r; i++) {
        for (long j = 0; j < c; j++) {
            m[i][j] = conv<zz_p>(data_ptr[i * c + j]);
        }
    }
}


void serialize_minbyte_mat_zz_p(uint8_t* v, const size_t s, const long r, const long c, const int nbits, const mat_zz_p& m) {
    long i, j, elem;
    int k, blk;
    size_t n;
    
    n = 0;
    blk = (nbits + 7) / 8;
    for(i = 0; i < r; i++) {
        for(j = 0; j < c; j++) {

            elem = NTL::conv<long>(m[i][j]);
            for(k=0; k < blk; ++k) {
                v[n++] = static_cast<uint8_t>((elem >> (8 * k)) & 0xFF);
            }
        }
    }
}


void deserialize_minbyte_mat_zz_p(mat_zz_p& m, const long r, const long c, const int nbits, const uint8_t* v, const size_t s) {
    long i, j, elem;
    int k, blk;
    size_t n;

    n = 0;
    blk = (nbits + 7) / 8;
    for(i = 0; i < r; ++i) {
        for(j = 0; j < c; ++j) {

            elem = 0;
            for(k = 0; k < blk; ++k) {
                if (n < s) {
                    elem |= (static_cast<long>(v[n++])) << (8 * k);
                }
            }
            m[i][j] = conv<zz_p>(elem);
        }
    }

}


void serialize_minbits_mat_zz_p(uint8_t* v, const size_t s, const long r, const long c, const int nbits, const mat_zz_p& m) {
    size_t l, b, blk, blk_b;
    size_t cnt = static_cast<size_t>(r) * static_cast<size_t>(c);
    long* vf = new long[cnt];
    long i, j, k;

    uint64_t n;
    int p;


    k = 0;
    for (i = 0; i < r; i++) {
        for (j = 0; j < c; j++) {
            vf[k++] = conv<long>(m[i][j]);
        }
    }
    for (l = 0; l < cnt; ++l) {
        b = l * nbits; blk = b / 8; blk_b = b % 8;
        n = static_cast<uint64_t>(vf[l]) & ((1ULL << nbits) - 1);
        for (p = 0; p < nbits; ++p) {
            if (n & (1ULL << p)) {
                v[blk + (blk_b + p) / 8] |= (1 << ((blk_b + p) % 8));
            }
        }
    }
    delete[] vf;
}


void deserialize_minbits_mat_zz_p(mat_zz_p& m, const long r, const long c, const int nbits, const uint8_t* v, const size_t s) {
    size_t l, b, blk, blk_b;
    size_t cnt = static_cast<size_t>(r) * static_cast<size_t>(c);
    long* vf = new long[cnt];
    long i, j, k;
    
    uint64_t n;
    int p;


    for (l = 0; l < cnt; ++l) {
        b = l * nbits; blk = b / 8; blk_b = b % 8;
        n = 0;
        for (p = 0; p < nbits; ++p) {
            if (v[blk + (blk_b + p) / 8] & (1 << ((blk_b + p) % 8))) {
                n |= (1ULL << p);
            }
        }
        vf[l] = static_cast<long>(n);
    }
    k = 0;
    for (i = 0; i < r; i++) {
        for (j = 0; j < c; j++) {
            m[i][j] = conv<zz_p>(vf[k++]);
        }
    }
    delete[] vf;
}
