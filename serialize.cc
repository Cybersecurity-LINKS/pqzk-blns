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


// compute size functions (size in bytes of vector v)
size_t calc_ser_size(const long r, const long c) {
    return ((r*c)*sizeof(long));
}

size_t calc_ser_size_minbyte(const long r, const long c, int nbits) {
    int b = static_cast<int>(ceil(nbits / 8.0));
    return (r*c*b);
}

size_t calc_ser_size_minbits(const long r, const long c, int nbits) {
    return static_cast<int>(ceil((r*c*nbits) / 8.0));;
}

size_t calc_ser_size_vec_poly(const long n, const long d) {
    return (n*d*sizeof(long));
}

size_t calc_ser_size_vec_ZZX(const long n, const long d) {
    return (n*d*sizeof(long));
}

size_t calc_ser_size_poly(long d) {
    return (d*sizeof(long));
}

size_t calc_ser_size_vec_zz_p(long l) {
    return (l*sizeof(long));
}

size_t calc_ser_size_vec_UL(long l) {
    return (l*sizeof(uint8_t));
}

size_t calc_ser_size_vec_ZZ(long l) {
    return(l*sizeof(long));
}

size_t calc_ser_size_ZZ(void) {
    return(sizeof(long));
}

size_t calc_ser_size_big_ZZ(long nbits) {
    return(static_cast<int>(ceil(nbits / 8.0)));
}

size_t calc_ser_size_vec_poly_minbyte(const long n, const long d, const int nbits) {
    int b = static_cast<int>(ceil(nbits / 8.0));
    return (n*d*b);
}

size_t calc_ser_size_poly_minbyte(const long d, const int nbits) {
    int b = static_cast<int>(ceil(nbits / 8.0));
    return (b*d);
}

size_t calc_ser_size_vec_zz_p_minbyte(const long l, const int nbits) {
    int b = static_cast<int>(ceil(nbits / 8.0));
    return (b*l);
}

size_t calc_ser_size_vec_ZZ_minbyte(const long l, const int nbits) {
    int b = static_cast<int>(ceil(nbits / 8.0));
    return (b*l);
}

size_t calc_ser_size_vec_poly_minbits(const long n, const long d, const int nbits) {
    return static_cast<int>(ceil((n*d*nbits) / 8.0));
}

size_t calc_ser_size_poly_minbits(const long d, const int nbits) {
    return static_cast<int>(ceil((d*nbits) / 8.0));
}

size_t calc_ser_size_vec_zz_p_minbits(const long l, const int nbits) {
    return static_cast<int>(ceil((l*nbits) / 8.0));
}


// serialize/deserialize functions for mat_zz_p, 64 bits per element (i.e. 1 long int) 
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
    m.SetDims(r, c);
    for (i = 0; i < r; i++) {
        for (j = 0; j < c; j++) {
            m[i][j] = data_ptr[i * c + j];
        }
    }
}


// serialize/deserialize functions for mat_zz_p, minimum number of bytes per element
void serialize_minbyte_mat_zz_p(uint8_t* v, const size_t s, const long r, const long c, const int nbits, const mat_zz_p& m) {
    long i, j, elem;
    int k, blk;
    size_t n;
    n = 0;
    blk = (nbits + 7) / 8;
    for(i = 0; i < r; i++) {
        for(j = 0; j < c; j++) {
            elem = conv<long>(m[i][j]);
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
    m.SetDims(r, c);
    for(i = 0; i < r; ++i) {
        for(j = 0; j < c; ++j) {
            elem = 0;
            for(k = 0; k < blk; ++k) {
                if (n < s) {
                    elem |= (static_cast<long>(v[n++])) << (8 * k);
                }
            }
            m[i][j] = elem;
        }
    }

}


// serialize/deserialize functions for mat_zz_p, minimum number of bits per element
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
    m.SetDims(r, c);
    for (i = 0; i < r; i++) {
        for (j = 0; j < c; j++) {
            m[i][j] = vf[k++];
        }
    }
    delete[] vf;
}


// serialize/deserialize functions for vec_zz_pX, 64 bits per coefficient (i.e. 1 long int) 
void serialize_vec_poly_zz_pX(uint8_t* v, const size_t s, const long n, const long d, const vec_zz_pX& p) {
    long i, j;
    long* data_ptr = reinterpret_cast<long*>(v);
    for (i = 0; i < n; i++) {
        for (j = 0; j < d; j++) {
            // data_ptr[i * d + j] = conv<long>(p[i][j]);
            data_ptr[i * d + j] = conv<long>(coeff(p[i], j));
        }
    }
}

void deserialize_vec_poly_zz_pX(vec_zz_pX& p, const long n, const long d, const uint8_t* v, const size_t s) {
    long i, j;
    const long* data_ptr = reinterpret_cast<const long*>(v);
    p.SetLength(n);
    for (i = 0; i < n; i++) {        
        p[i].SetLength(d);
        for (j = 0; j < d; j++) {
            p[i][j] = data_ptr[i * d + j];
            // SetCoeff(p[i], j, data_ptr[i * d + j]);
        }
        p[i].normalize();
    }
}


// serialize/deserialize functions for vec_ZZX, 64 bits per coefficient (i.e. 1 long int) 
void serialize_vec_ZZX(uint8_t* v, const size_t s, const long n, const long d, const vec_ZZX& p) {
    long i, j;
    long* data_ptr = reinterpret_cast<long*>(v);
    for (i = 0; i < n; i++) {
        for (j = 0; j < d; j++) {
            // data_ptr[i * d + j] = conv<long>(p[i][j]);
            data_ptr[i * d + j] = conv<long>(coeff(p[i], j));
        }
    }
}

void deserialize_vec_ZZX(vec_ZZX& p, const long n, const long d, const uint8_t* v, const size_t s) {
    long i, j;
    const long* data_ptr = reinterpret_cast<const long*>(v);
    p.SetLength(n);
    for (i = 0; i < n; i++) {        
        p[i].SetLength(d);
        for (j = 0; j < d; j++) {
            p[i][j] = data_ptr[i * d + j];
            // SetCoeff(p[i], j, data_ptr[i * d + j]);
        }
        p[i].normalize();
    }
}


// serialize/deserialize functions for zz_pX, 64 bits per coefficient (i.e. 1 long int) 
void serialize_poly_zz_pX(uint8_t* v, const size_t s, const long d, const zz_pX& p) {
    long i;
    long* data_ptr = reinterpret_cast<long*>(v);
    for (i = 0; i < d; i++) {
        data_ptr[i] = conv<long>(coeff(p, i));
    }
}

void deserialize_poly_zz_pX(zz_pX& p, const long d, const uint8_t* v, const size_t s) {
    long i;
    const long* data_ptr = reinterpret_cast<const long*>(v);
    p.SetLength(d);
    for (i = 0; i < d; i++) {
        p[i] = data_ptr[i];
    }
    p.normalize();   
}


// serialize/deserialize functions for vec_zz_p, 64 bits per element (i.e. 1 long int) 
void serialize_vec_zz_p(uint8_t* v, const size_t s, const long d, const vec_zz_p& p) {
    long i;
    long* data_ptr = reinterpret_cast<long*>(v);
    for (i = 0; i < d; i++) {
        data_ptr[i] = conv<long>(p[i]);
    }
}

void deserialize_vec_zz_p(vec_zz_p& p, const long d, const uint8_t* v, const size_t s) {
    long i;
    const long* data_ptr = reinterpret_cast<const long*>(v);
    p.SetLength(d);
    for (i = 0; i < d; i++) {
        p[i] = data_ptr[i];
    }
}


// serialize/deserialize functions for vec_UL, 8 bits per element (i.e. 1 uint8) 
void serialize_vec_UL(uint8_t* v, const size_t s, const long d, const vec_UL& p) {
    long i;
    for (i = 0; i < d; i++) {
        v[i] = (uint8_t)(p[i]);
    }
}

void deserialize_vec_UL(vec_UL& p, const long d, const uint8_t* v, const size_t s) {
    long i;
    p.SetLength(d);
    for (i = 0; i < d; i++) {
        p[i] = (long)(v[i]);
    }
}


// serialize/deserialize functions for vec_ZZ, 64 bits per element (i.e. 1 long int)
void serialize_vec_ZZ(uint8_t* v, const size_t s, const long d, const vec_ZZ& p) {
    long i;
    long* data_ptr = reinterpret_cast<long*>(v);
    for (i = 0; i < d; i++) {
        data_ptr[i] = conv<long>(p[i]);
    }
}

void deserialize_vec_ZZ(vec_ZZ& p, const long d, const uint8_t* v, const size_t s) {
    long i;
    const long* data_ptr = reinterpret_cast<const long*>(v);
    p.SetLength(d);
    for (i = 0; i < d; i++) {
        p[i] = data_ptr[i];
    }
}


// serialize/deserialize functions for "small" ZZ, max. 64 bits (i.e. 1 long int)
void serialize_ZZ(uint8_t* v, const size_t s, const ZZ& p) {
    long* data_ptr = reinterpret_cast<long*>(v);
    data_ptr[0] = conv<long>(p);
}

void deserialize_ZZ(ZZ& p, const uint8_t* v, const size_t s) {
    const long* data_ptr = reinterpret_cast<const long*>(v);
    p = data_ptr[0];
}


// serialize/deserialize functions for "big" ZZ, i.e. integers with more than 64 bits
void serialize_big_ZZ(uint8_t* v, const size_t s, const ZZ& p) {
    BytesFromZZ(reinterpret_cast<unsigned char*>(v), p, s);
}

void deserialize_big_ZZ(ZZ& p, const uint8_t* v, const size_t s) {
    ZZFromBytes(p, reinterpret_cast<const unsigned char*>(v), s);
}


// serialize/deserialize functions for vec_zz_pX, minimum number of bytes per element
void serialize_minbyte_vec_poly_zz_pX(uint8_t* v, const size_t s, const long n, const long d, const int nbits, const vec_zz_pX& p) {
    long i, j, c;
    int k, blk;
    size_t o;
    
    o = 0;
    blk = (nbits + 7) / 8;
    for(i = 0; i < n; i++) {
        for(j = 0; j < d; j++) {
            c = conv<long>(coeff(p[i], j));
            for(k = 0; k < blk; ++k) {
                v[o++] = static_cast<uint8_t>((c >> (8 * k)) & 0xFF);
            }
        }
    }
}

void deserialize_minbyte_vec_poly_zz_pX(vec_zz_pX& p, const long n, const long d, const int nbits, const uint8_t* v, const size_t s) {
    long i, j, elem;
    int k, blk;
    size_t o;

    o = 0;
    blk = (nbits + 7) / 8;
    p.SetLength(n);
    for(i = 0; i < n; ++i) { 
        p[i].SetLength(d);       
        for(j = 0; j < d; ++j) {
            elem = 0;
            for(k = 0; k < blk; ++k) {
                if (o < s) {
                    elem |= (static_cast<long>(v[o++])) << (8 * k);
                }
            }
            p[i][j] = elem;
        }
        p[i].normalize();
    }
}


// serialize/deserialize functions for zz_pX, minimum number of bytes per element
void serialize_minbyte_poly_zz_pX(uint8_t* v, const size_t s, const long d, const int nbits, const zz_pX& p) {
    long i, c;
    int k, blk;
    size_t o;
    
    o = 0;
    blk = (nbits + 7) / 8;
    for(i = 0; i < d; i++) {
        c = conv<long>(coeff(p, i));
        for(k = 0; k < blk; ++k) {
            v[o++] = static_cast<uint8_t>((c >> (8 * k)) & 0xFF);
        }
    }
}

void deserialize_minbyte_poly_zz_pX(zz_pX& p, const long d, const int nbits, const uint8_t* v, const size_t s) {
    long i, elem;
    int k, blk;
    size_t o;

    o = 0;
    blk = (nbits + 7) / 8;
    p.SetLength(d);
    for(i = 0; i < d; ++i) {
        elem = 0;
        for(k = 0; k < blk; ++k) {
            if (o < s) {
                elem |= (static_cast<long>(v[o++])) << (8 * k);
            }
        }
        p[i] = elem;
    }
    p.normalize();
}


// serialize/deserialize functions for vec_zz_p, minimum number of bytes per element
void serialize_minbyte_vec_zz_p(uint8_t* v, const size_t s, const long l, const int nbits, const vec_zz_p& p) {
    long i, elem;
    int k, blk;
    size_t n;
    
    n = 0;
    blk = (nbits + 7) / 8;
    for(i = 0; i < l; i++) {
        elem = conv<long>(p[i]);
        for(k = 0; k < blk; ++k) {
            v[n++] = static_cast<uint8_t>((elem >> (8 * k)) & 0xFF);
        }
    }
}

void deserialize_minbyte_vec_zz_p(vec_zz_p& p, const long l, const int nbits, const uint8_t* v, const size_t s) {
    long i, elem;
    int k, blk;
    size_t n;

    n = 0;
    blk = (nbits + 7) / 8;
    p.SetLength(l);
    for(i = 0; i < l; ++i) {
        elem = 0;
        for(k = 0; k < blk; ++k) {
            if (n < s) {
                elem |= (static_cast<long>(v[n++])) << (8 * k);
            }
        }
        p[i] = elem;
    }
}


// serialize/deserialize functions for vec_ZZ, minimum number of bytes per element
void serialize_minbyte_vec_ZZ(uint8_t* v, const size_t s, const long l, const int nbits, const vec_ZZ& p) {
    long i, elem;
    int k, blk;
    size_t n;
    
    n = 0;
    blk = (nbits + 7) / 8;
    for(i = 0; i < l; i++) {
        elem = conv<long>(p[i]);
        for(k = 0; k < blk; ++k) {
            v[n++] = static_cast<uint8_t>((elem >> (8 * k)) & 0xFF);
        }
    }
}

void deserialize_minbyte_vec_ZZ(vec_ZZ& p, const long l, const int nbits, const uint8_t* v, const size_t s) {
    long i, elem;
    int k, blk;
    size_t n;

    n = 0;
    blk = (nbits + 7) / 8;
    p.SetLength(l);
    for(i = 0; i < l; ++i) {
        elem = 0;
        for(k = 0; k < blk; ++k) {
            if (n < s) {
                elem |= (static_cast<long>(v[n++])) << (8 * k);
            }
        }
        p[i] = elem;
    }
}


// serialize/deserialize functions for vec_zz_pX, minimum number of bits per element
void serialize_minbits_vec_poly_zz_pX(uint8_t* v, const size_t s, const long n, const long d, const int nbits, const vec_zz_pX& p) {
    size_t l, b, blk, blk_b;
    size_t cnt = static_cast<size_t>(n) * static_cast<size_t>(d);
    long* vf = new long[cnt];
    long i, j, k;
    uint64_t e;
    int m;
    k = 0;
    for (i = 0; i < n; i++) {
        for (j = 0; j < d; j++) {
            vf[k++] = conv<long>(coeff(p[i], j));
        }
    }
    for (l = 0; l < cnt; ++l) {
        b = l * nbits; blk = b / 8; blk_b = b % 8;
        e = static_cast<uint64_t>(vf[l]) & ((1ULL << nbits) - 1);
        for (m = 0; m < nbits; ++m) {
            if (e & (1ULL << m)) {
                v[blk + (blk_b + m) / 8] |= (1 << ((blk_b + m) % 8));
            }
        }
    }
    delete[] vf;
}

void deserialize_minbits_vec_poly_zz_pX(vec_zz_pX& p, const long n, const long d, const int nbits, const uint8_t* v, const size_t s) {
    size_t l, b, blk, blk_b;
    size_t cnt = static_cast<size_t>(n) * static_cast<size_t>(d);
    long* vf = new long[cnt];
    long i, j, k;
    uint64_t e;
    int m;
    for (l = 0; l < cnt; ++l) {
        b = l * nbits; blk = b / 8; blk_b = b % 8;
        e = 0;
        for (m = 0; m < nbits; ++m) {
            if (v[blk + (blk_b + m) / 8] & (1 << ((blk_b + m) % 8))) {
                e |= (1ULL << m);
            }
        }
        vf[l] = static_cast<long>(e);
    }
    k = 0;
    p.SetLength(n);
    for (i = 0; i < n; i++) {
        p[i].SetLength(d);
        for (j = 0; j < d; j++) {
            p[i][j] = vf[k++];
        }
        p[i].normalize();
    }
    delete[] vf;
}


// serialize/deserialize functions for zz_pX, minimum number of bits per element
void serialize_minbits_poly_zz_pX(uint8_t* v, const size_t s, const long d, const int nbits, const zz_pX& p) {
    size_t l, b, blk, blk_b;
    size_t cnt = static_cast<size_t>(d);
    long* vf = new long[cnt];
    long i, k;
    uint64_t e;
    int m;
    k = 0;
    for (i = 0; i < d; i++) {
        vf[k++] = conv<long>(coeff(p, i));
    }
    for (l = 0; l < cnt; ++l) {
        b = l * nbits; blk = b / 8; blk_b = b % 8;
        e = static_cast<uint64_t>(vf[l]) & ((1ULL << nbits) - 1);
        for (m = 0; m < nbits; ++m) {
            if (e & (1ULL << m)) {
                v[blk + (blk_b + m) / 8] |= (1 << ((blk_b + m) % 8));
            }
        }
    }
    delete[] vf;
}

void deserialize_minbits_poly_zz_pX(zz_pX& p, const long d, const int nbits, const uint8_t* v, const size_t s) {
    size_t l, b, blk, blk_b;
    size_t cnt = static_cast<size_t>(d);
    long* vf = new long[d];
    long i, k;
    uint64_t e;
    int m;
    for (l = 0; l < cnt; ++l) {
        b = l * nbits; blk = b / 8; blk_b = b % 8;
        e = 0;
        for (m = 0; m < nbits; ++m) {
            if (v[blk + (blk_b + m) / 8] & (1 << ((blk_b + m) % 8))) {
                e |= (1ULL << m);
            }
        }
        vf[l] = static_cast<long>(e);
    }
    k = 0;
    p.SetLength(d);
    for (i = 0; i < d; i++) {
        p[i] = vf[k++];
    }
    p.normalize();
    delete[] vf;
}


// serialize/deserialize functions for vec_zz_p, minimum number of bits per element
void serialize_minbits_vec_zz_p(uint8_t* v, const size_t s, const long l, const int nbits, const vec_zz_p& p) {
    size_t c, b, blk, blk_b;
    size_t cnt = static_cast<size_t>(l);
    long* vf = new long[cnt];
    long i, k;
    uint64_t e;
    int m;
    k = 0;
    for (i = 0; i < l; i++) {
        vf[k++] = conv<long>(p[i]);
    }
    for (c = 0; c < cnt; ++c) {
        b = c * nbits; blk = b / 8; blk_b = b % 8;
        e = static_cast<uint64_t>(vf[c]) & ((1ULL << nbits) - 1);
        for (m = 0; m < nbits; ++m) {
            if (e & (1ULL << m)) {
                v[blk + (blk_b + m) / 8] |= (1 << ((blk_b + m) % 8));
            }
        }
    }
    delete[] vf;
}

void deserialize_minbits_vec_zz_p(vec_zz_p& p, const long l, const int nbits, const uint8_t* v, const size_t s) {
    size_t c, b, blk, blk_b;
    size_t cnt = static_cast<size_t>(l);
    long* vf = new long[cnt];
    long i, k;
    uint64_t e;
    int m;
    for (c = 0; c < cnt; ++c) {
        b = c * nbits; blk = b / 8; blk_b = b % 8;
        e = 0;
        for (m = 0; m < nbits; ++m) {
            if (v[blk + (blk_b + m) / 8] & (1 << ((blk_b + m) % 8))) {
                e |= (1ULL << m);
            }
        }
        vf[c] = static_cast<long>(e);
    }
    k = 0;
    p.SetLength(l);
    for (i = 0; i < l; i++) {
        p[i] = vf[k++];
    }
    delete[] vf;
}

