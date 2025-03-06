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

#ifndef BLNS_UTILS_H
#define BLNS_UTILS_H

#include "params.h"


ZZX         Phi();

ZZX         ModPhi(const ZZX& p);
zz_pX       ModPhi_q(const zz_pX& p);
ZZX         ModPhi_hat(const ZZX& p);
zz_pX       ModPhi_hat_q(const zz_pX& p);

void        OGS_Ortho( mat_D&  Bt, vec_D&  Norms2, const mat_L& B ); 

void        rot(         mat_L& M,   const ZZX& f );
void        rot_T(    mat_zz_p& M, const zz_pX& f );
void        rot_vect( mat_zz_p& R, const vec_zz_pX& v );

void        CoeffsX(vec_ZZ& coeffs_x, const vec_ZZX& x, const unsigned long& l);
void        CoeffsInv(vec_zz_pX& x, const vec_zz_p& c, const unsigned long& l);
void        CoeffsInvX(vec_ZZX& x, const vec_ZZ& c, const unsigned long& l);
void        CoeffsHat(vec_zz_p& coeffs_x, const vec_zz_pX& x, const unsigned long& l);
void        CoeffsInvHat(vec_zz_pX& x, const vec_zz_p& c, const unsigned long& l);

void        sigma_map(vec_zz_pX& N, const vec_zz_pX& M, const unsigned long& d);

zz_pX       poly_mult(     const vec_zz_pX& f, const vec_zz_pX& g );
zz_pX       poly_mult_hat( const vec_zz_pX& f, const vec_zz_pX& g );

zz_pX       Compute_f(     const mat_zz_p& B_f, const ZZ& x );

ZZ          Norm2(         const vec_ZZ&  v );
ZZ          Norm2m(        const vec_zz_p& v, const long& q);
ZZ          Norm2X(        const vec_ZZX& v, const long& d );
ZZ          Norm2Xm(       const vec_zz_pX& v, const long& d, const long& q);
double      Norm2D(        const vec_D&   v );

double      InnerProdD(    const vec_D& a, const vec_D& b );

#endif
