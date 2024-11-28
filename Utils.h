#ifndef BLNS_UTILS_H
#define BLNS_UTILS_H

#include "params.h"


ZZX         Phi();
ZZX         Phi_hat();

const ZZX   phi     =      Phi();
const ZZX   phi_hat =      Phi_hat();


void        GS_Ortho(      mat_RR& Bt, vec_RR& Norms2, const mat_L& B );
void        MGS_Ortho(     mat_RR& Bt, vec_RR& Norms2, const mat_L& B );
void        OGS_Ortho(     mat_D&  Bt, vec_D&  Norms2, const mat_L& B ); 

void        rot(           mat_L& M, const ZZX& f ); 
mat_ZZ_p    rot_T(         const ZZ_pX& f );    
mat_ZZ_p    rot_vect(      const vec_ZZ_pX& v ); 

vec_ZZ      Coeffs(        const vec_ZZ_pX x, const unsigned int l );
vec_ZZ      CoeffsX(       const vec_ZZX   x, const unsigned int l );
vec_ZZ_pX   CoeffsInv(     const vec_ZZ    c, const unsigned int l );
vec_ZZX     CoeffsInvX(    const vec_ZZ    c, const unsigned int l );
vec_ZZ      CoeffsHat(     const vec_ZZX   x, const unsigned int l );
vec_ZZ_pX   CoeffsInvHat(  const vec_ZZ_p  c, const unsigned int l );
vec_ZZX     CoeffsInvHatX( const vec_ZZ    c, const unsigned int l );

vec_ZZ_pX   sigma_map(     const vec_ZZ_pX& M, const unsigned int d );

ZZ_pX       poly_mult(     const vec_ZZ_pX& f, const vec_ZZ_pX& g );
ZZ_pX       poly_mult_hat( const vec_ZZ_pX& f, const vec_ZZ_pX& g );

ZZ_pX       Compute_f(     const mat_ZZ_p& B_f, const ZZ& x );

ZZ          Norm2(         const vec_ZZ&  v );
ZZ          Norm2X(        const vec_ZZX& v, const unsigned int d );
RR          Norm2R(        const vec_RR&  v );
double      Norm2D(        const vec_D&   v );

double      InnerProdD(    const vec_D& a, const vec_D& b );

#endif
