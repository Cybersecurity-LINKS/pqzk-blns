#ifndef BLNS_UTILS_H
#define BLNS_UTILS_H

#include "params.h"


ZZX         Phi();

ZZX         ModPhi(const ZZX& p);
zz_pX       ModPhi_q(const zz_pX& p);
ZZX         ModPhi_hat(const ZZX& p);
zz_pX       ModPhi_hat_q(const zz_pX& p);

void        OGS_Ortho(     mat_D&  Bt, vec_D&  Norms2, const mat_L& B ); 

void        rot(           mat_L& M, const ZZX& f ); 
mat_zz_p    rot_T(         const zz_pX& f );    
mat_zz_p    rot_vect(      const vec_zz_pX& v ); 

vec_ZZ      Coeffs(        const vec_zz_pX x, const unsigned int l );
vec_ZZ      CoeffsX(       const vec_ZZX   x, const unsigned int l );
vec_zz_pX   CoeffsInv(     const vec_ZZ    c, const unsigned int l );
vec_ZZX     CoeffsInvX(    const vec_ZZ    c, const unsigned int l );
vec_ZZ      CoeffsHat(     const vec_ZZX   x, const unsigned int l );
vec_zz_pX   CoeffsInvHat(  const vec_zz_p  c, const unsigned int l );
vec_ZZX     CoeffsInvHatX( const vec_ZZ    c, const unsigned int l );

vec_zz_pX   sigma_map(     const vec_zz_pX& M, const unsigned int d );

zz_pX       poly_mult(     const vec_zz_pX& f, const vec_zz_pX& g );
zz_pX       poly_mult_hat( const vec_zz_pX& f, const vec_zz_pX& g );

zz_pX       Compute_f(     const mat_zz_p& B_f, const ZZ& x );

ZZ          Norm2(         const vec_ZZ&  v );
ZZ          Norm2X(        const vec_ZZX& v, const unsigned int d );
RR          Norm2R(        const vec_RR&  v );
double      Norm2D(        const vec_D&   v );

double      InnerProdD(    const vec_D& a, const vec_D& b );

#endif
