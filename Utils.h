#ifndef BLNS_UTILS_H
#define BLNS_UTILS_H

#include "params.h"


ZZX         Phi();

ZZX         ModPhi(const ZZX& p);
zz_pX       ModPhi_q(const zz_pX& p);
ZZX         ModPhi_hat(const ZZX& p);
zz_pX       ModPhi_hat_q(const zz_pX& p);

void        OGS_Ortho(     mat_D&  Bt, vec_D&  Norms2, const mat_L& B ); 

void        rot(         mat_L& M,   const ZZX& f );
void        rot_T(    mat_zz_p& M, const zz_pX& f );
void        rot_vect( mat_zz_p& R, const vec_zz_pX& v );

void        CoeffsX(vec_ZZ& coeffs_x, const vec_ZZX& x, const unsigned int& l);
void        CoeffsInv(vec_zz_pX& x, const vec_ZZ& c, const unsigned int& l);
void        CoeffsInvX(vec_ZZX& x, const vec_ZZ& c, const unsigned int& l);
void        CoeffsHat(vec_ZZ& coeffs_x, const vec_ZZX& x, const unsigned int& l);
void        CoeffsInvHat(vec_zz_pX& x, const vec_zz_p& c, const unsigned int& l);
void        CoeffsInvHatX(vec_ZZX& x, const vec_ZZ& c, const unsigned int& l);

void        sigma_map(vec_zz_pX& N, const vec_zz_pX& M, const unsigned int& d);

zz_pX       poly_mult(     const vec_zz_pX& f, const vec_zz_pX& g );
zz_pX       poly_mult_hat( const vec_zz_pX& f, const vec_zz_pX& g );

zz_pX       Compute_f(     const mat_zz_p& B_f, const ZZ& x );

ZZ          Norm2(         const vec_ZZ&  v );
ZZ          Norm2X(        const vec_ZZX& v, const unsigned int& d );
double      Norm2D(        const vec_D&   v );

double      InnerProdD(    const vec_D& a, const vec_D& b );

#endif
