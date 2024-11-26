#ifndef BLNS_SQUARES_H
#define BLNS_SQUARES_H

#include "params.h"


void  sum_of_two_squares(ZZ& r0, ZZ& r1, const ZZ& p, const ZZ& t);
void  sum_of_four_squares_2_mod_4(ZZ& x, ZZ& y, ZZ& z, ZZ& w, const ZZ& n);
void  sum_of_four_squares(ZZ& x, ZZ& y, ZZ& z, ZZ& w, const ZZ& n);

void  fast_sum_of_squares(vec_ZZ& v, const ZZ& n);

#endif