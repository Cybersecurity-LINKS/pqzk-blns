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

#ifndef BLNS_SQUARES_H
#define BLNS_SQUARES_H

#include "params.h"


void  sum_of_two_squares(ZZ& r0, ZZ& r1, const ZZ& p, const ZZ& t);
void  sum_of_four_squares_2_mod_4(ZZ& x, ZZ& y, ZZ& z, ZZ& w, const ZZ& n);
void  sum_of_four_squares(ZZ& x, ZZ& y, ZZ& z, ZZ& w, const ZZ& n);

void  fast_sum_of_squares(vec_ZZ& v, const ZZ& n);

#endif