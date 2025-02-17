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

#ifndef BLNS_PARAMS_H
#define BLNS_PARAMS_H

#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <NTL/ZZX.h>
#include <NTL/lzz_pX.h>
#include <NTL/RR.h>

using namespace std;
using namespace NTL;

typedef Vec<double> vec_D;
typedef Mat<double> mat_D;
typedef Vec<long>   vec_L; 
typedef Mat<long>   mat_L; 


//============================================================================================
// Parameters for 128 bit security level, 
// derived from LaZer https://eprint.iacr.org/2024/1846.pdf

// Main parameters
#define lambda0 128             // λ, security level 
#define q0      12289           // q, modulus (positive integer)
#define d0      512             // d, degree of the cyclotomic ring R
#define sigma0  165.736617183   // \mathfrak{s}, preimages are sampled from D_(s)^{(m+2)d}
#define N0      (ZZ(1)<<t0)     // N,   domain size of the function f
#define h0      64              // h,   length of the string output by H_M
#define psi0    3               // ψ,   infinity norm of r and the output of H_M
#define lr0     2               // ℓr,  length of c1 (length of randomness r generated by the Holder)
#define lm0     1               // ℓm,  length of c0 (ℓm = ℓ·h / d) = (|idx_pub|+|idx_hid|)·h/d

// Proof systems parameters
#define d_hat   64              // \hat{d}, degree of the cyclotomic ring Rˆ (it is a power of two that divides λ)
#define q2_hat  1125899907006629 // \hat{q}_2 modulus of the proof system Π^ISIS_NIZK (must be divisible by q)
#define xi0     8               // ξ
#define nu0     140             // \nu
#define tau0    8               // τ,  natural number such that 1/p^τ_1 negl(λ)
#define n_ISIS  20              // n   for Π^ISIS_NIZK
#define n_Com   17              // n   for Π^Com_NIZK
#define m2_ISIS 64              // m_2 for Π^ISIS_NIZK, m2 ≥ n + ℓ
#define m2_Com  55              // m_2 for Π^Com_NIZK,  m2 ≥ n + ℓ
#define q1_hat  2199023288779   // \hat{q}_1 modulus of the proof system Π^Com_NIZK  (must be divisible by q)
#define p_bar   68527215        // \overline{p}, an odd integer co-prime to \hat{q}
#define n_i     28              // n_1 == n_2
#define eta_i   1               // η_1 == η_2

// Other parameters
#define t0      512             // t := log(N)      NOTE: t0 must be divisible by d_hat
#define l0      8               // ℓ,  number of attributes (d must divide ℓh) 
#define m0      3               // \textit{m}, length of a2
#define m1_Com  ((idx_hid*h0 + lr0*d0 + d_hat)/d_hat)                           // m_1 for Π^Com_NIZK
#define m1_ISIS (((m0+2)*d0 + (idx_hid*h0 + lr0*d0) + t0 + 2*d_hat) / d_hat)    // m_1 for Π^ISIS_NIZK
#define alpha_1     44          // α_1 ∈ O(√λ)
#define alpha_2     44          // α_2 ∈ O(√λ)
#define alpha_3     44          // α_3 ∈ O(√λ)
#define alpha_bar_1 44          // \overline{α}_1 ∈ O(√λ)
#define alpha_bar_2 44          // \overline{α}_2 ∈ O(√λ)
#define w_max       sqrt(337)   // ω_max(λ)
// NOTE: ω_max is a scalar defined for three specific values of the security level λ, in particular
//       ω_max(λ = 128) = √337, ω_max(λ = 192) = √530, ω_max(λ =256) = √767.
#define N1          128         // \textsf{N} ∈ O(λ)      NOTE: it is different from N0, domain size of the function f
#define idx_hid     4           // | \overline{\idx} |,   number of undisclosed attributes
#define idx_pub     4           // | idx |,               number of disclosed attributes 
// NOTE: idx_pub + idx_hid = l0 = ℓ, number of attributes (d must divide ℓh)

#define k0          32          // k, parameter related to c ∈ C (see BLNS Fig. 8, pag. 26)
//======================================================================================


const RR log2e_Const = (log(exp(RR(1))) / log(RR(2))); //log2(e) with e = EulerConstant

#endif
