#ifndef BLNS_PARAMS_H
#define BLNS_PARAMS_H

#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <NTL/vector.h>
#include <NTL/matrix.h>
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/lzz_pX.h>
#include <NTL/mat_ZZ.h>
#include <NTL/RR.h>
#include <NTL/mat_RR.h>

using namespace std;
using namespace NTL;

typedef Vec<double> vec_D;
typedef Mat<double> mat_D; 
typedef Mat<long>   mat_L; 

/*
//=====================================================================================
//// TOY parameters v1 (small, for TEST purposes ONLY!)
#define lambda0 128     // λ,   security level
#define q0      5       // q, modulus (positive integer)
// NOTE: q s.t. q1_hat / q is small
#define d0      16      // d, degree of the cyclotomic ring R
#define sigma0  6       // \mathfrak{s}, preimages are sampled from D_(s)^{(m+2)d}
#define N0     (1<<8)   // N,   domain size of the function f
#define h0      2       // h,   length of the string output by H_M
#define psi0    2       // ψ,   infinity norm of r and the output of H_M
#define lr0     2       // ℓr,  length of c1 (length of randomness r, must be divisible by \hat{d})
#define lm0     1       // ℓm,  length of c0 (ℓm = ℓ·h / d) = (|idx_pub|+|idx_hid|)·h/d

#define d_hat   8       // \hat{d}, degree of the cyclotomic ring Rˆ
#define q2_hat  15      // \hat{q}_2 modulus of the proof system Π^ISIS_NIZK (must be divisible by q)
#define xi0     17      // ξ
#define nu0     1       // \nu
#define tau0    2       // τ,  natural number such that 1/p^τ_1 negl(λ)
#define n_ISIS  3       // n   for Π^ISIS_NIZK
#define n_Com   3       // n   for Π^Com_NIZK
#define m2_ISIS 6       // m_2 for Π^ISIS_NIZK, m2 ≥ n + ℓ
#define m2_Com  5       // m_2 for Π^Com_NIZK
#define q1_hat  15      // \hat{q}_1 modulus of the proof system Π^Com_NIZK  (must be divisible by q)
#define p_bar   3       // \overline{p}, an odd integer co-prime to \hat{q}
#define n_i     3       // n_1 == n_2 
#define eta_i   1       // η_1 == η_2

#define t0      8       // t := log(N)
#define l0      8       // ℓ,   number of attributes (d must divide ℓh) 
#define m0      3       // \textit{m}, length of a2 (NOTE: different from m_prime!)
#define m_prime 6       // m,   an integer divisible by \hat{d} (e.g. m = d^)
#define m1_Com  ((idx_hid*h0 + lr0*d0 + d_hat)/d_hat)  // m_1 for Π^Com_NIZK
#define m1_ISIS (((m0+2)*d0 + (idx_hid*h0 + lr0*d0) + t0 + 2*d_hat) / d_hat)  // m_1 for Π^ISIS_NIZK

#define alpha_1     (3*11)      // α_1 ∈ O(√λ)
#define alpha_2     (3*11)      // α_2 ∈ O(√λ)
#define alpha_3     (3*11)      // α_3 ∈ O(√λ)
#define alpha_bar_1 (3*11)      // \overline{α}_1 ∈ O(√λ)
#define alpha_bar_2 (3*11)      // \overline{α}_2 ∈ O(√λ)
#define w_max       sqrt(337)     // ω_max(λ)
// NOTE: ω_max is a scalar defined for three specific values of the security level λ, in particular
//       ω_max(λ = 128) = √337, ω_max(λ = 192) = √530, ω_max(λ =256) = √767.
#define N1          128         // \textsf{N} ∈ O(λ)    NOTE: it is different from N0, domain size of the function f
#define idx_hid     4           // | \overline{\idx} |,   number of undisclosed attributes
#define idx_pub     4           // | idx |,               number of disclosed attributes 
// NOTE: idx_pub + idx_hid = l0 = ℓ, number of attributes (d must divide ℓh)

#define k0          32          // k, parameter related to c ∈ C (see BLNS Fig. 8, pag. 26)
*/

/*
//======================================================================================
//// TOY parameters v2 (from CrypTO, for TEST purposes ONLY!)
#define lambda0 128     // λ,   security level
#define q0      17      // q, modulus (positive integer) 
#define d0      16      // d, degree of the cyclotomic ring R
#define sigma0  11      // \mathfrak{s}, preimages are sampled from D_(s)^{(m+2)d}
#define N0     (1<<16)  // N,   domain size of the function f
#define h0      8       // h,   length of the string output by H_M
#define psi0    2       // ψ,   infinity norm of r and the output of H_M
#define lr0     2       // ℓr,  length of c1 (length of randomness r, must be divisible by \hat{d})
#define lm0     1       // ℓm,  length of c0 (ℓm = ℓ·h / d) = (|idx_pub|+|idx_hid|)·h/d

#define d_hat   8       // \hat{d}, degree of the cyclotomic ring Rˆ
#define q2_hat  85      // \hat{q}_2 modulus of the proof system Π^ISIS_NIZK (must be divisible by q)
#define xi0     17      // ξ
#define nu0     350     // \nu
#define tau0    10      // τ,  natural number such that 1/p^τ_1 negl(λ)
#define n_ISIS  30      // n   for Π^ISIS_NIZK
#define n_Com   18      // n   for Π^Com_NIZK
#define m2_ISIS 68      // m_2 for Π^ISIS_NIZK, m2 ≥ n + ℓ
#define m2_Com  51      // m_2 for Π^Com_NIZK,  m2 ≥ n + ℓ
#define p2_hat  17      // \hat{p}_2
#define q1_hat  85      // \hat{q}_1 modulus of the proof system Π^Com_NIZK  (must be divisible by q)
#define p_bar   68527215// \overline{p}, an odd integer co-prime to \hat{q}
#define n_i     28      // n_1 == n_2
#define eta_i   1       // η_1 == η_2

#define t0      16      // t := log(N)
#define l0      2       // ℓ,   number of attributes (d must divide ℓh)
#define m0      3       // \textit{m}, length of a2 (NOTE: different from m_prime!)
#define m_prime 8       // m,   an integer divisible by \hat{d} (e.g. m = d^)
#define m1_Com  ((idx_hid*h0 + lr0*d0 + d_hat)/d_hat)  // m_1 for Π^Com_NIZK
#define m1_ISIS (((m0+2)*d0 + (idx_hid*h0 + lr0*d0) + t0 + 2*d_hat) / d_hat)  // m_1 for Π^ISIS_NIZK

#define alpha_1     (3*11)      // α_1 ∈ O(√λ)
#define alpha_2     (3*11)      // α_2 ∈ O(√λ)
#define alpha_3     (3*11)      // α_3 ∈ O(√λ)
#define alpha_bar_1 (3*11)      // \overline{α}_1 ∈ O(√λ)
#define alpha_bar_2 (3*11)      // \overline{α}_2 ∈ O(√λ)
#define w_max       sqrt(337)   // ω_max(λ)
#define N1          128         // \textsf{N} ∈ O(λ)    NOTE: it is different from N0, domain size of the function f
#define idx_hid 1               // | \overline{\idx} |,   number of undisclosed attributes
#define idx_pub 1               // | idx |,               number of disclosed attributes 
// NOTE: idx_pub + idx_hid = l0 = ℓ, number of attributes (d must divide ℓh)

#define k0      32    // k, parameter related to c ∈ C (see BLNS Fig. 8, pag. 26)
*/


//============================================================================================
// DEMO parameters 

// Main parameters
#define lambda0 128     // λ, security level 
#define q0      12289   // q, modulus (positive integer)
#define d0      512     // d, degree of the cyclotomic ring R
#define sigma0  (1<<27) // \mathfrak{s}, preimages are sampled from D_(s)^{(m+2)d}
#define N0      (ZZ(1)<<t0) // N,   domain size of the function f
#define h0      64      // h,   length of the string output by H_M
#define psi0    2       // ψ,   infinity norm of r and the output of H_M
#define lr0     2       // ℓr,  length of c1 (length of randomness r, must be divisible by \hat{d})
#define lm0     1       // ℓm,  length of c0 (ℓm = ℓ·h / d) = (|idx_pub|+|idx_hid|)·h/d

// Proof systems parameters
#define d_hat   64      // \hat{d}, degree of the cyclotomic ring Rˆ (it is a power of two that divides λ)
#define q2_hat  12289   // \hat{q}_2 modulus of the proof system Π^ISIS_NIZK (must be divisible by q)
#define xi0     17      // ξ
#define nu0     350     // \nu
#define tau0    10      // τ,  natural number such that 1/p^τ_1 negl(λ)
#define n_ISIS  30      // n   for Π^ISIS_NIZK
#define n_Com   18      // n   for Π^Com_NIZK
#define m2_ISIS 68      // m_2 for Π^ISIS_NIZK, m2 ≥ n + ℓ
#define m2_Com  51      // m_2 for Π^Com_NIZK,  m2 ≥ n + ℓ
#define q1_hat  12289   // \hat{q}_1 modulus of the proof system Π^Com_NIZK  (must be divisible by q)
#define p_bar   68527215 // \overline{p}, an odd integer co-prime to \hat{q}
#define n_i     28      // n_1 == n_2
#define eta_i   1       // η_1 == η_2

// Other parameters
#define t0      d_hat   // t := log(N)      NOTE: t0 must be divisible by d_hat
#define l0      8       // ℓ,  number of attributes (d must divide ℓh) 
#define m0      3       // \textit{m}, length of a2 (NOTE: different from m_prime!)
#define m_prime 128     // m,  an integer divisible by \hat{d} (e.g. m = d^)
#define m1_Com  ((idx_hid*h0 + lr0*d0 + d_hat)/d_hat)  // m_1 for Π^Com_NIZK
#define m1_ISIS (((m0+2)*d0 + (idx_hid*h0 + lr0*d0) + t0 + 2*d_hat) / d_hat)  // m_1 for Π^ISIS_NIZK
#define alpha_1     1000        // α_1 ∈ O(√λ)
#define alpha_2     1000        // α_2 ∈ O(√λ)
#define alpha_3     1000        // α_3 ∈ O(√λ)
#define alpha_bar_1 1000        // \overline{α}_1 ∈ O(√λ)
#define alpha_bar_2 1000        // \overline{α}_2 ∈ O(√λ)
#define w_max       sqrt(337)   // ω_max(λ)
// NOTE: ω_max is a scalar defined for three specific values of the security level λ, in particular
//       ω_max(λ = 128) = √337, ω_max(λ = 192) = √530, ω_max(λ =256) = √767.
#define N1          128         // \textsf{N} ∈ O(λ)    NOTE: it is different from N0, domain size of the function f
#define idx_hid     4           // | \overline{\idx} |,   number of undisclosed attributes
#define idx_pub     4           // | idx |,               number of disclosed attributes 
// NOTE: idx_pub + idx_hid = l0 = ℓ, number of attributes (d must divide ℓh)

#define k0          32          // k, parameter related to c ∈ C (see BLNS Fig. 8, pag. 26)
//======================================================================================



const RR log2e_Const = (log(exp(RR(1))) / log(RR(2))); //log2(e) with e = EulerConstant

#endif
