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

#include "ISIS.h"


//==============================================================================
// Preprocessing_ISIS - Preprocessing function (PreprocessingProve^HISIS_ISIS). 
//      It takes as input (P, s, B_goth_s, C, r, B_goth_r)
//      with B_goth_s ≥ ||s|| and B_goth_r ≥ ||r||. 
//      It returns s and r such that B_goth_s = ||s|| and B_goth_r = ||r||, 
//      together with P and C filled with the appropriate number of zeros.
// 
// Inputs:
// - P:             matrix P  ∈ Z^[d × (m+2)·d]_(q_hat) 
// - s0:            vector s0 ∈ Z^((m+2)·d)_(q_hat)
// - B_goth_s2:     bound  B_goth_s^2 ∈ Z≥0 (it is a scalar)
// - C:             matrix C  ∈ Z^[d × (ℓm+ℓr)·d]_(q_hat)
// - r0:            vector r0 ∈ Z^(|idx_hid|·h + ℓr·d)_(q_hat)
// - B_goth_r2:     bound  B_goth_r^2 ∈ Z≥0 (it is a scalar)
//  
// Output:
// - P1:            matrix P1 ∈ Z^[d × (m+2)·d + d_hat]_(q_hat) 
// - s1:            vector s1 ∈ Z^((m+2)·d + d_hat)_(q_hat)
// - C1:            matrix C1 ∈ Z^[d × ((ℓm+ℓr)·d + d_hat)]_(q_hat)
// - r1:            vector r1 ∈ Z^(|idx_hid|·h + ℓr·d + d_hat)_(q_hat)
//==============================================================================
// NOTE: zero padding of P, C, already done in H_VerPres
void  Preprocessing_ISIS(vec_ZZ& s1, vec_ZZ& r1, const vec_ZZ& s0, const ZZ& B_goth_s2, const vec_ZZ& r0, const ZZ& B_goth_r2)
{    
    // NOTE: assuming that current modulus is q2_hat (not q0)
    long    i;
    ZZ      diff; //a1, a2, a3, a4;
    vec_ZZ  v;

    const long  m2d     = (m0 + 2)*d0;    // (m+2)·d
    // const long  lmlrd   = (lm0 + lr0)*d0; // (ℓm+ℓr)·d
    const long  idxhlrd = (idx_hid * h0) + (lr0 * d0); //|idx_hid|·h + ℓr·d


    // diff = B_goth_s^2 − ||s||^2
    diff = (B_goth_s2 - Norm2(s0));
        
    if (diff < 0)
    {
        cout << "ERROR! (B_goth_s^2 - ||s||^2) must be positive: " << diff << endl;
        assert(diff >= 0);
    }
    else if (diff == 0)
    {
        cout << "WARNING! (B_goth_s^2 - ||s||^2) == 0 " << endl;
    }

    // diff = diff % q2_hat;
    // // NOTE: diff mod q2_hat, to speed up sum_of_four_squares


    // // 1. (a1, a2, a3, a4) ← SumOfFourSquares(B_goth_s^2 − ||s||^2),   (a1, a2, a3, a4) ∈ Z^4
    // sum_of_four_squares(a1, a2, a3, a4, diff);


    // // 2. a ← (a1, a2, a3, a4, 0, ... , 0),  a ∈ Z^(d_hat) 
    // // NOTE: add d_hat − 4 zeros

    // 3. s1 ← (s, a),   s1 ∈ Z^((m+2)·d + d_hat)
    s1.SetLength(m2d + d_hat);
           
    for(i=0; i<m2d; i++)
    {
        s1[i] = s0[i];
    }

    // s1[m2d]   = a1;
    // s1[m2d+1] = a2;
    // s1[m2d+2] = a3;
    // s1[m2d+3] = a4;

    // NOTE: sum_of_four_squares (too slow) replaced with fast_sum_of_squares
    fast_sum_of_squares(v, diff);
    
    for(i=0; i<v.length(); i++)
    {
        s1[m2d+i] = v[i];
    }


    // 4. P1 ← [P,  0_(d × d_hat)],   P1 ∈ Z^[d × (m+2)·d + d_hat]_(q_hat)     
    // P1.SetDims(d0, (m2d + d_hat));

    // for(i=0; i<d0; i++)
    //     {   
    //     for(j=0; j<m2d; j++)
    //     {
    //         P1[i][j] = P[i][j];
    //     }
    // }


    // diff = B_goth_r^2 − ||r||^2
    diff = (B_goth_r2 - Norm2(r0));

    if (diff < 0)
    {
        cout << "ERROR! (B_goth_r^2 - ||r||^2) must be positive: " << diff << endl;
        assert(diff >= 0);
    }
     else if (diff == 0)
    {
        cout << "WARNING! (B_goth_r^2 - ||r||^2) == 0 " << endl;
    }

    // diff = diff % q2_hat;
    // // NOTE: diff mod q2_hat, to speed up sum_of_four_squares


    // // 5. (b1, b2, b3, b4) ← SumOfFourSquares(B_goth_r^2 − ||r||^2),   (b1, b2, b3, b4) ∈ Z^4
    // sum_of_four_squares(a1, a2, a3, a4, diff);

    // // 6. b ← (b1, b2, b3, b4, 0, ... , 0),  b ∈ Z^(d_hat) 
    // // NOTE: add d_hat − 4 zeros

    // 7. r1 ← (r, b),   r1 ∈ Z^(|idx_hid|·h + ℓr·d + d_hat)
    r1.SetLength(idxhlrd + d_hat);

    for(i=0; i<idxhlrd; i++)
    {
        r1[i]     = r0[i];
    }

    // r1[idxhlrd]   = a1;
    // r1[idxhlrd+1] = a2;
    // r1[idxhlrd+2] = a3;
    // r1[idxhlrd+3] = a4;

    // NOTE: sum_of_four_squares (too slow) replaced with fast_sum_of_squares
    fast_sum_of_squares(v, diff);

    for(i=0; i<v.length(); i++)
    {
        r1[idxhlrd+i] = v[i];
    }

    
    // 8. C1 ← [C,  0_(d × d_hat)],   C1 ∈ Z^[d × ((ℓm+ℓr)·d + d_hat)]_(q_hat)   
    // C1.SetDims(d0, (lmlrd + d_hat));

    // for(i=0; i<d0; i++)
    //     {   
    //     for(j=0; j<lmlrd; j++)
    //     {
    //         C1[i][j] = C[i][j];
    //     }
    // }

    // 9. return (P1, s1, C1, r1)
}


//==============================================================================
// Prove_ISIS   -   Computes the Commitment (Prove^HISIS_ISIS). 
//                  This function takes as input crs, x, and w.
//                  If all the checks pass, it returns the proof π.
// 
// Inputs:
// - inputStr:      string containing the initial seed for crs
// - crs:           structure crs_ISIS, generated by Hcrs from inputStr
// - ipk:           Issuer public key
// - (P, C, mex, B_f, 
//    Bounds, aux): they correspond to x
// -  w0:           it contains the vectors (s0, r0, u0)
//  
// Output:
// -  Pi:           proof (π) structure 
//==============================================================================
void Prove_ISIS(PROOF_I_t& Pi, const string& inputStr, const CRS_t& crs, const IPK_t& ipk, const mat_zz_p& P, const mat_zz_p& C, const vec_zz_p& mex, const mat_zz_p& B_f, const vec_ZZ& Bounds, const ZZ& aux, const Vec<vec_ZZ>& w0)
{
    // NOTE: assuming that current modulus is q2_hat (not q0)
   
    unsigned long        idx, i, j, k;
    long                 rst, b1, b2, b3;
    Mat<zz_pX>           B, D2; //D2_2_1
    vec_zz_pX            s_1_mod, s_2_mod, g; //t_B;
    vec_zz_pX            h_part1, h_part2;
    vec_zz_pX            r_j, p_j, Beta_j, c_r_j, mu, m, s_hat, tmp_vec, y;
    vec_zz_pX            sigma_s_1, d_1, acc_vec, D2_y;
    vec_ZZX              c_s1, c_s2;
    vec_ZZX              s, r, u, y_1, y_2, y_3, s_1, s_2;
    ZZX                  c;
    zz_pX                h_part3, h_part4, h_part5;
    zz_pX                acc, delta_1, delta_2, delta_3, f1; //d_0;
    Vec<vec_zz_pX>       e_, sigma_e_, sigma_p_, sigma_Beta_, sigma_c_r_;
    Vec<vec_zz_pX>       sigma_r_, sigma_r_s_, sigma_r_r_, sigma_r_u_;    
    stringstream         ss;
    mat_L                R_goth;
    vec_ZZ               s0, r0, u0, s1, r1, coeffs_s1, coeffs_y3, coeffs_R_goth_mult_s1;    
    mat_zz_p             gamma, C_m, C_r;
    vec_zz_p             ones;
    vec_zz_pX            s_mod, r_mod,  u_mod, coeffs_ones, u_m_ones, sigma_u_m_ones, sigma_ones;
    vec_zz_pX            sigma_s, sigma_r;
    vec_zz_p             e_tmp, m_C;   
    RR                   alpha_i,     B_goth;
    ZZ                   B_goth_s2,   B_goth_r2;
    zz_p                 B_goth_s2_p, B_goth_r2_p; //sums;

    // Initialise constants
    const unsigned long n           = n_ISIS;
    const unsigned long m1          = m1_ISIS;
    const unsigned long m2          = m2_ISIS;
    const unsigned long m2d         = (m0 + 2)*d0;    // (m+2)·d
    const unsigned long idxhlrd     = (idx_hid * h0) + (lr0 * d0); //|idx_hid|·h + ℓr·d    
    const unsigned long n256        = (256/d_hat);    
    const unsigned long t_d_hat     = (t0/d_hat);
    const unsigned long m1_n256_tau = 2*m1 + 2*(n256 + tau0);

    // Initialise the "goth" constants
    // M1 := exp(sqrt(2(λ+1)/log e) * 1/α_1 + 1/2α_1^2
    // M2 := exp(sqrt(2(λ+1)/log e) * 1/α_2 + 1/2α_2^2
    // M3 := exp(sqrt(2(λ+1)/log e) * 1/α_3 + 1/2α_3^2    
    alpha_i = RR(alpha_1);
    const RR  M_1 = exp( sqrt( RR(2*(lambda0 + 1)) / log2e_Const ) * 1/alpha_i + 1/(2*sqr(alpha_i)));
    alpha_i = RR(alpha_2);
    const RR  M_2 = exp( sqrt( RR(2*(lambda0 + 1)) / log2e_Const ) * 1/alpha_i + 1/(2*sqr(alpha_i)));
    alpha_i = RR(alpha_3);
    const RR  M_3 = exp( sqrt( RR(2*(lambda0 + 1)) / log2e_Const ) * 1/alpha_i + 1/(2*sqr(alpha_i)));
    
    
    // 1. Retrieve (A_1, A_2, B_y, B_g, b) ← crs_ISIS
    // A_1     = crs[0];    // ∈ R^^(n x m1)_(q_hat)
    // A_2     = crs[1];    // ∈ R^^(n x m2)_(q_hat)
    // B_y     = crs[2];    // ∈ R^^(256/d^ x m2)_(q_hat)
    // B_g     = crs[3];    // ∈ R^^(tau0 x m2)_(q_hat)
    // b       = crs[4][0]; // ∈ R^^(m2)_(q_hat)         NOTE: b in crs is a (1 x m_2) matrix
        
    // 2. (P, C, mex, B_f, Bounds, aux) ← x
    // NOTE: directly provided as inputs

    // 3. (B_goth_s^2, B_goth_r^2) ← Bounds,   (B_goth_s^2, B_goth_r^2) ∈ Z^2
    B_goth_s2 = Bounds[0];
    B_goth_r2 = Bounds[1];

    // B_goth = sqrt(B_goth_s^2 + B_goth_r^2 + t0)
    B_goth = sqrt( conv<RR>( B_goth_s2 + B_goth_r2 + t0 ) );

    B_goth_s2_p = conv<zz_p>(B_goth_s2);
    B_goth_r2_p = conv<zz_p>(B_goth_r2);

    // s1_goth = alpha_1*nu0*B_goth
    // s2_goth = alpha_2*nu0*sqrt(m2*d_hat)
    // s3_goth = alpha_3*w_max(lambda0)*B_goth    
    const RR  s1_goth = RR(alpha_1 * nu0) * B_goth;
    const RR  s2_goth = RR(alpha_2 * nu0) * sqrt( RR(m2*d_hat) );
    const RR  s3_goth = RR(alpha_3 * w_max) * B_goth;
    const double s1_goth_d = conv<double>(s1_goth);
    const double s2_goth_d = conv<double>(s2_goth);
    const double s3_goth_d = conv<double>(s3_goth);


    // 4. (s0, r0, u0) ← w0    
    s0 = w0[0]; // ∈ Z^((m+2)·d)
    r0 = w0[1]; // ∈ Z^(|idx_hid|·h + ℓr·d)
    u0 = w0[2]; // ∈ {0,1}^t
  
    // 5. (P, s, C, r) ← PreprocessingProve^HISIS_ISIS (P, s, B_goth_s, C, r, B_goth_r)
    // P ∈ Z^[d × (m+2)·d + d_hat]_(q_hat) 
    // s ∈ Z^((m+2)·d + d_hat)
    // C ∈ Z^[d × ((ℓm+ℓr)·d + d_hat)]_(q_hat)
    // r ∈ Z^(|idx_hid|·h + ℓr·d + d_hat)   
    // Preprocessing_ISIS(P, s1, C, r1, P1, s0, B_goth_s, C1, r0, B_goth_r);
    Preprocessing_ISIS(s1, r1, s0, B_goth_s2, r0, B_goth_r2);
    
    // 6. s ← Coeffs^−1(s1)    
    CoeffsInvHatX(s, s1, (m2d+d_hat)/d_hat); // s     ∈ R^^(((m+2)·d+d_hat)/d_hat)
    s_mod = conv<vec_zz_pX>( s );             // s_mod ∈ R^^(((m+2)·d+d_hat)/d_hat)_(q_hat)

    //    r ← Coeffs^−1(r1)
    CoeffsInvHatX(r, r1, (idxhlrd+d_hat)/d_hat); // r     ∈ R^^((|idx_hid|·h + ℓr·d + d_hat)/d_hat)
    r_mod = conv<vec_zz_pX>( r );                 // r_mod ∈ R^^((|idx_hid|·h + ℓr·d + d_hat)/d_hat)_(q_hat)

    //    u ← Coeffs^−1(u0) 
    CoeffsInvHatX(u, u0, t_d_hat); // u     ∈ R^^(t/d_hat)
    u_mod = conv<vec_zz_pX>( u );   // u_mod ∈ R^^(t/d_hat)_(q_hat) 
    
       
    // 7. Initialize rst ← 0,   rst ∈ Z, scalar
    rst     = 0;
    
    // 8. Initialize idx ← 0,   idx ∈ N, scalar
    idx     = 0;
    
    // 9. s_1 ← (s, r, u),   s_1 ∈ R^^(m1)
    s_1.SetLength(m1);
    // NOTE: m1 := ((m + 2)·d + d_hat      + (|idx_hid|·h + ℓr·d + d_hat)    + t) / d_hat
    //           = ((m0+2)*d0+d_hat)/d_hat + (idx_hid*h0+lr0*d0+d_hat)/d_hat + t0/d_hat
    //           = (m2d+d_hat) / d_hat     + (idxhlrd+d_hat) / d_hat         + t_d_hat
    
    k = 0;

    for(i=0; i<((m2d+d_hat)/d_hat); i++) 
    {
        s_1[k] = s[i]; 
        k++;       
    }
    
    for(i=0; i<((idxhlrd+d_hat)/d_hat); i++) 
    {
        s_1[k] = r[i];  
        k++;      
    }
   
    for(i=0; i<(t_d_hat); i++) 
    {
        s_1[k] = u[i];
        k++;      
    }

    s_1_mod = conv<vec_zz_pX>( s_1 );
    // NOTE: modulo q_hat on all coefficients

    // coeffs_s1 ← Coeffs(s_1),   coeffs_s1 ∈ Z^(m1*d_hat)   
    CoeffsHat(coeffs_s1, s_1, m1);
    // NOTE: coeffs_s1 contains the same coefficients listed in w0


    // Initialize e ∈ R^^(256 x 256/d_hat)_(q_hat)
    e_.SetLength(256);    
    // NOTE: defined as e_ to distinguish it from the Euler constant e 

    e_tmp.SetLength(256);

    for(k=0; k<256; k++)
    {
        e_tmp[k] = 0;
    }

    for(j=0; j<256; j++)
    {
        // Temporary coefficient vector to create e_j: it is a unit vector with its j-th coefficient being 1
        e_tmp[j] = 1;        

        // e_[j].SetLength(n256);
        CoeffsInvHat(e_[j], e_tmp, n256);

        // Reset the e_tmp coefficient vector
        e_tmp[j] = 0;
    }

    // Precompute σ(e_j), σ(p_j), σ(β_j), σ(s), σ(r), σ(s_1)
    sigma_e_.SetLength(256);
    sigma_p_.SetLength(d0);
    sigma_Beta_.SetLength(d0);
    sigma_s.SetLength((m2d+d_hat)/d_hat);
    sigma_r.SetLength((idxhlrd+d_hat)/d_hat);
    sigma_s_1.SetLength(m1);
    
    for(j=0; j<256; j++)        
    {
        sigma_map(sigma_e_[j], e_[j], d_hat);  
    }       

    for(j=0; j<d0; j++)        
    {
        CoeffsInvHat(p_j, P[j], (m2d+d_hat)/d_hat);
        CoeffsInvHat(Beta_j, B_f[j], t_d_hat); 
        sigma_map(sigma_p_[j], p_j, d_hat);        
        sigma_map(sigma_Beta_[j], Beta_j, d_hat);                  
    }

    sigma_map(sigma_s, s_mod, d_hat);
    sigma_map(sigma_r, r_mod, d_hat);
    sigma_map(sigma_s_1, s_1_mod, d_hat);

    // Create C_m and C_r
    C_m.SetDims(d0, (idx_pub*h0));
    C_r.SetDims(d0, (idxhlrd+d_hat)); 
    // NOTE: C = [C_m C_r] ∈ Z^[d × ((ℓm+ℓr)·d+d_hat)]_(q_hat), C_m has |idx_pub|·h columns

    for(i=0; i<d0; i++)
    {
        for(j=0; j<(idx_pub*h0); j++)        
        {
            C_m[i][j] = C[i][j];
        }

        for(j=0; j<(idxhlrd+d_hat); j++)        
        {
            C_r[i][j] = C[i][(idx_pub*h0)+j];
        }
    }    

    // Compute m_C := C_m * m
    m_C.SetLength(d0);
    m_C = ( C_m * mex );

    C_m.kill();

    // Precompute σ(C_r,j), h_part2
    sigma_c_r_.SetLength(d0);
    h_part2.SetLength(d0);

    for(j=0; j<d0; j++)        
    {
        CoeffsInvHat(c_r_j  , C_r[j], (idxhlrd+d_hat)/d_hat);      
        sigma_map(sigma_c_r_[j] , c_r_j, d_hat);
        h_part2[j] = poly_mult_hat(sigma_p_[j], s_mod) - poly_mult_hat(sigma_Beta_[j], u_mod) - m_C[j] - poly_mult_hat(sigma_c_r_[j], r_mod);
    }

    C_r.kill();

    // Precompute σ(u − Coeffs^−1(1^t)), h_part3, h_part4, h_part5
    ones.SetLength(t0);             // ∈ Z^^(t)_(q_hat)
    coeffs_ones.SetLength(t_d_hat); // ∈ R^^(t/d_hat)_(q_hat)
    u_m_ones.SetLength(t_d_hat);    // ∈ R^^(t/d_hat)_(q_hat)

    for(i=0; i<t0; i++) 
    {
        ones[i] = 1;
    }

    CoeffsInvHat(coeffs_ones, ones, t_d_hat);

    for(i=0; i<t_d_hat; i++) 
    {
        u_m_ones[i] = u_mod[i] - coeffs_ones[i]; 
    }

    sigma_map(sigma_ones    , coeffs_ones, d_hat);
    sigma_map(sigma_u_m_ones, u_m_ones, d_hat);

    h_part3 = poly_mult_hat(sigma_s, s_mod) - B_goth_s2_p;
    h_part4 = poly_mult_hat(sigma_r, r_mod) - B_goth_r2_p;
    h_part5 = poly_mult_hat(sigma_u_m_ones, u_mod);


    // 10. while (rst == 0 ∧ idx < N) do
    while((rst == 0) && (idx < N1))
    {
        b1 = 0; b2 = 0; b3 = 0;
       
        // 11. Increment idx
        idx = idx + 1;
        // cout << "idx = " << idx << endl;
        
        // 12. Random generation of s_2 ∈ R^^(m2)
        s_2.SetLength(m2);

        for(i=0; i<m2; i++)
        {
            s_2[i].SetLength(d_hat);

            for(j=0; j<d_hat; j++)
            {
                s_2[i][j] = RandomBnd(3) - 1;
                // NOTE: uniform distribution on ternary polynomials chi, that sample coeffs from {-1,0,1}
            }
        }

        s_2_mod = conv<vec_zz_pX>( s_2 ); 
        // NOTE: modulo q_hat on all coefficients


        // 13. t_A = A_1*s_1 + A_2*s_2,  tA ∈ R^^(n)_(q_hat)
        Pi.t_A.SetLength(n);
        acc.SetLength(d_hat);

        for(i=0; i<n; i++)
        {
            Pi.t_A[i].SetLength(d_hat);

            // acc = 0;
            clear(acc);     

            for(j=0; j<m1; j++)
            {
                acc += ModPhi_hat_q( crs[0][i][j] * s_1_mod[j] );
            }        

            for(j=0; j<m2; j++)
            {
                acc += ModPhi_hat_q( crs[1][i][j] * s_2_mod[j] ); 
            }            

            Pi.t_A[i] = acc;
        }

        
        // 14. Random generation of the y_1 ∈ R^^m1,  y_2 ∈ R^^m2,  y_3 ∈ R^^(256/d_hat)
        y_1.SetLength(m1);
        y_2.SetLength(m2);
        y_3.SetLength(n256);

        for(i=0; i<m1; i++)
        {
            y_1[i].SetLength(d_hat);

            for(j=0; j<d_hat; j++)
            {
                ZSampler(y_1[i][j], s1_goth_d, 0);
                // NOTE: implicitly sample the vector of coefficients and then convert it to a polynomial vector
            }
        }

        for(i=0; i<m2; i++)
        {
            y_2[i].SetLength(d_hat);

            for(j=0; j<d_hat; j++)
            {
                ZSampler(y_2[i][j], s2_goth_d, 0);
                // NOTE: implicitly sample the vector of coefficients and then convert it to a polynomial vector
            }
        }

        for(i=0; i<n256; i++)
        {
            y_3[i].SetLength(d_hat);

            for(j=0; j<d_hat; j++)
            {
                ZSampler(y_3[i][j], s3_goth_d, 0);
                // NOTE: implicitly sample the vector of coefficients and then convert it to a polynomial vector
            }
        }


        // 15. Random generation of g ∈ R^^(tau)_(q_hat)
        g.SetLength(tau0);

        for(i=0; i<tau0; i++)
        {
            g[i].SetLength(d_hat);
            g[i] = random_zz_pX(d_hat);            
            g[i][0] = 0;        
            // NOTE: the constant term of g (x^0) must be zero 
        }


        // 16. w = A1*y1 + A2*y2,  w ∈ R^^(n)_(q_hat)
        Pi.w.SetLength(n);
        // NOTE: it is different from the input w (= w0, from Prove_Init)
        acc.SetLength(d_hat);

        for(i=0; i<n; i++)
        {
            Pi.w[i].SetLength(d_hat);

            // acc = 0;
            clear(acc);

            for(j=0; j<m1; j++)
            {
                acc += ModPhi_hat_q( crs[0][i][j] * conv<zz_pX>( y_1[j] ) );
            }        

            for(j=0; j<m2; j++)
            {
                acc += ModPhi_hat_q( crs[1][i][j] * conv<zz_pX>( y_2[j] ) );
            }                

            Pi.w[i] = acc;
            // NOTE: modulo q_hat on all coefficients (zz_pX)
        }


        // 17. t_y = B_y*s2 + y3,  t_y ∈ R^^(256/d_hat)_(q_hat)
        Pi.t_y.SetLength(n256);
        acc.SetLength(d_hat);

        for(i=0; i<n256; i++)
        {
            Pi.t_y[i].SetLength(d_hat);

            // acc = 0;
            clear(acc);

            for(j=0; j<m2; j++)        
            {
                acc += ModPhi_hat_q( crs[2][i][j] * s_2_mod[j] );
            }
                
            acc += conv<zz_pX>( y_3[i] );
                    
            Pi.t_y[i] = acc;
            // NOTE: modulo q_hat on all coefficients (zz_pX)
        }


        // 18. t_g = B_g*s2 + g,  t_g ∈ R^^(tau)_(q_hat)
        Pi.t_g.SetLength(tau0);
        acc.SetLength(d_hat);

        for(i=0; i<tau0; i++)
        {
            Pi.t_g[i].SetLength(d_hat);

            // acc = 0;
            clear(acc);               

            for(j=0; j<m2; j++)        
            {
                acc += ModPhi_hat_q( crs[3][i][j] * s_2_mod[j] );
            }
                
            acc += g[i];
            
            Pi.t_g[i] = acc;
            // NOTE: modulo q_hat on all coefficients (zz_pX)
        }


        // 19. a_1 ← (t_A, t_y, t_g, w) 
        ss.str("");    ss.clear();
        // ss << crs << P << C << mex << B_f << Bounds << aux << t_A << t_y << t_g << w;
        ss << inputStr << ipk.a1 << ipk.a2 << ipk.c0 << ipk.c1 << mex << B_f << Bounds << aux << Pi.t_A << Pi.t_y << Pi.t_g << Pi.w;
        // NOTE: using inputStr, ipk, instead of crs, P, C to speedup Hash_Init

        // 20. (R_goth_0, R_goth_1) = H(1, crs, x, a_1)
        // 21. R_goth = R_goth_0 - R_goth_1
        HISIS1(R_goth, "1" + ss.str());
        // NOTE: R_goth ∈ {-1, 0, 1}^(256 x m_1*d_hat),
        //       equivalent to (R_goth_0 - R_goth_1) in BLNS
        
        // 22. coeffs_s1 ← Coeffs(s_1),   coeffs_s1 ∈ Z^(m1*d_hat)   
        // CoeffsHat(coeffs_s1, s_1, m1);
        // NOTE: precomputed after row 9
                
        // 23. coeffs_y3 ← Coeffs(y_3),   coeffs_y3 ∈ Z^(256)   
        CoeffsHat(coeffs_y3, y_3, n256);
        
        // 24.  z_3 = y_3 + R_goth*s_1,   z_3 ∈ Z^(256)   
        coeffs_R_goth_mult_s1.SetLength(256);
        Pi.z_3.SetLength(256);
        // NOTE: This equation is performed in Z not in polynomials, needing Coeffs() transformation of y_3, R_goth, and s_1

        for(i=0; i<256; i++)
        {
            coeffs_R_goth_mult_s1[i] = ( conv<vec_ZZ>(R_goth[i]) * coeffs_s1 );       
            // NOTE: this term corresponds to InnerProduct(result, coeffs_R_goth[i], coeffs_s1);

            Pi.z_3[i] = coeffs_y3[i] + coeffs_R_goth_mult_s1[i];
        }

        
        // 25. b3 ← Rej (z_3, R_goth * s_1, s3_goth, M_3),    b3 ∈ {0, 1}
        b3 = Rej_v_ZZ(Pi.z_3, coeffs_R_goth_mult_s1, s3_goth, M_3);
        
        // NOTE: if b3 == 0, continue the while loop (skip next rows until 53, then go to row 7)
        if (b3 == 0)
        {            
            rst = 0;
            continue;
        }

        
        // 26. a_2 ← z_3,   a2 ∈ Z^256    
        ss << Pi.z_3;

        // 27. gamma ← H(2, crs, x, a1, a2),   gamma ∈ Z^(tau0 x 256+d0+3)_q_hat
        HISIS2(gamma, "2" + ss.str());    
        // NOTE: gamma has 256+d0+3 columns in ISIS, while 256+d0+1 in Com

        
        // Precompute r_j, σ(r_j), σ(r_s,j), σ(r_r,j), σ(r_u,j), h_part1
        r_j.SetLength(m1);        
        sigma_r_.SetLength(256);
        sigma_r_s_.SetLength(256);
        sigma_r_r_.SetLength(256);
        sigma_r_u_.SetLength(256);
        h_part1.SetLength(256);

        for(j=0; j<256; j++)        
        {
            CoeffsInvHat(r_j,  conv<vec_zz_p>(R_goth[j]), m1 );
            sigma_map(sigma_r_[j], r_j, d_hat);

            // NOTE: (r_s,j , r_r,j , r_u,j ) ← r_j  at row 42, where:   
            //       r_s,j ∈ R^^(((m+2)d+d_hat)/d_hat)_(q_hat)
            //       r_r,j ∈ R^^((|idx_hid|·h+ℓr·d+d_hat)/d_hat)_(q_hat)
            //       r_u,j ∈ R^^(t/d_hat)_(q_hat)
            // NOTE: m1 = m1_ISIS = (((m+2)*d+d_hat)/d_hat) + (|idx_hid|·h + ℓr·d + d_hat)/d_hat + t/d_hat

            sigma_r_s_[j].SetLength((m2d+d_hat)/d_hat);
            sigma_r_r_[j].SetLength((idxhlrd+d_hat)/d_hat);
            sigma_r_u_[j].SetLength(t_d_hat); 

            for(k=0; k<((m2d+d_hat)/d_hat); k++)
            {
                sigma_r_s_[j][k] = sigma_r_[j][k];
            }

            for(k=0; k<((idxhlrd+d_hat)/d_hat); k++)
            {
                sigma_r_r_[j][k] = sigma_r_[j][k + (m2d+d_hat)/d_hat];
            }

            for(k=0; k<(t_d_hat); k++)
            {
                sigma_r_u_[j][k] = sigma_r_[j][k + (m2d+d_hat)/d_hat + (idxhlrd+d_hat)/d_hat];
            }

            h_part1[j] = poly_mult_hat(sigma_r_[j], s_1_mod) + poly_mult_hat(sigma_e_[j], conv<vec_zz_pX>( y_3 )) + (conv<zz_p>( - Pi.z_3[j] ));
        }

       
        // Initialize h ∈ R^^(tau)_(q_hat)
        Pi.h.SetLength(tau0);    
        acc.SetLength(d_hat);
        clear(acc);


        // 28. for i ∈ [τ] do
        for(i=0; i<tau0; i++) 
        {
            // 29. Compute h_i,   h_i ∈ R^_(q_hat)
            Pi.h[i].SetLength(d_hat);
            
            acc = g[i];

            for(j=0; j<256; j++)        
            {
                acc += gamma[i][j] * h_part1[j];
            }

            for(j=0; j<d0; j++)        
            {
                acc += gamma[i][256+j] * h_part2[j];
            }
            
            acc +=  gamma[i][256+d0]   * h_part3 +
                    gamma[i][256+d0+1] * h_part4 +  
                    gamma[i][256+d0+2] * h_part5;
                                
            // 30. h ← (h_1, . . . , h_τ),   h ∈ R^^(τ)_(q_hat)
            Pi.h[i] = acc;
        }


        // 31. a_3 ← h,   a_3 ∈ R^^(τ)_(q_hat)    
        ss << Pi.h;

        // 32. μ ← H(3, crs, x, a1, a2, a3),   μ ∈ R^^(τ)_(q_hat)
        HISIS3(mu, "3" + ss.str());

        // 33. B   ← [B_y; B_g],   B ∈ R^^((256/d_hat + tau) x m2)_(q_hat)
        B.SetDims((n256 + tau0), m2);

        // 34. t_B ← [t_y; t_g],   t_B ∈ R^^(256/d_hat + tau)_(q_hat)
        // t_B.SetLength(n256 + tau0);
        // NOTE: t_B is unused in Prove_ISIS

        // 35. m   ← [y_3; g],     m ∈ R^^(256/d_hat + tau)_(q_hat)
        m.SetLength(n256 + tau0);

        for(i=0; i<n256; i++) 
        {
            B[i]   = crs[2][i];
            // t_B[i] = t_y[i];
            m[i]   = conv<zz_pX>( y_3[i] );
        }

        for(i=n256; i<(n256 + tau0); i++) 
        {
            B[i]   = crs[3][i-n256];
            // t_B[i] = t_g[i-n256];
            m[i]   = g[i-n256];
        }

        // 36. s_hat ← (s1; σ(s1); m; σ(m)),     s_hat ∈ R^^(2*m1 + 2*(256/d_hat + tau))_(q_hat)
        s_hat.SetLength( m1_n256_tau );
        for(i=0; i<m1; i++) 
        {
            s_hat[i] = s_1_mod[i];        
        }

        k = 0;

        for(i=m1; i<(2*m1); i++) 
        {
            s_hat[i] = sigma_s_1[k];  
            k++;      
        }

        k = 0;

        for(i=(2*m1); i<(2*m1 + n256 + tau0); i++) 
        {
            s_hat[i] = m[k];
            k++; 
        }

        sigma_map(tmp_vec, m, d_hat);
        k = 0;

        for(i=(2*m1 + n256 + tau0); i<(m1_n256_tau); i++)  
        {
            s_hat[i] = tmp_vec[k];    
            k++;    
        }
        
        // 37. y ← (y1; σ(y1); −B*y2; σ(-B*y2)),     y ∈ R^^(2*m1 + 2*(256/d_hat + tau))_(q_hat)
        y.SetLength( m1_n256_tau );

        for(i=0; i<m1; i++) 
        {
            y[i] = conv<zz_pX>( y_1[i] );        
        }

        sigma_map(tmp_vec, conv<vec_zz_pX>(y_1), d_hat);
        k = 0;

        for(i=m1; i<(2*m1); i++) 
        {
            y[i] = tmp_vec[k];  
            k++;      
        }
            
        // Compute B*y2 in a temporary vector
        tmp_vec.SetLength(n256 + tau0);
        acc.SetLength(d_hat);

        for(i=0; i<(n256 + tau0); i++)
        {
            tmp_vec[i].SetLength(d_hat);

            // acc = 0;
            clear(acc);

            for(j=0; j<m2; j++)        
            {
                acc += ModPhi_hat_q( B[i][j] * conv<zz_pX>( y_2[j] ) );
            }           
            
            tmp_vec[i] = acc;
            // NOTE: modulo q_hat on all coefficients (zz_pX)
        }
        
        k = 0;

        for(i=(2*m1); i<(2*m1 + n256 + tau0); i++) 
        {
            tmp_vec[k] = -tmp_vec[k]; // -B*y2
            y[i] = tmp_vec[k];        // -B*y2
            k++; 
        }

        sigma_map(tmp_vec, tmp_vec, d_hat);
        k = 0;

        for(i=(2*m1 + n256 + tau0); i<( m1_n256_tau ); i++) 
        {
            y[i] = tmp_vec[k]; // σ(-B*y2)
            k++;    
        }
        
        // 38. δ_1 ← Sum_(i=1,τ){ μ_i · γ_(i,256+d+1) },   δ_1 ∈ R^^_(q_hat)
        // 39. δ_2 ← Sum_(i=1,τ){ μ_i · γ_(i,256+d+2) },   δ_2 ∈ R^^_(q_hat)
        // 40. δ_3 ← Sum_(i=1,τ){ μ_i · γ_(i,256+d+3) },   δ_3 ∈ R^^_(q_hat)
        delta_1.SetLength(d_hat);
        delta_2.SetLength(d_hat);
        delta_3.SetLength(d_hat);
        clear(delta_1);
        clear(delta_2);
        clear(delta_3);

        for(i=0; i<tau0; i++)        
        {
            delta_1 += mu[i]*gamma[i][256+d0];
            delta_2 += mu[i]*gamma[i][256+d0+1];
            delta_3 += mu[i]*gamma[i][256+d0+2];
        }

        // 41. Definition of D_2_(2,1) ∈ R^^(m1 x m1)_(q_hat)
        // NOTE: removed D2_2_1, D2 matrix directly filled using delta_1,2,3
        
        // 42. Construction of D_2 ∈ R^^((2*m1+2(256/d_hat+τ))×(2*m1+2(256/d_hat+τ)))_(q_hat)
        D2.SetDims(m1_n256_tau, m1_n256_tau);    
        
        for(i=0; i<m1_n256_tau; i++)
        {
            for(j=0; j<m1_n256_tau; j++)        
            {
                D2[i][j].SetLength(d_hat);
                clear(D2[i][j]);
            }
        }
        
        // NOTE: D2_2_1 in position (2,1), zeros in the rest
        // NOTE: m1 = m1_ISIS = (((m+2)*d+d_hat)/d_hat) + (|idx_hid|·h + ℓr·d + d_hat)/d_hat + t/d_hat
        k = m1;
        
        for(i=0; i<((m2d+d_hat)/d_hat); i++)
        {
            D2[k][i] = delta_1;
            k++;
        }

        for(i=((m2d+d_hat)/d_hat); i<((m2d+d_hat)/d_hat + (idxhlrd+d_hat)/d_hat); i++)
        {
            D2[k][i] = delta_2;
            k++;
        }

        for(i=((m2d+d_hat)/d_hat + (idxhlrd+d_hat)/d_hat); i<m1; i++)
        {
            D2[k][i] = delta_3;
            k++;
        }


        // 43.  (r_s,j , r_r,j , r_u,j ) ← r_j,   
        //       r_s,j ∈ R^^(((m+2)d+d_hat)/d_hat)_(q_hat)
        //       r_r,j ∈ R^^((|idx_hid|·h+ℓr·d+d_hat)/d_hat)_(q_hat)
        //       r_u,j ∈ R^^(t/d_hat)_(q_hat)
        //       c_r,j is the poly. vector with coeff. the j-th row of C_r
        // NOTE: σ(r_s,j), σ(r_r,j), σ(r_u,j), σ(c_r,j) already pre-computed
        
             
        // 44. Construction of d_1 ∈ R^^(2*m1+2(256/d_hat+τ))_(q_hat)
        d_1.SetLength(m1_n256_tau);

        for(i=0; i<m1_n256_tau; i++)
        {
            d_1[i].SetLength(d_hat);
            clear(d_1[i]);
        }    

        // 1st entry of d_1: ((m+2)d+d_hat)/d_hat polynomials
        acc_vec.SetLength((m2d+d_hat)/d_hat);
        
        for(i=0; i<tau0; i++)
        {
            // Reset acc_vec
            for(j=0; j<((m2d+d_hat)/d_hat); j++)
            {        
                acc_vec[j].SetLength(d_hat);
                
                // acc_vec[j] = 0;
                clear(acc_vec[j]);
            }               
                    
            for(j=0; j<256; j++)        
            {
                for(k=0; k<((m2d+d_hat)/d_hat); k++)        
                {
                    acc_vec[k] += gamma[i][j] * sigma_r_s_[j][k];  
                }
            }

            for(j=0; j<d0; j++)        
            {
                for(k=0; k<((m2d+d_hat)/d_hat); k++)        
                {                
                    acc_vec[k] += gamma[i][256+j] * sigma_p_[j][k];  
                }     
            }  
                        
            for(k=0; k<((m2d+d_hat)/d_hat); k++)        
            {
                // Fill d_1 by accumulating mu[i]*(...sums...) 
                d_1[k] += ModPhi_hat_q( mu[i] * acc_vec[k]); 
            }               
        }

        // 2nd entry of d_1: ((|idx_hid|·h + ℓr·d + d_hat)/d_hat) polynomials
        acc_vec.SetLength((idxhlrd+d_hat)/d_hat);
        
        for(i=0; i<tau0; i++)
        {
            // Reset acc_vec
            for(j=0; j<((idxhlrd+d_hat)/d_hat); j++)
            {        
                acc_vec[j].SetLength(d_hat);
                
                // acc_vec[j] = 0;
                clear(acc_vec[j]);
            }               
                    
            for(j=0; j<256; j++)        
            {
                for(k=0; k<((idxhlrd+d_hat)/d_hat); k++)        
                {
                    acc_vec[k] += gamma[i][j] * sigma_r_r_[j][k];  
                }
            }

            for(j=0; j<d0; j++)        
            {
                for(k=0; k<((idxhlrd+d_hat)/d_hat); k++)        
                {                
                    acc_vec[k] -= gamma[i][256+j] * sigma_c_r_[j][k];
                }     
            }  
                        
            for(k=0; k<((idxhlrd+d_hat)/d_hat); k++)        
            {
                // Fill d_1 by accumulating mu[i]*(...sums...) 
                d_1[k + ((m2d+d_hat)/d_hat)] += ModPhi_hat_q( mu[i] * acc_vec[k]); 
            }               
        }

        // 3rd entry of d_1: (t/d_hat) polynomials
        acc_vec.SetLength(t_d_hat);        
        
        for(i=0; i<tau0; i++)
        {
            // Reset acc_vec
            for(j=0; j<(t_d_hat); j++)
            {        
                acc_vec[j].SetLength(d_hat);
                
                // acc_vec[j] = 0;
                clear(acc_vec[j]);
            }               
                    
            for(j=0; j<256; j++)        
            {
                for(k=0; k<(t_d_hat); k++)        
                {
                    acc_vec[k] += gamma[i][j] * sigma_r_u_[j][k];  
                }
            }

            for(j=0; j<d0; j++)        
            {
                for(k=0; k<(t_d_hat); k++)        
                {                
                    acc_vec[k] -= gamma[i][256+j] * sigma_Beta_[j][k]; 
                }     
            }
            
            for(k=0; k<(t_d_hat); k++)        
            {                
                acc_vec[k] -= gamma[i][256+d0+2] * sigma_ones[k];
            }     
                        
            for(k=0; k<(t_d_hat); k++)        
            {
                // Fill d_1 by accumulating mu[i]*(...sums...) 
                d_1[k + ((m2d+d_hat)/d_hat + (idxhlrd+d_hat)/d_hat)] += ModPhi_hat_q( mu[i] * acc_vec[k]); 
            }               
        }

        // NOTE: skip 4th entry of d_1 (m1 zeros)

        // 5th entry of d_1 (256/d_hat polynomials)    
        acc_vec.SetLength(n256);
        
        for(i=0; i<tau0; i++)
        {
            // Reset acc_vec
            for(j=0; j<n256; j++)
            {        
                acc_vec[j].SetLength(d_hat); 

                // acc_vec[j] = 0;
                clear(acc_vec[j]);
            }                     
                    
            for(j=0; j<256; j++)        
            {
                for(k=0; k<n256; k++)        
                {
                    acc_vec[k] += gamma[i][j] * sigma_e_[j][k];
                }          
            }        
                
            for(k=0; k<n256; k++)        
            {
                // Fill d_1 by accumulating mu[i]*(...sums...) 
                d_1[k + (2*m1)] += ModPhi_hat_q( mu[i] * acc_vec[k]);             
            }
        } 
        
        // 6th entry of d_1 (tau0 polynomials)
        for(k=0; k<tau0; k++)
        {
            d_1[k + (2*m1 + n256)] = mu[k];
        }

        // NOTE: skip 7th entry of d_1 (tau0 + 256/d_hat zeros)
           
        
        // // 45. Definition of d_0 ∈ R^_(q_hat)    
        // d_0.SetLength(d_hat);
        // clear(d_0);        
        // // NOTE: d_0 (not d0 parameter) despite being in the pseudocode, it is not used in Prove_ISIS         
           
        // for(i=0; i<tau0; i++)
        // {            
        //     sums = 0;
            
        //     for(j=0; j<256; j++)
        //     {
        //         sums += gamma[i][j] * conv<zz_p>(Pi.z_3[j]);
        //     }
                        
        //     for(j=0; j<d0; j++)
        //     {
        //         sums += gamma[i][256+j] * m_C[j];  
        //     }

        //     sums += gamma[i][256+d0] * B_goth_s2_p + gamma[i][256+d0+1] * B_goth_r2_p;
        
        //     d_0 = d_0 - ModPhi_hat_q( mu[i] * ( sums + Pi.h[i] ));
        // }

        
        // 46. Definition of f1 ∈ R^_(q_hat)
        f1.SetLength(d_hat);
        clear(f1);

        // 1st addend of f1, (s_hat^T * D2 * y)
        acc_vec.SetLength(m1_n256_tau);

        // Compute  (D2 * y),  (2*m1 + 2*(256/d_hat + tau0)) polynomials
        for(i=0; i<m1_n256_tau; i++)
        {
            acc_vec[i].SetLength(d_hat); 

            // acc_vec[i] = 0;
            clear(acc_vec[i]);

            acc_vec[i] = poly_mult_hat(D2[i], y);
        }

        // Accumulate  s_hat^T * (D2 * y)
        D2_y = acc_vec;
        f1 += poly_mult_hat(s_hat, acc_vec);

        // 2nd addend of f1,  (y^T * D2 * s_hat)
        acc_vec.SetLength(m1_n256_tau);

        // Compute  (D2 * s_hat),  (2*m1 + 2*(256/d_hat + tau0)) polynomials
        for(i=0; i<m1_n256_tau; i++)    
        {
            acc_vec[i].SetLength(d_hat); 

            // acc_vec[i] = 0;
            clear(acc_vec[i]);     

            acc_vec[i] = poly_mult_hat(D2[i], s_hat);
        }

        // Accumulate  y^T * (D2 * s_hat)    
        f1 += poly_mult_hat(y, acc_vec);

        // 3rd addend of f1,   (d_1^T * y)
        f1 += poly_mult_hat(d_1, y);


        // 47. Definition of f0 ∈ R^_(q_hat)
        Pi.f0.SetLength(d_hat);
        clear(Pi.f0);
        
        Pi.f0 = poly_mult_hat(y, D2_y) + poly_mult_hat(crs[4][0], conv<vec_zz_pX>(y_2));
        // NOTE: D2_y = (D2 * y) was already computed in row 45 (1st addend of f1) 
    
        
        // 48. Definition of t ∈ R^_(q_hat)
        Pi.t.SetLength(d_hat);
        clear(Pi.t);

        Pi.t = poly_mult_hat(crs[4][0], s_2_mod) + f1;


        // 49. a_4 ← (t, f0),   a_4 ∈ R^_(q_hat) x R^_(q_hat)  
        ss << Pi.t << Pi.f0;


        // 50. c ← H(4, crs, x, a1, a2, a3, a4),   c ∈ C ⊂ R^
        HISIS4(c, "4" + ss.str());


        // 51. for i ∈ {1, 2} do
        // NOTE: for simplicity, next operations are duplicated with suffixes _1 and _2

        // 52. z_i ← y_i + c*s_i,   z_i ∈ R^^(m_i)    
        Pi.z_1.SetLength(m1);
        Pi.z_2.SetLength(m2);

        c_s1.SetLength(m1);
        c_s2.SetLength(m2);
        
        for(i=0; i<m1; i++)
        {
            Pi.z_1[i].SetLength(d_hat);
            c_s1[i] = ModPhi_hat( c * s_1[i] );
            Pi.z_1[i]  = y_1[i] + c_s1[i]; 
        }
            
        for(i=0; i<m2; i++)
        {
            Pi.z_2[i].SetLength(d_hat);
            c_s2[i] = ModPhi_hat( c * s_2[i] );
            Pi.z_2[i]  = y_2[i] + c_s2[i]; 
        }
            

        // 53. b_i ← Rej(z_i, c*s_i, s_i_goth, M_i),   b_i ∈ {0, 1} 
        b1 = Rej_v_ZZX(Pi.z_1, c_s1, s1_goth, M_1);

        // NOTE: if b1 == 0, continue the while loop (skip next rows until 53, then go to row 7)
        if (b1 == 0)
        {            
            rst = 0;
            continue;
        }

        b2 = Rej_v_ZZX(Pi.z_2, c_s2, s2_goth, M_2); 

        // NOTE: if b2 == 0, continue the while loop (skip next rows until 53, then go to row 7)
        if (b2 == 0)
        {            
            rst = 0;
            continue;
        }

        
        // 55. rst ← b1*b2*b3
        rst = b1*b2*b3;        
    
    } // End of while loop (row 9)

    // 56. if rst = 1 then return π
    if (rst == 1)
    {
        // 54. π ← (t_A, t_y, t_g, w, z_3, h, t, f0, z_1, z_2)
        // Pi.t_A   = t_A;
        // Pi.t_y   = t_y;
        // Pi.t_g   = t_g;
        // Pi.w     = w;
        // Pi.z_3   = z_3;
        // Pi.h     = h;
        // Pi.t     = t;
        // Pi.f0    = f0;
        // Pi.z_1   = z_1;
        // Pi.z_2   = z_2;
        // NOTE: additional flag, to identify a valid proof
        Pi.valid = 1;
    }
    // 57. else return ⊥
    else // (rst == 0)      
    {
        // NOTE: invalid proof, other data fields in PROOF_I_t structure are empty
        Pi.valid = 0;        
    }
    
    // return Pi;
}


//==============================================================================
// Verify_ISIS  -   Verify the Commitment (Verify^HISIS_ISIS). 
//                  This function takes as input crs, x, and the proof π.
//                  It returns as output “reject” or “accept”.
// 
// Inputs:
// -  inputStr:     string containing the initial seed for crs
// -  crs:          structure crs_ISIS, generated by Hcrs from inputStr
// -  ipk:          Issuer public key
// - (P, C, mex, B_f, 
//    Bounds, aux): they correspond to x
// -  Pi:           proof (π) structure 
//  
// Output:
// -  0 or 1:       reject or accept 
//==============================================================================
long Verify_ISIS(const string& inputStr, const CRS_t& crs, const IPK_t& ipk, const mat_zz_p& P, const mat_zz_p& C, const vec_zz_p& mex, const mat_zz_p& B_f, const vec_ZZ& Bounds, const ZZ& aux, const PROOF_I_t& Pi)
{
    // NOTE: assuming that current modulus is q2_hat (not q0) 
   
    unsigned long        i, j, k;
    Mat<zz_pX>           B, D2; //D2_2_1
    vec_zz_pX            t_B, z, d_1, acc_vec, coeffs_ones, sigma_ones;
    vec_zz_pX            r_j, p_j, Beta_j, c_r_j, mu, tmp_vec, tmp_vec2, z_1_mod, z_2_mod;
    ZZX                  c0;
    zz_pX                acc, delta_1, delta_2, delta_3, c, d_0;
    Vec<vec_zz_pX>       e_, sigma_e_, sigma_p_, sigma_Beta_, sigma_c_r_;
    Vec<vec_zz_pX>       sigma_r_, sigma_r_s_, sigma_r_r_, sigma_r_u_;    
    stringstream         ss;
    mat_L                R_goth;
    mat_zz_p             gamma, C_m, C_r;
    vec_zz_p             ones, e_tmp, m_C;
    RR                   B_goth;    
    ZZ                   B_goth_s2,   B_goth_r2;
    zz_p                 B_goth_s2_p, B_goth_r2_p, sums;
    ZZ                   norm2_z1, norm2_z2, norm2_z3;  
    RR                   norm_z1,  norm_z2,  norm_z3;

    // Initialise constants
    const unsigned long n           = n_ISIS;
    const unsigned long m1          = m1_ISIS;
    const unsigned long m2          = m2_ISIS;
    const unsigned long m2d         = (m0 + 2)*d0;    // (m+2)·d
    const unsigned long idxhlrd     = (idx_hid * h0) + (lr0 * d0); //|idx_hid|·h + ℓr·d    
    const unsigned long n256        = (256/d_hat);
    const unsigned long t_d_hat     = (t0/d_hat);
    const unsigned long m1_n256_tau = 2*m1 + 2*(n256 + tau0);
    
    // Initialise the "goth" constants
    // (B_goth_s^2, B_goth_r^2) ← Bounds,   (B_goth_s^2, B_goth_r^2) ∈ Z^2
    B_goth_s2 = Bounds[0];
    B_goth_r2 = Bounds[1];

    // B_goth = sqrt(B_goth_s^2 + B_goth_r^2 + t0)
    B_goth = sqrt( conv<RR>( B_goth_s2 + B_goth_r2 + t0 ) );

    B_goth_s2_p = conv<zz_p>(B_goth_s2);
    B_goth_r2_p = conv<zz_p>(B_goth_r2);

    // s1_goth = alpha_1*nu0*B_goth
    // s2_goth = alpha_2*nu0*sqrt(m2*d_hat)
    // s3_goth = alpha_3*w_max(lambda0)*B_goth
    const RR  s1_goth  = conv<RR>( alpha_1 * nu0   ) * B_goth;
    const RR  s2_goth  = conv<RR>( alpha_2 * nu0   ) * sqrt( m2*d_hat );
    const RR  s3_goth  = conv<RR>( alpha_3 * w_max ) * B_goth;
    const RR  B_goth_1 = s1_goth * sqrt( 2*m1*d_hat );
    const RR  B_goth_2 = s2_goth * sqrt( 2*m2*d_hat );
    const RR  B_goth_3 = 1.7 * s3_goth * sqrt( 256 );


    // 1. Retrieve (A_1, A_2, B_y, B_g, b) ← crs_ISIS
    // A_1     = crs[0];    // ∈ R^^(n x m1)_(q_hat)
    // A_2     = crs[1];    // ∈ R^^(n x m2)_(q_hat)
    // B_y     = crs[2];    // ∈ R^^(256/d^ x m2)_(q_hat)
    // B_g     = crs[3];    // ∈ R^^(tau0 x m2)_(q_hat)
    // b       = crs[4][0]; // ∈ R^^(m2)_(q_hat)         NOTE: b in crs is a (1 x m_2) matrix

    // 2. (P, C, mex, B_f, Bounds, aux) ← x
    // NOTE: directly provided as inputs
            
    // 3. P ← [P0,  0_(d × d_hat)],   P ∈ Z^[d × (m+2)·d + d_hat]_(q_hat)     
    // NOTE: zero padding of P already done in V_Verify
 
    // 4. C ← [C0,  0_(d × d_hat)],   C ∈ Z^[d × ((ℓm+ℓr)·d + d_hat)]_(q_hat)
    // NOTE: zero padding of C already done in V_Verify

    // 5. (t_A, t_y, t_g, w, z_3, h, t, f0, z_1, z_2) ← π
    // Check if Pi contains a valid proof   
    if (Pi.valid != 1)
    {
        cout << "ERROR! Pi does not contain a valid proof" << endl;
        return 0;
    }
    // NOTE: to save memory, proof values will be directly accessed as Pi.{name},
    //       apart z_1 and z_2, that are often used modulo q2_hat      
    z_1_mod = conv<vec_zz_pX>( Pi.z_1 );
    z_2_mod = conv<vec_zz_pX>( Pi.z_2 );
    
    // 6. a_1 ← (t_A, t_y, t_g, w)     
    // a_1 << Pi.t_A << Pi.t_y << Pi.t_g << Pi.w;
    
    // 7. a_2 ← z_3,   a2 ∈ Z^256    
    // a_2 << Pi.z_3;
    
    // 8. a_3 ← h,   a_3 ∈ R^^(tau)_(q_hat)    
    // a_3 << Pi.h;
    
    // 9. a_4 ← (t, f0),   a_4 ∈ R^_(q_hat) x R^_(q_hat)  
    // a_4 << Pi.t << Pi.f0;
    
    // 10. (R_goth_0, R_goth_1) = H(1, crs, x, a_1)
    // ss << crs << P << C << mex << B_f << Bounds << aux << Pi.t_A << Pi.t_y << Pi.t_g << Pi.w;
    ss << inputStr << ipk.a1 << ipk.a2 << ipk.c0 << ipk.c1 << mex << B_f << Bounds << aux << Pi.t_A << Pi.t_y << Pi.t_g << Pi.w;
    // NOTE: using inputStr, ipk, instead of crs, P, C to speedup Hash_Init
    HISIS1(R_goth, "1" + ss.str());
    // NOTE: R_goth ∈ {-1, 0, 1}^(256 x m_1*d_hat),
    //       equivalent to (R_goth_0 - R_goth_1) in BLNS

    // 11. gamma ← H(2, crs, x, a1, a2),   gamma ∈ Z^(tau0 x 256+d0+3)_q_hat
    ss << Pi.z_3;
    HISIS2(gamma, "2" + ss.str());    
    // NOTE: gamma has 256+d0+3 columns in ISIS, while 256+d0+1 in Com

    // 12. μ ← H(3, crs, x, a1, a2, a3),   μ ∈ R^^(τ)_(q_hat)
    ss << Pi.h;
    HISIS3(mu, "3" + ss.str());

    // 13. c ← H(4, crs, x, a1, a2, a3, a4),   c ∈ C ⊂ R^_(q_hat)
    ss << Pi.t << Pi.f0;
    HISIS4(c0, "4" + ss.str());
    c = conv<zz_pX>( c0 );
    // NOTE: Verify_ISIS only uses c mod q_hat 

    // 14. B   ← [B_y; B_g],   B ∈ R^^((256/d_hat + tau) x m2)_(q_hat)
    B.SetDims((n256 + tau0), m2);

    // 15. t_B ← [t_y; t_g],   t_B ∈ R^^(256/d_hat + tau)_(q_hat)
    t_B.SetLength(n256 + tau0);    

    for(i=0; i<n256; i++) 
    {
        B[i]   = crs[2][i];
        t_B[i] = Pi.t_y[i];
    }

    for(i=n256; i<(n256 + tau0); i++) 
    {
        B[i]   = crs[3][i-n256];
        t_B[i] = Pi.t_g[i-n256];
    }


    // 16. z ← (z_1; σ(z_1); (c*t_B − B*z_2); σ(c*t_B − B*z_2)),     z ∈ R^^(2*m1 + 2*(256/d_hat + tau))_(q_hat)
    z.SetLength( m1_n256_tau ); 

    for(i=0; i<m1; i++) 
    {
        z[i] = z_1_mod[i]; // z_1
    }

    sigma_map(tmp_vec, z_1_mod, d_hat);
    k = 0;

    for(i=m1; i<(2*m1); i++) 
    {
        z[i] = tmp_vec[k]; // σ(z_1)
        k++;      
    }

    // Compute (c*t_B − B*z_2) in a temporary vector
    tmp_vec.SetLength(n256 + tau0);
    acc.SetLength(d_hat);

    for(i=0; i<(n256 + tau0); i++)
    {
        tmp_vec[i].SetLength(d_hat);

        // acc = 0;
        clear(acc);

        // c*t_B
        acc = ModPhi_hat_q( c * t_B[i] );

        // − B*z_2
        for(j=0; j<m2; j++)        
        {
            acc += ModPhi_hat_q(- B[i][j] * z_2_mod[j] );
        }           
        
        tmp_vec[i] = acc;
        // NOTE: modulo q_hat on all coefficients (zz_pX)
    }
    
    k = 0;

    for(i=(2*m1); i<(2*m1 + n256 + tau0); i++) 
    {
        z[i] = tmp_vec[k]; // (c*t_B − B*z_2)
        k++; 
    }

    // Compute σ(c*t_B − B*z_2) in a temporary vector
    sigma_map(tmp_vec2, tmp_vec, d_hat); 
    k = 0;

    for(i=(2*m1 + n256 + tau0); i<(m1_n256_tau); i++)  
    {
        z[i] = tmp_vec2[k]; // σ(c*t_B − B*z_2)   
        k++;    
    }


    // 17. δ_1 ← Sum_(i=1,τ){ μ_i · γ_(i,256+d+1) },   δ_1 ∈ R^^_(q_hat)
    // 18. δ_2 ← Sum_(i=1,τ){ μ_i · γ_(i,256+d+2) },   δ_2 ∈ R^^_(q_hat)
    // 19. δ_3 ← Sum_(i=1,τ){ μ_i · γ_(i,256+d+3) },   δ_3 ∈ R^^_(q_hat)
    delta_1.SetLength(d_hat);
    delta_2.SetLength(d_hat);
    delta_3.SetLength(d_hat);
    clear(delta_1);
    clear(delta_2);
    clear(delta_3);

    for(i=0; i<tau0; i++)        
    {
        delta_1 += mu[i]*gamma[i][256+d0];
        delta_2 += mu[i]*gamma[i][256+d0+1];
        delta_3 += mu[i]*gamma[i][256+d0+2];
    }

    // 20. Definition of D_2_(2,1) ∈ R^^(m1 x m1)_(q_hat)
    // NOTE: removed D2_2_1, D2 matrix directly filled using delta_1,2,3
    
    // 21. Construction of D_2 ∈ R^^((2*m1+2(256/d_hat+τ))×(2*m1+2(256/d_hat+τ)))_(q_hat)
    D2.SetDims(m1_n256_tau, m1_n256_tau);    
    
    for(i=0; i<m1_n256_tau; i++)
    {
        for(j=0; j<m1_n256_tau; j++)        
        {
            D2[i][j].SetLength(d_hat);
            clear(D2[i][j]);
        }
    }
    
    // NOTE: D2_2_1 in position (2,1), zeros in the rest
    // NOTE: m1 = m1_ISIS = (((m+2)*d+d_hat)/d_hat) + (|idx_hid|·h + ℓr·d + d_hat)/d_hat + t/d_hat
    k = m1;

    for(i=0; i<((m2d+d_hat)/d_hat); i++)
    {
        D2[k][i] = delta_1;
        k++;
    }

    for(i=((m2d+d_hat)/d_hat); i<((m2d+d_hat)/d_hat + (idxhlrd+d_hat)/d_hat); i++)
    {
        D2[k][i] = delta_2;
        k++;
    }

    for(i=((m2d+d_hat)/d_hat + (idxhlrd+d_hat)/d_hat); i<m1; i++)
    {
        D2[k][i] = delta_3;
        k++;
    }


    // 22.  (r_s,j , r_r,j , r_u,j ) ← r_j,
    //       r_s,j ∈ R^^(((m+2)d+d_hat)/d_hat)_(q_hat)
    //       r_r,j ∈ R^^((|idx_hid|·h+ℓr·d+d_hat)/d_hat)_(q_hat)
    //       r_u,j ∈ R^^(t/d_hat)_(q_hat)
    //       c_r,j is the poly. vector with coeff. the j-th row of C_r
    
    // Precompute r_j, σ(r_j), σ(r_s,j), σ(r_r,j), σ(r_u,j)      
    r_j.SetLength(m1);        
    sigma_r_.SetLength(256);
    sigma_r_s_.SetLength(256);
    sigma_r_r_.SetLength(256);
    sigma_r_u_.SetLength(256);

    for(j=0; j<256; j++)        
    {
        CoeffsInvHat(r_j,  conv<vec_zz_p>(R_goth[j]), m1 );
        sigma_map(sigma_r_[j], r_j, d_hat);
        
        // NOTE: m1 = m1_ISIS = (((m+2)*d+d_hat)/d_hat) + (|idx_hid|·h + ℓr·d + d_hat)/d_hat + t/d_hat

        sigma_r_s_[j].SetLength((m2d+d_hat)/d_hat);
        sigma_r_r_[j].SetLength((idxhlrd+d_hat)/d_hat);
        sigma_r_u_[j].SetLength(t_d_hat);

        for(k=0; k<((m2d+d_hat)/d_hat); k++)
        {
            sigma_r_s_[j][k] = sigma_r_[j][k];
        }

        for(k=0; k<((idxhlrd+d_hat)/d_hat); k++)
        {
            sigma_r_r_[j][k] = sigma_r_[j][k + (m2d+d_hat)/d_hat];
        }

        for(k=0; k<(t_d_hat); k++)
        {
            sigma_r_u_[j][k] = sigma_r_[j][k + (m2d+d_hat)/d_hat + (idxhlrd+d_hat)/d_hat];
        }
    }
    
    // Initialize e ∈ R^^(256 x 256/d_hat)_(q_hat)
    e_.SetLength(256);    
    // NOTE: defined as e_ to distinguish it from the Euler constant e 

    e_tmp.SetLength(256);

    for(k=0; k<256; k++)
    {
        e_tmp[k] = 0;
    }

    for(j=0; j<256; j++)
    {
        // Temporary coefficient vector to create e_j: it is a unit vector with its j-th coefficient being 1
        e_tmp[j] = 1;        

        // e_[j].SetLength(n256);
        CoeffsInvHat(e_[j], e_tmp, n256);

        // Reset the e_tmp coefficient vector
        e_tmp[j] = 0;
    }
    
    // Precompute σ(e_j), σ(p_j), σ(β_j)
    sigma_e_.SetLength(256);
    sigma_p_.SetLength(d0);
    sigma_Beta_.SetLength(d0);
        
    for(j=0; j<256; j++)        
    {
        sigma_map(sigma_e_[j], e_[j], d_hat);  
    }       

    for(j=0; j<d0; j++)        
    {
        CoeffsInvHat(p_j    , P[j], (m2d+d_hat)/d_hat);
        CoeffsInvHat(Beta_j , B_f[j], t_d_hat); 
        sigma_map(sigma_p_[j]    , p_j, d_hat);        
        sigma_map(sigma_Beta_[j] , Beta_j, d_hat);                  
    }

    // Create C_m and C_r
    C_m.SetDims(d0, (idx_pub*h0));
    C_r.SetDims(d0, (idxhlrd+d_hat)); 
    // NOTE: C = [C_m C_r] ∈ Z^[d × ((ℓm+ℓr)·d+d_hat)]_(q_hat), C_m has |idx_pub|·h columns
    
    for(i=0; i<d0; i++)
    {
        for(j=0; j<(idx_pub*h0); j++)        
        {
            C_m[i][j] = C[i][j];
        }

        for(j=0; j<(idxhlrd+d_hat); j++)          
        {
            C_r[i][j] = C[i][(idx_pub*h0)+j];
        }
    }

    // Compute m_C := C_m * m
    m_C.SetLength(d0);
    m_C = ( C_m * mex );

    C_m.kill(); 

    // Precompute σ(C_r,j)    
    sigma_c_r_.SetLength(d0);

    for(j=0; j<d0; j++)        
    {
        CoeffsInvHat(c_r_j  , C_r[j], (idxhlrd+d_hat)/d_hat);
        sigma_map(sigma_c_r_[j] , c_r_j, d_hat);                  
    }

    C_r.kill(); 

    // Precompute Coeffs^−1(1^t)
    ones.SetLength(t0);             // ∈ Z^^(t)_(q_hat)
    coeffs_ones.SetLength(t_d_hat); // ∈ R^^(t/d_hat)_(q_hat)
    
    for(i=0; i<t0; i++) 
    {
        ones[i] = 1;
    }

    CoeffsInvHat(coeffs_ones, ones, t_d_hat);
    sigma_map(sigma_ones , coeffs_ones, d_hat);


    // 23. Construction of d_1 ∈ R^^(2*m1+2(256/d_hat+τ))_(q_hat)
    d_1.SetLength(m1_n256_tau);

    for(i=0; i<m1_n256_tau; i++)
    {
        d_1[i].SetLength(d_hat);
        clear(d_1[i]);
    }    

    // 1st entry of d_1: ((m+2)d+d_hat)/d_hat polynomials
    acc_vec.SetLength((m2d+d_hat)/d_hat);
    
    for(i=0; i<tau0; i++)
    {
        // Reset acc_vec
        for(j=0; j<((m2d+d_hat)/d_hat); j++)
        {        
            acc_vec[j].SetLength(d_hat);
            
            // acc_vec[j] = 0;
            clear(acc_vec[j]);
        }               
                
        for(j=0; j<256; j++)        
        {
            for(k=0; k<((m2d+d_hat)/d_hat); k++)        
            {
                acc_vec[k] += gamma[i][j] * sigma_r_s_[j][k];  
            }
        }

        for(j=0; j<d0; j++)        
        {
            for(k=0; k<((m2d+d_hat)/d_hat); k++)        
            {                
                acc_vec[k] += gamma[i][256+j] * sigma_p_[j][k];  
            }     
        }  
                    
        for(k=0; k<((m2d+d_hat)/d_hat); k++)        
        {
            // Fill d_1 by accumulating mu[i]*(...sums...) 
            d_1[k] += ModPhi_hat_q( mu[i] * acc_vec[k]); 
        }               
    }

    // 2nd entry of d_1: ((|idx_hid|·h + ℓr·d + d_hat)/d_hat) polynomials
    acc_vec.SetLength((idxhlrd+d_hat)/d_hat);
    
    for(i=0; i<tau0; i++)
    {
        // Reset acc_vec
        for(j=0; j<((idxhlrd+d_hat)/d_hat); j++)
        {        
            acc_vec[j].SetLength(d_hat);
            
            // acc_vec[j] = 0;
            clear(acc_vec[j]);
        }               
                
        for(j=0; j<256; j++)        
        {
            for(k=0; k<((idxhlrd+d_hat)/d_hat); k++)        
            {
                acc_vec[k] += gamma[i][j] * sigma_r_r_[j][k];  
            }
        }

        for(j=0; j<d0; j++)        
        {
            for(k=0; k<((idxhlrd+d_hat)/d_hat); k++)        
            {                
                acc_vec[k] -= gamma[i][256+j] * sigma_c_r_[j][k];
            }     
        }  
                    
        for(k=0; k<((idxhlrd+d_hat)/d_hat); k++)        
        {
            // Fill d_1 by accumulating mu[i]*(...sums...) 
            d_1[k + ((m2d+d_hat)/d_hat)] += ModPhi_hat_q( mu[i] * acc_vec[k]); 
        }               
    }

    // 3rd entry of d_1: (t/d_hat) polynomials
    acc_vec.SetLength(t_d_hat);        
    
    for(i=0; i<tau0; i++)
    {
        // Reset acc_vec
        for(j=0; j<(t_d_hat); j++)
        {        
            acc_vec[j].SetLength(d_hat);
            
            // acc_vec[j] = 0;
            clear(acc_vec[j]);
        }               
                
        for(j=0; j<256; j++)        
        {
            for(k=0; k<(t_d_hat); k++)        
            {
                acc_vec[k] += gamma[i][j] * sigma_r_u_[j][k];  
            }
        }

        for(j=0; j<d0; j++)        
        {
            for(k=0; k<(t_d_hat); k++)        
            {                
                acc_vec[k] -= gamma[i][256+j] * sigma_Beta_[j][k]; 
            }     
        }
        
        for(k=0; k<(t_d_hat); k++)        
        {                
            acc_vec[k] -= gamma[i][256+d0+2] * sigma_ones[k];
        }     
                    
        for(k=0; k<(t_d_hat); k++)        
        {
            // Fill d_1 by accumulating mu[i]*(...sums...) 
            d_1[k + ((m2d+d_hat)/d_hat + (idxhlrd+d_hat)/d_hat)] += ModPhi_hat_q( mu[i] * acc_vec[k]); 
        }               
    }

    // NOTE: skip 4th entry of d_1 (m1 zeros)

    // 5th entry of d_1 (256/d_hat polynomials)    
    acc_vec.SetLength(n256);
    
    for(i=0; i<tau0; i++)
    {
        // Reset acc_vec
        for(j=0; j<n256; j++)
        {        
            acc_vec[j].SetLength(d_hat); 

            // acc_vec[j] = 0;
            clear(acc_vec[j]);
        }                     
                
        for(j=0; j<256; j++)        
        {
            for(k=0; k<n256; k++)        
            {
                acc_vec[k] += gamma[i][j] * sigma_e_[j][k];
            }          
        }        
            
        for(k=0; k<n256; k++)        
        {
            // Fill d_1 by accumulating mu[i]*(...sums...) 
            d_1[k + (2*m1)] += ModPhi_hat_q( mu[i] * acc_vec[k]);             
        }
    } 
    
    // 6th entry of d_1 (tau0 polynomials)
    for(k=0; k<tau0; k++)
    {
        d_1[k + (2*m1 + n256)] = mu[k];
    }

    // NOTE: skip 7th entry of d_1 (tau0 + 256/d_hat zeros)


    // 24. Definition of d_0 ∈ R^_(q_hat)    
    d_0.SetLength(d_hat);
    clear(d_0);        
    // NOTE: d_0 (not d0 parameter) 
        
    for(i=0; i<tau0; i++)
    {            
        sums = 0;
        
        for(j=0; j<256; j++)
        {
            sums += gamma[i][j] * conv<zz_p>(Pi.z_3[j]);
        }
                    
        for(j=0; j<d0; j++)
        {
            sums += gamma[i][256+j] * m_C[j];  
        }

        sums += gamma[i][256+d0] * B_goth_s2_p + gamma[i][256+d0+1] * B_goth_r2_p;
    
        d_0 = d_0 - ModPhi_hat_q( mu[i] * ( sums + Pi.h[i] ));
    }


    // 25.  if one of the 4 conditions below does not hold, then return 0
       
    // Compute ||z_i||, Euclidean norm of each z_i
    norm_z1 = sqrt( conv<RR>( Norm2X(Pi.z_1, d_hat ) ) ); 
    norm_z2 = sqrt( conv<RR>( Norm2X(Pi.z_2, d_hat ) ) ); 
    norm_z3 = sqrt( conv<RR>( Norm2( Pi.z_3) ) );

    
    // 25.1 First condition: ||z_1|| ≤ B_goth_1, ||z_2|| ≤ B_goth_2, ||z_3|| ≤ B_goth_3
    // NOTE: equations in RR  
    if ( norm_z1 > B_goth_1)
    { 
        cout << "First condition failed - Invalid z_1 norm!" << endl; 
        return 0;
    }

    if ( norm_z2 > B_goth_2)
    { 
        cout << "First condition failed - Invalid z_2 norm!" << endl; 
        return 0;
    }

    if ( norm_z3 > B_goth_3)
    { 
        cout << "First condition failed - Invalid z_3 norm!" << endl; 
        return 0;
    }


    // 25.2 Second condition: h˜_i == 0 for i ∈ [τ] 
    // NOTE: equations in R^^_(q_hat)
    for(i=0; i<tau0; i++)
    {
        if ( coeff(Pi.h[i], 0) != 0 )
        {
            cout << "Second condition failed! \n  h = " << Pi.h << endl;
            return 0;
        }
    }


    // 25.3 Third condition: A_1*z_1 + A_2*z_2 == w + c*t_A 
    // NOTE: equations in R^^(n)_(q_hat)    
    tmp_vec.SetLength(n);
    tmp_vec2.SetLength(n);
    acc.SetLength(d_hat);

    for(i=0; i<n; i++)
    {
        tmp_vec[i].SetLength(d_hat);
        tmp_vec2[i].SetLength(d_hat);

        // acc = 0;
        clear(acc);  
        
        for(j=0; j<m1; j++)
        {
            acc += ModPhi_hat_q( crs[0][i][j] * z_1_mod[j] ); 
        }        

        for(j=0; j<m2; j++)
        {
            acc += ModPhi_hat_q( crs[1][i][j] * z_2_mod[j] ); 
        }  

        // A_1*z_1 + A_2*z_2
        tmp_vec[i] = acc;

        // w + c*t_A 
        tmp_vec2[i] = Pi.w[i] + ModPhi_hat_q( c * Pi.t_A[i] );
    }

    if (tmp_vec != tmp_vec2)
    {
        cout << "Third condition failed!" << endl; 
        return 0;
    }


    // 25.4 Fourth condition: z^T*D2*z + c*d_1^T*z + c^2*d_0 − (c*t − b^T*z_2) == f0 
    // NOTE: equations in R^^_(q_hat)
    acc.SetLength(d_hat);
    // acc = 0;
    clear(acc);  

    // 1st addend (z^T * D2 * z)  
    // Compute  (D2 * z),  (2*m1 + 2*(256/d_hat + tau0)) polynomials
    acc_vec.SetLength(m1_n256_tau);
    
    for(i=0; i<m1_n256_tau; i++)    
    {
        acc_vec[i].SetLength(d_hat); 

        // acc_vec[i] = 0;
        clear(acc_vec[i]);

        acc_vec[i] = poly_mult_hat(D2[i], z);
    }

    // Accumulate  z^T * (D2 * z)    
    acc += poly_mult_hat(z, acc_vec); 
    
    // 2nd addend (c * d_1^T * z)   
    // Compute  (c * d_1^T),  (2*m1 + 2*(256/d_hat + tau0)) polynomials
    for(i=0; i<m1_n256_tau; i++)    
    {
        // acc_vec[i] = 0;
        clear(acc_vec[i]);

        acc_vec[i] = ModPhi_hat_q( c * d_1[i] );
    }

    // Accumulate  (c * d_1^T) * z 
    acc += poly_mult_hat(acc_vec, z);   
    
    // 3rd addend (c^2 * d_0)
    acc += ModPhi_hat_q( ModPhi_hat_q( sqr(c) ) * d_0 );
      
    // 4rd addend −(c*t − b^T * z_2)
    acc -= ( ModPhi_hat_q( c * Pi.t ) - poly_mult_hat(crs[4][0], z_2_mod) );
    
    if (acc != Pi.f0)
    {
        cout << "Fourth condition failed!" << endl; 
        cout << acc << " != " << Pi.f0 << endl; 
        return 0;
    }
    
    // cout << "# Verify_ISIS: OK!" << endl;

    // 26. else, return 1
    return 1;
}
