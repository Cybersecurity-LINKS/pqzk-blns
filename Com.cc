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

#include "Com.h"


//==============================================================================
// Preprocessing_Com -  Preprocessing function (PreprocessingProve^HCom_Com). 
//                      It takes as input (P, s, B_goth) with B_goth ≥ ||s||. 
//                      It returns s such that B_goth = ||s|| 
//                      and P filled with the appropriate number of zeros.
// 
// Inputs:
// - P:             matrix P  ∈ Z^[d x (|idx_hid|·h + ℓr·d)]_(q_hat)
// - s:             vector s  ∈ Z^(|idx_hid|·h + ℓr·d)_(q_hat)
// - B_goth2:       bound  B_goth^2 ∈ Z≥0 (it is a scalar)
//  
// Output:
// - P1:            matrix P1 ∈ Z^[d x (|idx_hid|·h + ℓr·d + d_hat]_(q_hat)
// - s0:            vector s0 ∈ Z^(|idx_hid|·h + ℓr·d + d_hat)_(q_hat)
//==============================================================================
// NOTE: zero padding of P, s already done in H_VerCred1
void  Preprocessing_Com(vec_ZZ& s0, const ZZ& B_goth2)
{    
    // NOTE: assuming that current modulus is q1_hat (not q0)
    long    i;//j;
    ZZ      diff; //a1, a2, a3, a4;
    vec_ZZ  a;
    
    const long  idxhlrd = (idx_hid * h0) + (lr0 * d0); //|idx_hid|·h + ℓr·d
    
    // diff = B_goth^2 − ||s||^2
    diff = (B_goth2 - Norm2(s0));
    
    if (diff < 0)
    {
        cout << "ERROR! (B_goth^2 - ||s||^2) must be positive: " << diff << endl;
        assert(diff >= 0);
    }
    else if (diff == 0)
    {
        cout << "WARNING! (B_goth^2 - ||s||^2) == 0 " << endl;
    }

    // diff = diff % q1_hat;
    // // NOTE: diff mod q1_hat, to speed up sum_of_four_squares


    // // 1. (a1, a2, a3, a4) ← SumOfFourSquares(B_goth^2 − ||s||^2),   (a1, a2, a3, a4) ∈ Z^4
    // sum_of_four_squares(a1, a2, a3, a4, diff);

    // // 2. a ← (a1, a2, a3, a4, 0, ... , 0),  a ∈ Z^(d_hat)
    // // NOTE: add d_hat − 4 zeros

    // 3. s0 ← (s, a),   s0 ∈ Z^(|idx_hid|·h + ℓr·d + d_hat)
    // s0.SetLength(idxhlrd + d_hat);

    // for(i=0; i<idxhlrd; i++)
    // {
    //     s0[i] = s[i];
    // }

    // s0[idxhlrd]   = a1;
    // s0[idxhlrd+1] = a2;
    // s0[idxhlrd+2] = a3;
    // s0[idxhlrd+3] = a4;

    // NOTE: sum_of_four_squares (too slow) replaced with fast_sum_of_squares
    fast_sum_of_squares(a, diff);
    
    for(i=0; i<a.length(); i++)
    {
        s0[idxhlrd+i] = a[i];
    }


    // 4. P1 ← [P,  0_(d × d_hat)],   P1 ∈ Z^[d x (|idx_hid|·h + ℓr·d + d_hat]_(q_hat)
    // P1.SetDims(d0, (idxhlrd + d_hat));

    // for(i=0; i<d0; i++)
    //     {   
    //     for(j=0; j<idxhlrd; j++)
    //     {
    //         P1[i][j] = P[i][j];
    //     }
    // }  

    // 5. return (P1, s0)
}


//==============================================================================
// Prove_Com   -    Computes the Commitment (Prove^HCom_Com). 
//                  This function takes as input crs, (P, u, B_goth), and w.
//                  If all the checks pass, it returns the proof π.
// 
// Inputs:
// - inputStr:      string containing the initial seed for crs
// - crs:           structure crs_Com, generated by Hcrs from inputStr
// - ipk:           Issuer public key
// - (P,u,B_goth2): this triplet corresponds to x, with:
//                  P ∈ Z^[d x (m1*d_hat)]_(q_hat)
//                  u ∈ Z^(d)_(q_hat)
//                  B_goth^2 ∈ Z≥0 (it is a scalar)
// - w0:            it contains the vector s ∈ Z^(m1*d_hat)
//                  NOTE: (|idx_hid|*h + l_r*d + d_hat) == (m1 * d_hat)
//  
// Output:
// - Pi:            proof (π) structure 
//==============================================================================
void Prove_Com(PROOF_C_t& Pi, const string& inputStr, const CRS_t& crs, const IPK_t& ipk, const mat_zz_p& P, const vec_zz_p& u0, const ZZ& B_goth2, const vec_ZZ& w0)
{
    // NOTE: assuming that current modulus is q1_hat (not q0)

    unsigned long       i, j, k, idx;
    long                rst, b1, b2, b3, bbar1, bbar2;
    vec_ZZ              s0;
    Mat<zz_pX>          B, D2_2_1;
    vec_zz_pX           u, g;
    vec_zz_pX           h_part1, h_part2;
    vec_zz_pX           r_j, p_j, mu, tmp_vec, y;
    vec_zz_pX           d_1, acc_vec, sigma_s_1, D2_y, sigma_y_1;
    vec_zz_pX           s_1, s_2, y_1, y_2, y_3;
    vec_zz_pX           c_s1, c_s2;
    zz_pX               c, h_part3, acc, sum, f1;
    RR                  alpha_i;
    Mat<zz_pX>          e_, e_prime;
    Mat<zz_pX>          sigma_r_, sigma_p_, sigma_e_, sigma_e_prime_;
    LHC_ST_t            st_1, st_2;
    stringstream        ss;
    mat_zz_p            R_goth, gamma;
    vec_zz_p            coeffs_s0, e_tmp, coeffs_R_goth_mult_s1, coeffs_y3;
    HASH_STATE_t       *state0, *state;

    // Initialise constants    
    const unsigned long n           = n_Com;
    const unsigned long m1          = m1_Com;
    const unsigned long m2          = m2_Com;
    const unsigned long d_d_hat     = (d0/d_hat);
    const unsigned long n256        = (256/d_hat);
    const unsigned long m1_n256_tau = 2*m1 + 2*(n256 + tau0);

    // Initialise the "goth" constants
    const RR  B_goth = sqrt(conv<RR>(B_goth2));
    // s1_goth = alpha_1*nu0*B_goth
    // s2_goth = alpha_2*nu0*sqrt(m2*d_hat)
    // s3_goth = alpha_3*w_max(lambda0)*B_goth    
    const RR  s1_goth = RR(alpha_1 * nu0) * B_goth;
    const RR  s2_goth = RR(alpha_2 * nu0) * sqrt( RR(m2*d_hat) );
    const RR  s3_goth = RR(alpha_3 * w_max) * B_goth;
    const double s1_goth_d = conv<double>(s1_goth);
    const double s2_goth_d = conv<double>(s2_goth);
    const double s3_goth_d = conv<double>(s3_goth);
    
    // M1 := exp(sqrt(2(λ+1)/log e) * 1/α_1 + 1/2α_1^2
    // M2 := exp(sqrt(2(λ+1)/log e) * 1/α_2 + 1/2α_2^2
    // M3 := exp(sqrt(2(λ+1)/log e) * 1/α_3 + 1/2α_3^2    
    alpha_i = RR(alpha_1);
    const RR  M_1 = exp( sqrt( RR(2*(lambda0 + 1)) / log2e_Const ) * 1/alpha_i + 1/(2*sqr(alpha_i)));
    alpha_i = RR(alpha_2);
    const RR  M_2 = exp( sqrt( RR(2*(lambda0 + 1)) / log2e_Const ) * 1/alpha_i + 1/(2*sqr(alpha_i)));
    alpha_i = RR(alpha_3);
    const RR  M_3 = exp( sqrt( RR(2*(lambda0 + 1)) / log2e_Const ) * 1/alpha_i + 1/(2*sqr(alpha_i)));


    // 1. Retrieve from crs_Com
    // A_1     = crs[0];    // ∈ (R^_q^)^(n x m1)
    // A_2     = crs[1];    // ∈ (R^_q^)^(n x m2)
    // B_y     = crs[2];    // ∈ (R^_q^)^(256/d^ x m2)
    // B_g     = crs[3];    // ∈ (R^_q^)^(tau0 x m2)
    // b       = crs[4][0]; // ∈ (R^_q^)^(m2)         NOTE: b in crs is a (1 x m_2) matrix
    // Abar_1  = crs[5];    // ∈ (R^_q^)^(m1 x n1)
    // Abar_2  = crs[6];    // ∈ (R^_q^)^(m2 x n2)
    // Bbar_1  = crs[7];    // ∈ (R^_q^)^(m1 x n1)
    // Bbar_2  = crs[8];    // ∈ (R^_q^)^(m2 x n2)

    // 2. (P, u, B_goth) ← x
    // NOTE: directly provided as inputs
    
    // 3. s ← w
    s0 = w0;    
    
    // 4. (P, s) ← PreprocessingProve^HCom_Com (P, s, B_goth)
    // P ∈ Z^[d x (|idx_hid|·h + ℓr·d + d_hat]_(q_hat)
    // s ∈ Z^(|idx_hid|·h + ℓr·d + d_hat)_(q_hat)
    Preprocessing_Com(s0, B_goth2);
    coeffs_s0 = conv<vec_zz_p>(s0);
    
    // 5. Initialize rst ∈ Z, scalar
    rst     = 0;

    // 6. Initialize idx ∈ N, scalar
    idx     = 0;

    // 7.1 Convert vector s0 into polynomial vector s_1 ∈ R^^(m1)_(q_hat)
    CoeffsInvHat(s_1, coeffs_s0, m1);

    // 7.2 Convert vector u0 into polynomial vector u ∈ R^^(d/d_hat)_(q_hat)
    CoeffsInvHat(u, u0, d_d_hat);


    // Initialize e ∈ R^^(256 x 256/d_hat)_(q_hat)
    e_.SetDims(256, n256);
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

    // Initialize e_prime ∈ R^^(d0 x d0/d_hat)_(q_hat)
    e_prime.SetDims(d0, d_d_hat);
    e_tmp.SetLength(d0);

    for(k=0; k<d0; k++)
    {
        e_tmp[k] = 0;
    }

    for(j=0; j<d0; j++)
    {
        // Temporary coefficient vector to create e_j: it is a unit vector with its j-th coefficient being 1
        e_tmp[j] = 1;        

        // e_prime[j].SetLength(d_d_hat);
        CoeffsInvHat(e_prime[j], e_tmp, d_d_hat);

        // Reset the e_tmp coefficient vector
        e_tmp[j] = 0;
    }
   
    // Precompute σ(e_j), σ(p_j), σ(e′_j), σ(s_1), h_part2, h_part3
    sigma_e_.SetDims(256, n256);
    sigma_p_.SetDims(d0, m1);
    sigma_e_prime_.SetDims(d0, d_d_hat);
    sigma_s_1.SetLength(m1);
    h_part2.SetLength(d0);
    
    for(j=0; j<256; j++)        
    {
        sigma_map(sigma_e_[j], e_[j], d_hat);  
    }       

    for(j=0; j<d0; j++)        
    {
        CoeffsInvHat(p_j, P[j], m1);
        sigma_map(sigma_p_[j], p_j, d_hat);
        sigma_map(sigma_e_prime_[j], e_prime[j], d_hat);   
        h_part2[j] = poly_mult_hat(sigma_p_[j], s_1) - poly_mult_hat(sigma_e_prime_[j], u);
    }  

    sigma_map(sigma_s_1, s_1, d_hat);
    h_part3   = poly_mult_hat(sigma_s_1, s_1) + conv<zz_p>(-B_goth2);

    // Initialize the custom Hash function
    // ss << crs << P << u0 << B_goth2;
    ss << inputStr << ipk.c0 << ipk.c1 << idx_hid << u0 << B_goth2;
    // NOTE: using inputStr, ipk.c0, ipk.c1, idx_hid, instead of crs, P to speedup Hash_Init
    state0 = Hash_Init(ss.str());

    
    // 8. while (rst == 0 ∧ idx < N) do
    while((rst == 0) && (idx < N1))
    {
        b1 = 0; b2 = 0; b3 = 0;
        bbar1 = 0;   bbar2 = 0;
    
        // 9. Increment idx
        idx = idx + 1;
        // cout << "idx = " << idx << endl;
        
        // 10. Random generation of s_2 ∈ R^^(m2)_(q_hat) 
        s_2.SetLength(m2);

        for(i=0; i<m2; i++)
        {
            s_2[i].SetLength(d_hat);

            for(j=0; j<d_hat; j++)
            {
                // s_2[i][j] = conv<zz_p>( RandomBnd(3) - 1 );
                SetCoeff( s_2[i], j, conv<zz_p>( RandomBnd(3) - 1 ) );
                // NOTE: uniform distribution on ternary polynomials chi, that sample coeffs from {-1,0,1} mod q1_hat
            }
        }


        // 11. t_A = A_1*s_1 + A_2*s_2,  t_A ∈ R^^(n)_(q_hat)
        Pi.t_A.SetLength(n);
        
        for(i=0; i<n; i++)
        {
            Pi.t_A[i] = poly_mult_hat(crs[0][i], s_1) + poly_mult_hat(crs[1][i], s_2);
        }
        
        // 12. Random generation of the y_1 ∈ R^^m1_(q_hat),  y_2 ∈ R^^m2_(q_hat),  y_3 ∈ R^^(256/d_hat)_(q_hat)
        y_1.SetLength(m1);
        y_2.SetLength(m2);
        y_3.SetLength(n256);

        for(i=0; i<m1; i++)
        {
            polySampler_hat(y_1[i], s1_goth_d);
        }

        for(i=0; i<m2; i++)
        {
            polySampler_hat(y_2[i], s2_goth_d);
        }

        for(i=0; i<n256; i++)
        {
            polySampler_hat(y_3[i], s3_goth_d);
        }

        // 13. Random generation of g ∈ R^^(tau)_(q_hat)
        g.SetLength(tau0);

        for(i=0; i<tau0; i++)
        {
            g[i] = random_zz_pX(d_hat);            
            // g[i][0] = 0;
            SetCoeff(g[i], 0, 0);     
            // NOTE: the constant term of g (x^0) must be zero 
        }
    
        // 14. w = A1*y1 + A2*y2,  w ∈ R^^(n)_(q_hat)
        Pi.w.SetLength(n);
        // NOTE: it is different from the input w (= w0, from Prove_Init)
        
        for(i=0; i<n; i++)
        {
            Pi.w[i] = poly_mult_hat(crs[0][i], y_1) + poly_mult_hat(crs[1][i], y_2);
        }

        // 15. t_y = B_y*s2 + y3,  t_y ∈ R^^(256/d_hat)_(q_hat)
        Pi.t_y.SetLength(n256);
        
        for(i=0; i<n256; i++)
        {
            Pi.t_y[i] = poly_mult_hat(crs[2][i], s_2) + y_3[i];
        }

        // 16. t_g = B_g*s2 + g,  t_g ∈ R^^(tau)_(q_hat)
        Pi.t_g.SetLength(tau0);

        for(i=0; i<tau0; i++)
        {
            Pi.t_g[i] = poly_mult_hat(crs[3][i], s_2) + g[i];
        }
                                            
        // 17. (com_1, st_1) = LHC_Com(1, crs_LHC1, s_1, y_1)    
        LHC_Com(Pi.com_1, st_1, 1, crs[5], crs[7], s_1, y_1);

        // 18. (com_2, st_2) = LHC_Com(2, crs_LHC2, s_2, y_2)
        LHC_Com(Pi.com_2, st_2, 2, crs[6], crs[8], s_2, y_2);

        // 19. a1 ← (t_A, t_y, t_g, w, com_1, com_2) 
        ss.str("");
        ss << Pi.t_A << Pi.t_y << Pi.t_g << Pi.w << Pi.com_1.t_1 << Pi.com_1.t_2 << Pi.com_1.w_1 << Pi.com_1.w_2 << Pi.com_2.t_1 << Pi.com_2.t_2 << Pi.com_2.w_1 << Pi.com_2.w_2;
        // NOTE: copy the initial status structure, already initialized with (crs, x) before row 8        
        state = Hash_Copy(state0);
        Hash_Update(state, ss.str());

        // 20. (R_goth_0, R_goth_1) = H(1, crs, x, a_1)
        // 21. R_goth = R_goth_0 - R_goth_1
        HCom1(R_goth, state, "1");
        // NOTE: R_goth ∈ {-1, 0, 1}^(256 x m_1*d_hat) ⊂ Z^(256 x m_1*d_hat)_(q_hat)
        //       equivalent to (R_goth_0 - R_goth_1) in BLNS

        // 22. coeffs_y3 ← Coeffs(y_3),   coeffs_y3 ∈ Z^(256)_(q_hat)   
        CoeffsHat(coeffs_y3, y_3, n256);
        
        // 23.  z_3 = y_3 + R_goth*s,   z_3 ∈ Z^(256)_(q_hat)
        coeffs_R_goth_mult_s1.SetLength(256);
        Pi.z_3.SetLength(256);
        // NOTE: This equation is performed in Z not in polynomials, needing Coeffs() transformation of y_3, R_goth, and s_1

        for(i=0; i<256; i++)
        {
            coeffs_R_goth_mult_s1[i] = R_goth[i] * coeffs_s0;
            // NOTE: this term corresponds to InnerProduct(result, R_goth[i], coeffs_s0);

            Pi.z_3[i] = coeffs_y3[i] + coeffs_R_goth_mult_s1[i];            
        }
        
        // 24. b3 ← Rej (z_3, R_goth * s, s3_goth, M_3),    b3 ∈ {0, 1}
        b3 = Rej_v_zzp(Pi.z_3, coeffs_R_goth_mult_s1, q1_hat, s3_goth, M_3);
        
        // NOTE: if b3 == 0, continue the while loop (skip next rows until 49, then go to row 8)
        if (b3 == 0)
        {            
            rst = 0;
            continue;
        }
        
        // 25. a2 ← z_3,   a2 ∈ Z^256_(q_hat)
        ss.str("");
        ss << Pi.z_3;
        Hash_Update(state, ss.str());
        
        // 26. gamma ← H(2, crs, x, a1, a2),   gamma ∈ Z^(tau0 x 256+d0+1)_(q_hat)
        HCom2(gamma, state, "2");


        // Initialize h ∈ R^^(tau)_(q_hat)
        Pi.h.SetLength(tau0);

        // Precompute σ(r_j), h_part1
        sigma_r_.SetDims(256, m1);
        h_part1.SetLength(256); 

        for(j=0; j<256; j++)        
        {
            CoeffsInvHat(r_j, R_goth[j], m1); 
            sigma_map(sigma_r_[j], r_j, d_hat);  
            h_part1[j] = poly_mult_hat(sigma_r_[j], s_1) + poly_mult_hat(sigma_e_[j], y_3) - Pi.z_3[j];
        }


        // 27. for i ∈ [τ] do
        for(i=0; i<tau0; i++) 
        {
            // 28. Compute h_i,   h_i ∈ R^_(q_hat)
            acc = g[i];

            for(j=0; j<256; j++)        
            {
                acc += gamma[i][j] * h_part1[j];
            }
                        
            for(j=0; j<d0; j++)        
            {
                acc += gamma[i][256+j] * h_part2[j];
            }
            
            acc += gamma[i][256+d0] * h_part3;
                    
            // 29. h ← (h_1, . . . , h_τ),   h ∈ R^^(tau)_(q_hat)
            Pi.h[i] = acc;
        }

        // 30. a_3 ← h,   a_3 ∈ R^^(tau)_(q_hat)
        ss.str("");
        ss << Pi.h;
        Hash_Update(state, ss.str());

        // 31. μ ← H(3, crs, x, a1, a2, a3),   μ ∈ R^^(tau)_(q_hat)        
        HCom3(mu, state, "3");

        // 32. B   ← [B_y; B_g],   B ∈ R^^((256/d_hat + tau) x m2)_(q_hat)
        B.SetDims((n256 + tau0), m2);
        
        for(i=0; i<n256; i++) 
        {
            B[i]   = crs[2][i];
        }

        for(i=n256; i<(n256 + tau0); i++) 
        {
            B[i]   = crs[3][i-n256];
        }


        // 33. y ← (y1; σ(y1); −B*y2; σ(-B*y2)),     y ∈ R^^(2*m1 + 2*(256/d_hat + tau))_(q_hat)
        y.SetLength( m1_n256_tau );

        for(i=0; i<m1; i++) 
        {
            y[i] = y_1[i];        
        }

        sigma_map(sigma_y_1, y_1, d_hat);
        k = 0;

        for(i=m1; i<(2*m1); i++) 
        {
            y[i] = sigma_y_1[k];  
            k++;      
        }
            
        // Compute -B*y2 in a temporary vector
        tmp_vec.SetLength(n256 + tau0);

        for(i=0; i<(n256 + tau0); i++)
        {
            tmp_vec[i] = -poly_mult_hat(B[i], y_2);
        }
        
        k = 0;

        for(i=(2*m1); i<(2*m1 + n256 + tau0); i++) 
        {
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


        // 34. Definition of D_2_(2,1) ∈ R^^(m1 x m1)_(q_hat) 
        D2_2_1.SetDims(m1, m1);        
        // sum = 0;
        clear(sum);
                
        for(i=0; i<tau0; i++)
        {
            sum += ( mu[i] * gamma[i][256+d0] );           
        }
        
        for(i=0; i<m1; i++)
        {
            D2_2_1[i][i] = sum;
        }
        // NOTE: sum in the diagonal of D2_2_1, zeros in the rest


        // 35. Construction of d_1 ∈ R^^(2*m1+2(256/d_hat+τ))_(q_hat)
        d_1.SetLength(m1_n256_tau);

        for(i=0; i<m1_n256_tau; i++)
        {
            clear(d_1[i]);
        }

        // 1st entry of d_1 (m1 polynomials)
        acc_vec.SetLength(m1);
        
        for(i=0; i<tau0; i++)
        {
            // Reset acc_vec
            for(j=0; j<m1; j++)
            {        
                // acc_vec[j] = 0;
                clear(acc_vec[j]);
            }               
                    
            for(j=0; j<256; j++)        
            {
                for(k=0; k<m1; k++)        
                {
                    acc_vec[k] += gamma[i][j] * sigma_r_[j][k];  
                }
            }

            for(j=0; j<d0; j++)        
            {
                for(k=0; k<m1; k++)        
                {                
                    acc_vec[k] += gamma[i][256+j] * sigma_p_[j][k];  
                }     
            }  
                        
            for(k=0; k<m1; k++)        
            {
                // Fill d_1 (first m1 polynomials) by accumulating mu[i]*(...sums...) 
                d_1[k] += ModPhi_hat_q( mu[i] * acc_vec[k]); 
            }               
        }

        // NOTE: skip 2nd entry of d_1 (m1 zeros)

        // 3rd entry of d_1 (256/d_hat polynomials)    
        acc_vec.SetLength(n256);
        
        for(i=0; i<tau0; i++)
        {
            // Reset acc_vec
            for(j=0; j<n256; j++)
            {        
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
                // Fill d_1 (256/d_hat polynomials) by accumulating mu[i]*(...sums...) 
                d_1[(2*m1)+k] += ModPhi_hat_q( mu[i] * acc_vec[k]);             
            }
        } 
        
        // 4th entry of d_1 (tau0 polynomials)    
        k = 0;

        for(i=(2*m1 + n256); i<(2*m1 + n256 + tau0); i++) 
        {
            d_1[i] = mu[k];
            k++;    
        }
        // NOTE: skip 5th entry of d_1 (256/d_hat + tau0 zeros)
        
                        
        // 36. Definition of f1 ∈ R^_(q_hat)
        // 1st addend of f1, (σ(s_1)^T * D2_2_1 * y_1)
        D2_y.SetLength(m1);

        // Compute  (D2_2_1 * y_1),  m1 polynomials
        for(i=0; i<m1; i++)    
        {
            D2_y[i] = poly_mult_hat(D2_2_1[i], y_1);
        }

        // Accumulate  σ(s_1)^T * (D2_2_1 * y_1)
        f1 = poly_mult_hat(sigma_s_1, D2_y);

        // 2nd addend of f1,  (σ(y_1)^T * D2_2_1 * s_1)
        acc_vec.SetLength(m1);

        // Compute  (D2_2_1 * s_1),  m1 polynomials
        for(i=0; i<m1; i++)    
        {
            acc_vec[i] = poly_mult_hat(D2_2_1[i], s_1);
        }

        // Accumulate  σ(y_1)^T * (D2_2_1 * s_1)    
        f1 += poly_mult_hat(sigma_y_1, acc_vec);

        // 3rd addend of f1,   (d_1^T * y)
        f1 += poly_mult_hat(d_1, y);


        // 37. Definition of f0 ∈ R^_(q_hat)
        Pi.f0 = poly_mult_hat(sigma_y_1, D2_y) + poly_mult_hat(crs[4][0], y_2);
        // NOTE: D2_y = (D2_2_1 * y_1) was already computed in row 37 (1st addend of f1) 
    
        
        // 38. Definition of t ∈ R^_(q_hat)
        Pi.t = poly_mult_hat(crs[4][0], s_2) + f1;

        // 39. a_4 ← (t, f0),   a_4 ∈ R^^2_(q_hat)
        ss.str("");
        ss << Pi.t << Pi.f0;
        Hash_Update(state, ss.str());

        // 40. c ← H(4, crs, x, a1, a2, a3, a4),   c ∈ C ⊂ R^_(q_hat)
        HCom4(c, state, "4");


        // 41. for i ∈ {1, 2} do
        // NOTE: for simplicity, next operations are duplicated with suffixes _1 and _2

        // 42. z_i ← y_i + c*s_i,   z_i ∈ R^^(m_i)_(q_hat)
        Pi.z_1.SetLength(m1);
        Pi.z_2.SetLength(m2);

        c_s1.SetLength(m1);
        c_s2.SetLength(m2);
        
        for(i=0; i<m1; i++)
        {
            c_s1[i] = ModPhi_hat_q( c * s_1[i] );
            Pi.z_1[i] = y_1[i] + c_s1[i];
        }
            
        for(i=0; i<m2; i++)
        {
            c_s2[i] = ModPhi_hat_q( c * s_2[i] );
            Pi.z_2[i] = y_2[i] + c_s2[i];
        }
            

        // 43. b_i ← Rej(z_i, c*s_i, s_i_goth, M_i),   b_i ∈ {0, 1} 
        b1 = Rej_v_zzpX(Pi.z_1, c_s1, q1_hat, s1_goth, M_1);

        // NOTE: if b1 == 0, continue the while loop (skip next rows until 49, then go to row 8)
        if (b1 == 0)
        {            
            rst = 0;
            continue;
        }

        b2 = Rej_v_zzpX(Pi.z_2, c_s2, q1_hat, s2_goth, M_2);

        // NOTE: if b2 == 0, continue the while loop (skip next rows until 49, then go to row 8)
        if (b2 == 0)
        {            
            rst = 0;
            continue;
        }


        // 44. op_i ← LHC.Open(i, c, st_i),    op_i ∈ {⊥} ∪ R^^(n_i)_(q_hat) × R^^(m_i)_(q_hat) × R^^(m_i)_(q_hat)
        LHC_Open(Pi.op_1, 1, c, st_1);        

        // 45. if op_i = ⊥ then b_bar_i = 0
        if (Pi.op_1.valid == 0)
        {
            bbar1 = 0;
            // NOTE: if bbar1 == 0, continue the while loop (skip next rows until 49, then go to row 8)
            rst = 0;
            continue;
        }
        // 46. else b_bar_i = 1
        else
        {
            bbar1 = 1;
        }

        LHC_Open(Pi.op_2, 2, c, st_2);

        if (Pi.op_2.valid == 0)
        {
            bbar2 = 0;
        }
        else
        {
            bbar2 = 1;
        }

        // 47. π ← (t_A, t_y, t_g, w, com_1, com_2, z_3, h, t, f0, z_1, z_2, op_1, op_2)    
        // Pi.t_A   = t_A;
        // Pi.t_y   = t_y;
        // Pi.t_g   = t_g;
        // Pi.w     = w;
        // Pi.com_1 = com_1;
        // Pi.com_2 = com_2;
        // Pi.z_3   = z_3;
        // Pi.h     = h;
        // Pi.t     = t;
        // Pi.f0    = f0;
        // Pi.z_1   = z_1;
        // Pi.z_2   = z_2;
        // Pi.op_1  = op_1;
        // Pi.op_2  = op_2;        

        // 48. rst ← b1*b2*b3*b_bar_1*b_bar_2
        rst = b1*b2*b3*bbar1*bbar2;        
    
    } // End of while loop (row 8)

    delete  state0;
    delete  state;

    // 49. if rst = 1 then return π
    if (rst == 1)
    {
        // NOTE: additional flag, to identify a valid proof
        Pi.valid = 1;
    }
    // 50. else return ⊥
    else // (rst == 0)      
    {
        // NOTE: invalid proof
        Pi.valid = 0;        
    }
    
    // return Pi;
}


//==============================================================================
// Verify_Com   -   Verify the Commitment (Verify^HCom_Com). 
//                  This function takes as input crs, (P, u, B_goth), and the proof π.
//                  It returns as output “reject” or “accept”.
// 
// Inputs:
// - inputStr:      string containing the initial seed for crs
// - crs:           structure crs_Com, generated by Hcrs from inputStr
// - ipk:           Issuer public key
// - (P,u,B_goth2): this triplet corresponds to x, with:
//                  P ∈ Z^[d x (m1*d_hat)]_(q_hat)
//                  u ∈ Z^(d)_(q_hat)
//                  B_goth^2 ∈ Z≥0 (it is a scalar)
// - Pi:            proof (π) structure 
//  
// Output:
// - 0 or 1:        reject or accept 
//==============================================================================
long Verify_Com(const string& inputStr, const CRS_t& crs, const IPK_t& ipk, const mat_zz_p& P, const vec_zz_p& u0, const ZZ& B_goth2, const PROOF_C_t& Pi)
{
    // NOTE: assuming that current modulus is q1_hat (not q0)

    unsigned long       i, j, k;
    long                b1, b2;
    Mat<zz_pX>          B, D2_2_1;
    vec_zz_pX           u, t_B, mu, z, tmp_vec, tmp_vec2, r_j, p_j;
    vec_zz_pX           d_1, acc_vec, sigma_z_1;
    zz_pX               c, acc, sum, d_0, sum_sigma_e_u;
    Mat<zz_pX>          e_, e_prime;
    Mat<zz_pX>          sigma_r_, sigma_p_, sigma_e_, sigma_e_prime_;
    stringstream        ss;
    mat_zz_p            R_goth, gamma;
    vec_zz_p            e_tmp;
    zz_p                sum_z3, B_goth_p;
    ZZ                  norm2_z1, norm2_z2, norm2_z3;
    HASH_STATE_t*       state;

    // Initialise constants and variables
    const unsigned long n           = n_Com;
    const unsigned long m1          = m1_Com;
    const unsigned long m2          = m2_Com;
    const unsigned long d_d_hat     = (d0/d_hat);
    const unsigned long n256        = (256/d_hat);
    const unsigned long m1_n256_tau = 2*m1 + 2*(n256 + tau0);

    // Initialise the "goth" constants
    // NOTE: all values are squared, to avoid sqrt and floating points
    const ZZ  s1_goth2  = sqr(ZZ(alpha_1)  * ZZ(nu0))   * B_goth2;
    const ZZ  s2_goth2  = sqr(ZZ(alpha_2)  * ZZ(nu0))   * ( m2 * ZZ(d_hat) );
    const ZZ  s3_goth2  = sqr(ZZ(alpha_3)) * ZZ(w_max2) * B_goth2;
    const ZZ  B_goth2_1 = s1_goth2 * ( 2*m1*ZZ(d_hat) );
    const ZZ  B_goth2_2 = s2_goth2 * ( 2*m2*ZZ(d_hat) );
    const ZZ  B_goth2_3 = RoundToZZ( sqr(RR(1.7)) * conv<RR>(s3_goth2) * RR(256) );
    // NOTE:  B_goth2_3 is a floating point (due to 1.7 factor), rounded to nearest integer


    // 1. Retrieve from crs_Com
    // A_1     = crs[0];    // ∈ (R^_q^)^(n x m1)
    // A_2     = crs[1];    // ∈ (R^_q^)^(n x m2)
    // B_y     = crs[2];    // ∈ (R^_q^)^(256/d^ x m2)
    // B_g     = crs[3];    // ∈ (R^_q^)^(tau0 x m2)
    // b       = crs[4][0]; // ∈ (R^_q^)^(m2)         NOTE: b in crs is a (1 x m_2) matrix
    // Abar_1  = crs[5];    // ∈ (R^_q^)^(m1 x n1)
    // Abar_2  = crs[6];    // ∈ (R^_q^)^(m2 x n2)
    // Bbar_1  = crs[7];    // ∈ (R^_q^)^(m1 x n1)
    // Bbar_2  = crs[8];    // ∈ (R^_q^)^(m2 x n2)

    // 2. Retrieve (P, u, B_goth) ← x
    // NOTE: (P, u0, B_goth) already provided as inputs, so we just need to 
    //       convert u0 ∈ Z^(d)_(q_hat)  to  u ∈ R^^(d/d_hat)_(q_hat), needed at row 19
    CoeffsInvHat(u, u0, d_d_hat);

    // 3. P ← [P1,  0_(d × d_hat)],   P ∈ Z^[d x (|idx_hid|·h + ℓr·d + d_hat]_(q_hat)    
    // NOTE: zero padding of P already done in I_VerCred

    // 4. (t_A, t_y, t_g, w, com1, com2, z_3, h, t, f0, z_1, z_2, op1, op2) ← π       
    // Check if Pi contains a valid proof   
    if (Pi.valid != 1)
    {
        cout << "ERROR! Pi does not contain a valid proof" << endl;
        return 0;
    }
    // NOTE: to save memory, proof values will be directly accessed as Pi.{name}

    // Initialize the custom Hash function
    // ss << crs << P << u0 << B_goth2;
    ss << inputStr << ipk.c0 << ipk.c1 << idx_hid << u0 << B_goth2;
    // NOTE: using inputStr, ipk.c0, ipk.c1, idx_hid, instead of crs, P to speedup Hash_Init
    state = Hash_Init(ss.str());

    // 5. a1 ← (t_A, t_y, t_g, w, com_1, com_2) 
    // a_1 << Pi.t_A << Pi.t_y << Pi.t_g << Pi.w << Pi.com_1 << Pi.com_2;
    ss.str("");
    ss << Pi.t_A << Pi.t_y << Pi.t_g << Pi.w << Pi.com_1.t_1 << Pi.com_1.t_2 << Pi.com_1.w_1 << Pi.com_1.w_2 << Pi.com_2.t_1 << Pi.com_2.t_2 << Pi.com_2.w_1 << Pi.com_2.w_2;
    Hash_Update(state, ss.str());

    // 6. a2 ← z_3,   a2 ∈ Z^256_(q_hat)
    // a_2 << Pi.z_3;
    
    // 7. a_3 ← h,   a_3 ∈ R^^(tau)_(q_hat)    
    // a_3 << Pi.h;
    
    // 8. a_4 ← (t, f0),   a_4 ∈ R^^2_(q_hat)
    // a_4 << Pi.t << Pi.f0;

    // 9. (R_goth_0, R_goth_1) = H(1, crs, x, a_1)    
    // 10. R_goth = R_goth_0 - R_goth_1
    HCom1(R_goth, state, "1");
    // NOTE: R_goth ∈ {-1, 0, 1}^(256 x m_1*d_hat) ⊂ Z^(256 x m_1*d_hat)_(q_hat)
    //       equivalent to (R_goth_0 - R_goth_1) in BLNS

    // 11. gamma ← H(2, crs, x, a1, a2),   gamma ∈ Z^(tau0 x 256+d0+1)_(q_hat)     
    ss.str("");
    ss << Pi.z_3;
    Hash_Update(state, ss.str());
    HCom2(gamma, state, "2");

    // 12. μ ← H(3, crs, x, a1, a2, a3),   μ ∈ R^^(tau)_(q_hat)        
    ss.str("");
    ss << Pi.h;
    Hash_Update(state, ss.str());
    HCom3(mu, state, "3");

    // 13. c ← H(4, crs, x, a1, a2, a3, a4),   c ∈ C ⊂ R^_(q_hat)
    ss.str("");
    ss << Pi.t << Pi.f0;
    Hash_Update(state, ss.str());
    HCom4(c, state, "4");

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
        z[i] = Pi.z_1[i]; // z_1
    }

    sigma_map(sigma_z_1, Pi.z_1, d_hat);
    k = 0;

    for(i=m1; i<(2*m1); i++) 
    {
        z[i] = sigma_z_1[k]; // σ(z_1)
        k++;      
    }

    // Compute (c*t_B − B*z_2) in a temporary vector
    tmp_vec.SetLength(n256 + tau0);
    
    for(i=0; i<(n256 + tau0); i++)
    {
        tmp_vec[i] = ModPhi_hat_q( c * t_B[i] ) - poly_mult_hat(B[i], Pi.z_2);
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


    // 17. Definition of D_2_(2,1) ∈ R^^(m1 x m1)_(q_hat) 
    D2_2_1.SetDims(m1, m1);        
    // sum = 0;
    clear(sum);
            
    for(i=0; i<tau0; i++)
    {
        sum += ( mu[i] * gamma[i][256+d0] );           
    }

    for(i=0; i<m1; i++)
    {
        D2_2_1[i][i] = sum;
    }
    // NOTE: sum in the diagonal of D2_2_1, zeros in the rest

    
    // Initialize e ∈ R^^(256 x 256/d_hat)_(q_hat)
    e_.SetDims(256, n256);
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

    // Initialize e_prime ∈ R^^(d0 x d0/d_hat)_(q_hat)
    e_prime.SetDims(d0, d_d_hat);
    e_tmp.SetLength(d0);

    for(k=0; k<d0; k++)
    {
        e_tmp[k] = 0;
    }

    for(j=0; j<d0; j++)
    {
        // Temporary coefficient vector to create e_j: it is a unit vector with its j-th coefficient being 1
        e_tmp[j] = 1;        

        // e_prime[j].SetLength(d_d_hat);
        CoeffsInvHat(e_prime[j], e_tmp, d_d_hat);

        // Reset the e_tmp coefficient vector
        e_tmp[j] = 0;
    }
    
    // Precompute σ(r_j), σ(p_j), σ(e_j), σ(e′_j) 
    sigma_r_.SetDims(256, m1);          
    sigma_e_.SetDims(256, n256);
    sigma_p_.SetDims(d0, m1);
    sigma_e_prime_.SetDims(d0, d_d_hat);
    
    for(j=0; j<256; j++)        
    {
        CoeffsInvHat(r_j,  R_goth[j], m1);
        sigma_map(sigma_r_[j], r_j, d_hat);
        sigma_map(sigma_e_[j], e_[j], d_hat);  
    }

    for(j=0; j<d0; j++)        
    {
        CoeffsInvHat(p_j, P[j], m1);
        sigma_map(sigma_p_[j], p_j, d_hat);
        sigma_map(sigma_e_prime_[j], e_prime[j], d_hat);            
    }

    
    // 18. Construction of d_1 ∈ R^^(2*m1+2(256/d_hat+τ))_(q_hat)
    d_1.SetLength(m1_n256_tau);

    for(i=0; i<m1_n256_tau; i++)
    {
        clear(d_1[i]);
    }    

    // 1st entry of d_1 (m1 polynomials)
    acc_vec.SetLength(m1);
    
    for(i=0; i<tau0; i++)
    {
        // Reset acc_vec
        for(j=0; j<m1; j++)
        {        
            // acc_vec[j] = 0;
            clear(acc_vec[j]);
        }               
                
        for(j=0; j<256; j++)        
        {
            for(k=0; k<m1; k++)        
            {
                acc_vec[k] += gamma[i][j] * sigma_r_[j][k];  
            }
        }

        for(j=0; j<d0; j++)        
        {
            for(k=0; k<m1; k++)        
            {                
                acc_vec[k] += gamma[i][256+j] * sigma_p_[j][k];  
            }     
        }  
                    
        for(k=0; k<m1; k++)        
        {
            // Fill d_1 (first m1 polynomials) by accumulating mu[i]*(...sums...) 
            d_1[k] += ModPhi_hat_q( mu[i] * acc_vec[k]); 
        }               
    }

    // NOTE: skip 2nd entry of d_1 (m1 zeros)

    // 3rd entry of d_1 (256/d_hat polynomials)    
    acc_vec.SetLength(n256);
  
    for(i=0; i<tau0; i++)
    {
        // Reset acc_vec
        for(j=0; j<n256; j++)
        {        
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
            // Fill d_1 (256/d_hat polynomials) by accumulating mu[i]*(...sums...) 
            d_1[(2*m1)+k] += ModPhi_hat_q( mu[i] * acc_vec[k]);             
        }
    }
    
    // 4th entry of d_1 (tau0 polynomials)    
    k = 0;

    for(i=(2*m1 + n256); i<(2*m1 + n256 + tau0); i++) 
    {
        d_1[i] = mu[k];
        k++;    
    }
    // NOTE: skip 5th entry of d_1 (256/d_hat + tau0 zeros)


    // 19. Definition of d_0 ∈ R^_(q_hat)
    clear(d_0);
    // NOTE: d_0, not d0 parameter

    B_goth_p = conv<zz_p>(B_goth2);

    for(i=0; i<tau0; i++)
    {            
        sum_z3 = 0;
        
        for(j=0; j<256; j++)
        {
            sum_z3 += gamma[i][j] * Pi.z_3[j];
        }
            
        clear(sum_sigma_e_u);
        
        for(j=0; j<d0; j++)
        {
            sum_sigma_e_u += gamma[i][256+j] * poly_mult_hat(sigma_e_prime_[j], u);        
        }
    
        d_0 = d_0 - ModPhi_hat_q( mu[i] * ( sum_z3 + sum_sigma_e_u + gamma[i][256+d0] * B_goth_p + Pi.h[i] )); 
    }


    // 20.  if one of the 4 conditions below does not hold, then return 0
       
    // Compute ||z_i||^2, squared Euclidean norm of each z_i
    norm2_z1 = Norm2Xm(Pi.z_1, d_hat, q1_hat);
    norm2_z2 = Norm2Xm(Pi.z_2, d_hat, q1_hat);
    norm2_z3 = Norm2m( Pi.z_3, q1_hat );
    // NOTE: norms computed using values in {-(q-1)/2, ..., (q-1)/2}


    // 20.1 First condition: ||z_1|| ≤ B_goth_1, ||z_2|| ≤ B_goth_2, ||z_3|| ≤ B_goth_3
    // NOTE: equations in ZZ, with squared norms and thresholds
    if ( norm2_z1 > B_goth2_1)
    { 
        cout << "First condition failed - Invalid z_1 norm!" << endl; 
        return 0;
    }

    if ( norm2_z2 > B_goth2_2)
    { 
        cout << "First condition failed - Invalid z_2 norm!" << endl; 
        return 0;
    }

    if ( norm2_z3 > B_goth2_3)
    { 
        cout << "First condition failed - Invalid z_3 norm!" << endl; 
        return 0;
    }


    // 20.2 Second condition: h˜_i == 0 for i ∈ [τ] 
    // NOTE: equations in R^^_(q_hat)
    for(i=0; i<tau0; i++)
    {
        if ( coeff(Pi.h[i], 0) != 0 )
        {
            cout << "Second condition failed! \n  h = " << Pi.h << endl;
            return 0;
        }
    }


    // 20.3 Third condition: A_1*z_1 + A_2*z_2 == w + c*t_A 
    // NOTE: equations in R^^(n)_(q_hat)    
    tmp_vec.SetLength(n);
    tmp_vec2.SetLength(n);
    
    for(i=0; i<n; i++)
    {
        // A_1*z_1 + A_2*z_2
        tmp_vec[i] = poly_mult_hat(crs[0][i], Pi.z_1) + poly_mult_hat(crs[1][i], Pi.z_2);

        // w + c*t_A 
        tmp_vec2[i] = Pi.w[i] + ModPhi_hat_q( c * Pi.t_A[i] );
    }

    if (tmp_vec != tmp_vec2)
    {
        cout << "Third condition failed!" << endl; 
        return 0;
    }


    // 20.4 Fourth condition: σ(z_1)^T*D2_2_1*z_1 + c*d_1^T*z + c^2*d_0 − (c*t − b^T*z_2) == f0 
    // NOTE: equations in R^^_(q_hat)
    
    // 1st addend σ(z_1)^T * (D2_2_1 * z_1)
    // Compute  (D2_2_1 * z_1),  m1 polynomials
    acc_vec.SetLength(m1);
    
    for(i=0; i<m1; i++)    
    {
        acc_vec[i] = poly_mult_hat(D2_2_1[i], Pi.z_1);
    }

    // Accumulate  σ(z_1)^T * (D2_2_1 * z_1)    
    acc = poly_mult_hat(sigma_z_1, acc_vec);
    
    // 2nd addend (c * d_1^T * z)
    acc_vec.SetLength(m1_n256_tau);
    // Compute  (c * d_1^T),  (2*m1 + 2*(256/d_hat + tau0)) polynomials
    for(i=0; i<m1_n256_tau; i++)    
    {
        acc_vec[i] = ModPhi_hat_q( c * d_1[i] );
    }

    // Accumulate  (c * d_1^T) * z 
    acc += poly_mult_hat(acc_vec, z);   
    
    // 3rd addend (c^2 * d_0)
    acc += ModPhi_hat_q( ModPhi_hat_q( sqr(c) ) * d_0 );
      
    // 4rd addend −(c*t − b^T * z_2)
    acc -= ( ModPhi_hat_q( c * Pi.t ) - poly_mult_hat(crs[4][0], Pi.z_2) );
    
    if (acc != Pi.f0)
    {
        cout << "Fourth condition failed!" << endl; 
        return 0;
    }
    

    // 20.5 Fifth condition: for i ∈ {1, 2}, LHC.Verify_i((A_i, B_i), (com_i, c), (z_i, op_i)) == 1
    b1 = LHC_Verify(1, crs[5], crs[7], Pi.com_1, c, Pi.z_1, Pi.op_1);
    b2 = LHC_Verify(2, crs[6], crs[8], Pi.com_2, c, Pi.z_2, Pi.op_2);
    
    if ((b1==0) || (b2==0))
    {
        cout << "Fifth condition failed!" << endl; 
        return 0; 
    }

    // cout << "# Verify_Com: OK!" << endl;


    // 21. else, return 1
    return 1;
}