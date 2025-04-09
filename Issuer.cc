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

#include "Issuer.h"


//==============================================================================
// I_KeyGen - Issuer.KeyGen function. It generates Issuer Public Key and 
//            Issuer Secret Key from parameters d and q.
//
// Inputs:
// - None
//
// Outputs:
// - ipk:   Issuer Public key, containing a1, a2, c0, c1 vectors of polynomials
// - isk:   Issuer Secret Key, i.e. matrix of integers B (size 2d * 2d)
//
// NOTE: Issuer.KeyGen in BLNS pseudocode, corresponding to
//       Fig. 18: AnonCreds.Init, pag. 52 in [BLNS23]
//==============================================================================
void I_KeyGen(IPK_t& ipk, ZZX& f, ZZX& g, ZZX& F, ZZX& G)
{
    unsigned long i;

#ifdef ENABLE_FALCON
    // Keygen algorithm from the Falcon reference implementation 
    Falcon_keygen(ipk.a1, f, g, F, G);
#else
    // NTRU.TrapGen(q, d) algorithm in [BLNS23]
    NTRU_TrapGen(ipk.a1, f, g, F, G);    
#endif
    // NOTE: a1 is a Polynomial with d coefficients modulo q (i.e. h in [DLP14])
    
    ipk.a2.SetLength(m0);   
    // NOTE: a2 is a vector of m Polynomials with d coefficients modulo q

    for(i=0; i<m0; i++)
    {
        ipk.a2[i] = random_zz_pX(d0);
    }

    ipk.c0.SetLength(lm0);  
    // NOTE: c0 is a vector of l_m Polynomials with d coefficients modulo q

    for(i=0; i<lm0; i++)
    {
        ipk.c0[i] = random_zz_pX(d0);
    }

    ipk.c1.SetLength(lr0);  
    // NOTE: c1 is a vector of l_r Polynomials with d coefficients modulo q

    for(i=0; i<lr0; i++)
    {
        ipk.c1[i] = random_zz_pX(d0);
    }
          
    // Output Issuer Public Key and Issuer Secret Key
    // ipk ← (a1, a2, c0, c1)
    // isk ← B
}


//==============================================================================
// I_VerCred    -   Issuer.VerCred function
// 
// Inputs:
// - crs_seed:      initial public seed for crs structure
// - crs:           structure with the pair (crs_ISIS, crs_Com), generated by Hcrs
// - B_f:           public random matrix B_f ∈ Z^(nd×t)_q
// - ipk:           Issuer Public key, containing a1, a2, c0, c1 vectors of polynomials
// - isk:           Issuer Secret Key, i.e. matrix of integers B (size 2d * 2d)
// - attrs_prime:   disclosed attributes (attrs′)
// - idx_pub:       |idx|, number of disclosed attributes
// - u, Pi_ptr:     commitment u and proof π, corresponding to the structure ρ_1 
// 
// Outputs:
// - s_0:           short vector (output of GSampler),       s_0 ∈ Z^(2d)
// - w:             polynomial vector (output of GSampler),  w ∈ R^m
// - x:             random integer, uniformly sampled from the set [N]
//                  NOTE: (s_0, w, x) correspond to the structure ρ_2
//==============================================================================
void I_VerCred(vec_ZZ& s_0, vec_ZZX& w, ZZ& x, const unsigned char* crs_seed, const CRS2_t& crs, const mat_zz_p& B_f, const IPK_t& ipk, const ZZX& f, const ZZX& g, const ZZX& F, const ZZX& G, const Vec<string>& attrs_prime, const zz_pX& u, uint8_t** Pi_ptr)
{    
    // NOTE: assuming that current modulus is q0 (not q_hat)
    unsigned long   i, j, k, result;
    zz_pX           a1, fx_u;
    vec_zz_pX       a2, c0, c1; //mex_prime;
    vec_ZZ          m_i, coeffs_m;
    mat_zz_p        P0, P1, P; 
    vec_zz_p        u_vect, prod;   
    ZZ              B_goth2;
    long            mul;
    mat_L           A, B;

    
    const unsigned long idxhlrd = (idx_hid * h0) + (lr0 * d0); //|idx_hid|·h + ℓr·d

    
    // 1. (a'_1, ... , a'_k) ← attrs',  a'_i ∈ {0, 1}∗
    // NOTE: l0 = idx_hid + idx_pub = len(attrs),  d0 must divide l0*h0
    // NOTE: for every variable of l0 elements, the first are the idx_hid elements, the last are the idx_pub elements
    
    // 2. (a1, a2, c0, c1) ← ipk,   ipk ∈ R_q × R^m_q × R^ℓm_q × R^ℓr_q
    a1 = ipk.a1;
    a2 = ipk.a2;
    c0 = ipk.c0;
    c1 = ipk.c1;

    // 3. (u, π) ← ρ1
    // NOTE: u, Pi provided as separate inputs

    // 4. B ← isk,   B ∈ Z^(2d×2d)    
    // NOTE: B replaced with isk in next rows

    // 5. m′← Coeffs^−1(H_M(a′_1), . . . , H_M(a′_k )) ∈ R^ℓm_q    
    coeffs_m.SetLength(l0 * h0);
    k = 0;

    for(i=0; i<l0; i++)
    {                  
        // a_i =  attrs_prime[i];        
        HM(m_i, attrs_prime[i] );        

        for(j=0; j<h0; j++)     
        {
            coeffs_m[k] = m_i[j];
            k++;
        }
    }    

    // mex_prime = CoeffsInv(coeffs_m, lm0);
    // NOTE: coeffs_m is directly used instead of mex_prime


    // 6. P ← [rot(c0^T)_(idx_hid) | rot(c1^T)],    P ∈ Z_q^(d × (|idx_hid|·h + ℓr·d))     
    P.SetDims(d0, (idxhlrd + d_hat));
    // NOTE: zero padding of P (d_hat columns) anticipated here, from Verify_Com

    P0.SetDims(d0, lm0*d0);
    rot_vect(P0, c0);

    // NOTE: only first idx_hid*h0 columns of P0 (corresponding to undisclosed attributes) 
    //       are copied into P, while P1 is fully copied into P.
    k = 0;

    for(j=0; j<(idx_hid*h0); j++)
    {
        for(i=0; i<d0; i++)
        {   
            P[i][k] = P0[i][j];
        }
        k++;
    }    

    P1.SetDims(d0, lr0*d0);
    rot_vect(P1, c1);

    for(j=0; j<(lr0*d0); j++)
    {
        for(i=0; i<d0; i++)
        {   
            P[i][k] = P1[i][j];
        }
        k++;
    }

    P1.kill();

   
    // 7. u ← Coeffs(u) − rot(c0^T)_idx * Coeffs(m')_idx ∈ Z_q^d
    u_vect.SetLength(d0); 
    prod.SetLength(d0);

    // NOTE: only last idx_pub*h0 columns of P0 and coeffs_m (corresponding to disclosed attributes) 
    //       are considered in the product rot(c0^T)_idx * Coeffs(m')_idx 
    for(j=0; j<d0; j++)
    {
       for(i=0; i<(idx_pub*h0); i++)
        {     
            k = idx_hid*h0 + i;       
            prod[j] += P0[j][k] * conv<zz_p>( coeffs_m[k] );
        }
    }

    P0.kill();

    for(i=0; i<d0; i++)
    {
        // u_vect[i] = u[i] - prod[i];
        u_vect[i] = coeff(u, i) - prod[i];
    }


    // 8. if  Verify_Com(crs_Com, (q1_hat/q*P, q1_hat/q*u_vect, ψ*sqrt(h*|idx|+ℓr*d)), π) == 0:
    if (not(divide( ZZ(q1_hat), q0)))
    {
        cout << " ERROR: q1_hat must be divisible by q! " << endl;
    }

    mul = long(q1_hat) / long(q0);
    
    // B_goth = psi0 * sqrt(conv<RR>( idxhlrd ));    
    B_goth2 = sqr( ZZ(psi0)) * ZZ( idxhlrd );
    
    {
        zz_pPush push(q1_hat); 
        // NOTE: backup current modulus q0, temporarily set to q1_hat (i.e., zz_p::init(q1_hat))

        result = Verify_Com(crs_seed, crs[1], ipk, (mul * P), (mul * u_vect), B_goth2, Pi_ptr);
        // NOTE: P, u are converted from modulo q0 to q1_hat
    }

    P.kill();
        
    if (result == 0)
    {
        // 9. return ⊥       
        cout << "\n Invalid proof Pi! (return ⊥ )" << endl;
        s_0.SetLength(0);
        x = -1;
        return;
    }


    // 10. x ← [N],   x ∈ {1, 2, ... , N}
    x = RandomBnd(N0) + 1;

    
    // 11. (s_0, w) ← GSampler(a1, a2, B, s_goth, f(x) + u),   (s_0, w) ∈ Z^(2d) × R^m_q
    s_0.SetLength(2*d0);    

    // Compute  f(x) + u
    fx_u = Compute_f(B_f, x) + u;
    
    

        // B ← [rot(g) −rot(f) 
    //      rot(G) −rot(F)] ∈ Z^(2d×2d)
    B.SetDims(2*d0, 2*d0);
    A.SetDims(d0, d0);

    rot(A, g); 

    for(i=0; i<d0; i++)
    {
        for(j=0; j<d0; j++)
        {
            B[i][j] = A[i][j];
        }
    }
    
    rot(A, -f); 

    for(i=0; i<d0; i++)
    {
        for(j=0; j<d0; j++)
        {
            B[i][j+d0] = A[i][j];
        }
    }

    rot(A, G); 

    for(i=0; i<d0; i++)
    {
        for(j=0; j<d0; j++)
        {
            B[i+d0][j] = A[i][j];
        }
    }

    rot(A, -F); 

    for(i=0; i<d0; i++)
    {
        for(j=0; j<d0; j++)
        {
            B[i+d0][j+d0] = A[i][j];
        }
    }

    A.kill();

    // Gaussian sampling  
    GSampler(s_0, w, ipk.a1, ipk.a2, B, sigma0, fx_u);
        

    // 12. ρ_2 ← (s_0, w, x),   ρ_2 ∈ Z^(2d) × R^m × N

    // 13. return ρ_2    
    // NOTE: (s_0, w, x) provided as separate outputs         
}
