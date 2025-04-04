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

#include "Holder.h"


//==============================================================================
// H_Init  -    Initialization of common random string and matrices 
// 
// Input:
// - inputStr:  string containing the input message (initial seed)
//
// Outputs:
// - crs:       structure with the pair (crs_ISIS, crs_Com), generated by Hcrs
// - attrs:     attributes
// - idx_hid:   | \overline{\idx} |, number of undisclosed attributes
// - idx_pub:   | idx |, number of disclosed attributes 
//==============================================================================
void H_Init(CRS2_t& crs, Vec<string>& attrs, const string& inputStr)
{    
    // NOTE: assuming that current modulus is q0 (not q_hat)
    unsigned long   i;
    
    // Generation of crs structure, using H_crs custom Hash function     
    Hcrs(crs, inputStr);
    // NOTE: crs contains 3D uniformly random matrices mod q_hat

    // Initialize attributes
    attrs.SetLength(l0);

    for(i=0; i<l0; i++)
    {
        attrs[i] = to_string(i+1) + "-" + inputStr; 
        // NOTE: dummy attributes (l0 = 8)
    }

    // return(crs, attrs, idx_pub, idx_hid)        
}


//==============================================================================
// H_VerCred1   -   Holder.VerCred1 function
// 
// Inputs:
// - inputStr:      string containing the initial seed for crs
// - crs:           structure with the pair (crs_ISIS, crs_Com), generated by Hcrs
// - ipk:           Issuer public key
// - attrs:         attributes
// - idx_hid:       | \overline{\idx} |, number of undisclosed attributes
// - idx_pub:       | idx |, number of disclosed attributes 
// 
// Outputs:
// - u, Pi_ptr:     commitment u and proof π, corresponding to the structure ρ_1 
// - state:         structure that contains the polynomial vectors m and r
//==============================================================================
void H_VerCred1(zz_pX& u, uint8_t** Pi_ptr, STATE_t& state, const string& inputStr, const CRS2_t& crs, const IPK_t& ipk, const Vec<string>& attrs)
{
    // NOTE: assuming that current modulus is q0 (not q_hat)
    unsigned long   i, j, k;
    vec_zz_pX       c0, c1;
    vec_ZZX         mex, r;
    vec_ZZ          m_i, coeffs_m, coeffs_r, s;
    mat_zz_p        P0, P1, P; 
    vec_zz_p        u_vect, prod;
    ZZ              range, B_goth2;
    long            mul;

    const unsigned long idxhlrd = (idx_hid * h0) + (lr0 * d0); //|idx_hid|·h + ℓr·d

       
    // 1. (a_1, ... , a_l) ← attrs,  a_i ∈ {0, 1}∗
    // NOTE: l0 = idx_hid + idx_pub = len(attrs),  d0 must divide l0*h0
    // NOTE: for every variable of l0 elements, the first are the idx_hid elements, the last are the idx_pub elements
    
    // 2. (c0, c1) ← ipk,   (c0, c1) ∈ R^ℓm_q × R^ℓr_q
    c0 = ipk.c0;
    c1 = ipk.c1;

    // 3. m ← Coeffs^−1( H_M(a1), ... , H_M(a_l) ) ∈ R^ℓm
    mex.SetLength(lm0);    
    coeffs_m.SetLength(l0 * h0);
    k = 0;

    for(i=0; i<l0; i++)
    {                  
        // a_i = attrs[i];        
        HM(m_i, attrs[i] );        

        for(j=0; j<h0; j++)     
        {
            coeffs_m[k] = m_i[j];
            k++;
        }
    }    

    CoeffsInvX(mex, coeffs_m, lm0);


    // 4. r ← S^ℓr_ψ,   r ∈ R^ℓr
    r.SetLength(lr0);   
    range = 2*psi0 + 1;
    
    for(i=0; i<lr0; i++)
    {
        r[i].SetLength(d0);
        
        for(j=0; j<d0; j++)     
        {
            r[i][j] = RandomBnd(range) - psi0;
            // NOTE: each coefficient is in the range [−psi0, psi0];
        }
    }


    // 5. u ← c0^T * m + c1^T * r ∈ R_q
    u.SetLength(d0);
    u = poly_mult(c0, conv<vec_zz_pX>(mex)) + poly_mult(c1, conv<vec_zz_pX>(r));
    

    // 6. P ← [rot(c0^T)_(idx_hid) | rot(c1^T)],    P ∈ Z_q^(d × (|idx_hid|·h + ℓr·d))   
    P.SetDims(d0, (idxhlrd + d_hat));
    // NOTE: zero padding of P (d_hat columns) anticipated here, from Preprocessing_Com

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
    
                
    // 7. s ← (Coeffs(m)_(idx_hid), Coeffs(r)),   s ∈ Z^(|idx_hid|·h + ℓr·d)
    s.SetLength(idxhlrd + d_hat);
    // NOTE: zero padding of s (d_hat values) anticipated here, from Preprocessing_Com
    
    coeffs_r.SetLength(lr0*d0);
    CoeffsX(coeffs_r, r, lr0);
    
    // NOTE: only first idx_hid*h0 coeffs of m (corresponding to undisclosed attributes) 
    //       are copied into s, while coeffs of r is fully copied into s.
    k = 0;

    for(j=0; j<(idx_hid*h0); j++)
    {
        s[k] = coeffs_m[j];
        k++;
    }    

    for(j=0; j<(lr0*d0); j++)
    {
        s[k] = coeffs_r[j];
        k++;
    }


    // 8. u ← Coeffs(u) − rot(c0^T)_idx * Coeffs(m)_idx ∈ Z_q^d
    u_vect.SetLength(d0); 
    prod.SetLength(d0);
    
    // NOTE: only last idx_pub*h0 columns of P0 and coeffs_m (corresponding to disclosed attributes) 
    //       are considered in the product rot(c0^T)_idx * Coeffs(m)_idx 
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


    // 9. π ← Prove_Com^HCom (crs_Com, (q1_hat/q·P, q1_hat/q·u_vect, ψ·sqrt(h·|idx_hid| + ℓr·d, s)
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
   
        Prove_Com(Pi_ptr, inputStr, crs[1], ipk, (mul * P), (mul * u_vect), B_goth2, s);
        // NOTE: P, u are converted from modulo q0 to q1_hat
    }

    P.kill();
   

    // 10. ρ_1 ← (u, π)
    // NOTE: u & Pi are kept separate in the output, for simplicity
          
    // 11. state ← (m, r) state ∈ R^ℓm × R^ℓr
    state.m = mex;
    state.r = r;
          
    // 12. return (ρ_1, state) 
}


//==============================================================================
// H_VerCred2   -    Holder.VerCred2 function
// 
// Inputs:
// - ipk:            Issuer public key
// - B_f:            public random matrix B_f ∈ Z^(nd×t)_q
// - s_0:            short vector (output of GSampler),  s_0 ∈ Z^(2d)
// - w:              polynomial vector (output of GSampler),  w ∈ R^m
// - x:              random integer, uniformly sampled from the set [N]
// - state:          structure that contains the polynomial vectors m and r
//                   NOTE: (s_0, w, x) correspond to the structure ρ_2
// 
// Outputs:
// - cred = (s,r,x): triple that corresponds to the credential
//==============================================================================
void H_VerCred2(CRED_t& cred, const IPK_t& ipk, const mat_zz_p& B_f, const vec_ZZ& s_0, const vec_ZZX& w, const ZZ& x, const STATE_t& state)
{
    // NOTE: assuming that current modulus is q0 (not q_hat)
    unsigned long   i, j;
    zz_pX           a1, left, right;
    vec_zz_pX       a2, c0, c1, a;
    vec_ZZX         m, r, s;
    ZZ              acc;
    RR              norm_s, norm_r, th_s, th_r;

    // 1. (m, r) ← state,   state ∈ R^(ℓm) × R^(ℓr)
    m = state.m;
    r = state.r;
    
    // 2. (a1, a2, c0, c1) ← ipk,   ipk ∈ R_q × R^m_q × R^(ℓm)_q × R^(ℓr)_q
    a1 = ipk.a1;
    a2 = ipk.a2;
    c0 = ipk.c0;
    c1 = ipk.c1;
    
    // 3. (s_0, w, x) ← ρ2,   ρ_2 ∈ Z^(2d) × R^m × N
    // NOTE: (s_0, w, x) provided as separate inputs  

    // 4. s ← [Coeffs^(−1)(s_0) | w],   s ∈ R^(m+2)
    // NOTE: s_0 and s are different from s in Holder.VerCred1
    s.SetLength(m0+2);

    for(i=0; i<2; i++)
    {
        s[i].SetLength(d0);

        for(j=0; j<d0; j++)     
        {
            // s[i][j] =  s_0[d0*i + j] );
            SetCoeff(s[i], j, s_0[d0*i + j] );
        }        
    }     

    for(i=0; i<m0; i++)
    {
        s[i+2] = w[i];
    }
  
   
    // 5. if {∥s∥ > sigma0·√((m + 2)d))} ∨ {∥r∥ > ψ·√(ℓr·d)} ∨ {[1|a1|a2^T]*s != f(x) + c0^T*m + c1^T*r}
    norm_s = sqrt( conv<RR>( Norm2X(s, d0) ) );
    th_s   = conv<RR>(sigma0) * sqrt(conv<RR>( (m0+2)*d0) );

    norm_r = sqrt( conv<RR>( Norm2X(r, d0) ) );
    th_r   = conv<RR>(psi0) * sqrt(conv<RR>( lr0*d0 ) );

    // a ← [1|a1|a2^T]
    a.SetLength(m0+2);
    
    a[0].SetLength(d0);
    a[0] = zz_pX(1); 

    a[1] = a1;

    for(i=0; i<m0; i++)
    {
        a[2+i] = a2[i];
    }

    // left ← [1|a1|a2^T] * s = a * s,   left ∈ R_q
    left = poly_mult(a, conv<vec_zz_pX>(s) );
    left.normalize();

    // right ← f(x) + c0^T * m + c1^T * r,   right ∈ R_q    
    right = Compute_f(B_f, x) + poly_mult(c0, conv<vec_zz_pX>(m)) + poly_mult(c1, conv<vec_zz_pX>(r) );
    right.normalize();
    
    cred.valid = 0;

    // 5.1 if {∥s∥ > sigma0·√((m + 2)d))} 
    if (norm_s > th_s)
    {
        // 6.    return ⊥
        cout << "First  condition failed - Invalid s norm!" << endl;
        // cout << " norm_s = " << norm_s << " > " << th_s << endl;
        return;
    }
    
    // 5.2 ... or {∥r∥ > ψ·√(ℓr·d)} 
    else if (norm_r > th_r)
    {
        // 6.    return ⊥
        cout << "Second condition failed - Invalid r norm!" << endl;
        // cout << " norm_r = " << norm_r << " > " << th_r << endl;        
        return;               
    }

    // 5.3 ... or {[1|a1|a2^T]*s != f(x) + c0^T*m + c1^T*r}
    else if (left != right)
    {
        // 6.    return ⊥
        cout << "Third  condition failed - left != right!" << endl;
        // cout << " " << left << " != " << right << endl;        
        return; 
    }
    
    // 7. else
    else
    {
        // 8. cred ← (s, r, x),   (s, r, x) ∈ R^(m+2)_q × R^(ℓr) × N
        cred.s = s;
        cred.r = r;
        cred.x = x;
        cred.valid = 1;
    }

    // 9. return cred    
}


//==============================================================================
// H_VerPres    -    Holder.VerPres function
// 
// Inputs:
// - cred = (s,r,x): triple that corresponds to the credential
// - inputStr:       string containing the initial seed for crs
// - crs:            structure with the pair (crs_ISIS, crs_Com), generated by Hcrs
// - ipk:            Issuer public key
// - B_f:            public random matrix B_f ∈ Z^(nd×t)_q
// - attrs:          attributes
// - idx_hid:        | \overline{\idx} |, number of undisclosed attributes
// - idx_pub:        | idx |, number of disclosed attributes 
// 
// Output:
// - VP:             structure for the Verifiable Presentation
//==============================================================================
void H_VerPres(VP_t& VP, const CRED_t& cred, const string& inputStr, const CRS2_t& crs, const IPK_t& ipk, const mat_zz_p& B_f, const Vec<string>& attrs)
{
    // NOTE: assuming that current modulus is q0 (not q_hat)
    unsigned long   i, j, k;
    zz_pX           a1;
    vec_zz_pX       a2, c0, c1, a; 
    vec_ZZX         s, r; //mex
    ZZ              x;   
    vec_ZZ          m_i, coeffs_m, coeffs_s, coeffs_r, r_vect, coeffs_u;
    vec_zz_p        coeffs_m_idx;
    mat_zz_p        P, C0, C1, C; 
    vec_ZZ          Bounds;
    Vec<vec_ZZ>     sig;
    long            mul, Pi_valid;

    const unsigned long m2d     = (m0 + 2)*d0;    // (m+2)·d
    const unsigned long lmlrd   = (lm0 + lr0)*d0; // (ℓm+ℓr)·d
    const unsigned long idxhlrd = (idx_hid * h0) + (lr0 * d0); //|idx_hid|·h + ℓr·d

       
    // 1. (a_1, ... , a_l) ← attrs,   a_i ∈ {0, 1}∗
    // NOTE: l0 = idx_hid + idx_pub = len(attrs),  d0 must divide l0*h0
    // NOTE: for every variable of l0 elements, the first are the idx_hid elements, the last are the idx_pub elements

    // 2. (a1, a2, c0, c1) ← ipk,   ipk ∈ R_q × R^m_q × R^(ℓm)_q × R^(ℓr)_q
    a1 = ipk.a1;
    a2 = ipk.a2;
    c0 = ipk.c0;
    c1 = ipk.c1;

    // 3. (s, r, x) ← cred,   cred ∈ R^(m+2) × R^(ℓr) × N
    if (cred.valid)
    {
        s = cred.s;        
        r = cred.r;
        x = cred.x;
    }
    else
    {
        cout << "\n Invalid credential!" << endl;
        return;
    }

      
    // 4. m ← Coeffs^−1( H_M(a1), ... , H_M(a_l) ) ∈ R^ℓm
    // mex.SetLength(lm0);    
    coeffs_m.SetLength(l0 * h0);
    k = 0;

    for(i=0; i<l0; i++)
    {                  
        // a_i = attrs[i];        
        HM(m_i, attrs[i] );

        for(j=0; j<h0; j++)     
        {
            coeffs_m[k] = m_i[j];
            k++;
        }
    }    

    // mex = CoeffsInvX(coeffs_m, lm0);
    // NOTE: coeffs_m is directly used instead of mex

    
    // a ← [1|a1|a2^T]  
    a.SetLength(m0+2);
    
    a[0].SetLength(d0);
    a[0] = zz_pX(1);  
       
    a[1] = a1;

    for(i=0; i<m0; i++)
    {
        a[2+i] = a2[i];
    } 
          
          
    // 5. P ← rot([1|a1|a2^T]) = rot(a),   P ∈ Z^(d×(m+2)d)_q      
    P.SetDims(d0, (m2d + d_hat)); 
    // NOTE: zero padding of P (d_hat columns) anticipated here, from Preprocessing_ISIS   
    rot_vect(P, a);
    

    // 6. C ← [rot(c0^T)_(idx_pub) | rot(c0^T)_(idx_hid) | rot(c1^T)],   C ∈ Z_q^(d × (ℓm+ℓr)d)
    C.SetDims(d0, (lmlrd + d_hat));
    // NOTE: zero padding of C (d_hat columns) anticipated here, from Preprocessing_ISIS

    C0.SetDims(d0, lm0*d0);
    rot_vect(C0, c0);

    // NOTE: first copy in C the columns for disclosed attributes, then those for undisclosed attributes
    // NOTE: lm0*d0 = l0*h0 = (idx_pub + idx_hid) * h0       
    k = 0;

    for(j=0; j<(idx_pub*h0); j++)
    {
        for(i=0; i<d0; i++)
        {   
            C[i][k] = C0[i][idx_hid*h0+j];
        }
        k++;
    }    
   
    for(j=0; j<(idx_hid*h0); j++)
    {
        for(i=0; i<d0; i++)
        {   
            C[i][k] = C0[i][j];
        }
        k++;
    } 

    C0.kill();   

    C1.SetDims(d0, lr0*d0);  
    rot_vect(C1, c1);    

    for(j=0; j<(lr0*d0); j++)
    {
        for(i=0; i<d0; i++)
        {   
            C[i][k] = C1[i][j];
        }
        k++;
    }

    C1.kill(); 


    // 7. m ← Coeffs(m)_(idx_pub),   m ∈ Z_q^(idx_pub·h)
    coeffs_m_idx.SetLength(idx_pub*h0);

    for(i=0; i<(idx_pub*h0); i++)
    {
        coeffs_m_idx[i] = conv<zz_p>(coeffs_m[(idx_hid*h0)+i]);
    }


    // 8. s ← Coeffs(s),   s ∈ Z^((m+2)d)
    coeffs_s.SetLength(m2d + d_hat);
    // NOTE: zero padding of coeffs_s (d_hat values) anticipated here, from Preprocessing_ISIS
    CoeffsX(coeffs_s, s, (m0+2));
    
    
    // 9. r ← (Coeffs(m)_(idx_hid), Coeffs(r)),   r ∈ Z^(|idx_hid|·h + ℓr·d)
    r_vect.SetLength(idxhlrd + d_hat);
    // NOTE: zero padding of r_vect (d_hat values) anticipated here, from Preprocessing_ISIS
    
    coeffs_r.SetLength(lr0*d0);    
    CoeffsX(coeffs_r, r, lr0);
    
    // NOTE: only first idx_hid*h0 coeffs of m (corresponding to undisclosed attributes) 
    //       are copied into r_vect, while coeffs of r is fully copied into r_vect.
    k = 0;

    for(j=0; j<(idx_hid*h0); j++)
    {
        r_vect[k] = coeffs_m[j];
        k++;
    }    

    for(j=0; j<(lr0*d0); j++)
    {
        r_vect[k] = coeffs_r[j];
        k++;
    }


    // 10. u ← enc(x) ∈ {0, 1}^t  
    // Compute enc(x), the binary decomposition of (x−1) 
    coeffs_u.SetLength(t0);

    for(i=0; i<t0; i++)
    {
        coeffs_u[i] = bit(x-1, i);
    }
    

    // 11. Bounds ← ( sigma0·√((m + 2)d), ψ·√(h·|idx| + ℓr·d) ),    Bounds ∈ Z^2
    Bounds.SetLength(2);
    // Bounds[0] = conv<RR>(sigma0) * sqrt(conv<RR>( (m0+2)*d0) );
    // Bounds[1] = psi0 * sqrt(conv<RR>( idxhlrd ));
    Bounds[0] = sqr( ZZ(sigma0) ) * ZZ( (m0+2)*d0 );
    Bounds[1] = sqr( ZZ(psi0)   ) * ZZ( idxhlrd   );

      
    // 12. VP ← emptyVP()
    // NOTE: see VP_t

    // 13. VP.cp ← VC.cp
    // VP.cp = VC.cp;

    // 14. VP.ipk ← VC.ipk
    VP.ipk = ipk;

    // 15. VP.attrs′ ← (attrs_{idx_pub} | {0}_{idx_hid})
    VP.attrs_prime = attrs; 
    // NOTE: select only disclosed attributes and fill with zeros for each i ∈ idx_hid 
    for(i=0; i<idx_hid; i++)
    {
        VP.attrs_prime[i] = "0"; // Zero padding
    }
    // cout << "  attrs  = " << attrs << endl;
    // cout << "  attrs' = " << VP.attrs_prime << endl;


    // 16. VP.idx ← idx
    // VP.idx = idx_pub; 
         
        
    // 17. VP.pi ← Prove^HISIS_ISIS( crs_ISIS, ˆq2/q · P, ˆq2/q · C, m, ˆq2/q · B_f, Bounds, idx), (s, r, u) )
    if (not(divide( ZZ(q2_hat), q0)))
    {
        cout << " ERROR: q2_hat must be divisible by q! " << endl;
    }  

    mul = long(q2_hat) / long(q0);

    sig.SetLength(3);
    sig[0] = coeffs_s;
    sig[1] = r_vect;
    sig[2] = coeffs_u;

    {
        zz_pPush push(q2_hat);
        // NOTE: backup current modulus q0, temporarily set to q2_hat (i.e., zz_p::init(q2_hat)) 
    
        Prove_ISIS(&(VP.Pi), inputStr, crs[0], ipk, (mul * P), (mul * C), coeffs_m_idx, (mul * B_f), Bounds, idx_pub, sig );
        // NOTE: P, C, B_f are converted from modulo q0 to q2_hat
    }

    P.kill();
    C.kill();
    
    Pi_valid = VP.Pi[0];

    // 18. return VP  
    if (Pi_valid)
    {
        VP.valid = 1;
    }
    else
    {
        VP.valid = 0;
    }
}