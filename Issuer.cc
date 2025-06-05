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
#include "serialize.h"


//==============================================================================
// I_KeyGen - Issuer.KeyGen function: generate the Issuer Public Key and 
//            Issuer Secret Key from parameters d and q.
//
// Input:
// - None
//
// Outputs:
// - ipk_prt:    Issuer Public Key, pointer to serialized bytes (a1 + seed_ipk)
// - f, g, F, G: Issuer Secret Key (i.e. base polynomials to build matrix B)
//
// NOTE: Issuer.KeyGen in BLNS pseudocode, corresponding to
//       Fig. 18: AnonCreds.Init, pag. 52 in [BLNS23]
//==============================================================================
void I_KeyGen(uint8_t** ipk_prt, ZZX& f, ZZX& g, ZZX& F, ZZX& G)
{
    zz_pX           a1;
    size_t          len_a1, len_ipk;
    uint8_t         *ipk_bytes;
    unsigned char   *seed_bytes;

    const int   nbits   = ceil(log2(conv<double>(q0-1)));


#ifdef ENABLE_FALCON
    // Keygen algorithm from the Falcon reference implementation 
    Falcon_keygen(a1, f, g, F, G);
#else
    // NTRU.TrapGen(q, d) algorithm in [BLNS23]
    NTRU_TrapGen(a1, f, g, F, G);    
#endif
    // NOTE: a1 is a Polynomial with d coefficients modulo q (i.e. h in [DLP14])


    // Allocate a vector of bytes to store ipk (a1 + seed_ipk)
    len_a1    = calc_ser_size_poly_minbyte(d0, nbits); // zz_pX
    len_ipk   = len_a1 + SEED_LEN;
    *ipk_prt  = new uint8_t[len_ipk];
    ipk_bytes = *ipk_prt;
    cout << "  Size ipk: " << (len_ipk/1024.0) << " KiB" << endl; // 1 KiB kibibyte = 1024 bytes
    
    // Serialize a1 (first bytes in ipk_bytes)
    serialize_minbyte_poly_zz_pX(ipk_bytes, len_a1, d0, nbits, a1);   

    // Initialize a 32 byte (256 bit) public seed for completing ipk (i.e. a2, c0, c1),
    // using the cryptographically strong pseudo-random number generator from NTL.
    // Store seed_ipk in the last bytes (SEED_LEN) of ipk_bytes
    seed_bytes = reinterpret_cast<unsigned char*>(ipk_bytes + len_a1);
    RandomStream& RS = GetCurrentRandomStream();    
    RS.get(seed_bytes, SEED_LEN);
    // for(long i=0; i<SEED_LEN; i++)
    // {
    //     printf("%0x", seed_bytes[i]);
    // }
    // printf("\n");
              
    // Output Issuer Public Key and Issuer Secret Key (i.e. B)
    // CompleteIPK(ipk, ipk_bytes);
    // ipk ← (a1, a2, c0, c1)
    // isk ← (f, g, F, G)
}


//==============================================================================
// CompleteIPK - CompleteIPK function: deserialize the Issuer Public Key
//               and complete it by generating a2, c0, c1 from a public seed.
//
// Input:
// - ipk_bytes: Serialized Issuer Public Key (a1 + seed_ipk)
//
// Output:
// - ipk:       Deserialized IPK (a1, a2, c0, c1 vectors of polynomials)
//==============================================================================
void CompleteIPK(IPK_t& ipk, const uint8_t* ipk_bytes)
{   
    // NOTE: assuming that current modulus is q0

    unsigned long   i;
    HASH_STATE_t    *state;  
    size_t          len_a1;

    // Compute the minimum number of bits and bytes to represent each coefficient of a1, a2, c0, c1
    const int    nbits    = ceil(log2(conv<double>(q0-1)));
    const size_t b_coeffs = ceil(log2(conv<double>(q0-1)) / 8.0);
    

    // Deserialize a1 (first bytes in ipk_bytes)
    len_a1 = calc_ser_size_poly_minbyte(d0, nbits); // zz_pX
    deserialize_minbyte_poly_zz_pX(ipk.a1, d0, nbits, ipk_bytes, len_a1);
    // NOTE: a1 is a Polynomial with d coefficients modulo q (i.e. h in [DLP14])

    // Retrieve seed_ipk (last bytes in ipk_bytes) to generate a2, c0, c1
    for(i=0; i<SEED_LEN; i++)
    {
        ipk.seed_ipk[i] = ipk_bytes[len_a1 + i];
        // printf("%0x", ipk.seed_ipk[i]);
    }
    // printf("\n");
        
    // Initialize the Hash function with seed_ipk, to generate a2, c0, c1
    state = Hash_Init(ipk.seed_ipk, SEED_LEN);

    ipk.a2.SetLength(m0);   
    // NOTE: a2 is a vector of m Polynomials with d coefficients modulo q

    for(i=0; i<m0; i++)
    {
        Hash_zz_pX(ipk.a2[i], state, d0, b_coeffs);
    }

    ipk.c0.SetLength(lm0);  
    // NOTE: c0 is a vector of l_m Polynomials with d coefficients modulo q

    for(i=0; i<lm0; i++)
    {
        Hash_zz_pX(ipk.c0[i], state, d0, b_coeffs);
    }

    ipk.c1.SetLength(lr0);  
    // NOTE: c1 is a vector of l_r Polynomials with d coefficients modulo q

    for(i=0; i<lr0; i++)
    {
        Hash_zz_pX(ipk.c1[i], state, d0, b_coeffs);
    }

    delete state;

    // Return ipk ← (a1, a2, c0, c1)
}



//==============================================================================
// I_VerCred    -   Issuer.VerCred function
// 
// Inputs:
// - seed_crs:      initial public seed for crs structure
// - crs:           structure with the pair (crs_ISIS, crs_Com), generated by Hcrs
// - B_f:           public random matrix B_f ∈ Z^(nd×t)_q
// - ipk_bytes:     serialized Issuer Public Key
// - f, g, F, G:    base polynomials used to build Issuer Secret Key (i.e. matrix of integers B)
// - attrs_prime:   disclosed attributes (attrs′)
// - idx_pub:       |idx|, number of disclosed attributes
// - Rho1:          structure ρ_1 that contains the commitment u and proof π
// 
// Output:
// - Rho2_ptr:      pointer to the structure ρ_2 = (s_0, w, x) where:
//      * s_0:      short vector (output of GSampler),       s_0 ∈ Z^(2d)
//      * w:        polynomial vector (output of GSampler),  w ∈ R^m
//      * x:        random integer, uniformly sampled from the set [N]
//==============================================================================
void I_VerCred(uint8_t** Rho2_ptr, const unsigned char* seed_crs, const CRS2_t& crs, const mat_zz_p& B_f, const uint8_t* ipk_bytes, const ZZX& f, const ZZX& g, const ZZX& F, const ZZX& G, const Vec<string>& attrs_prime, RHO1_t& Rho1)
{    
    // NOTE: assuming that current modulus is q0 (not q_hat)
    unsigned long   i, j, k, result;
    IPK_t           ipk;
    vec_ZZ          s_0, m_i, coeffs_m;
    vec_ZZX         w;
    ZZ              x, B_goth2;
    zz_pX           u, fx_u;
    mat_zz_p        P0, P1, P; 
    vec_zz_p        u_vect, prod;
    long            mul;
    size_t          len_u, len_s0, len_w, len_x, len_Rho2;
    uint8_t        *Rho2_bytes;
    
    const unsigned long idxhlrd = (idx_hid * h0) + (lr0 * d0); //|idx_hid|·h + ℓr·d
    const int           nbits   = ceil(log2(conv<double>(q0-1)));

    
    // 1. (a'_1, ... , a'_k) ← attrs',  a'_i ∈ {0, 1}∗
    // NOTE: l0 = idx_hid + idx_pub = len(attrs),  d0 must divide l0*h0
    // NOTE: for every variable of l0 elements, the first are the idx_hid elements, the last are the idx_pub elements
    
    // 2. (a1, a2, c0, c1) ← ipk,   ipk ∈ R_q × R^m_q × R^ℓm_q × R^ℓr_q
    CompleteIPK(ipk, ipk_bytes);

    // 3. (u, π) ← ρ1   
    len_u = calc_ser_size_poly_minbyte(d0, nbits);
    // cout << "  Size u: " << (len_u/1024.0) << " KiB" << endl; // 1 KiB kibibyte = 1024 bytes

    // Deserialize u
    deserialize_minbyte_poly_zz_pX(u, d0, nbits, Rho1.u, len_u);
    // Rho1.u += len_u;
    
    // Free the vector with serialized u
    delete[] Rho1.u;


    // 4. B ← isk,   B ∈ Z^(2d×2d)
    // NOTE: using (f, g, F, G) instead of B

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

    // vec_zz_pX mex_prime = CoeffsInv(coeffs_m, lm0);
    // NOTE: coeffs_m is directly used instead of mex_prime


    // 6. P ← [rot(c0^T)_(idx_hid) | rot(c1^T)],    P ∈ Z_q^(d × (|idx_hid|·h + ℓr·d))     
    P.SetDims(d0, (idxhlrd + d_hat));
    // NOTE: zero padding of P (d_hat columns) anticipated here, from Verify_Com

    P0.SetDims(d0, lm0*d0);
    rot_vect(P0, ipk.c0);

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
    rot_vect(P1, ipk.c1);

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
    B_goth2 = sqr(ZZ(psi0)) * ZZ(idxhlrd); // B_goth^2
    
    {
        zz_pPush push(q1_hat); 
        // NOTE: backup current modulus q0, temporarily set to q1_hat (i.e., zz_p::init(q1_hat))

        result = Verify_Com(seed_crs, crs[1], ipk, (mul * P), (mul * u_vect), B_goth2, &(Rho1.Pi));
        // NOTE: P, u_vect are converted from modulo q0 to q1_hat
        // NOTE: Verify_Com deserializes the proof π in Rho1.Pi
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

    
    // 11. (s_0, w) ← GSampler(a1, a2, B, s_goth, f(x) + u),   (s_0, w) ∈ Z^(2d) × R^m
    s_0.SetLength(2*d0);    

    // Compute  f(x) + u
    fx_u = Compute_f(B_f, x) + u;
    
    
    #ifdef ENABLE_FALCON

        Falcon_GSampler(s_0, w, ipk.a1, ipk.a2, f, g, F, G, fx_u);

    #else

        mat_L           A, B;

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
        GSampler(s_0, w, ipk.a1, ipk.a2, B, fx_u);

        B.kill();
            
    #endif
    
    
    // 12. ρ_2 ← (s_0, w, x),   ρ_2 ∈ Z^(2d) × R^m × N
    
    // Compute the number of bytes for each component of the structure ρ_2
    len_s0 = calc_ser_size_vec_ZZ(2*d0);    // vec_ZZ (2*d0*long)     -  8192 bytes
    len_w  = calc_ser_size_vec_ZZX(m0, d0); // vec_ZZX (m0*d0*long)   - 12288 bytes
    len_x  = calc_ser_size_big_ZZ(t0);      // big ZZ (t0 = 512 bits) -    64 bytes    
    len_Rho2 = len_s0 + len_w + len_x;      //                          20544 bytes
    cout << "  Size Rho2: " << (len_Rho2/1024.0) << " KiB" << endl; // 1 KiB kibibyte = 1024 bytes
  
    // Allocate a vector of bytes to store the structure ρ_2
    *Rho2_ptr = new uint8_t[len_Rho2];
    Rho2_bytes = *Rho2_ptr;

    // Serialize (s_0, w, x) in ρ_2
    serialize_vec_ZZ(Rho2_bytes, len_s0, 2*d0, s_0);
    Rho2_bytes+= len_s0;
    serialize_vec_ZZX(Rho2_bytes, len_w, m0, d0, w);
    Rho2_bytes+= len_w;
    serialize_big_ZZ(Rho2_bytes, len_x, x-1);    
    // Rho2_bytes+= len_x;


    // 13. return ρ_2
}
