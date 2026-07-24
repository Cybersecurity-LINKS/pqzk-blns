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

#include "Hash.h"


//==============================================================================
// Hash_Init  - Initialize the Custom Hash function, implemented using SHAKE128
// 
// Inputs:
// - v:         vector of bytes containing the input message (initial seed)
// - len:       length of v (number of bytes)
//
// Output:
// - state:     status structure
//==============================================================================
HASH_STATE_t* Hash_Init(const uint8_t* v, const size_t len)
{
    HASH_STATE_t *state = new HASH_STATE_t();

    _shake128_init(state);
    _shake128_absorb(state, v, len); 

    return state;
}

void Hash_Init_To(HASH_STATE_t& state, const uint8_t* v, const size_t len)
{
    _shake128_init(&state);
    _shake128_absorb(&state, v, len);
}


//==============================================================================
// Hash_Update - Update the Custom Hash function with a new input message
// 
// Inputs:
// - state:      status structure
// - v:          vector of bytes containing a new input message
// - len:        length of v (number of bytes)
//
// Output:
// - state:      updated status structure
//==============================================================================
void Hash_Update(HASH_STATE_t *state, const uint8_t* v, const size_t len)
{
    _shake128_absorb(state, v, len);

    // return state;
}


//==============================================================================
// Hash_Copy  - Copy the status structure of the Custom Hash function
// 
// Inputs:
// - state0:    initial status structure
//
// Output:
// - state:     copy of the status structure
//==============================================================================
HASH_STATE_t* Hash_Copy(const HASH_STATE_t *state0)
{
    HASH_STATE_t *state = new HASH_STATE_t();

    copy(state0->s, state0->s + 25, state->s);
    // state->s[25] = state0->s[25];
    state->pos      = state0->pos;
    state->final    = state0->final;
     
    return state;
}


void Hash_Copy_To(HASH_STATE_t& dst, const HASH_STATE_t& src)
{
    copy(src.s, src.s + 25, dst.s);
    dst.pos   = src.pos;
    dst.final = src.final;
}


//==============================================================================
// Hash_zz_pX - Generate a random polynomial using the Custom Hash function
// 
// Inputs:
// - state:     status structure
// - n_coeffs:  number of coefficients of the random polynomial (i.e. d_hat)
// - b_coeffs:  number of bytes for each coefficient (i.e. |q_hat-1|)
//
// Output:
// - out_poly:  random polynomial with n_coeffs coefficients (mod q_hat)
// - state:     updated status structure
//==============================================================================
void Hash_zz_pX(zz_pX& out_poly, HASH_STATE_t *state, const long& n_coeffs, const size_t& b_coeffs)
{    
    // NOTE: the current modulus (q_hat or q0) must already be set by the calling function

    long i;
    unsigned char* y_arr = new unsigned char[b_coeffs];   
           
    out_poly.SetLength(n_coeffs);

    for(i=0; i < n_coeffs; i++)
    {
        _shake128_squeeze(state, y_arr, b_coeffs);
        
        out_poly[i] = conv<zz_p>(ZZFromBytes(y_arr, b_coeffs)); 
    }
    
    out_poly.normalize();
    
    delete[] y_arr;   

    // return out_poly;
}

void Hash_zz_pX_buf(
    zz_pX& out_poly,
    HASH_STATE_t* state,
    const long n_coeffs,
    const size_t b_coeffs,
    unsigned char* buf)
{
    long i;
    out_poly.SetLength(n_coeffs);

    for (i = 0; i < n_coeffs; i++)
    {
        _shake128_squeeze(state, buf, b_coeffs);
        out_poly[i] = conv<zz_p>(ZZFromBytes(buf, b_coeffs));
    }

    out_poly.normalize();
}


//==============================================================================
// Hash_v_zz_p - Generate a random vector (mod q_hat) using Custom Hash function
// 
// Inputs:
// - state:     status structure
// - n_elems:   number of elements of the random vector (i.e. 256+n+1)
// - b_num:     number of bytes for each random number (i.e. |q_hat-1|)
//
// Output:
// - out_vec:   vector of random numbers (modulo q_hat)
// - state:     updated status structure
//==============================================================================
void Hash_v_zz_p(vec_zz_p& out_vec, HASH_STATE_t *state, const long& n_elems, const size_t& b_num)
{    
    // NOTE: the current modulus (q_hat or q0) must already be set by the calling function
    long      i;
    unsigned char* y_arr = new unsigned char[b_num];    

    out_vec.SetLength(n_elems);
    
    for(i=0; i < n_elems; i++)
    {
        _shake128_squeeze(state, y_arr, b_num);

        out_vec[i] = conv<zz_p>( ZZFromBytes(y_arr, b_num) );
    }
    
    delete[] y_arr;   

    // return out_vec;
}

void Hash_v_zz_p_buf(
    vec_zz_p& out_vec,
    HASH_STATE_t* state,
    const long n_elems,
    const size_t b_num,
    unsigned char* buf)
{
    long i;
    out_vec.SetLength(n_elems);

    for (i = 0; i < n_elems; i++)
    {
        _shake128_squeeze(state, buf, b_num);
        out_vec[i] = conv<zz_p>(ZZFromBytes(buf, b_num));
    }
}


//==============================================================================
// Hash_R_goth - Generate a random vector for R_goth using the Custom Hash function
// 
// Inputs:
// - state:     status structure
// - n_elems:   number of elements of the random vector (i.e. m1*d_hat)
//
// Output:
// - out:       random vector with n_elems elements in {-1, 0, 1} mod q_hat,
//              equivalent to the pair (R_goth_0 - R_goth_1) in BLNS
// - state:     updated status structure
//==============================================================================
void Hash_R_goth(vec_zz_p& out, HASH_STATE_t *state, const long& n_elems)
{    
    long            i, j, k, curr_byte, R_goth_0, R_goth_1;
    unsigned char*  y_arr;
    
    // Compute the minimum number of bytes needed to fill the vector  
    
    //const long n_bytes = ceil(2*n_elems / 8.0);
    const long n_bytes = (n_elems + 3) >> 2;

    y_arr = new unsigned char[n_bytes];

    _shake128_squeeze(state, y_arr, n_bytes);
           
    // out.SetLength(n_elems);
    k = 0;

    for(i=0; i < n_bytes; i++)
    {
        curr_byte = (long)(y_arr[i]);
        
        // NOTE: each byte will fill 4 elements, 2 bits per element (R_goth_0, R_goth_1)
        for(j=0; j < 4; j++)
        {
            if(k < n_elems)
            {
                R_goth_0 = ( curr_byte & 1 );
                curr_byte = curr_byte >> 1;
                R_goth_1 = ( curr_byte & 1 );
                curr_byte = curr_byte >> 1;
                out[k] = conv<zz_p>(R_goth_0 - R_goth_1);
                // NOTE: each element is in {-1, 0, 1} mod q_hat
            }
            k++;
        }
    } 
    
    delete[] y_arr;

    // return out;
}



//==============================================================================
// Hash_ZZ_xi0 - Generate a random integer modulo (xi0+1) using Custom Hash function
// 
// Inputs:
// - state:      status structure
// - b_num:      number of bytes of the random integer
//
// Output:
// - out:        random integer modulo (xi0+1), i.e. from 0 to xi0
// - state:      updated status structure
//==============================================================================
void Hash_ZZ_xi0(ZZ& out, HASH_STATE_t *state, const size_t& b_num)
{    
    unsigned char* y_arr = new unsigned char[b_num];    

    _shake128_squeeze(state, y_arr, b_num);

    out = (ZZFromBytes(y_arr, b_num)) % (xi0+1);
        
    delete[] y_arr;   

    // return out;
}

void Hash_ZZ_xi0_buf(
    ZZ& out,
    HASH_STATE_t* state,
    const size_t b_num,
    unsigned char* buf)
{
    _shake128_squeeze(state, buf, b_num);

    out = ZZFromBytes(buf, b_num) % (xi0 + 1);
}


//==============================================================================
// Hcrs    -    H_crs, custom Hash function needed in BLNS for crs. 
//              It generates the pair of common random string (crs_ISIS, crs_Com)
//              and the random matrix B_f, from the seed seed_crs.
// 
// Input:
// - seed_crs:      initial public seed for crs structure
// - num_idx_hid:   number of undisclosed attributes (hidden)
//
// Output:
// - crs:           structure with the pair (crs_ISIS, crs_Com)
// - B_f:           random matrix B_f ∈ Z^(d×t)_q
//==============================================================================
void Hcrs(CRS2_t& crs, mat_zz_p& B_f, const uint8_t* seed_crs, const long &num_idx_hid)
{
    long            i, j, n, m1, m2, n256;
    HASH_STATE_t    state;
    size_t          b_coeffs;
       
    Hash_Init_To(state, seed_crs, SEED_LEN);

    // Create the crs structure  
    crs.SetLength(2); 
    crs[0].SetLength(5); // crs_ISIS
    crs[1].SetLength(5); // crs_Com

    // ###########################  crs_ISIS  #########################################
    {
        // NOTE: elements of all matrices in crs_ISIS are mod q2_hat
        zz_pPush push(q2_hat); 
        // NOTE: backup current modulus q0, temporarily set to q2_hat (i.e., zz_p::init(q2_hat))    

        // Compute the minimum number of bytes to represent each coefficient
        b_coeffs = ceil(log2( conv<double>(q2_hat-1) ) / 8.0);    
        
        n    = n_ISIS;
        m1   = (((m0+2)*d0 + (num_idx_hid*h0 + lr0*d0) + t0 + 2*d_hat) / d_hat); // m_1 for Π^ISIS_NIZK
        m2   = m2_ISIS;
        n256 = 256/d_hat;

        if ( (256 % d_hat) != 0)
        {
            cout << "ERROR! 256 must be divisible by d_hat" << endl;
            assert((256 % d_hat) == 0);
        }

        // Create the crs_ISIS structure, i.e. crs[0]
        // crs[0][0] = A_1;
        // crs[0][1] = A_2;
        // crs[0][2] = B_y;
        // crs[0][3] = B_g;
        // crs[0][4] = b;
            
        // Random generation of A_1 ∈ R^(n x m_1)_(q_hat)   
        //                      A_2 ∈ R^(n x m_2)_(q_hat)      
        crs[0][0].SetDims(n, m1);    
        crs[0][1].SetDims(n, m2);

        vector<unsigned char> buf(b_coeffs);

        for(i=0; i<n; i++)
        {
            for(j=0; j<m1; j++)
            {
                Hash_zz_pX_buf(crs[0][0][i][j], &state, d_hat, b_coeffs, buf.data());
            }
            for(j=0; j<m2; j++)
            {
                Hash_zz_pX_buf(crs[0][1][i][j], &state, d_hat, b_coeffs, buf.data());
            }
        }       

        // Random generation of B_y ∈ R^(256/d_hat x m_2)_(q_hat)  
        crs[0][2].SetDims(n256, m2);

        for(i=0; i<n256; i++)
        {
            for(j=0; j<m2; j++)
            {
                Hash_zz_pX_buf(crs[0][2][i][j], &state, d_hat, b_coeffs, buf.data());
            }
        }

        // Random generation of B_g ∈ R^(tau_ISIS^ x m_2)_(q_hat)
        crs[0][3].SetDims(tau_ISIS, m2);

        for(i=0; i<tau_ISIS; i++)
        {
            for(j=0; j<m2; j++)
            {
                Hash_zz_pX_buf(crs[0][3][i][j], &state, d_hat, b_coeffs, buf.data());
            }
        }

        // Random generation of b ∈ R^(m_2)_(q_hat)
        crs[0][4].SetDims(1, m2);
        // NOTE: b is (1 x m_2) matrix, not a vector!

        for(i=0; i<m2; i++)
        {
            Hash_zz_pX_buf(crs[0][4][0][i], &state, d_hat, b_coeffs, buf.data());
        }        
    }

      
    // ###########################  crs_Com  #########################################    
    {
        // NOTE: elements of all matrices in crs_Com are mod q1_hat
        zz_pPush push(q1_hat); 
        // NOTE: backup current modulus q0, temporarily set to q1_hat (i.e., zz_p::init(q1_hat))

        // Compute the minimum number of bytes to represent each coefficient
        b_coeffs = ceil(log2( conv<double>(q1_hat-1) ) / 8.0);

        n    = n_Com;
        m1   = ((num_idx_hid*h0 + lr0*d0 + d_hat)/d_hat); // m_1 for Π^Com_NIZK
        m2   = m2_Com;
        // n256 = 256/d_hat;

        // Create the crs_Com structure, i.e. crs[1]        
        // crs[1][0] = A_1;
        // crs[1][1] = A_2;
        // crs[1][2] = B_y;
        // crs[1][3] = B_g;
        // crs[1][4] = b;
            
        // Random generation of A_1 ∈ R^(n x m_1)_(q_hat)   
        //                      A_2 ∈ R^(n x m_2)_(q_hat)      
        crs[1][0].SetDims(n, m1);    
        crs[1][1].SetDims(n, m2); 
        
        vector<unsigned char> buf(b_coeffs);

        for(i=0; i<n; i++)
        {
            for(j=0; j<m1; j++)
            {
                Hash_zz_pX_buf(crs[1][0][i][j], &state, d_hat, b_coeffs, buf.data());
            }
            for(j=0; j<m2; j++)
            {
                Hash_zz_pX_buf(crs[1][1][i][j], &state, d_hat, b_coeffs, buf.data());
            }
        }       

        // Random generation of B_y ∈ R^(256/d_hat x m_2)_(q_hat)  
        crs[1][2].SetDims(n256, m2);

        for(i=0; i<n256; i++)
        {
            for(j=0; j<m2; j++)
            {
                Hash_zz_pX_buf(crs[1][2][i][j], &state, d_hat, b_coeffs, buf.data());
            }
        }

        // Random generation of B_g ∈ R^(tau_Com^ x m_2)_(q_hat)
        crs[1][3].SetDims(tau_Com, m2);

        for(i=0; i<tau_Com; i++)
        {
            for(j=0; j<m2; j++)
            {
                Hash_zz_pX_buf(crs[1][3][i][j], &state, d_hat, b_coeffs, buf.data());
            }
        }

        // Random generation of b ∈ R^(m_2)_(q_hat)
        crs[1][4].SetDims(1, m2);
        // NOTE: b is (1 x m_2) matrix, not a vector!

        for(i=0; i<m2; i++)
        {
            // for(j=0; j<1; j++)
            {
                Hash_zz_pX_buf(crs[1][4][0][i], &state, d_hat, b_coeffs, buf.data());
            }
        }
    }


    // ###########################  B_f  ######################################### 
    // NOTE: assuming that current modulus is q0 (not q1/q2_hat)

    // Compute the minimum number of bytes to represent each coefficient
    b_coeffs = ceil(log2( conv<double>(q0-1) ) / 8.0);

    // Initialize a random matrix B_f ∈ Z^(d×t)_q
    B_f.SetDims(d0, t0);

    vector<unsigned char> buf(b_coeffs);

    for(i=0; i<d0; i++)
    {
        Hash_v_zz_p_buf(B_f[i], &state, t0, b_coeffs, buf.data());
    }
    
    // return crs, B_f;    
}


//==============================================================================
// Hash_Sigma_R_goth_encoded
//
// Expand pseudorandom bytes from SHAKE128 into one sparse sigma-encoded row of
// R_goth.
//
// Conceptually, the row is a concatenation of m1 polynomial blocks, each of
// length d_hat. Every raw coefficient belongs to {-1, 0, +1} and is derived
// from 2 bits according to:
//
//   00 ->  0
//   01 -> +1
//   10 -> -1
//   11 ->  0
//
// The coefficient is then stored after applying the sigma mapping:
//   raw position 0 -> sigma position 0
//   raw position c -> sigma position d_hat - c, for c != 0
//
// For c != 0, the sign is flipped.
//
// Each non-zero sigma-encoded coefficient is packed into one byte:
//
//   bit 7     : sign   (1 => +1, 0 => -1)
//   bits 0..6 : coefficient position inside the current polynomial block
//
// out.offsets partitions out.entries into m1 contiguous segments, one per
// polynomial block.
//
// Assumption:
//   d_hat <= 128, so positions 0..d_hat-1 fit into 7 bits.
//==============================================================================
void Hash_Sigma_R_goth_encoded(
    TernaryCoeffStructure& out,
    HASH_STATE_t* state,
    const ulong m1,
    unsigned char* y_arr)
{
    assert(d_hat > 0 && d_hat <= 128);
    assert((m1 * d_hat) <= 65535);

    const ulong n_elems = m1 * d_hat;
    const ulong n_bytes = (n_elems + 3) >> 2;

    _shake128_squeeze(state, y_arr, n_bytes);

    ternary_blocks_init(out, m1, n_elems / 2 + 128);

    ulong blk = 0;
    ulong c   = 0;

    for (ulong i = 0; i < n_bytes && blk < m1; i++)
    {
        unsigned char curr = y_arr[i];

        for (int step = 0; step < 4 && blk < m1; step++)
        {
            // Decode one coefficient from the 2 least significant bits.
            // Mapping:
            //   00 ->  0
            //   01 -> +1
            //   10 -> -1
            //   11 ->  0
            const int val = (curr & 1) - ((curr >> 1) & 1);
            curr >>= 2;

            // Append this coefficient after applying sigma.
            ternary_append_sigma(
                out,
                static_cast<uint8_t>(c),
                val,
                d_hat
            );

            c++;

            if (c == d_hat)
            {
                c = 0;
                blk++;

                if (blk <= m1)
                {
                    ternary_finish_block(out, blk);
                }
            }
        }
    }
}


//==============================================================================
// HCom1   -    H_Com, custom Hash function needed in BLNS for commitment. 
//              It generates the 1st challenge used in the NIZK proof system.
// 
// Input:
// - state0:    initial status structure
// - m1:        m1_Com parameter
//
// Output:
// - R_goth:    matrix of {-1, 0, 1} mod q1_hat values values,
//              equivalent to (R_goth_0 - R_goth_1) in BLNS
//==============================================================================
void HCom1(mat_zz_p& R_goth, const HASH_STATE_t *state0, const ulong &m1)
{
    long         i;
    HASH_STATE_t state;

    Hash_Copy_To(state, *state0);

    const uint8_t v[1] = {1};
    Hash_Update(&state, v, 1);
    
    // Create the R_goth matrix  
    R_goth.SetDims(256, m1*d_hat);   
    
    //Random generation of R_goth ∈ {-1, 0, 1}^(256 x m_1*d_hat) mod q1_hat
    for(i=0; i<256; i++)
    { 
        Hash_R_goth(R_goth[i], &state, m1*d_hat);
    }

    // return R_goth;
}


//==============================================================================
// HCom2   -    H_Com, custom Hash function needed in BLNS for commitment. 
//              It generates the 2nd challenge used in the NIZK proof system.
// 
// Input:
// - state0:    initial status structure
//
// Output:
// - gamma:     matrix of integers modulo q1_hat
//==============================================================================
void HCom2(mat_zz_p& gamma, const HASH_STATE_t *state0)
{
    // NOTE: assuming that current modulus is q1_hat (not q0)
    long         i, n257;
    HASH_STATE_t state; 
    
    Hash_Copy_To(state, *state0);

    const uint8_t v[1] = {2};
    Hash_Update(&state, v, 1);

    // Compute the minimum number of bytes to represent each coefficient
    const size_t b_coeffs = ceil(log2( conv<double>(q1_hat-1) ) / 8.0);    
    
    n257 = 256 + d0 + 1; 
    // NOTE: gamma has 256+d+1 columns in Com, while 256+d+3 in ISIS   

    // Random generation of gamma ∈ Z^(tau_Com x 256+d0+1)_q1_hat
    gamma.SetDims(tau_Com, n257);

    vector<unsigned char> buf(b_coeffs);

    for(i=0; i<tau_Com; i++)
    {
        Hash_v_zz_p_buf(gamma[i], &state, n257, b_coeffs, buf.data());
    }

    // return gamma;
}


//==============================================================================
// HCom3   -    H_Com, custom Hash function needed in BLNS for commitment. 
//              It generates the 3rd challenge used in the NIZK proof system.
// 
// Input:
// - state0:    initial status structure
//
// Output:
// - mu:        vector with tau_Com polynomials with d_hat coefficients modulo q1_hat
//==============================================================================
void HCom3(vec_zz_pX& mu, const HASH_STATE_t *state0)
{
    // NOTE: assuming that current modulus is q1_hat (not q0)
    long         i;
    HASH_STATE_t state;

    // Compute the minimum number of bytes to represent each coefficient
    const size_t b_coeffs = ceil(log2( conv<double>(q1_hat-1) ) / 8.0);   

    Hash_Copy_To(state, *state0);

    const uint8_t v[1] = {3};
    Hash_Update(&state, v, 1);

    // Random generation of mu ∈ R^(tau_Com)_q1_hat
    mu.SetLength(tau_Com);

    vector<unsigned char> buf(b_coeffs);

    for (i = 0; i < tau_Com; i++)
    {
        Hash_zz_pX_buf(mu[i], &state, d_hat, b_coeffs, buf.data());
    }
        
    // return mu;
}


//==============================================================================
// HCom4   -    H_Com, custom Hash function needed in BLNS for commitment. 
//              It generates the 4th challenge used in the NIZK proof system.
// 
// Input:
// - state0:    initial status structure
//
// Output:
// - c:         polynomial with d_hat coefficients, c ∈ C ⊂ R^_(q1_hat)
//==============================================================================
void HCom4(zz_pX& c, const HASH_STATE_t *state0)
{
    long         i;
    HASH_STATE_t state;    
    ZZ           norm1_c, c_i;    
    ZZX          c0, c_2k;
  
        
    // Compute the minimum number of bytes to represent each coefficient
    const size_t b_coeffs = ceil(log2(xi0) / 8.0);
    
    // Compute (nu0)^(2*k0)
    const ZZ    nu0_2k = power(conv<ZZ>(nu0), 2*k0);
      
    // Initialize the variable norm1_c = ||c^(2k)||_1
    norm1_c = 2*nu0_2k;
    
    Hash_Copy_To(state, *state0);

    const uint8_t v[1] = {4};
    Hash_Update(&state, v, 1);

    c0.SetLength(d_hat);

    vector<unsigned char> buf(b_coeffs);
    // Loop to ensure that (2k)√(||c^(2k)||_1 ≤ nu0,  
    // i.e.  ||c^(2k)||_1 ≤ (nu0)^(2k)
    while(norm1_c > nu0_2k)
    {
        // Random generation of c ∈ R^_(xi0+1)
        Hash_ZZ_xi0_buf(c_i, &state, b_coeffs, buf.data());
        // NOTE: generate each coefficient c[i] ∈ [0, xi0], to ensure ||c||∞ ≤ ξ
        
        // c[0] = c_i;
        SetCoeff(c0, 0, c_i);       
                
        for(i=1; i<(d_hat/2); i++)
        {
            Hash_ZZ_xi0(c_i, &state, b_coeffs);
            
            // c[i] = c_i;
            SetCoeff(c0, i, c_i);

            // c[d_hat-i] = -c[i];
            SetCoeff(c0, (d_hat-i), -c_i);
            // NOTE: this ensures that σ(c) = c
        }
        c0.normalize();

        c = conv<zz_pX>(c0);
        
        // NOTE: avoid (rare) cases with c == 0
        if (IsZero(c))
        {
            continue;
        }
        
        // c_2k = power(c, (2*k0));
        PowerOfTwo_ModPhi_hat(c_2k, c0, 2*k0);

        // Compute ||c^(2k)||_1
        norm1_c = 0;

        for(i=0; i<=deg(c_2k); i++)
        {
            // norm1_c = norm1_c + c_2k[i];
            norm1_c += abs(coeff(c_2k, i)); 
        }
    }
         
    // return c;
}


//==============================================================================
// HISIS1   -   H_ISIS, custom Hash function needed in BLNS for ISIS. 
//              It generates the 1st challenge used in the NIZK proof system.
// 
// Input:
// - state0:    initial status structure
// - m1:        m1_ISIS parameter
//
// Output:
// - R_goth:    matrix of {-1, 0, 1} mod q2_hat values, 
//              equivalent to (R_goth_0 - R_goth_1) in BLNS
//==============================================================================
// NOTE: HISIS1 is identical to HCom1, apart m1
void HISIS1(mat_zz_p& R_goth, const HASH_STATE_t *state0, const ulong &m1)
{
    long         i;
    HASH_STATE_t state;

    Hash_Copy_To(state, *state0);

    const uint8_t v[1] = {1};
    Hash_Update(&state, v, 1);
    
    // Create the R_goth matrix  
    R_goth.SetDims(256, m1*d_hat);   
    
    // Random generation of R_goth ∈ {-1, 0, 1}^(256 x m_1*d_hat) mod q2_hat
    for(i=0; i<256; i++)
    { 
        Hash_R_goth(R_goth[i], &state, m1*d_hat);
    }

    // return R_goth;
}

//==============================================================================
// HISIS1   -   H_ISIS, custom Hash function needed in BLNS for ISIS. 
//              It generates the 1st challenge used in the NIZK proof system.
//
// Build the 256 sparse packed rows of R_goth used by the ISIS hash layer.
//
// Inputs:
// - state0 : initial hash state
// - m1     : number of polynomial blocks in each row of R_goth
//
// Output:
// - R_goth : vector of 256 sparse packed rows
//==============================================================================

void HISIS1_optimized(std::vector<TernaryCoeffStructure>& R_goth,
                  const HASH_STATE_t* state0,
                  const ulong& m1)
{
    HASH_STATE_t state;
    Hash_Copy_To(state, *state0);

    // Domain-separation byte for the HISIS1 expansion.
    const uint8_t v[1] = {1};
    Hash_Update(&state, v, 1);

    // Total number of coefficients in one row of R_goth.
    const ulong n_elems = m1 * d_hat;

    // Each byte encodes 4 coefficients, since each coefficient uses 2 bits.
    const ulong n_bytes = (n_elems + 3) >> 2;

    // Temporary buffer used to receive squeezed bytes from SHAKE128.
    std::vector<unsigned char> buf(n_bytes);

    // Ensure the output vector contains exactly 256 rows.
    //R_goth.resize(256);

    // Build each row independently.
    for (ulong r = 0; r < 256; ++r)
    {
        Hash_Sigma_R_goth_encoded(R_goth[r], &state, m1, buf.data());
    }
}




//==============================================================================
// HISIS2   -   H_ISIS, custom Hash function needed in BLNS for ISIS. 
//              It generates the 2nd challenge used in the NIZK proof system.
// 
// Input:
// - state0:    initial status structure
//
// Output:
// - gamma:     matrix of integers modulo q2_hat
//==============================================================================
void HISIS2(mat_zz_p& gamma, const HASH_STATE_t *state0)
{
    // NOTE: assuming that current modulus is q2_hat (not q0)
    long         i, n259;       
    HASH_STATE_t state; 
    
    Hash_Copy_To(state, *state0);

    const uint8_t v[1] = {2};
    Hash_Update(&state, v, 1);

    // Compute the minimum number of bytes to represent each coefficient
    const size_t b_coeffs = ceil(log2( conv<double>(q2_hat-1) ) / 8.0);    
    
    n259 = 256 + d0 + 3;
    // NOTE: gamma has 256+d+3 columns in ISIS, while 256+d+1 in Com 

    // Random generation of gamma ∈ R^(tau_ISIS x 256+d+3)_q2_hat
    gamma.SetDims(tau_ISIS, n259);

    vector<unsigned char> buf(b_coeffs);

    for (i = 0; i < tau_ISIS; i++)
    {
        Hash_v_zz_p_buf(gamma[i], &state, n259, b_coeffs, buf.data());
    }
    // return gamma;
}


//==============================================================================
// HISIS3   -   H_ISIS, custom Hash function needed in BLNS for ISIS. 
//              It generates the 3rd challenge used in the NIZK proof system.
// 
// Input:
// - state0:    initial status structure
//
// Output:
// - mu:        vector with tau_ISIS polynomials with d_hat coefficients modulo q2_hat
//==============================================================================
void HISIS3(vec_zz_pX& mu, const HASH_STATE_t *state0)
// NOTE: HISIS3 is identical to HCom3, apart the modulo  
{     
    // NOTE: assuming that current modulus is q2_hat (not q0)
    long         i;
    HASH_STATE_t state;

    // Compute the minimum number of bytes to represent each coefficient
    const size_t b_coeffs = ceil(log2( conv<double>(q2_hat-1) ) / 8.0);   

    Hash_Copy_To(state, *state0);

    const uint8_t v[1] = {3};
    Hash_Update(&state, v, 1);

    // Random generation of mu ∈ R^(tau_ISIS)_q2_hat
    mu.SetLength(tau_ISIS);

    vector<unsigned char> buf(b_coeffs);

    for (i = 0; i < tau_ISIS; i++)
    {
        Hash_zz_pX_buf(mu[i], &state, d_hat, b_coeffs, buf.data());
    }
        
    // return mu;
}


//==============================================================================
// HISIS4   -   H_ISIS, custom Hash function needed in BLNS for ISIS. 
//              It generates the 4th challenge used in the NIZK proof system.
// 
// Input:
// - state0:    initial status structure
//
// Output:
// - c:         polynomial with d_hat coefficients, c ∈ C ⊂ R^_(q2_hat)
//==============================================================================
void HISIS4(zz_pX& c, const HASH_STATE_t *state0)
// NOTE: HISIS4 is identical to HCom4, apart the modulo
{
    long         i;
    HASH_STATE_t state;    
    ZZ           norm1_c, c_i;
    ZZX          c0, c_2k;

    //c_2k.SetMaxLength(d_hat);
        
    // Compute the minimum number of bytes to represent each coefficient
    const size_t b_coeffs = ceil(log2(xi0) / 8.0);
    
    // Compute (nu0)^(2*k0)
    const ZZ    nu0_2k = power(conv<ZZ>(nu0), 2*k0);
      
    // Initialize the variable norm1_c = ||c^(2k)||_1
    norm1_c = 2*nu0_2k;
    
    Hash_Copy_To(state, *state0);

    const uint8_t v[1] = {4};
    Hash_Update(&state, v, 1);

    c0.SetLength(d_hat);

    vector<unsigned char> buf(b_coeffs);
    // Loop to ensure that (2k)√(||c^(2k)||_1 ≤ nu0,  
    // i.e.  ||c^(2k)||_1 ≤ (nu0)^(2k)
    while(norm1_c > nu0_2k)
    {
        
        // Random generation of c ∈ R^_(xi0+1)
        Hash_ZZ_xi0_buf(c_i, &state, b_coeffs, buf.data());
        // NOTE: generate each coefficient c[i] ∈ [0, xi0], to ensure ||c||∞ ≤ ξ
        
        // c[0] = c_i;
        SetCoeff(c0, 0, c_i);       
                
        for(i=1; i<(d_hat/2); i++)
        {
            Hash_ZZ_xi0(c_i, &state, b_coeffs);
            
            // c[i] = c_i;
            SetCoeff(c0, i, c_i);

            // c[d_hat-i] = -c[i];
            SetCoeff(c0, (d_hat-i), -c_i);
            // NOTE: this ensures that σ(c) = c
        }
        c0.normalize();

        c = conv<zz_pX>(c0);

        if (IsZero(c0))
        {
            continue;
        }
        
        // c_2k = power(c, (2*k0)); 
        PowerOfTwo_ModPhi_hat(c_2k, c0, 2*k0);
        // Compute ||c^(2k)||_1
        norm1_c = 0;

        for(i=0; i<=deg(c_2k); i++)
        {
            // norm1_c = norm1_c + c_2k[i];
            norm1_c += abs(coeff(c_2k, i)); 
        }
    }
         
    // return c;
}


//==============================================================================
// HM     -     H_M, custom Hash function needed in BLNS for hashing attributes. 
//              It hashes an attribute a_i into a vector of length h0, 
//              with coefficients in the range [−ψ, ψ].
// 
// Input:
// - a_i:       attribute, string of bits of arbitrary length a_i ∈ {0, 1}∗
//
// Output:
// - m_i:       vector with h0 coefficients in the range [−psi0, psi0]
//==============================================================================
void HM(vec_ZZ& m_i, const string& a_i)
{
    long         k, range;
    HASH_STATE_t state;
    
    // Compute the numerical range of each coefficient
    range = 2*psi0 + 1;

    zz_pPush push(range);
    // NOTE: backup current modulus q0, temporarily set to range (i.e., zz_p::init(range))
    vec_zz_p    tmp;
    
    // Compute the minimum number of bytes to represent each coefficient
    const size_t b_coeffs = ceil(log2(range-1) / 8.0);

    Hash_Init_To(state, reinterpret_cast<const uint8_t*>(&a_i[0]), a_i.length());

    vector<unsigned char> buf(b_coeffs);

    // Random generation of m_i (modulo range)
    Hash_v_zz_p_buf(tmp, &state, h0, b_coeffs, buf.data());
    m_i = conv<vec_ZZ>( tmp );

    for(k=0; k<h0; k++)
    {
        m_i[k] = m_i[k] - psi0;
        // NOTE: now each coefficient is in the range [−psi0, psi0]
    }

    // return m_i;
}
