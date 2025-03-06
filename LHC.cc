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

#include "LHC.h"


//==============================================================================
// LHC_Com  -   Computes the i-th LHC Commitment. 
//              This function takes as inputs index, A_i, B_i, s and y.
//              It returns as outputs com and st.
// 
// Inputs:
// - index:     i (1 or 2)
// - A_i, B_i:  uniformly random matrices (i.e., structure crs_LHC) 
//              A_i, B_i ∈ R_hat^(m_i x n_i)_(q_hat)
// - s, y:      vectors of polynomials, s, y ∈ R_hat^m_i  
//  
// Outputs:
// - com:       commitment structure    
// - st:        status structure
//==============================================================================
void LHC_Com(LHC_COM_t& com, LHC_ST_t& st, const long& index, const Mat<zz_pX>& A_i, const Mat<zz_pX>& B_i, const vec_zz_pX& s, const vec_zz_pX& y) 
{
    // NOTE: assuming that current modulus is q1_hat (not q0)
    long        i, m;
    double      alpha_i;
    
    // Manage the invocation with index 1 or 2
    const long n   = n_i;
    const long eta = eta_i; 
 
    if (index == 1) 
    {
        alpha_i = double(alpha_bar_1);
        m       = m1_Com;
    }
    else if (index == 2)
    {
        alpha_i = double(alpha_bar_2);
        m       = m2_Com;
    }
    else
    {
        cout << "ERROR! index must be 1 or 2" << endl;
        assert((index == 1) || (index == 2));
        m = 0;
    }     
    // \overline{\mathfrak{s}}_1 (or 2)    
    const double s_goth = alpha_i * double(eta_i * nu0) * sqrt(double((n + 2*m) * d0));
    
    
    // 1. Retrieve A_i, B_i from crs
    // NOTE: A_i, B_i directly provided as inputs
        
    // 2. Random generation of e_i
    {
        zz_pPush push(eta+1); 
        // NOTE: backup current modulus q1_hat, temporarily set to eta+1 (i.e., zz_p::init(eta+1))

        st.e_1.SetLength(n);
        st.e_2.SetLength(m);
        st.e_3.SetLength(m);      

        for(i=0; i<n; i++)
        {
            st.e_1[i] = random_zz_pX(d_hat);
        }
        for(i=0; i<m; i++)
        {
            st.e_2[i] = random_zz_pX(d_hat);
            st.e_3[i] = random_zz_pX(d_hat);
        }
    }
    // NOTE: e_1, e_2, e_3 coefficients are in the range [0, eta]

    // 3. / 4. Compute t_i
    com.t_1.SetLength(m);
    com.t_2.SetLength(m);
       
    for(i=0; i<m; i++)
    {
        com.t_1[i] = (poly_mult_hat(A_i[i], st.e_1) + st.e_2[i]) * p_bar;
        com.t_2[i] = (poly_mult_hat(B_i[i], st.e_1) + st.e_3[i]) * p_bar + s[i];
    }
   
    // 5. Initialization of f_i
    st.f_1.SetLength(n);
    st.f_2.SetLength(m);
    st.f_3.SetLength(m);  
    
    // Discrete gaussian random generation (using Zsampler) for f_i   
    for(i=0; i<n; i++)
    {
        polySampler_hat(st.f_1[i], s_goth);
    }
    for(i=0; i<m; i++)
    {
        polySampler_hat(st.f_2[i], s_goth);
        polySampler_hat(st.f_3[i], s_goth);
    }


    // 6. / 7. Compute w_i
    com.w_1.SetLength(m);
    com.w_2.SetLength(m);

    for(i=0; i<m; i++)
    {
        com.w_1[i] = (poly_mult_hat(A_i[i], st.f_1) + st.f_2[i]) * p_bar;
        com.w_2[i] = (poly_mult_hat(B_i[i], st.f_1) + st.f_3[i]) * p_bar + y[i];
    }
        
    // 8. Store the results in com
    // com.t_1 = t_1;
    // com.t_2 = t_2;
    // com.w_1 = w_1;
    // com.w_2 = w_2;   

    // 9. Store the results in st
    // st.e_1 = e_1;
    // st.e_2 = e_2;
    // st.e_3 = e_3;
    // st.f_1 = f_1;
    // st.f_2 = f_2;
    // st.f_3 = f_3;
  
    // 10. return(com, st)  
}


//==============================================================================
// Rej_v_zzp  - Rejection function, version for vec_zz_p. 
//              It takes as inputs 2 vectors of the same length (z, v), 
//              their modulo q and 2 scalars (s, M). 
//              It returns as output “reject” (0) or “accept” (1).
//
// Inputs:
// - z, v:  vectors of scalars
// - q:     modulo for the scalars
// - s:     scalar, it is the standard deviation
// - M:     scalar, it is a coefficient
//
// Output:
// - 0 | 1: reject or accept
//==============================================================================
long Rej_v_zzp(const vec_zz_p& z, const vec_zz_p& v, const long& q, const RR& s, const RR& M)
{
    long    i, len;
    RR      u, mul, den, eq;
    ZZ      dot_prod, norm2, z_i, v_i, thresh;

    thresh = q/2; 
    // NOTE: thresh = floor(q/2);

    len = z.length();
   
    if (len != v.length())
    {
        cout << "ERROR! Two input vectors must have the same dimensions" << endl;
        assert(len == v.length());
    }
            
    // u <--[0,1), uniformly distributed
    u = random_RR();

    // <z, v> : Dot product between z and v    
    dot_prod = 0;
    // ||v||^2: Square of Euclidean norm of v 
    norm2 = 0;

    for(i=0; i<len; i++)
    {
        z_i = conv<ZZ>( z[i] );

        if (z_i > thresh)
        {
            z_i -= q;
        }

        v_i = conv<ZZ>( v[i] );

        if (v_i > thresh)
        {
            v_i -= q;
        }
        
        // dot_prod += z[i] * v[i];
        dot_prod    += z_i  * v_i;
        // norm2    += v[i] * v[i];
        norm2       += sqr(v_i);
    }    

    mul = 1.0 / M;
    den = 2*sqr(s);
    eq = mul * exp(conv<RR>(-2*dot_prod + norm2) / den);

    // Condition for accepting or rejecting
    if (u > eq)
    {
        return(0); // reject
    }
    else
    {
        return(1); // accept
    }   
}


//==============================================================================
// Rej_v_zzpX - Rejection function, version for vec_zz_pX. 
//              It takes as inputs 2 vectors of the same length (z, v), 
//              their modulo q and 2 scalars (s, M). 
//              It returns as output “reject” (0) or “accept” (1).
//
// Inputs:
// - z, v:  vectors of (m1 or m2) polynomials of length d_hat
// - q:     modulo for the coefficients of the polynomials
// - s:     scalar, it is the standard deviation
// - M:     scalar, it is a coefficient
//
// Output:
// - 0 | 1: reject or accept
//==============================================================================
long Rej_v_zzpX(const vec_zz_pX& z, const vec_zz_pX& v, const long& q, const RR& s, const RR& M)
{
    long    i, j, len;
    RR      u, mul, den, eq;
    ZZ      dot_prod, norm2, z_ij, v_ij, thresh;
    
    thresh = q/2; 
    // NOTE: thresh = floor(q/2);

    len = z.length();
   
    if (len != v.length())
    {
        cout << "ERROR! Two input vectors must have the same dimensions" << endl;
        assert(len == v.length());
    }
            
    // u <--[0,1), uniformly distributed
    u = random_RR();
   
    // <z, v> : Dot product between z and v    
    dot_prod = 0;
    // ||v||^2: Square of Euclidean norm of v 
    norm2 = 0;

    for(i=0; i<len; i++)
    {
        for(j=0; j<d_hat; j++)
        {
            z_ij = conv<ZZ>(coeff(z[i], j)); // z[i][j]

            if (z_ij > thresh)
            {
                z_ij -= q;
            }

            v_ij = conv<ZZ>(coeff(v[i], j)); // v[i][j]

            if (v_ij > thresh)
            {
                v_ij -= q;
            }
            
            // dot_prod += z[i][j] * v[i][j];             
            dot_prod    += z_ij    * v_ij;
            // norm2    += v[i][j] * v[i][j];
            norm2       += sqr(v_ij);  
        }
    }
    
    mul = 1.0 / M;
    den = 2*sqr(s);
    eq = mul * exp(conv<RR>(-2*dot_prod + norm2) / den);

    // Condition for accepting or rejecting
    if (u > eq)
    {
        return(0); // reject
    }
    else
    {
        return(1); // accept
    }   
}


//==============================================================================
// LHC_Open - This function takes as inputs index, c and st. 
//            If the check passes, it returns as output op.
// 
// Inputs:
// - index:   i (1 or 2)
// - c:       polynomial, generated using H(4,...) in Prove_Com
// - st:      status structure
// 
// Output:  
// - op:      list of (n + m + m) polynomials of d_hat length, if accept, 
//            otherwise op = [] (i.e. op = ⊥, reject)
//==============================================================================
void LHC_Open(LHC_OP_t& op, const long& index, const zz_pX& c, const LHC_ST_t& st)
{
    long        i, m, b;
    RR          alpha_i;
    vec_zz_pX   v_1, v_2, v_3, z, v;
    
    // Manage the invocation with index 1 or 2  
    const long n   = n_i;    

    if (index == 1) 
    {
        alpha_i = RR(alpha_bar_1);
        m       = m1_Com;
    }
    else if (index == 2)
    {
        alpha_i = RR(alpha_bar_2);
        m       = m2_Com;
    }
    else
    {
        cout << "ERROR! index must be 1 or 2" << endl;
        assert((index == 1) || (index == 2));
        m = 0;
    }

    const long n2m = n + 2*m;
    
    // \overline{\mathfrak{s}}_1 (or 2)
    const RR s_goth = alpha_i * RR(eta_i * nu0) * sqrt(RR(n2m * d0));

    // \overline{M}_1 (or 2)
    const RR M_bar = exp( sqrt( RR(2*(lambda0 + 1)) / log2e_Const ) * 1/alpha_i + 1/(2*sqr(alpha_i)));
   
    // Initialize z_i, v_i for i = 1, 2, 3
    op.z_1.SetLength(n);
    op.z_2.SetLength(m);
    op.z_3.SetLength(m);
    v_1.SetLength(n);
    v_2.SetLength(m);
    v_3.SetLength(m);  
    

    // 1. Retrieve e_i and f_i from st
    // e_1 = st.e_1;
    // e_2 = st.e_2;
    // e_3 = st.e_3;
    // f_1 = st.f_1;
    // f_2 = st.f_2;
    // f_3 = st.f_3;
    
    // 2. / 3. Compute z_i
    for(i=0; i<n; i++)
    {
        v_1[i] = ModPhi_hat_q(c * st.e_1[i]);
        op.z_1[i] = st.f_1[i] + v_1[i];
    }
    for(i=0; i<m; i++)
    {
        v_2[i] = ModPhi_hat_q(c * st.e_2[i]);
        op.z_2[i] = st.f_2[i] + v_2[i];

        v_3[i] = ModPhi_hat_q(c * st.e_3[i]);
        op.z_3[i] = st.f_3[i] + v_3[i];
    }
    
    // Initialize v, z, to be passed to Rej  
    z.SetLength(n2m);  
    v.SetLength(n2m);     
    
    // Construction of z=[z_1||z_2||z_3] and v=[v_1||v_2||v_3]
    //                     n    m    m           n    m    m
    for(i=0; i<(n2m); i++)
    {
        if (i < n)  // filling first n positions
        {
            z[i] = op.z_1[i];
            v[i] = v_1[i];
        }
        else if (i < (n + m)) // filling m positions starting from the n-th 
        {
            z[i] = op.z_2[i-n];
            v[i] = v_2[i-n];
        }
        else // filling last m positions
        {
            z[i] = op.z_3[i-n-m];
            v[i] = v_3[i-n-m];
        }        
    }
    
    // 4. Call Rej function
    b = Rej_v_zzpX(z, v, q1_hat, s_goth, M_bar);
    
    // 5. Reject or accept
    if (b == 0)
    {
        // op = [] (i.e. ⊥, reject)
        op.z_1.kill();
        op.z_2.kill();
        op.z_3.kill();
        // NOTE: additional flag, to identify an invalid op        
        op.valid = 0;
    }
    else // (accept)
    {        
        // NOTE: valid op
        op.valid = 1;
    }
    // 6. return(op);
}


//==============================================================================
// LHC_Verify - This function takes as inputs index, A_i, B_i, com, c, z, op. 
//              It returns as output “accept” (1) or “reject” (0).
// 
// Inputs:
// - index:     i (1 or 2)
// - A_i, B_i:  uniformly random matrices (i.e., structure crs_LHC) 
//              A_i, B_i ∈ R_hat^(m_i x n_i)_(q_hat)
// - com:       commitment structure
// - c:         polynomial, generated using H(4,...) in Prove_Com
// - z:         vector of polynomials, z ∈ R_hat^(m_i)_(q_hat) (z = c*s+y leads to result=1)
// - op:        list of (n + m + m) polynomials of d_hat length
// 
// Output:
// - 0 or 1:    reject or accept 
//==============================================================================
long LHC_Verify(const long& index, const Mat<zz_pX>& A_i, const Mat<zz_pX>& B_i, const LHC_COM_t& com, const zz_pX& c, const vec_zz_pX& z, const LHC_OP_t& op)
{
    // NOTE: assuming that current modulus is q1_hat (not q0)
    long        i, m, n, flag;
    ZZ          alpha_i, s_goth2, thres, norm2_z1, norm2_z2, norm2_z3;
    vec_zz_pX   z_a, z_b;
    
    if (op.valid == 0)
    {
        // NOTE: reject because op = [] (i.e. ⊥)
        cout << "\n Reject because op = []" << endl;
        return 0;
    }

    // Manage the invocation with index 1 or 2  
    n = n_i;
    if (index == 1) 
    {
        alpha_i = alpha_bar_1;
        m       = m1_Com;
    }
    else if (index == 2)
    {
        alpha_i = alpha_bar_2;
        m       = m2_Com;
    }
    else
    {
        cout << "\n ERROR! index must be 1 or 2" << endl;
        assert((index == 1) || (index == 2));
        m = 0;
    }     
    // Compute the square of \overline{\mathfrak{s}}_1 (or 2) 
    s_goth2 = sqr(alpha_i * eta_i * nu0) * ((n + 2*m) * d0);

        
    // 1. Retrieve t_i, w_i from com   
    // t_1 = com.t_1;
    // t_2 = com.t_2;
    // w_1 = com.w_1;
    // w_2 = com.w_2;

    // 2. Retrieve z_i from op 
    // z_1 = op.z_1;
    // z_2 = op.z_2;
    // z_3 = op.z_3;
   
    // Compute ||z_i||^2: Square of Euclidean norm of each z_i
    norm2_z1 = Norm2Xm(op.z_1, d_hat, q1_hat);
    norm2_z2 = Norm2Xm(op.z_2, d_hat, q1_hat);
    norm2_z3 = Norm2Xm(op.z_3, d_hat, q1_hat);
       
    // 3. Check z_1 norm
    thres = s_goth2 * (2 * n * d0);

    if ( norm2_z1 > thres)
    { 
        cout << "\n Invalid z_1 norm!" << endl; 
        return 0;
    }

    // 4. Check z_2 norm
    thres = s_goth2 * (2 * m * d_hat);

    if ( norm2_z2 > thres)
    { 
        cout << "\n Invalid z_2 norm!" << endl; 
        return 0;
    }

    // 5. Check z_3 norm
    if ( norm2_z3 > thres)
    { 
        cout << "\n Invalid z_3 norm!" << endl; 
        return 0;
    }


    // 6. / 7. Compute z_a, z_b
    z_a.SetLength(m);
    z_b.SetLength(m);
           
    for(i=0; i<m; i++)
    {
        z_a[i] = ModPhi_hat_q( c * com.t_1[i] ) + com.w_1[i] - (poly_mult_hat(A_i[i], op.z_1) + op.z_2[i]) * p_bar;
        z_b[i] = ModPhi_hat_q( c * com.t_2[i] ) + com.w_2[i] - (poly_mult_hat(B_i[i], op.z_1) + op.z_3[i]) * p_bar;
    }
    

    // 8. Check final conditions to reject (0) or accept (1)
    flag = 0;

    for(i=0; i<m; i++)
    {
        //  Check if (z_a[i] != 0)
        z_a[i].normalize();
        flag += (z_a[i] != 0);
        
        //  Check if (z_b[i] != z[i])
        flag += (z_b[i] != z[i]);
    }    

    // 9. Return reject (0) or accept (1)
    if (flag != 0) // ((z_a != 0) | (z_b != z))
    {
        cout << "\n REJECT!" << endl;
        return 0;
    }
    else
    {
        return 1;
    }
}