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
void LHC_Com(Vec<vec_zz_pX>& com, Vec<vec_ZZX>& st, const long& index, const Mat<zz_pX>& A_i, const Mat<zz_pX>& B_i, const vec_ZZX& s, const vec_ZZX& y) 
{
    // NOTE: assuming that current modulus is q1_hat (not q0)
    long        i, j, m, n, eta;
    RR          alpha_i, s_goth;    
    vec_ZZX     e_1, e_2, e_3, f_1, f_2, f_3; 
    vec_zz_pX   t_1, t_2, w_1, w_2;
    ZZX         acc_1, acc_2;
    
    // Manage the invocation with index 1 or 2
    n   = n_i;
    eta = eta_i; 
 
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
    // \overline{\mathfrak{s}}_1 (or 2)    
    s_goth = alpha_i * RR(eta_i * nu0) * sqrt(RR((n + 2*m) * d0));
               
    // Random generation of e_i        
    e_1.SetLength(n);
    e_2.SetLength(m);
    e_3.SetLength(m);      

    for(i=0; i<n; i++)
    {
        e_1[i].SetLength(d_hat);

        for(j=0; j<d_hat; j++)
        {
            e_1[i][j] = RandomBnd(eta+1);            
        }
    }
     for(i=0; i<m; i++)
    {
        e_2[i].SetLength(d_hat);
        e_3[i].SetLength(d_hat);

        for(j=0; j<d_hat; j++)
        {
            e_2[i][j] = RandomBnd(eta+1);  
            e_3[i][j] = RandomBnd(eta+1);          
        }
    }    

    // Initialization of t_i
    t_1.SetLength(m);
    t_2.SetLength(m);

    for(i=0; i<m; i++)
    {
        t_1[i].SetLength(d_hat);
        t_2[i].SetLength(d_hat);
    }   
    
    // Other variables useful for computations
    acc_1.SetLength(d_hat);
    acc_2.SetLength(d_hat);   
       
    for(i=0; i<m; i++)
    {
        // acc_1 = 0;        
        // acc_2 = 0;          
        clear(acc_1);
        clear(acc_2);
              
        for(j=0; j<n; j++)
        {
            // acc_1 += A_i[i,j] * e_i,1[j]
            acc_1 += ModPhi_hat(conv<ZZX>(A_i[i][j]) * e_1[j]); 
            // acc_2 += B_i[i,j] * e_i,1[j]
            acc_2 += ModPhi_hat(conv<ZZX>(B_i[i][j]) * e_1[j]);
        }      
        
        t_1[i] = conv<zz_pX>( p_bar * (acc_1 + e_2[i]) );
        t_2[i] = conv<zz_pX>( p_bar * (acc_2 + e_3[i]) + s[i] );         
        // NOTE: modulo q_hat on all coefficients (zz_pX)        
    }
   
    // Initialization of f_i
    f_1.SetLength(n);
    f_2.SetLength(m);
    f_3.SetLength(m);  
    
    // Discrete gaussian random generation (using Zsampler) for f_i   
    for(i=0; i<n; i++)
    {
        f_1[i].SetLength(d_hat);

        for(j=0; j<d_hat; j++)
        {
            ZSampler(f_1[i][j], s_goth, RR(0));
        }
    }
    for(i=0; i<m; i++)
    {
        f_2[i].SetLength(d_hat);
        f_3[i].SetLength(d_hat);

        for(j=0; j<d_hat; j++)
        {
            ZSampler(f_2[i][j], s_goth, RR(0));
            ZSampler(f_3[i][j], s_goth, RR(0));
        }
    }

    // Initialization of w_i    
    w_1.SetLength(m);
    w_2.SetLength(m);

    for(i=0; i<m; i++)
    {
        w_1[i].SetLength(d_hat);
        w_2[i].SetLength(d_hat);
    
        // acc_1 = 0;
        // acc_2 = 0; 
        clear(acc_1);
        clear(acc_2);  
              
        for(j=0; j<n; j++)
        {
            // acc_1 += A_i[i,j] * f_i,1[j]
            acc_1 += ModPhi_hat(conv<ZZX>(A_i[i][j]) * f_1[j]); 
            // acc_2 += B_i[i,j] * f_i,1[j]
            acc_2 += ModPhi_hat(conv<ZZX>(B_i[i][j]) * f_1[j]);
        }      
        
        w_1[i] = conv<zz_pX>( p_bar * (acc_1 + f_2[i]) );
        w_2[i] = conv<zz_pX>( p_bar * (acc_2 + f_3[i]) + y[i] );       
        // NOTE: modulo q_hat on all coefficients (zz_pX)
    }
        
    // Store the results in com
    com.SetLength(4);
    com[0] = t_1;
    com[1] = t_2;
    com[2] = w_1;
    com[3] = w_2;   

    // Store the results in st
    st.SetLength(6);
    st[0] = e_1;
    st[1] = e_2;
    st[2] = e_3;
    st[3] = f_1;
    st[4] = f_2;
    st[5] = f_3;
  
    // return(com, st)  
}


//==============================================================================
// Rej  -   Rejection function. It takes as inputs the index (i),  
//          2 vectors of the same length (z, v), and 2 scalars (s, M). 
//          It returns as output “reject” (0) or “accept” (1).
//
// Inputs:
// - index: i (1 or 2)
// - z, v:  vectors of n+2*m polynomials of length d_hat
// - s:     scalar, it is the standard deviation
// - M:     scalar, it is a coefficient
//
// Output:
// - 0 | 1: reject or accept
//==============================================================================
long Rej(const long& index, const vec_ZZX& z, const vec_ZZX& v, const RR& s, const RR& M)
{
    long i, j, m;
    RR u, mul, den, eq;
    ZZ dot_prod, norm2, v_ij;    

    // Manage the invocation with index 1 or 2        
    if (index == 1) 
    {        
        m = m1_Com;
    }
    else if (index == 2)
    {
        m = m2_Com;
    }
    else
    {
        cout << "ERROR! index must be 1 or 2" << endl;
        assert((index == 1) || (index == 2));
        m = 0;
    }     
            
    // u <--[0,1), uniformly distributed
    u = random_RR();

    // <z, v> : Dot product between z and v    
    dot_prod = 0;
    // ||v||^2: Square of Euclidean norm of v 
    norm2 = 0;

    for(i=0; i<(n_i+2*m); i++)
    {
        for(j=0; j<d_hat; j++)
        {
            v_ij = coeff(v[i], j); // v[i][j]
            
            // dot_prod += z[i][j] * v[i][j];
            dot_prod    += coeff(z[i], j) * v_ij;
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
// Rej_v_ZZ  -  Rejection function, modified version for vec_ZZ. 
//              It takes as inputs 2 vectors of the same length (z, v), and 2 scalars (s, M). 
//              It returns as output “reject” (0) or “accept” (1).
//
// Inputs:
// - z, v:  vectors of scalars
// - s:     scalar, it is the standard deviation
// - M:     scalar, it is a coefficient
//
// Output:
// - 0 | 1: reject or accept
//==============================================================================
long Rej_v_ZZ(const vec_ZZ& z, const vec_ZZ& v, const RR& s, const RR& M)
{
    long i, len;
    RR u, mul, den, eq;
    ZZ dot_prod, norm2;    

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
        dot_prod += z[i] * v[i];
        // norm2 += v[i] * v[i]; 
        norm2    += sqr(v[i]);   
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
// Rej_v_ZZX  - Rejection function, modified version for vec_ZZX. 
//              It takes as inputs 2 vectors of the same length (z, v), and 2 scalars (s, M). 
//              It returns as output “reject” (0) or “accept” (1).
//
// Inputs:
// - z, v:  vectors of (m1 or m2) polynomials of length d_hat
// - s:     scalar, it is the standard deviation
// - M:     scalar, it is a coefficient
//
// Output:
// - 0 | 1: reject or accept
//==============================================================================
long Rej_v_ZZX(const vec_ZZX& z, const vec_ZZX& v, const RR& s, const RR& M)
{
    long i, j, len;
    RR u, mul, den, eq;
    ZZ dot_prod, norm2, z_Z, v_Z;    
    
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
            z_Z = coeff(z[i], j); // z[i][j]            
            v_Z = coeff(v[i], j); // v[i][j]            
            
            // dot_prod += z[i][j] * v[i][j];             
            dot_prod    += z_Z     * v_Z;
            // norm2    += v[i][j] * v[i][j];
            norm2       += sqr(v_Z);  
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
void LHC_Open(Vec<vec_ZZX>& op, const long& index, const ZZX& c, const Vec<vec_ZZX>& st)
{
    long i, m, n, b;
    RR alpha_i, s_goth, M_bar;
    vec_ZZX  e_1, e_2, e_3, f_1, f_2, f_3, z_1, z_2, z_3, z, v;
        
    // Manage the invocation with index 1 or 2  
    n = n_i;

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
    // \overline{\mathfrak{s}}_1 (or 2)
    s_goth = alpha_i * RR(eta_i * nu0) * sqrt(RR((n + 2*m) * d0));
    // # \overline{M}_1 (or 2)
    M_bar = exp( sqrt( RR(2*(lambda0 + 1)) / log2e_Const ) * 1/alpha_i + 1/(2*sqr(alpha_i)));

    const long n2m = n + 2*m;
   
    // Initialize e_i, f_i, z_i, for i = 1, 2, 3
    e_1.SetLength(n);
    e_2.SetLength(m);
    e_3.SetLength(m);   
    f_1.SetLength(n);
    f_2.SetLength(m);
    f_3.SetLength(m);    
    z_1.SetLength(n);
    z_2.SetLength(m);
    z_3.SetLength(m);  
    for(i=0; i<n; i++)
    {
        e_1[i].SetLength(d_hat);
        f_1[i].SetLength(d_hat);
        z_1[i].SetLength(d_hat);
    }
     for(i=0; i<m; i++)
    {
        e_2[i].SetLength(d_hat);
        e_3[i].SetLength(d_hat);
        f_2[i].SetLength(d_hat);
        f_3[i].SetLength(d_hat);
        z_2[i].SetLength(d_hat);
        z_3[i].SetLength(d_hat);
    }

    // Retrieve e_i and f_i from st
    e_1 = st[0];
    e_2 = st[1];
    e_3 = st[2];
    f_1 = st[3];
    f_2 = st[4];
    f_3 = st[5];
    
    // Compute z_i
    for(i=0; i<n; i++)
    {
        z_1[i] = f_1[i] + ModPhi_hat(c * e_1[i]);
    }
     for(i=0; i<m; i++)
    {
        z_2[i] = f_2[i] + ModPhi_hat(c * e_2[i]);
        z_3[i] = f_3[i] + ModPhi_hat(c * e_3[i]);
    }
    
    // Initialize v, z, to be passed to Rej  
    z.SetLength(n2m);  
    v.SetLength(n2m);     
    
    // Construction of z=[z_1||z_2||z_3] and v=[v_1||v_2||v_3]
    //                     n    m    m           n    m    m
    for(i=0; i<(n2m); i++)
    {
        z[i].SetLength(d_hat);
        v[i].SetLength(d_hat);
        
        if (i < n)  // filling first n positions
        {
            z[i] = z_1[i];
            v[i] = ModPhi_hat(c * e_1[i]);
        }
        else if (i < (n + m)) // filling m positions starting from the n-th 
        {
            z[i] = z_2[i-n];
            v[i] = ModPhi_hat(c * e_2[i-n]);
        }
        else // filling last m positions
        {
        z[i] = z_3[i-n-m];
        v[i] = ModPhi_hat(c * e_3[i-n-m]);
        }        
    }
    
    // Call Rej function to accept or reject
    b = Rej(index, z, v, s_goth, M_bar);
        
    if (b == 0)
    {
        //     op = [] (i.e. ⊥, reject)
        op.kill();
    }
    else //if (b == 1) // (accept)
    {        
        op.SetLength(3);
        op[0] = z_1;
        op[1] = z_2;
        op[2] = z_3;
    }
    // return(op);
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
long LHC_Verify(const long& index, const Mat<zz_pX>& A_i, const Mat<zz_pX>& B_i, const Vec<vec_zz_pX>& com, const ZZX& c, const vec_ZZX& z, const Vec<vec_ZZX>& op)
{
    // NOTE: assuming that current modulus is q1_hat (not q0)
    long        i, j, m, n, flag;
    RR          alpha_i, s_goth, thres;
    vec_zz_pX   t_1, t_2, w_1, w_2, z_a, z_b, zi_mod; 
    vec_ZZX     z_1, z_2, z_3;
    ZZ          norm2_z1, norm2_z2, norm2_z3;
    ZZX         acc_1, acc_2;
    zz_pX       c_mod;
    
    if (op.length() == 0)
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
    // \overline{\mathfrak{s}}_1 (or 2)    
    s_goth = ((alpha_i * eta_i * nu0) * sqrt(conv<RR>((n + 2*m) * d0)));
       
    // Initialize t_i, w_i, z_i
    t_1.SetLength(m);
    t_2.SetLength(m);
    w_1.SetLength(m);
    w_2.SetLength(m);    
    z_1.SetLength(n);
    z_2.SetLength(m);
    z_3.SetLength(m);  

    for(i=0; i<n; i++)
    {
        z_1[i].SetLength(d_hat);
    }
    for(i=0; i<m; i++)
    {
        t_1[i].SetLength(d_hat);
        t_2[i].SetLength(d_hat);
        w_1[i].SetLength(d_hat);
        w_2[i].SetLength(d_hat);
        z_2[i].SetLength(d_hat);
        z_3[i].SetLength(d_hat);
    }   

    // Retrieve t_i, w_i, z_i from com and op   
    t_1 = com[0];
    t_2 = com[1];
    w_1 = com[2];
    w_2 = com[3];
    z_1 = op[0];
    z_2 = op[1];
    z_3 = op[2];
   
    // Compute ||z_i||^2: Square of Euclidean norm of each z_i
    norm2_z1 = Norm2X(z_1, d_hat);
    norm2_z2 = Norm2X(z_2, d_hat);
    norm2_z3 = Norm2X(z_3, d_hat);
       
    // Check first conditions to reject (0)
    thres = s_goth * sqrt(conv<RR>(2 * n * d0));

    if ( sqrt(conv<RR>(norm2_z1)) > thres)
    { 
        cout << "\n Invalid z_1 norm!" << endl; 
        return 0;
    }

    thres = s_goth * sqrt(conv<RR>(2 * m * d_hat));

    if ( sqrt(conv<RR>(norm2_z2)) > thres)
    { 
        cout << "\n Invalid z_2 norm!" << endl; 
        return 0;
    }

    if ( sqrt(conv<RR>(norm2_z3)) > thres)
    { 
        cout << "\n Invalid z_3 norm!" << endl; 
        return 0;
    }

    // Initialize z_a, z_b
    z_a.SetLength(m);
    z_b.SetLength(m);
    for(i=0; i<m; i++)
    {
        z_a[i].SetLength(d_hat);
        z_b[i].SetLength(d_hat);
    }

    // Other variables useful for computations
    acc_1.SetLength(d_hat);
    acc_2.SetLength(d_hat); 
    c_mod = conv<zz_pX>(c);
       
    for(i=0; i<m; i++)
    {
        // acc_1 = 0;
        // acc_2 = 0; 
        clear(acc_1);
        clear(acc_2);
              
        for(j=0; j<n; j++)
        {
            // acc_1 += A_i[i,j] * z_i,1[j]
            acc_1 += ModPhi_hat( conv<ZZX>(A_i[i][j]) * z_1[j] ); 
            // acc_2 += B_i[i,j] * z_i,1[j]
            acc_2 += ModPhi_hat( conv<ZZX>(B_i[i][j]) * z_1[j] );
        }      
        
        z_a[i] = ModPhi_hat_q( c_mod * t_1[i] ) + w_1[i] - conv<zz_pX>(p_bar * (acc_1 + z_2[i]));
        z_b[i] = ModPhi_hat_q( c_mod * t_2[i] ) + w_2[i] - conv<zz_pX>(p_bar * (acc_2 + z_3[i]));     
        // NOTE: modulo q_hat on all coefficients (zz_pX)
    }
    
    // Check final conditions to reject (0) or accept (1)
    flag = 0;
    zi_mod.SetLength(m);

    for(i=0; i<m; i++)
    {
        //  Check if (z_a[i] != 0)
        z_a[i].normalize();
        flag += (z_a[i] != 0);
        
        //  Check if (z_b[i] != z[i])  
        zi_mod[i].SetLength(d_hat);      
        zi_mod[i] = conv<zz_pX>(z[i]);        
        flag += (z_b[i] != zi_mod[i]);
    }    

    // Return reject (0) or accept (1)
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
