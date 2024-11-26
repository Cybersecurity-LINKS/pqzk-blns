#include "Utils.h"


ZZX Phi()
{
    ZZX phi;

    phi.SetLength(d0+1);
    // phi[0]  = 1;
    SetCoeff(phi, 0, 1);
    // phi[d0] = 1;
    SetCoeff(phi, d0, 1);

    return phi;
}

ZZX Phi_hat()
{
    ZZX phi_hat;

    phi_hat.SetLength(d_hat+1);
    // phi_hat[0]     = 1;
    SetCoeff(phi_hat, 0, 1);
    // phi_hat[d_hat] = 1;
    SetCoeff(phi_hat, d_hat, 1);

    return phi_hat;
}


//==============================================================================
// GS_Ortho  - (Classical) GramSchmidt.Orthogonalization function.
//             For an input matrix              B      ∈ Z^(2d×2d) 
//             it returns the matrix            Bt     ∈ R^(2d×2d)
//             and the vector of squared norms  Norms2 ∈ R^(2d)
//==============================================================================
// NOTE: Classical Gram-Schmidt algorithm is numerically UNSTABLE, 
//       see https://en.wikipedia.org/wiki/Gram%E2%80%93Schmidt_process#Numerical_stability
void GS_Ortho(mat_RR& Bt, vec_RR& Norms2, const mat_L& B)
{
    long    j, s;
    vec_RR  Sum;

    Bt.SetDims(2*d0, 2*d0);
    Sum.SetLength(2*d0);
    Norms2.SetLength(2*d0);
        
    // 1. b˜_1 ← b_1
    Bt[0] = conv<vec_RR>(B[0]);

    Norms2[0] = Norm2R(Bt[0]);

    // 2. for j = 2 : n
    for(j=1; j<(2*d0); j++)
    {
        // 3. b˜_j ← b_j − Sum_{s=1}^{j−1} ⟨b_j, b˜_s / ∥b˜_s∥^2 ⟩ * b˜_s
        Bt[j] = conv<vec_RR>(B[j]);
        
        clear(Sum); // = 0

        for(s=0; s<j; s++)
        {
            Sum +=  ((Bt[j] * Bt[s]) / Norms2[s]) * Bt[s];
        }

        Bt[j] -= Sum;

        Norms2[j] = Norm2R(Bt[j]);
    }

    // 4. return Bt ← (b˜_1, ..., b˜_n)
}


//==============================================================================
// MGS_Ortho - Modified Gram-Schmidt Orthogonalization function.
//             For an input matrix              B      ∈ Z^(2d×2d) 
//             it returns the matrix            Bt     ∈ R^(2d×2d)
//             and the vector of squared norms  Norms2 ∈ R^(2d)
//==============================================================================
// NOTE: implemented as in https://ocw.mit.edu/courses/18-335j-introduction-to-numerical-methods-spring-2019/be0cdadd9de56ff8d20d9a0c6d6d9206_MIT18_335JS19_lec9_reading.pdf
void MGS_Ortho(mat_RR& Bt, vec_RR& Norms2, const mat_L& B)
{
    long    i, j;    
    RR      rii, rij;
    vec_RR  qi;

    Bt.SetDims(2*d0, 2*d0);
    qi.SetLength(2*d0);
    Norms2.SetLength(2*d0);

    for(i=0; i<(2*d0); i++)
    {
        Bt[i] = conv<vec_RR>(B[i]);
    }
        
    for(i=0; i<(2*d0); i++)
    {
        Norms2[i] = Norm2R(Bt[i]);
        rii = sqrt(Norms2[i]);

        for(j=0; j<(2*d0); j++)
        {
            qi[j] = Bt[i][j] / rii;
        }

        for(j=(i+1); j<(2*d0); j++)
        {   
            rij    = qi * Bt[j];
            Bt[j] -= rij * qi;
        }
    }

    // return Bt ← (b˜_1, ..., b˜_n)
}


//==============================================================================
// rot(f) - Returns anticircular matrix associated to polynomial f and integer d
// 
// NOTE: A_N(f) matrix as in Definition 1, page 28 of [DLP], 
//       it corresponds to transpose(rot(f)) with respect to [BLNS]
//==============================================================================
void rot(mat_L& M, const ZZX& f)
{
    unsigned int    i, j, dfu;
    int             df;
    
    M.SetDims(d0, d0);
    df = deg(f);
       
    if(df!=-1)
    {
        dfu = ((unsigned) df);

        if(dfu>=d0)
        {
            assert(dfu<d0);
        }
            
        for(i=0; i<d0; i++)
        {
            for(j=0; j<i; j++)     
            {
                // M[i][j] = -f[d0-i+j];
                M[i][j] = conv<long>( -coeff(f, d0-i+j) );            
            }
            for(j=i; j<d0; j++)    
            {
                // M[i][j] = f[j-i];
                M[i][j] = conv<long>( coeff(f, j-i) );
            }
        }  
    }

    // return M; 
}


//==============================================================================
// rot_T(f) - rot function for a polynomial f modulo q
// 
// NOTE: multiplication matrix rot(f) as in page 8 of [BLNS], 
//       it corresponds to the transpose of the A_N(f) matrix in [DLP]
//==============================================================================
mat_ZZ_p rot_T(const ZZ_pX& f)
{
    unsigned int    i, j, dfu;
    int             df;
    mat_ZZ_p        M;  

    M.SetDims(d0, d0);
    df = deg(f);

    if(df==-1)
    {
        return M;
    }
    
    dfu = ((unsigned) df);

    if(dfu>=d0)
    {
        assert(dfu<d0);
    }
        
    for(i=0; i<d0; i++)
    {
        for(j=0; j<=i; j++)     
        {
            // M[i][j] = f[i-j];
            M[i][j] = coeff(f, i-j);
        }
        for(j=i+1; j<d0; j++)    
        {
            // M[i][j] = -f[d0-j+i];
            M[i][j] = -coeff(f, d0-j+i);
        }
    }    

    return M;
}


//==============================================================================
// rot_vect(v) - same as rot_T function, but suitable to a polynomial vector v
// 
// NOTE: it applies the transpose(v) operation by using rot_T() 
//==============================================================================
mat_ZZ_p  rot_vect( const vec_ZZ_pX& v )
{
    unsigned int    i, j, k, r, len;
    mat_ZZ_p        M, R;
    
    len = v.length();
    M.SetDims(d0, d0);
    R.SetDims(d0, len*d0);
    r = 0;

    for(i=0; i<len; i++)
    {
        M = rot_T( v[i] );

        for(j=0; j<d0; j++)
        {
            for(k=0; k<d0; k++)
            {
                R[k][r] = M[k][j];
                // NOTE: subsequent M = rot(v[i]) appended in R as a row (not column!) of matrixes
            }
            r++;
        }
    }

    return R;
}


//==============================================================================
// Coeffs(x) - For an input polynomial vector x ∈ R^l_(q), 
//             it returns the coefficient vector of x, Coeffs(x) ∈ Z^(l*d)
//==============================================================================
vec_ZZ Coeffs(const vec_ZZ_pX  x, const unsigned int l)
{
    unsigned int    i, j, ld;
    vec_ZZ          coeffs_x;
    
    ld = l * d0;   
    coeffs_x.SetLength(ld);
   
    for(i=0; i<l; i++)
    {
        for(j=0; j<d0; j++)     
        {
            // coeffs_x[d0*i + j] = conv<ZZ>( x[i][j] );
            coeffs_x[d0*i + j] = conv<ZZ>( coeff(x[i], j) );
        }        
    }      

    return coeffs_x;
}


//==============================================================================
// CoeffsX(x) - For an input polynomial vector x ∈ R^l, 
//              it returns the coefficient vector of x, Coeffs(x) ∈ Z^(l*d)
//==============================================================================
vec_ZZ CoeffsX(const vec_ZZX  x, const unsigned int l)
{
    unsigned int    i, j, ld;
    vec_ZZ          coeffs_x;
    
    ld = l * d0;   
    coeffs_x.SetLength(ld);
   
    for(i=0; i<l; i++)
    {
        for(j=0; j<d0; j++)     
        {
            // coeffs_x[d0*i + j] = x[i][j];
            coeffs_x[d0*i + j] = coeff(x[i], j);
        }        
    }      

    return coeffs_x;
}


//==============================================================================
// CoeffsInv(c) - For an input vector of coefficients c ∈ Z^(l*d),
//                it returns the polynomial vector x = Coeffs^{−1}(c) ∈ R^l_q
//==============================================================================
vec_ZZ_pX  CoeffsInv(const vec_ZZ c, const unsigned int l)
{
    unsigned int    i, j;
    vec_ZZ_pX       x;
    
    x.SetLength(l);
   
    for(i=0; i<l; i++)
    {
        x[i].SetLength(d0);

        for(j=0; j<d0; j++)     
        {
            // x[i][j] = conv<ZZ_p>( c[d0*i + j] );
            SetCoeff(x[i], j, conv<ZZ_p>( c[d0*i + j] ) );
        }        
    }      

    return x;
}


//==============================================================================
// CoeffsInvX(c) - For an input vector of coefficients c ∈ Z^(l*d),
//                 it returns the polynomial vector x = Coeffs^{−1}(c) ∈ R^l
//==============================================================================
vec_ZZX  CoeffsInvX(const vec_ZZ c, const unsigned int l)
{
    unsigned int    i, j;
    vec_ZZX         x;
    
    x.SetLength(l);
   
    for(i=0; i<l; i++)
    {
        x[i].SetLength(d0);

        for(j=0; j<d0; j++)     
        {
            // x[i][j] = c[d0*i + j];
            SetCoeff( x[i], j, c[d0*i + j] );
        }        
    }      

    return x;
}


//==============================================================================
// CoeffsHat(x) - For an input polynomial vector x ∈ R_hat^l, 
//                it returns the coefficient vector of x, Coeffs(x) ∈ Z^(l*d_hat)
//==============================================================================
vec_ZZ CoeffsHat(const vec_ZZX  x, const unsigned int l)
{
    unsigned int    i, j, ld;
    vec_ZZ          coeffs_x;
    
    ld = l * d_hat;   
    coeffs_x.SetLength(ld);    
   
    for(i=0; i<l; i++)
    {
        for(j=0; j<d_hat; j++)     
        {
            // coeffs_x[d_hat*i + j] = x[i][j];
            coeffs_x[d_hat*i + j] = coeff(x[i], j);
        }        
    }      

    return coeffs_x;
}


//==============================================================================
// CoeffsInvHat(c) - For an input vector of coefficients c ∈ Z^(l*d_hat)_(q_hat),
//                   it returns the polynomial vector x = Coeffs^{−1}(c) ∈ R_hat^l_(q_hat)
//==============================================================================
vec_ZZ_pX  CoeffsInvHat(const vec_ZZ_p c, const unsigned int l)
{
    unsigned int    i, j;
    vec_ZZ_pX       x;
    
    x.SetLength(l);
   
    for(i=0; i<l; i++)
    {
        x[i].SetLength(d_hat);

        for(j=0; j<d_hat; j++)     
        {
            // x[i][j] = c[d_hat*i + j];
            SetCoeff( x[i], j, c[d_hat*i + j] );
        }        
    }      

    return x;
}


//==============================================================================
// CoeffsInvHatX(c) - For an input vector of coefficients c ∈ Z^(l*d_hat),
//                    it returns the polynomial vector x = Coeffs^{−1}(c) ∈ R_hat^l
//==============================================================================
vec_ZZX  CoeffsInvHatX(const vec_ZZ c, const unsigned int l)
{
    unsigned int    i, j;
    vec_ZZX         x;
    
    x.SetLength(l);
   
    for(i=0; i<l; i++)
    {
        x[i].SetLength(d_hat);

        for(j=0; j<d_hat; j++)     
        {
            // x[i][j] = c[d_hat*i + j];
            SetCoeff( x[i], j, c[d_hat*i + j] );
        }        
    }      

    return x;
}


//==============================================================================
// sigma_map(M, d) - This is the sigma automorphism that maps X --> X^(d-1),
//                   for example sigma(2X^2 + 3X + 5) = -2X^{d-2} -3X^{d-1} + 5.   
//                   Note that the result is mod d and mod phi_hat2 = (X^d + 1).
//                   It takes as input a polynomial vector and its degree.
//                   It outputs the result of the automorphism.
//==============================================================================
vec_ZZ_pX  sigma_map(const vec_ZZ_pX& M, const unsigned int d)
{    
    unsigned int i, j, len;  
    vec_ZZ_pX    N;
  
    len = M.length();  
    N.SetLength(len);

    for(i=0; i<len; i++)
    {       
        N[i].SetLength(d);

        // N[i][0] = M[i][0];
        SetCoeff( N[i], 0, coeff(M[i], 0) );

        for(j=1; j<d; j++) // NOTE: j starts from 1 (not 0)
        {            
            // N[i][d - j] = -M[i][j];
            SetCoeff( N[i], (d - j), -coeff(M[i], j) );
        }        
    }

    return N;
}


//=====================================================================================
// poly_mult  -  scalar product between two vectors of polynomials of length d0.
//=====================================================================================
ZZ_pX  poly_mult(const vec_ZZ_pX& f, const vec_ZZ_pX& g)
{
    const ZZ_pX     phi_2 = conv<ZZ_pX>(phi);
    
    unsigned int    i, len;
    ZZ_pX           h;
    
    len = f.length();

    if (len != g.length())
    {
        cout << "ERROR! Two input vectors must have the same dimensions" << endl;
        assert(len == g.length());
    }

    h.SetLength(d0);

    for(i=0; i<len; i++)
    {
        h = h + (f[i] * g[i]) % (phi_2);        
    }

    return h;   
}


//=====================================================================================
// poly_mult_hat  -  scalar product between two vectors of polynomials of length d_hat.
//=====================================================================================
ZZ_pX  poly_mult_hat(const vec_ZZ_pX& f, const vec_ZZ_pX& g)
{
    const ZZ_pX     phi_hat2 = conv<ZZ_pX>(phi_hat);
    
    unsigned int    i, len;
    ZZ_pX           h;
    
    len = f.length();

    if (len != g.length())
    {
        cout << "ERROR! Two input vectors must have the same dimensions" << endl;
        assert(len == g.length());
    }

    h.SetLength(d_hat);

    for(i=0; i<len; i++)
    {
        h = h + (f[i] * g[i]) % (phi_hat2);        
    }

    return h;   
}


//=====================================================================================
// Compute_f -  Compute the function f(x) associated with the ISIS_f problem
//                    f(x) := Coeffs^(−1)(B_f · enc(x)) ∈ R^n_q   
//              where B_f ∈ Z^(nd×t)_q  is a randomly chosen matrix.
//=====================================================================================
ZZ_pX   Compute_f(const mat_ZZ_p& B_f, const ZZ& x)
{    
    const unsigned int n = 1;    
    // NOTE: assuming  n = 1, thus B_f ∈ Z^(d×t)_q  and  f(x) ∈ R_q 

    unsigned int    i;    
    vec_ZZ_p        enc_x;
    vec_ZZ_pX       vec_f;
    ZZ_pX           f_x;
    
    // Compute enc(x) ∈ {0, 1}^t, the binary decomposition of (x−1)   
    enc_x.SetLength(t0);

    for(i=0; i<t0; i++)
    {
        enc_x[i] = bit(x-1, i);
    }

    // Compute f(x) := Coeffs^(−1)(B_f · enc(x)) 
    vec_f = CoeffsInv(conv<vec_ZZ>(B_f*enc_x), n);

    f_x = vec_f[0];

    return f_x;
}


//==============================================================================
// Computes the squared norm of a vector v of integers, without modulo.
//==============================================================================
ZZ  Norm2(const vec_ZZ& v)
{
    long    i;
    ZZ      norm2;

    norm2 = 0;

    for(i=0; i<(v.length()); i++)
    {
        // norm2 = norm2 + v[i]*v[i];
        norm2 += sqr( v[i] );
    }   

    return norm2;
}


//==============================================================================
// Computes the squared norm of a vector v of polynomials with d coefficients.
//==============================================================================
ZZ  Norm2X(const vec_ZZX& v, const unsigned int d)
{
    long    i, j;
    ZZ      norm2;

    norm2 = 0;

    for(i=0; i<(v.length()); i++)
    {
        for(j=0; j<d; j++)
        { 
            // norm2 = norm2 + v[i][j] * v[i][j];
            norm2 += sqr( coeff(v[i], j) );
        }
    }  

    return norm2;
}


//=================================================================================
// Computes the squared norm of a vector v of arbitrary-precision floating points.
//=================================================================================
RR  Norm2R(const vec_RR& v)
{
    long    i;
    RR      norm2;

    norm2 = 0;

    for(i=0; i<(v.length()); i++)
    {
        // norm2 = norm2 + v[i]*v[i];
        norm2 += sqr( v[i] );
    }   

    return norm2;
}
