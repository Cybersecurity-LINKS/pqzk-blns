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


//==============================================================================
// ModPhi(_hat)(_q) - Fast algorithms to compute p % phi (or  p % phi_hat), 
//                    without or with modulo q (or q_hat) on all coefficients
//==============================================================================
ZZX ModPhi(const ZZX& p)
{
    return (trunc(p, d0) - RightShift(p, d0));
}

zz_pX ModPhi_q(const zz_pX& p)
{
    return (trunc(p, d0) - RightShift(p, d0));
}

ZZX ModPhi_hat(const ZZX& p)
{
    return (trunc(p, d_hat) - RightShift(p, d_hat));
}

zz_pX ModPhi_hat_q(const zz_pX& p)
{
    return (trunc(p, d_hat) - RightShift(p, d_hat));
}


//==============================================================================
// OGS_Ortho - Optimized Gram-Schmidt Orthogonalization function.
//             For an input matrix              B      ∈ Z^(2d×2d) 
//             it returns the matrix            Bt     ∈ R^(2d×2d)
//             and the vector of squared norms  Norms2 ∈ R^(2d)
//==============================================================================
// NOTE: optimized version with double instead of RR,
//       based on Chapter 4 in https://tprest.github.io/pdf/pub/thesis-thomas-prest.pdf
void OGS_Ortho(mat_D& Bt, vec_D& Norms2, const mat_L& B)
{
    long    i, k;  
    double  Ck, Dk, CDk;
    vec_D   v, v1;

    Bt.SetDims(2*d0, 2*d0);
    v.SetLength(2*d0);
    v1.SetLength(2*d0);
    Norms2.SetLength(2*d0);
    
    Bt[0] = conv<vec_D>(B[0]);
   
    for(i=0; i<d0-1; i++)
    {    
        v[i]    = Bt[0][i+1];
        v[i+d0] = Bt[0][i+1+d0];
    }

    v[d0-1]   = -Bt[0][0];
    v[2*d0-1] = -Bt[0][d0];

    v1 = v;

    Ck = InnerProdD(v1, Bt[0]);
    Dk = Norm2D(v1);
    Norms2[0] = Dk;

    for(k=1; k<d0; k++)
    {
        CDk = Ck/Dk;
        Bt[k][0]  = -Bt[k-1][d0-1]   + CDk*v[d0-1];
        Bt[k][d0] = -Bt[k-1][2*d0-1] + CDk*v[2*d0-1];
        
        for(i=1; i<d0; i++)
        {
            Bt[k][i]    = Bt[k-1][i-1]    - CDk*v[i-1];
            Bt[k][i+d0] = Bt[k-1][i+d0-1] - CDk*v[i+d0-1];
        }

        for(i=0; i<2*d0; i++)
        {
            v[i] -= CDk*Bt[k-1][i];
        }
        
        Dk = Dk - Ck*CDk;
        Ck = InnerProdD(v1, Bt[k]); 
        Norms2[k] = Dk;
    }

    for(i=0; i<d0; i++)
    {    
        Bt[d0][d0+i] =  Bt[d0-1][d0-1-i]*q0/Dk;
        Bt[d0][i]    = -Bt[d0-1][2*d0-1-i]*q0/Dk;
    }

    for(i=0; i<d0-1; i++)
    {
        v[i]    = Bt[d0][i+1];
        v[i+d0] = Bt[d0][i+1+d0];
    }

    v[d0-1]   = -Bt[d0][0];
    v[2*d0-1] = -Bt[d0][d0];
    
    v1 = v;
   
    Ck = InnerProdD(v1, Bt[d0]);
    Dk = Norm2D(Bt[d0]);
    Norms2[d0] = Dk;

    for(k=d0+1; k<2*d0; k++)
    {
        CDk = Ck/Dk;
        Bt[k][0]  = -Bt[k-1][d0-1]   + CDk*v[d0-1];
        Bt[k][d0] = -Bt[k-1][2*d0-1] + CDk*v[2*d0-1];
        
        for(i=1; i<d0; i++)
        {
            Bt[k][i]    = Bt[k-1][i-1]    - CDk*v[i-1];
            Bt[k][i+d0] = Bt[k-1][i+d0-1] - CDk*v[i+d0-1];
        }
        
        for(i=0; i<2*d0; i++)
        {
            v[i] -= CDk*Bt[k-1][i];
        }
        
        Dk = Dk - Ck*CDk;
        Ck = InnerProdD(v1, Bt[k]); 
        Norms2[k] = Dk;  
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
    unsigned long   i, j, dfu;
    long            df;
    
    // M.SetDims(d0, d0);
    // NOTE: the size of the output matrix M must be (d0, d0)
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
void rot_T(mat_zz_p& M, const zz_pX& f)
{
    unsigned long   i, j, dfu;
    long            df;

    // M.SetDims(d0, d0);
    // NOTE: the size of the output matrix M must be (d0, d0)
    df = deg(f);

    if(df==-1)
    {
        M.kill();
        return;
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

    // return M;
}


//==============================================================================
// rot_vect(v) - same as rot_T function, but suitable to a polynomial vector v
// 
// NOTE: it applies the transpose(v) operation by using rot_T() 
//==============================================================================
void rot_vect( mat_zz_p& R, const vec_zz_pX& v )
{
    unsigned long   i, j, k, r, len;
    mat_zz_p        M;
    
    len = v.length();
    M.SetDims(d0, d0);
    // R.SetDims(d0, len*d0);
    // NOTE: the size of the output matrix R must be at least (d0, len*d0)

    r = 0;

    for(i=0; i<len; i++)
    {
        rot_T( M, v[i] );

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

    // return R;
}


//==============================================================================
// CoeffsX(x) - For an input polynomial vector x ∈ R^l, 
//              it returns the coefficient vector of x, Coeffs(x) ∈ Z^(l*d)
//==============================================================================
void CoeffsX(vec_ZZ& coeffs_x, const vec_ZZX& x, const unsigned long& l)
{
    unsigned long i, j;
    
    // unsigned long ld = l * d0;
    // coeffs_x.SetLength(ld);
    // NOTE: the size of the output vector coeffs_x must be at least (l * d0)
   
    for(i=0; i<l; i++)
    {
        for(j=0; j<d0; j++)     
        {
            // coeffs_x[d0*i + j] = x[i][j];
            coeffs_x[d0*i + j] = coeff(x[i], j);
        }        
    }      

    // return coeffs_x;
}


//==============================================================================
// CoeffsInv(c) - For an input vector of coefficients c ∈ Z^(l*d),
//                it returns the polynomial vector x = Coeffs^{−1}(c) ∈ R^l_q
//==============================================================================
void CoeffsInv(vec_zz_pX& x, const vec_ZZ& c, const unsigned long& l)
{
    unsigned long   i, j;
    
    x.SetLength(l);
   
    for(i=0; i<l; i++)
    {
        x[i].SetLength(d0);

        for(j=0; j<d0; j++)     
        {
            // x[i][j] = conv<zz_p>( c[d0*i + j] );
            SetCoeff(x[i], j, conv<zz_p>( c[d0*i + j] ) );
        }        
    }      

    // return x;
}


//==============================================================================
// CoeffsInvX(c) - For an input vector of coefficients c ∈ Z^(l*d),
//                 it returns the polynomial vector x = Coeffs^{−1}(c) ∈ R^l
//==============================================================================
void CoeffsInvX(vec_ZZX& x, const vec_ZZ& c, const unsigned long& l)
{
    unsigned long   i, j;
    
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

    // return x;
}


//==============================================================================
// CoeffsHat(x) - For an input polynomial vector x ∈ R_hat^l, 
//                it returns the coefficient vector of x, Coeffs(x) ∈ Z^(l*d_hat)
//==============================================================================
void CoeffsHat(vec_ZZ& coeffs_x, const vec_ZZX& x, const unsigned long& l)
{
    unsigned long   i, j, ld;
    
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

    // return coeffs_x;
}


//==============================================================================
// CoeffsHat_q(x) - For an input polynomial vector x ∈ R_hat^l_(q_hat), 
//                  it returns the coefficient vector of x, Coeffs(x) ∈ Z^(l*d_hat)_(q_hat)
//==============================================================================
void CoeffsHat_q(vec_zz_p& coeffs_x, const vec_zz_pX& x, const unsigned long& l)
{
    unsigned long   i, j, ld;
    
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

    // return coeffs_x;
}


//==============================================================================
// CoeffsInvHat(c) - For an input vector of coefficients c ∈ Z^(l*d_hat)_(q_hat),
//                   it returns the polynomial vector x = Coeffs^{−1}(c) ∈ R_hat^l_(q_hat)
//==============================================================================
void CoeffsInvHat(vec_zz_pX& x, const vec_zz_p& c, const unsigned long& l)
{
    unsigned long   i, j;
    
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

    // return x;
}


//==============================================================================
// CoeffsInvHatX(c) - For an input vector of coefficients c ∈ Z^(l*d_hat),
//                    it returns the polynomial vector x = Coeffs^{−1}(c) ∈ R_hat^l
//==============================================================================
void CoeffsInvHatX(vec_ZZX& x, const vec_ZZ& c, const unsigned long& l)
{
    unsigned long   i, j;
    
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

    // return x;
}


//==============================================================================
// sigma_map(M, d) - This is the sigma automorphism that maps X --> X^(d-1),
//                   for example sigma(2X^2 + 3X + 5) = -2X^{d-2} -3X^{d-1} + 5.   
//                   Note that the result is mod d and mod phi = (X^d + 1).
//                   It takes as input a polynomial vector and its degree.
//                   It outputs the result of the automorphism.
//==============================================================================
void sigma_map(vec_zz_pX& N, const vec_zz_pX& M, const unsigned long& d)
{    
    unsigned long i, j, len;
  
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

    // return N;
}


//=====================================================================================
// poly_mult  -  scalar product between two vectors of polynomials of length d0.
//=====================================================================================
zz_pX  poly_mult(const vec_zz_pX& f, const vec_zz_pX& g)
{
    long    i, len;
    zz_pX   h;
    
    len = f.length();

    if (len != g.length())
    {
        cout << "ERROR! Two input vectors must have the same dimensions" << endl;
        assert(len == g.length());
    }

    h.SetLength(d0);

    for(i=0; i<len; i++)
    {
        h += ModPhi_q( f[i] * g[i]);
    }

    return h;   
}


//=====================================================================================
// poly_mult_hat  -  scalar product between two vectors of polynomials of length d_hat.
//=====================================================================================
zz_pX  poly_mult_hat(const vec_zz_pX& f, const vec_zz_pX& g)
{
    long    i, len;
    zz_pX   h;
    
    len = f.length();

    if (len != g.length())
    {
        cout << "ERROR! Two input vectors must have the same dimensions" << endl;
        assert(len == g.length());
    }

    h.SetLength(d_hat);

    for(i=0; i<len; i++)
    {
        h += ModPhi_hat_q( f[i] * g[i] );
    }

    return h;   
}


//=====================================================================================
// Compute_f -  Compute the function f(x) associated with the ISIS_f problem
//                    f(x) := Coeffs^(−1)(B_f · enc(x)) ∈ R^n_q   
//              where B_f ∈ Z^(nd×t)_q  is a randomly chosen matrix.
//=====================================================================================
zz_pX   Compute_f(const mat_zz_p& B_f, const ZZ& x)
{    
    const unsigned long n = 1;    
    // NOTE: assuming  n = 1, thus B_f ∈ Z^(d×t)_q  and  f(x) ∈ R_q 

    unsigned long   i;    
    vec_zz_p        enc_x;
    vec_zz_pX       vec_f;
    zz_pX           f_x;
    
    // Compute enc(x) ∈ {0, 1}^t, the binary decomposition of (x−1)   
    enc_x.SetLength(t0);

    for(i=0; i<t0; i++)
    {
        enc_x[i] = bit(x-1, i);
    }

    // Compute f(x) := Coeffs^(−1)(B_f · enc(x)) 
    CoeffsInv(vec_f, conv<vec_ZZ>(B_f*enc_x), n);

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
// Computes the squared norm of a vector v of integers with modulo q.
//==============================================================================
ZZ  Norm2m(const vec_zz_p& v, const long& q)
{
    long    i;
    ZZ      norm2, thresh, v_i;

    thresh = q/2; 
    // NOTE: thresh = floor(q/2);

    norm2 = 0;

    for(i=0; i<(v.length()); i++)
    {
        v_i = conv<ZZ>( v[i] );

        if (v_i > thresh)
        {
            v_i -= q;
        }
        
        // norm2 = norm2 + v[i]*v[i];
        norm2 += sqr( v_i );
    }   

    return norm2;
}


//==============================================================================
// Computes the squared norm of a vector v of polynomials with d coefficients.
//==============================================================================
ZZ  Norm2X(const vec_ZZX& v, const long& d)
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


//==============================================================================
// Computes the squared norm of a vector v of polynomials 
// with d coefficients with modulo q.
//==============================================================================
ZZ  Norm2Xm(const vec_zz_pX& v, const long& d, const long& q)
{
    long    i, j;
    ZZ      norm2, thresh, v_ij;

    thresh = q/2; 
    // NOTE: thresh = floor(q/2);

    norm2 = 0;

    for(i=0; i<(v.length()); i++)
    {
        for(j=0; j<d; j++)
        { 
            v_ij = conv<ZZ>( coeff(v[i], j) );

            if (v_ij > thresh)
            {
                v_ij -= q;
            }
            
            // norm2 = norm2 + v[i][j] * v[i][j];
            norm2 += sqr( v_ij );
        }
    }  

    return norm2;
}


//=================================================================================
// Computes the squared norm of a vector v of doubles.
//=================================================================================
double  Norm2D(const vec_D& v)
{
    long    i;
    double  norm2;

    norm2 = 0;

    for(i=0; i<(v.length()); i++)
    {
        norm2 += v[i] * v[i];
    }   

    return norm2;
}


//=================================================================================
// Computes the inner product of two vectors a & b of doubles.
//=================================================================================
double  InnerProdD(const vec_D& a, const vec_D& b)
{
    long    i, len;
    double      prod;

    len = a.length();
    assert(len == b.length());

    prod = 0;
    
    for(i=0; i<len; i++)
    {
        prod += a[i] * b[i];
    }   

    return prod;
}
