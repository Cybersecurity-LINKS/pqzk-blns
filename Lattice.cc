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

#include "Lattice.h"


//==============================================================================
// NTRU_TrapGen - Generates a1 and B from parameters d and q.
//
// Inputs:
// - None
//
// Outputs:
// - a1 : polynomial 
// - B  : matrix with size (2d * 2d), basis of the 2d-dimensional lattice
//
// NOTE: Algorithm 2 Master Keygen(N, q) at pag. 15 in [DLP14]
//       equivalent to algorithm NTRU.TrapGen(q, d) in [BLNS23]  
//==============================================================================
void NTRU_TrapGen(zz_pX& a1, mat_L& B)
{
    long    i, j, valid;
    ZZX     f, g, F, G, fr, gr, num, den, inv_den, iphi, a, b, rho_f, rho_g, k;
    zz_pX   inv_f;
    ZZ      acc, res, out, R_f, R_g, u, v;
    RR      gamma, gamma2;
    mat_L   A;

    const ZZ    q          = conv<ZZ>(q0);
    const ZZX   phi         = Phi();
    const RR    max_gamma   = RR(1.17) * sqrt( RR(q0) );
    
    // 1. σ_f ← 1.17*√(q/2d),  σ_f ∈ R
    const double sigma_f    = (double)(1.17) * sqrt( (double)(q0) / (double)(2*d0) );
    
    
    // Loop to find a valid basis (see steps 10, 13, 15)
    valid = 0;
    
    while(valid == 0) 
    {
        // 2. for i = 0 : d−1
        // 3.   f_i ← ZSampler(σ_f, 0),  f_i ∈ Z
        // 4.   g_i ← ZSampler(σ_f, 0),  g_i ∈ Z
        
        // 5. f ← Sum_{i=0}^{d−1}(f_i * x^i),  f ∈ R
        f.SetLength(d0);
        polySampler(f, sigma_f);

        // 6. g ← Sum_{i=0}^{d−1}(g_i * x^i),  g ∈ R
        g.SetLength(d0);
        polySampler(g, sigma_f);

        // 7. fr ← f_0 − Sum_{i=0}^{d−1}(f_(d−i) * x^i),  fr ∈ R
        fr.SetLength(d0);
        // fr[0] = f[0];
        SetCoeff(fr, 0, coeff(f, 0));

        // 8. gr ← g_0 − Sum_{i=0}^{d−1}(g_(d−i) * x^i),  gr ∈ R
        gr.SetLength(d0);
        // gr[0] = g[0];
        SetCoeff(gr, 0, coeff(g, 0));

        for(i=1; i<d0; i++)
        {
            // fr[i] = -f[d0-i];
            SetCoeff(fr, i, -coeff(f, d0-i));

            // gr[i] = -g[d0-i];
            SetCoeff(gr, i, -coeff(g, d0-i));
        }


        // 9. γ ← max{ ∥(g, −f)∥, ∥(q*fr /(f*fr+g*gr)), (q*gr /(f*fr+g*gr))∥ },  γ ∈ R
        acc = 0; 

        // Compute gamma = ∥(g, −f)∥
        for(i=0; i<d0; i++)
        {
            // acc = acc + g[i]^2 + f[i]^2;
            acc += sqr( coeff(g, i) ) + sqr( coeff(f, i) );
        }  
        gamma = sqrt(conv<RR>(acc));
        
        // Compute gamma2 = ∥(q*fr /(f*fr+g*gr)), (q*gr /(f*fr+g*gr))∥
        den.SetLength(d0);
        den = ModPhi(f * fr + g * gr);
        
        // inv_den = inv(den) % phi
        XGCD(res, inv_den, iphi, den, phi, 0);
        a = ModPhi(fr * inv_den);
        b = ModPhi(gr * inv_den);
        acc = 0;
        
        for(i=0; i<d0; i++)
        {        
            // acc = acc + (a[i] * a[i] + b[i] * b[i]);
            acc += sqr( coeff(a, i) ) + sqr( coeff(b, i) );
        }
        
        gamma2 = RR(q0) * sqrt( conv<RR>(acc) ) / conv<RR>(res);
        
        if(gamma2 > gamma)
        {
            gamma = gamma2;
        }
        

        // 10. if γ > 1.17*√q: go to step 2
        if (gamma > max_gamma)
        {
            valid = 0;
            continue; // skip next steps, go to 2
        }
        else
        {
            valid = 1;
        }

        
        // 11.  Using XGCD algorithm, compute ρ_f, ρ_g ∈ R and R_f, R_g ∈ Z
        // 12.1 such that:   −ρ_f*f = R_f mod (x^d + 1)         
        XGCD(R_f, rho_f, iphi, -f, phi, 0); 
        
        // 13.2 if GCD(R_f, q) ≠ 1,  go to step 2
        if(GCD(R_f, q)!=1)
        {
            valid = 0;
            continue; // skip next steps, go to 2
        }
        else
        {
            valid = 1;
        }
        
        // 12.2              −ρ_g*g = R_g mod (x^d + 1)
        XGCD(R_g, rho_g, iphi, -g, phi, 0); 
        
        // 13.1 if GCD(R_f, R_g) ≠ 1,  go to step 2        
        // NOTE: this check is combined with steps 14, 15

        // 14. Using XGCD algorithm, compute u, v ∈ Z        
        XGCD(out, u, v, R_f, R_g);

        // 15. such that u*R_f + v*R_g = 1
        if(out!=1)
        {
            valid = 0;
            continue; // skip next steps, go to 2
        }
        else
        {
            valid = 1;
        }
              
    } // End of while loop (valid f, g pair)


    // 16. F ← q·v·ρ_g F ∈ R
    F.SetLength(d0);
    F =  q * v * rho_g;

    // 17. G ← −q·u·ρ_f G ∈ R
    G.SetLength(d0);  
    G = -q * u * rho_f;

    k = 1;

    // NOTE: loop to properly reduce F, G
    while( deg(k) >= 0 )
    {
        // 18. k = [(F*fr + G*gr)/(f*fr + g*gr)] ∈ R
        num = ModPhi(F * fr + G * gr);
        // den = (f * fr + g * gr) % (phi);    
        // XGCD(res, inv_den, iphi, den, phi, 0);    
        k = ModPhi(num * inv_den);
                
        for(i=0; i<d0; i++)
        {
            // NOTE: k[i] = floor( k[i] / res )            
            SetCoeff(k, i, ( coeff(k, i) / res ) );
        }
        
        // 19. F ← F − k*f,  F ∈ R
        F = ModPhi(F - k*f);

        // 20. G ← G − k*g,  G ∈ R
        G = ModPhi(G - k*g);
    }
          
    // 21. a1 ← g*f^(−1) mod q ∈ R_q   
    a1.SetLength(d0);
    inv_f = InvMod( conv<zz_pX>(f), conv<zz_pX>(phi) );
    a1    = ModPhi_q( conv<zz_pX>(g) * inv_f );

    
    // 22. B ← [rot(g) −rot(f); 
    //          rot(G) −rot(F)] ∈ Z^(2d×2d)
    B.SetDims(2*d0, 2*d0);

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
     

    // 23. return (a1, B)
}


//==============================================================================
// preGSampler - Computes the Gaussian sampling 
//
// Inputs:
// - B:     matrix B ∈ Z^(2d×2d)
// - sigma: standard deviation sigma > 0
// - c:     center, vector of integers c ∈ Z^(2d)
//
// Output:
// - v:     vector of integers, v ∈ Z^(2d), 
//          such that s := c − v is a short vector
//
// NOTE: Algorithm 1 Gaussian Sampler(B, sigma, c) at pag. 13 in [DLP14]
//       Dimension n in [DLP14] corresponds to 2*d0 in [BLNS23]
// NOTE: Optimized version, with OGS_Ortho & all double instead of RR.
//==============================================================================
void preGSampler(vec_ZZ& v, const mat_L& B, const double& sigma, const vec_ZZ& c)
{
    long    i, valid;
    mat_D   Bt;
    vec_D   Norms2;    
    vec_ZZ  ci, zibi;
    double  cpi, spi;
    ZZ      zi;

    const ZZ thres_s = sqr( ZZ(sigma0) ) * ZZ(2*d0);

    zibi.SetLength(2*d0);

    // 1. (b˜_1, ..., b˜_(2d)) ← GramSchmidt.Orthogonalization(B),   b˜_i ∈ R^(2d)
    OGS_Ortho(Bt, Norms2, B);


    // NOTE: loop to find a valid s = c - v  (i.e. small norm, see Holder.VerCred2, row 5)
    valid = 0;
    v.SetLength(2*d0);
    ci.SetLength(2*d0);
        
    while(valid == 0) 
    {
        // 2. v_(2d) ← 0,   v_(2d) ∈ Z^(2d)    
        clear(v);

        // 3. c_(2d) ← c,  c_(2d) ∈ Z^(2d)    
        ci = c;

        // 4. for i ← 2d, ..., 1 do
        for(i=(2*d0-1); i>=0; i--)    
        {
            // 5. c′_i ← ⟨c_i, b˜_i⟩ / ∥b˜_i∥^2,   c′_i ∈ R
            cpi   = InnerProdD( conv<vec_D>(ci), Bt[i] ) / Norms2[i];

            // 6. σ′_i ← σ / ∥b˜_i∥,   σ′_i ∈ R
            spi = sigma / sqrt(Norms2[i]);

            // 7. z_i ← ZSampler(σ′_i, c′_i),   z_i ∈ Z
            ZSampler(zi, spi, cpi); 

            // 8. c_(i−1) ← c_i − z_i*b_i,   c_(i−1) ∈ Z^(2d)
            zibi = zi * conv<vec_ZZ>( B[i] );
            ci  -= zibi;

            // 9. v_(i−1) ← v_i + z_i*b_i,   v_(i−1) ∈ Z^(2d)
            // v += zibi;
        }

        // NOTE: step 9 moved outside of the loop, for efficiency
        v = c - conv<vec_ZZ>( ci );

        if ( Norm2(c - v) <=  thres_s )
        {
            valid = 1;
        }
    }

    // 10. return v_0
}


//==============================================================================
// GSampler - Computes the Gaussian sampling 
// Inputs:
// - h:     part of public key, i.e. polynomial        a_1 ∈ R_q
// - a:     part of public key, i.e. polynomial vector a_2 ∈ R_q^m
// - B:     secret key,         i.e. matrix            B   ∈ Z^(2d×2d)
// - sigma: standard deviation  sigma > 0
// - d:     center (= f(x)+u ), i.e. polynomial  d ∈ R_q
//
// Outputs:  
// - s:     short vector,       s ∈ Z^(2d)
// - w:     polynomial vector,  w ∈ R^m
//==============================================================================
void GSampler(vec_ZZ& s, vec_ZZX& w, const zz_pX& h, const vec_zz_pX& a, const mat_L& B, const double& sigma, const zz_pX& d)
{
    long    i, j, valid;  
    ZZX     u;
    mat_ZZ  A, Id;
    mat_L   R;        
    vec_ZZ  c, d_u, v;  

    const ZZ thres_w = sqr( ZZ(sigma0) ) * ZZ(d0*m0);             
 
    // NOTE: loop to find a valid w (i.e. small norm, see Holder.VerCred2, row 5)
    valid = 0;
    w.SetLength(m0);
    u.SetLength(d0); // u ∈ R
    
    while(valid == 0) 
    {
        // 1. u ← 0    
        clear(u);

        // 2. a ← (a_1, ... , a_m) ∈ R_q^m,   a_i ∈ R_q
        
        // 3. for i ← 1, ... , m  do 
        for(i=0; i<m0; i++)
        {
            // 4. w_i ← polySampler(σ, 0),  w_i ∈ R     
            polySampler(w[i], sigma);

            // 5. u ← u + w_i * a_i,   u ∈ R            
            u += ModPhi( w[i] * conv<ZZX>( a[i]) );
        }
        // 6. w ← (w_1, ... , w_m),   w ∈ R^m 
        
        if ( Norm2X(w, d0)  <=  thres_w )
        {
            valid = 1;
        }        
    }

    
    // 7. A ← [  I_d   ]
    //        [ rot(h) ] ∈ Z^(2d×d)
    A.SetDims(2*d0, d0);

    Id.SetDims(d0, d0);    
    Id = ident_mat_ZZ(d0);
 
    for(i=0; i<d0; i++)
    {
        for(j=0; j<d0; j++)
        {
            A[i][j] = Id[i][j];
        }
    }

    Id.kill();
    
    // R = rot(h) 
    rot(R, conv<ZZX>(h)); 

    for(i=0; i<d0; i++)
    {
        for(j=0; j<d0; j++)
        {
            A[i+d0][j] = conv<ZZ>( R[i][j] );
        }
    }

    R.kill();
    
    
    // 8. c ← LinearSolve(c*A = d − Coeffs(u)),   c ∈ Z^(2d)
    c.SetLength(2*d0);

    // Compute (d − Coeffs(u))
    d_u.SetLength(d0);

    for(i=0; i<d0; i++)
    {
        d_u[i] = conv<ZZ>( coeff(d, i) ) - coeff(u, i); 
    }

    // NOTE: LinearSolve is not strictly needed, set the straightforward solution 
    //       c = [(d_u), 0 ... 0] ∈ Z^(2d)
    for(i=0; i<d0; i++)
    {
        c[i]    = d_u[i]; 
        c[i+d0] = 0;
    }
              
           
    // 9. v ← preGSampler(B, σ, c),   v ∈ Z^(2d)    
    preGSampler(v, B, sigma, c); 

    // 10. s ← c − v,   s ∈ Z^(2d)    
    s.SetLength(2*d0);    
    s = c - v;    

    // 11. return (s, w)
}


// ==============================================================================
// ZSampler - Generates a discrete Gaussian sample  
//
// Inputs:
// - sigma: standard deviation sigma > 0
// - c:     center
//
// Output:
// - x:     integer sampled from the discrete Gaussian distribution
// ==============================================================================
void ZSampler(ZZ& x, const double& sigma, const double& c) 
{
    long    b;
    ZZ      left, right, range, c_int;    
    double  delta, iden, val, p, u, c_floor, c_frac;

    if (sigma<=0)
    {
        cout << "ERROR! Sigma must be > 0 (sigma = " << (sigma) << ")" << endl;
        assert(sigma>0);
    }
    
    // NOTE: split c, to reduce numerical approximations
    c_floor = floor(c);
    c_frac  = c - c_floor;
    c_int   = conv<ZZ>(c_floor);
               
    b = 0;
    delta = sigma * log2(lambda0);
    left  = conv<ZZ>(ceil ( c_frac - delta ));
    right = conv<ZZ>(floor( c_frac + delta ));
    range = right - left + 1;
    
    if (range < 2)
    {
        cout << "WARNING! ZSampler: small sigma (sigma = "   << sigma << ", range = " << range << ")" << endl;
        
        // NOTE: workaround with inaccurate sigma, but at least correct center c
        u = conv<double>(random_RR());

        if (u >= c_frac)
        {
            x = c_int;
            return;
        }
        else // (u < c_frac)
        {
            x = c_int+1;
            return;
        }
    }

    iden   = 1 / (2 * sigma * sigma);

    while(b == 0)
    {     
        // x = RandomBnd(range) + left + c_int; 
        x = RandomBnd(range) + left;  
        // NOTE: add c_int at the end, to reduce numerical approximations
        
        // p ← ρ_{σ,c}(x)        
        val = conv<double>(x) - c_frac;
        val = -(val * val) * iden;
        
        if (val > 0)
        {
            p = 1;
        }
        else if (val < -1000)
        {
            // NOTE: exp(-1000) = 5,076e-435
            p = 0;
            continue;
        }
        else
        {        
            p = exp(val); 
        }        

        // b ← Bernoulli(p) 
        u = conv<double>(random_RR());

        if (u < p)
        {
            b = 1;                
        }
    }
    
    // return (x + c_int);
    x += c_int;
}


//==============================================================================
// polySampler - Generates a random polynomial with d elements 
//               and fixed standard deviation
//
// Input:
// - sigma: standard deviation sigma > 0
//
// Output:
// - s:     polynomial, s ∈ R
//==============================================================================
void polySampler(ZZX& s, const double& sigma)
{
    long    i;
    ZZ      x;

    s.SetLength(d0);
   
    for(i=0; i<d0; i++)
    {
        ZSampler(x, sigma, 0);
        SetCoeff(s, i, x);
    }       

    // NOTE: ensure that deg(s) == d0-1, to avoid errors in NTRU_TrapGen    
    if (coeff(s, d0-1) == 0)
    {
        SetCoeff(s, d0-1, 1);
    }
      
    // return s; // s ∈ R
}
