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

#include "Rej.h"


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
// - z, v:  vectors of m1 (or m2) polynomials of length d_hat
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
// BimodalRej_v_zzp - Bimodal rejection sampling (LNP22, Lemma 2.14(3)).
//
// Accepts with probability
//      p = 2 / ( e^(y+x) + e^(y-x) ),
//      x = <z,v> / s^2,
//      y = log(M) - ||v||^2 / (2 s^2).
//
// Input:
// - z, v:   vectors mod q (centered representatives are used internally)
// - q:      modulus
// - s:      standard deviation of the Gaussian mask
// - log_M:  log(M), precomputed by the caller.
//           For M = exp(1/(2*alpha^2)) it is exactly 1/(2*alpha^2).
//
// Output:
// - 1 if the sample is accepted, 0 if it is rejected.
//
// NOTE: <z,v> and ||v||^2 are computed EXACTLY with 128-bit integers.
//       The acceptance probability is still computed in RR (150-bit precision),
//       so there is no precision loss with respect to a pure RR/ZZ version.
//==============================================================================
long BimodalRej_v_zzp(const vec_zz_p& z, const vec_zz_p& v, const long& q,
                      const RR& s, const RR& log_M)
{
    if (z.length() != v.length())
    {
        throw std::invalid_argument("BimodalRej_v_zzp: vector length mismatch");
    }

    if (s <= 0)
    {
        throw std::invalid_argument("BimodalRej_v_zzp: s must be positive");
    }

    if (log_M < 0)
    {
        throw std::invalid_argument("BimodalRej_v_zzp: M must be at least one");
    }

    const long half = q / 2;

    // Temporaries declared once, outside the loop, to avoid
    // a memory allocation for every coefficient.
    ZZ dot(0), norm(0), z_zz, v_zz, tmp;

    for (long i = 0; i < z.length(); ++i)
    {
        long z_i = rep(z[i]);      // representative in [0, q)
        long v_i = rep(v[i]);

        // Center the representatives in (-q/2, q/2]
        if (z_i > half) z_i -= q;
        if (v_i > half) v_i -= q;

        conv(z_zz, z_i);
        conv(v_zz, v_i);

        mul(tmp, z_zz, v_zz);   add(dot,  dot,  tmp);   // dot  += z_i * v_i
        sqr(tmp, v_zz);         add(norm, norm, tmp);   // norm += v_i^2
    }

    // RR (150-bit precision)
    const RR s_squared = sqr(s);

    const RR x = conv<RR>(dot) / s_squared;

    const RR y = log_M - conv<RR>(norm) / (RR(2) * s_squared);

    const RR acceptance_probability = RR(2) / (exp(y + x) + exp(y - x));

    if (acceptance_probability < 0 || acceptance_probability > 1)
    {
        throw std::logic_error("BimodalRej_v_zzp: invalid acceptance probability");
    }

    const RR uniform = random_RR();

    return uniform <= acceptance_probability ? 1 : 0;
}
