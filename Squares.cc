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

#include "Squares.h"


//==============================================================================
// sum_of_two_squares(p, t)
// 
// Inputs:
// - p:         nonnegative integer p 
// - t:         integer t ∈ (1, p) such that t^2 ≡ −1 (mod p)
//
// Outputs:
// - r0, r1:    two nonnegative integers such that r0^2 + r1^2 = p,
//              or (-1, -1) if it fails (i.e. ⊥)
// 
// NOTE:        it never fails if p is a prime number equal to 1 modulo 4.
//==============================================================================
void  sum_of_two_squares(ZZ& r0, ZZ& r1, const ZZ& p, const ZZ& t)
{
    ZZ  tmp;  
        
    if (t > (p / 2))
    {
        r1 = p - t;
    }
    else
    {
        r1 = t;
    }

    r0 = p;

    while (sqr(r0) >= p)
    {
        // r0, r1 ← r1, r0 mod r1
        tmp = r0; 
        r0  = r1; 
        r1  = tmp % r1;               
    }

    // if (sqr(r0) + sqr(r1) == p)
    //     return r0, r1        
    // else return ⊥
    if (sqr(r0) + sqr(r1) != p)
    {
        r0 = -1;
        r1 = -1;
    }
}


//==============================================================================
// sum_of_four_squares_2_mod_4(n)
// 
// Input:
// - n:         nonnegative integer n with n ≡ 2 (mod 4)
//
// Outputs:
// - x,y,z,w:   four nonnegative integers such that x^2 + y^2 + z^2 + w^2 = n
//==============================================================================
void  sum_of_four_squares_2_mod_4(ZZ& x, ZZ& y, ZZ& z, ZZ& w, const ZZ& n)
{
    ZZ  s, p, t;

    while (1)
    {
        z = 2 * RandomBnd(n/2 + 1);
        w = 2 * RandomBnd(n/2 + 1) + 1;
        s = sqr(z) + sqr(w);

        while (s > n)
        {
            z = 2 * RandomBnd(n/2 + 1);
            w = 2 * RandomBnd(n/2 + 1) + 1;
            s = sqr(z) + sqr(w);
        }

        p = n - s;

        if (p == 1)
        {
            // return 0, 1, z, w
            x = 0; y = 1;
            break;
        }

        s = RandomBnd(p-1) + 1;       
        t = PowerMod(s, (p - 1)/4, p);

        if ((sqr(t) + 1) % p != 0)
        {
            // Go back to the beginning of the while loop
            continue; 
        }
        
        sum_of_two_squares(x, y, p, t);

        if ((x == -1) && (y == -1))
        {
            // Go back to the beginning of the while loop
            continue;
        }

        // return x, y, z, w
        break;
    }
}


//==============================================================================
// sum_of_four_squares(n)
// 
// Input:
// - n:         nonnegative integer n
//
// Outputs:
// - x,y,z,w:   four nonnegative integers such that x^2 + y^2 + z^2 + w^2 = n
//==============================================================================
void  sum_of_four_squares(ZZ& x, ZZ& y, ZZ& z, ZZ& w, const ZZ& n)
{
    long v;
    ZZ   m, x1, y1, z1, w1, pow2;

    // if (n == 0)
    //     return 0, 0, 0, 0
    x = 0; y = 0; z = 0; w = 0;

    // else
    if (n != 0)
    {
        // Compute the unique integers v, m such that n = 2^v * m
        // where v ≥ 0 and m is an odd integer        
        m = n;
        v = MakeOdd(m);
        
        // x1, y1, z1, w1 = sum_of_four_squares_2_mod_4(2*m)
        sum_of_four_squares_2_mod_4(x1, y1, z1, w1, 2*m);

        if (v % 2 == 0) // v is even
        {
            if ((x1 - z1) % 2 == 0)
            {
                // x1, y1, z1, w1 = x1, z1, y1, w1
                swap(y1, z1);
            }
            else if ((x1 - w1) % 2 == 0)
            {
                // x1, y1, z1, w1 = x1, w1, y1, z1
                swap(y1, w1);
                swap(w1, z1);
            }
            // assert((x1 - y1) % 2 == 0);
            // assert((z1 - w1) % 2 == 0);
           
            if (v == 0)
            {
                // 2^(v/2 - 1) == 2^(-1) == 1/2
                x = abs(x1 - y1) / 2;
                y = abs(x1 + y1) / 2;
                z = abs(z1 - w1) / 2;
                w = abs(z1 + w1) / 2;
            }
            else // v > 0
            {
                pow2 = power2_ZZ(v/2 - 1);
                x = pow2 * abs(x1 - y1);
                y = pow2 * abs(x1 + y1);
                z = pow2 * abs(z1 - w1);
                w = pow2 * abs(z1 + w1);
            }
        }
        else // v is odd
        {
            pow2 = power2_ZZ((v - 1)/2);
            x = pow2 * x1;
            y = pow2 * y1;
            z = pow2 * z1;
            w = pow2 * w1;
        }      
    }
    // return x, y, z, w
}


//==============================================================================
// fast_sum_of_squares(n)
// 
// Input:
// - n:         nonnegative integer n
//
// Outputs:
// - v:         vector with length <= d_hat, containing nonnegative integers 
//              such that  Sum_{i=0}^{d_hat}( v[i]^2 ) = n
//==============================================================================
void  fast_sum_of_squares(vec_ZZ& v, const ZZ& n)
{
    long i = 0;
    ZZ   diff, r;
    
    v.SetLength(d_hat);
    diff = n;

    while (diff > 0)
    {
        r     = FloorToZZ( sqrt( conv<RR>( diff ) ) ) ;
        diff -= sqr(r);

        assert(i < d_hat);
        // NOTE: length of v must not exceed d_hat

        v[i]  = r;
        i++;        
    }

    v.SetLength(i);
    
    // return v
}