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

#include "params.h"
#include "Issuer.h"
#include "Holder.h"
#include "Verifier.h"


//=========================================================================================
// Main - Implementation of the framework for Post-Quantum Anonymous Verifiable Credentials 
//        defined by Bootle, Lyubashevsky, Nguyen, and Sorniotti (BLNS) in:
//        https://eprint.iacr.org/2023/560.pdf
//=========================================================================================
int main()
{
    zz_p::init(q0); // Initialize modulus q
    
    mat_L           isk;
    IPK_t           ipk;
    unsigned char   crs_seed[SEED_LEN];
    Vec<string>     attrs, attrs_prime;
    mat_zz_p        B_f;
    zz_pX           u;
    vec_ZZ          s;
    vec_ZZX         w;
    ZZ              x;
    CRS2_t          crs;
    uint8_t*        Pi; 
    STATE_t         state;
    CRED_t          cred;
    VP_t            VP;
    long            iter, N, i, valid;
    double          t1, t2, t3, ta, tb;  

    N = 10;  // Number of iterations, for demonstration purposes
    
    for(iter=1; iter<=N; iter++)
    {
        cout << "\n#####################################################################" << endl;
        cout << "  ITERATION: " << iter << " of " << N << endl;
        cout << "#####################################################################" << endl;
        
        cout << "\n- Issuer.KeyGen    (key generation)" << endl;
        t1 = GetWallTime();
        I_KeyGen(ipk, isk);
        t2 = GetWallTime();
        cout << "  CPU time: " << (t2 - t1) << " s" << endl;
        
        cout << "\n- Holder.Init      (initialize common random string and matrices)" << endl;
        ta = GetWallTime();

        // Initialize a 32 byte (256 bit) public seed for common random string (crs) structure,
        // using the cryptographically strong pseudo-random number generator from NTL
        RandomStream& RS = GetCurrentRandomStream();
        RS.get(crs_seed, SEED_LEN);
        // for(i=0; i<SEED_LEN; i++)
        // {
        //     printf("%0x", crs_seed[i]);
        // }
        // printf("\n");

        H_Init(crs, attrs, crs_seed);

        // Initialize a random matrix B_f ∈ Z^(d×t)_q
        B_f = random_mat_zz_p(d0, t0);
        tb = GetWallTime();        
        cout << "  CPU time: " << (tb - ta) << " s" << endl;

        
        cout << "\n=====================================================================" << endl;
        cout << "  ISSUING PROTOCOL" << endl;
        cout << "=====================================================================" << endl;
        ta = GetWallTime();     
        cout << "\n- Holder.VerCred1  (prove knowledge of undisclosed attributes)" << endl;
        H_VerCred1(u, &Pi, state, crs_seed, crs, ipk, attrs);
        tb = GetWallTime();        
        cout << "  CPU time: " << (tb - ta) << " s" << endl;

        // Select disclosed attributes, fill with zeros hidden attributes 
        attrs_prime = attrs;
        
        for(i=0; i<idx_hid; i++)
        {
            attrs_prime[i] = "0"; // Zero padding
        }
        // cout << "  attrs  = " << attrs << endl;
        // cout << "  attrs' = " << attrs_prime << endl;
   
        ta = GetWallTime();
        cout << "\n- Issuer.VerCred   (verify proof and compute blind signature)" << endl;
        I_VerCred(s, w, x, crs_seed, crs, B_f, ipk, isk, attrs_prime, u, &Pi);
        tb = GetWallTime();        
        cout << "  CPU time: " << (tb - ta) << " s" << endl;

        cout << "\n- Holder.VerCred2  (unblind signature and store credential)" << endl;
        ta = GetWallTime();        
        H_VerCred2(cred, ipk, B_f, s, w, x, state);
        tb = GetWallTime();        
        cout << "  CPU time: " << (tb - ta) << " s" << endl;
        assert(cred.valid);
        
        
        cout << "\n=====================================================================" << endl;
        cout << "  PRESENTATION PROTOCOL" << endl;
        cout << "=====================================================================" << endl;
        ta = GetWallTime();
        cout << "\n- Holder.VerPres   (prove knowledge of signature and attributes)" << endl;
        H_VerPres(VP, cred, crs_seed, crs, ipk, B_f, attrs);
        tb = GetWallTime();        
        cout << "  CPU time: " << (tb - ta) << " s" << endl;

        ta = GetWallTime();  
        cout << "\n- Verifier.Verify  (verify proof and authorize)" << endl;
        valid = V_Verify(VP, crs_seed, crs, B_f);
        tb = GetWallTime();        
        cout << "  CPU time: " << (tb - ta) << " s" << endl;

        if (valid)
        {
            cout << "  OK!" << endl;
        }   
        assert(valid == 1);
        
        t3 = GetWallTime();
        cout << "\n=====================================================================\n";
        cout << "  TOT time: " << (t3 - t1) << " s (" << (t3 - t2) << " s)" << endl;
    }

    return 0;
}