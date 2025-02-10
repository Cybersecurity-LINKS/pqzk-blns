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


//=====================================================================
// Main - Proof-of-Concept of BLNS Framework for Anonymous Credentials
//=====================================================================
int main()
{
    zz_p::init(q0); // Initialize modulus q
    
    mat_L           isk;
    IPK_t           ipk;
    string          randomSeed;
    Vec<string>     attrs, attrs_prime;
    mat_zz_p        B_f;
    zz_pX           u;
    vec_ZZ          s;
    vec_ZZX         w;
    ZZ              x;
    CRS2_t          crs;
    PROOF_C_t       Pi;
    STATE_t         state;
    CRED_t          cred;
    VP_t            VP;
    long            iter, N, i, valid;
    double          t1, t2, t3;  
        

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
        
        randomSeed = to_string( RandomBnd(12345678) ); 
        cout << "  seed = " << randomSeed << endl; 

        H_Init(crs, attrs, randomSeed);

        // Initialize a random matrix B_f ∈ Z^(d×t)_q
        B_f = random_mat_zz_p(d0, t0);

        
        cout << "\n=====================================================================" << endl;
        cout << "  ISSUING PROTOCOL" << endl;
        cout << "=====================================================================" << endl;
     
        cout << "\n- Holder.VerCred1  (prove knowledge of undisclosed attributes)" << endl;
        H_VerCred1(u, Pi, state, randomSeed, crs, ipk, attrs);

        // Select disclosed attributes, fill with zeros hidden attributes 
        attrs_prime = attrs;
        
        for(i=0; i<idx_hid; i++)
        {
            attrs_prime[i] = "0"; // Zero padding
        }
        // cout << "  attrs  = " << attrs << endl;
        // cout << "  attrs' = " << attrs_prime << endl;
   

        cout << "\n- Issuer.VerCred   (verify proof and compute blind signature)" << endl;
        I_VerCred(s, w, x, randomSeed, crs, B_f, ipk, isk, attrs_prime, u, Pi);
 

        cout << "\n- Holder.VerCred2  (unblind signature and store credential)" << endl;
        H_VerCred2(cred, ipk, B_f, s, w, x, state);

        assert(cred.valid);
        
        
        cout << "\n=====================================================================" << endl;
        cout << "  PRESENTATION PROTOCOL" << endl;
        cout << "=====================================================================" << endl;

        cout << "\n- Holder.VerPres   (prove knowledge of signature and attributes)" << endl;
        H_VerPres(VP, cred, randomSeed, crs, ipk, B_f, attrs);


        cout << "\n- Verifier.Verify  (verify proof and authorize)" << endl;
        valid = V_Verify(VP, randomSeed, crs, B_f);

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