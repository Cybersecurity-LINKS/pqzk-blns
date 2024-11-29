// Copyright 2024 Fondazione LINKS

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
    ZZ_p::init(q1);
    
    mat_L           isk;
    IPK_STRUCT      ipk;
    string          randomSeed;
    Vec<string>     attrs, attrs_prime;
    mat_ZZ_p        B_f;
    ZZ_pX           u;
    vec_ZZ          s;
    vec_ZZX         w;
    ZZ              x;
    CRS_Data2       crs;
    PROOF_Com       Pi;
    STATE_STRUCT    state;
    CRED_STRUCT     cred;
    VP_STRUCT       VP;
    int             iter, N, i, valid;
        

    N = 10;  // Number of iterations, for demonstration purposes
    
    for(iter=1; iter<=N; iter++)
    {
        cout << "\n#####################################################################" << endl;
        cout << "  ITERATION: " << iter << " of " << N << endl;
                        

        cout << "\n- Issuer.KeyGen    (key generation)" << endl;        
        I_KeyGen(ipk, isk);            


        cout << "\n- Holder.Init      (initialize common random string and matrices)" << endl;
        
        randomSeed = to_string( RandomBnd(12345678) ); 
        cout << "  seed = " << randomSeed << endl; 

        H_Init(crs, attrs, randomSeed);

        // Initialize a random matrix B_f ∈ Z^(d×t)_q
        B_f = random_mat_ZZ_p(d0, t0);

        
        cout << "\n=====================================================================" << endl;
        cout << "  ISSUING PROTOCOL" << endl;
        cout << "=====================================================================" << endl;
     
        cout << "\n- Holder.VerCred1  (prove knowledge of undisclosed attributes)" << endl;
        H_VerCred1(u, Pi, state, crs, ipk, attrs);

        // Select disclosed attributes, fill with zeros hidden attributes 
        attrs_prime = attrs;
        
        for(i=0; i<idx_hid; i++)
        {
            attrs_prime[i] = "0"; // Zero padding
        }
        // cout << "  attrs  = " << attrs << endl;
        // cout << "  attrs' = " << attrs_prime << endl;
   

        cout << "\n- Issuer.VerCred   (verify proof and compute blind signature)" << endl;
        I_VerCred(s, w, x, crs, B_f, ipk, isk, attrs_prime, u, Pi);
 

        cout << "\n- Holder.VerCred2  (unblind signature and store credential)" << endl;
        H_VerCred2(cred, ipk, B_f, s, w, x, state);

        assert(cred.valid);
        
        
        cout << "\n=====================================================================" << endl;
        cout << "  PRESENTATION PROTOCOL" << endl;
        cout << "=====================================================================" << endl;

        cout << "\n- Holder.VerPres   (prove knowledge of signature and attributes)" << endl;
        H_VerPres(VP, cred, crs, ipk, B_f, attrs);


        cout << "\n- Verifier.Verify  (verify proof and authorize)" << endl;
        valid = V_Verify(VP, crs, B_f);
                
        if (valid)
        {
            cout << "  OK!" << endl;
        }   
        assert(valid == 1);
    }

    return 0;
}