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
#include "Utils.h"
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
    
    vec_UL          idx_pub, idx_hid;
    ISK_t           isk;
    uint8_t        *ipk;
    uint8_t         seed_crs[SEED_LEN];
    Vec<string>     attrs, attrs_prime;
    mat_zz_p        B_f;
    CRS2_t          crs;
    RHO1_t          Rho1;
    uint8_t        *Rho2;
    STATE_t         state;
    CRED_t          cred;
    VP_t            VP;
    long            iter, N, valid;
    double          t1, t2, t3, ta, tb;  

    idx_pub = conv<vec_UL>("[0 1 2 3]");    // Indexes of disclosed attributes (revealed, i.e. idx)
    idx_hid = Compute_idx_hid(idx_pub);     // Indexes of undisclosed attributes (hidden, i.e. \overline{\idx})
    // NOTE: both are vectors of non-negative integers in ascending order (one could be the empty array)
    // NOTE: in principle, Holder can use different indexes during Issuing and Presentation protocols
 
    N = 10;  // Number of iterations, for demonstration purposes
    
    for(iter=1; iter<=N; iter++)
    {
        cout << "\n#####################################################################" << endl;
        cout << "  ITERATION: " << iter << " of " << N << endl;
        cout << "#####################################################################" << endl;
        
        cout << "\n- Issuer.KeyGen    (key generation)" << endl;
        t1 = GetWallTime();
        I_KeyGen(&ipk, isk);
        t2 = GetWallTime();
        cout << "  CPU time: " << (t2 - t1) << " s" << endl;
        
        cout << "\n- Holder.Init      (initialize common random string and matrices)" << endl;
        ta = GetWallTime();
        H_Init(crs, B_f, seed_crs, attrs, idx_hid.length());
        tb = GetWallTime();
        cout << "  CPU time: " << (tb - ta) << " s" << endl;

        
        cout << "\n=====================================================================" << endl;
        cout << "  ISSUING PROTOCOL" << endl;
        cout << "=====================================================================" << endl;
        ta = GetWallTime();     
        cout << "\n- Holder.VerCred1  (prove knowledge of undisclosed attributes)" << endl;
        H_VerCred1(Rho1, state, seed_crs, crs, ipk, attrs, idx_pub);
        tb = GetWallTime();        
        cout << "  CPU time: " << (tb - ta) << " s" << endl;

        // Select disclosed attributes, fill with zeros hidden attributes 
        attrs_prime = attrs;
        
        for(auto &i: idx_hid)
        {
            attrs_prime[i] = "0"; // Zero padding
        }
        // cout << "  attrs  = " << attrs << endl;
        // cout << "  attrs' = " << attrs_prime << endl;

        ta = GetWallTime();
        cout << "\n- Issuer.VerCred   (verify proof and compute blind signature)" << endl;
        I_VerCred(&Rho2, seed_crs, crs, B_f, ipk, isk, attrs_prime, idx_pub, Rho1);
        tb = GetWallTime();        
        cout << "  CPU time: " << (tb - ta) << " s" << endl;

        cout << "\n- Holder.VerCred2  (unblind signature and store credential)" << endl;
        ta = GetWallTime();        
        H_VerCred2(cred, ipk, B_f, &Rho2, state);
        tb = GetWallTime();        
        cout << "  CPU time: " << (tb - ta) << " s" << endl;
        assert(cred.valid);

        
        cout << "\n=====================================================================" << endl;
        cout << "  PRESENTATION PROTOCOL" << endl;
        cout << "=====================================================================" << endl;
        ta = GetWallTime();
        cout << "\n- Holder.VerPres   (prove knowledge of signature and attributes)" << endl;
        H_VerPres(VP, cred, seed_crs, crs, ipk, B_f, attrs, idx_pub);
        tb = GetWallTime();        
        cout << "  CPU time: " << (tb - ta) << " s" << endl;

        ta = GetWallTime();  
        cout << "\n- Verifier.Verify  (verify proof and authorize)" << endl;
        valid = V_Verify(VP, seed_crs, crs, B_f, idx_pub);
        tb = GetWallTime();        
        cout << "  CPU time: " << (tb - ta) << " s" << endl;

        if (valid)
        {
            cout << "  OK!" << endl;
        }   
        assert(valid == 1);

        // Free up memory
        delete[] ipk;
        
        t3 = GetWallTime();
        cout << "\n=====================================================================\n";
        cout << "  TOT time: " << (t3 - t1) << " s (" << (t3 - t2) << " s)" << endl;
    }

    return 0;
}