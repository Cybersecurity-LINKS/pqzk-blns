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
    
    long            i;
    mat_L           isk; 
    IPK_STRUCT      ipk;
    float           diff;
    clock_t         tstart, t1, t2;  
    string          inputStr;     
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
    int             result; 
    
    {             
        cout << "\n=====================================================================\n";
        cout << " ISSUING PROTOCOL";
        cout << "\n=====================================================================\n";
        
        cout << "\n- Issuer.KeyGen \n"; 
               
        // Issuer.KeyGen (i.e. AnonCreds.Init in [BLNS23])
        tstart = clock();
        t1 = tstart;
        I_KeyGen(ipk, isk);                
        t2 = clock();

        diff = ((float)t2 - (float)t1)/(float)CLOCKS_PER_SEC;    
        cout << "CPU time: " << diff << " seconds" << endl;

        //============================================================================== 
        cout << "\n- Holder.Init \n";      
        
        inputStr = to_string( RandomBnd(12345) ); 
        cout << "inputStr = " << inputStr << endl;   

        // Initialize the public randomly chosen matrix B_f ∈ Z^(nd×t)_q, with n = 1
        B_f = random_mat_ZZ_p(d0, t0);
        // cout << " B_f = \n" << B_f << endl;  
       
        t1 = clock();
        H_Init(crs, attrs, inputStr);                        
        t2 = clock();
    
        diff = ((float)t2 - (float)t1)/(float)CLOCKS_PER_SEC;    
        cout << "CPU time: " << diff << " seconds" << endl;

        //============================================================================== 
        cout << "\n- Holder.VerCred1 \n";
    
        t1 = clock();
        H_VerCred1(u, Pi, state, crs, ipk, attrs);                        
        t2 = clock();
   
        diff = ((float)t2 - (float)t1)/(float)CLOCKS_PER_SEC;    
        cout << "CPU time: " << diff << " seconds" << endl;

        //==============================================================================
        cout << "\n- Issuer.VerCred \n";

        // attrs′ ← (attrs_{idx_pub} | {0}_{idx_hid})
        attrs_prime = attrs;

        // NOTE: select disclosed attributes, fill with zeros for each i ∈ idx_hid 
        for(i=0; i<idx_hid; i++)
        {
            attrs_prime[i] = "0"; // Zero padding
        }
        // cout << "  attrs  = " << attrs << endl;
        // cout << "  attrs' = " << attrs_prime << endl;
        
        t1 = clock();
        I_VerCred(s, w, x, crs, B_f, ipk, isk, attrs_prime, u, Pi);
        t2 = clock();
        
        diff = ((float)t2 - (float)t1)/(float)CLOCKS_PER_SEC;    
        cout << "CPU time: " << diff << " seconds" << endl;

        //==============================================================================
        cout << "\n- Holder.VerCred2 \n";

        t1 = clock();
        H_VerCred2(cred, ipk, B_f, s, w, x, state);
        t2 = clock();

        assert(cred.valid);

        diff = ((float)t2 - (float)t1)/(float)CLOCKS_PER_SEC;    
        cout << "CPU time: " << diff << " seconds" << endl;
        
        
        cout << "\n=====================================================================\n";
        cout << " PRESENTATION PROTOCOL";
        cout << "\n=====================================================================\n";

        cout << "\n- Holder.VerPres \n";

        t1 = clock();
        H_VerPres(VP, cred, crs, ipk, B_f, attrs);
        t2 = clock();

        diff = ((float)t2 - (float)t1)/(float)CLOCKS_PER_SEC;    
        cout << "CPU time: " << diff << " seconds" << endl;

        //==============================================================================        
        cout << "\n- Verifier.Verify \n";

        t1 = clock();
        result = V_Verify(VP, crs, B_f);
        t2 = clock();
        
        assert(result == 1);
   
        diff = ((float)t2 - (float)t1)/(float)CLOCKS_PER_SEC;    
        cout << "CPU time: " << diff << " seconds" << endl;
        //==============================================================================

        cout << "\n=====================================================================\n";
        diff = ((float)t2 - (float)tstart)/(float)CLOCKS_PER_SEC;    
        cout << "TOT time: " << diff << " seconds" << endl;
    }

    return 0;
}