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

#include <fstream>
#include <numeric> 

extern long idx_Com, idx_ISIS; // Global variables, for benchmarking purposes
mat_D       Perfo;


//=================================================================================
// stats - Compute statistics on an input vector, return the results as a string.
//=================================================================================
string stats(const vec_D v)
{
    double  min, max, avg, sigma, sum, diff;
    long    i, len;

    // Find minimum and maximum values
    len = v.length();
    min = v[0];
    max = v[0];

    for (i = 1; i < len; ++i) 
    {
        if (v[i] < min)
        {
           min = v[i];
        }
        if (v[i] > max)
        {
           max = v[i];
        }
    }

    // Compute average    
    avg = accumulate(v.begin(), v.end(), 0.0) / len; 
    
    // Compute standard deviation
    sum = 0.0;    

    for (i = 0; i < len; ++i) 
    { 
        diff = v[i] - avg;
        sum += diff * diff; 
    }
    sigma = sqrt(sum / len);

    // Return a string with results
    char buf[100];
    sprintf(buf, "min = %-8.5f  max = %-8.5f  avg = %-8.5f  std = %-8.5f", min, max, avg, sigma);

    return string(buf);
}



//=========================================================================================
// bench - Benchmarks for the framework for Post-Quantum Anonymous Verifiable Credentials 
//         defined by Bootle, Lyubashevsky, Nguyen, and Sorniotti (BLNS) in:
//         https://eprint.iacr.org/2023/560.pdf
//=========================================================================================
int main()
{
    zz_p::init(q0); // Initialize modulus q
    
    vec_UL          idx_pub, idx_hid;
    ISK_t           isk;
    uint8_t        *ipk;
    uint8_t         seed_crs[SEED_LEN];
    Vec<string>     attrs;
    mat_zz_p        B_f;
    CRS2_t          crs;
    CRED_t          cred;
    VP_t            VP;
    long            iter, N, valid;
    double          t1, t2, ta, tb;  

    idx_pub = conv<vec_UL>("[4 5 6 7]");    // Indexes of disclosed attributes (revealed, i.e. idx)
    idx_hid = Compute_idx_hid(idx_pub);     // Indexes of undisclosed attributes (hidden, i.e. \overline{\idx})
    // NOTE: both are vectors of non-negative integers in ascending order (one could be the empty array)
    // NOTE: in principle, Holder can use different indexes during Issuing and Presentation protocols
 
    N = 100; //1100; // Number of iterations, for benchmarking purposes

    Perfo.SetDims(10, N);
    
    for(iter=0; iter<N; iter++)
    {
        #ifdef VERBOSE
        cout << "\n#####################################################################" << endl;
        cout << "  ITERATION: " << iter+1 << " of " << N << endl;
        cout << "#####################################################################" << endl;
        
        cout << "\n- Issuer.KeyGen         (key generation)" << endl;
        #endif
        t1 = GetWallTime();
        I_KeyGen(&ipk, isk);
        t2 = GetWallTime();
        #ifdef VERBOSE
        cout << "  CPU time: " << (t2 - t1) << " s" << endl;
        #endif
        Perfo[0][iter] = t2 - t1;

        #ifdef VERBOSE
        cout << "\n- Holder.Init           (init common random string and matrices)" << endl;
        #endif
        ta = GetWallTime();
        H_Init(crs, B_f, seed_crs, attrs, idx_hid.length());
        tb = GetWallTime();
        #ifdef VERBOSE
        cout << "  CPU time: " << (tb - ta) << " s" << endl;
        #endif
        Perfo[1][iter] = tb - ta;

        #ifdef USE_ISSUER_SIGNATURE // Issuer Signature on Plaintext VC

            uint8_t        *Rho;
        
            #ifdef VERBOSE
            cout << "\n=====================================================================" << endl;
            cout << "  ISSUING PROTOCOL  --  Issuer Signature" << endl;
            cout << "=====================================================================" << endl;

            cout << "\n- Issuer.VerCred_Plain  (sign plaintext attributes)" << endl;
            #endif
            ta = GetWallTime();            
            I_VerCred_Plain(&Rho, B_f, ipk, isk, attrs);
            tb = GetWallTime();
            #ifdef VERBOSE
            cout << "  CPU time: " << (tb - ta) << " s" << endl;
            #endif
            // cout << "  attrs  = " << attrs << endl;
            Perfo[3][iter] = tb - ta;

            #ifdef VERBOSE
            cout << "\n- Holder.VerCred_Plain  (verify signature and store VC)" << endl;
            #endif
            ta = GetWallTime();        
            H_VerCred_Plain(cred, ipk, B_f, &Rho, attrs);
            tb = GetWallTime();
            #ifdef VERBOSE       
            cout << "  CPU time: " << (tb - ta) << " s" << endl;
            #endif
            Perfo[4][iter] = tb - ta;
            assert(cred.valid);
            
        #endif
        // #else
        #ifdef USE_ISSUER_BLIND_SIGNATURE

            Vec<string>     attrs_prime;
            RHO1_t          Rho1;
            uint8_t        *Rho2;
            STATE_t         state;
            
            #ifdef VERBOSE
            cout << "\n=====================================================================" << endl;
            cout << "  ISSUING PROTOCOL  --  Blind Signature" << endl;
            cout << "=====================================================================" << endl;

            cout << "\n- Holder.VerCred1       (prove knowledge of undisclosed attributes)" << endl;
            #endif
            ta = GetWallTime();
            H_VerCred1(Rho1, state, seed_crs, crs, ipk, attrs, idx_pub);
            tb = GetWallTime();        
            #ifdef VERBOSE
            cout << "  CPU time: " << (tb - ta) << " s" << endl;
            cout << "  Prove_Com  trials: " << idx_Com << endl;
            #endif
            Perfo[8][iter] = idx_Com;
            Perfo[2][iter] = tb - ta;

            // Select disclosed attributes, fill with zeros hidden attributes 
            attrs_prime = attrs;
            
            for(auto &i: idx_hid)
            {
                attrs_prime[i] = "0"; // Zero padding
            }
            // cout << "  attrs  = " << attrs << endl;
            // cout << "  attrs' = " << attrs_prime << endl;

            #ifdef VERBOSE
            cout << "\n- Issuer.VerCred        (verify proof and compute blind signature)" << endl;
            #endif
            ta = GetWallTime();
            I_VerCred(&Rho2, seed_crs, crs, B_f, ipk, isk, attrs_prime, idx_pub, Rho1);
            tb = GetWallTime();
            #ifdef VERBOSE
            cout << "  CPU time: " << (tb - ta) << " s" << endl;
            #endif
            Perfo[3][iter] = tb - ta;

            #ifdef VERBOSE
            cout << "\n- Holder.VerCred2       (unblind signature and store credential)" << endl;
            #endif
            ta = GetWallTime();        
            H_VerCred2(cred, ipk, B_f, &Rho2, state);
            tb = GetWallTime();
            #ifdef VERBOSE       
            cout << "  CPU time: " << (tb - ta) << " s" << endl;
            #endif
            Perfo[4][iter] = tb - ta;
            assert(cred.valid);

        #endif
        
        #ifdef VERBOSE
        cout << "\n=====================================================================" << endl;
        cout << "  PRESENTATION PROTOCOL" << endl;
        cout << "=====================================================================" << endl;
        
        cout << "\n- Holder.VerPres        (prove knowledge of signature and attributes)" << endl;
        #endif
        ta = GetWallTime();        
        H_VerPres(VP, cred, seed_crs, crs, ipk, B_f, attrs, idx_pub);
        tb = GetWallTime();
        #ifdef VERBOSE    
        cout << "  CPU time: " << (tb - ta) << " s" << endl;
        cout << "  Prove_ISIS  trials: " << idx_ISIS << endl;
        #endif
        Perfo[9][iter] = idx_ISIS;
        Perfo[5][iter] = tb - ta;
        
        #ifdef VERBOSE
        cout << "\n- Verifier.Verify       (verify proof and authorize)" << endl;
        #endif
        ta = GetWallTime();
        valid = V_Verify(VP, seed_crs, crs, B_f, idx_pub);
        tb = GetWallTime();
        Perfo[6][iter] = tb - ta;
        #ifdef VERBOSE    
        cout << "  CPU time: " << (tb - ta) << " s" << endl;
        
        if (valid)
        {
            cout << "  OK!" << endl;
        }
        #endif
        assert(valid == 1);

        double t3 = GetWallTime();
        Perfo[7][iter] = t3 - t1;

        #ifdef VERBOSE            
            cout << "\n=====================================================================\n";
            cout << "  TOT time: " << (t3 - t1) << " s  (" << (t3 - t2) << " s)" << endl;    
        #else
            cout << "  ITERATION: " << iter+1 << " of " << N << " - TOT time: " << (t3 - t1) << " s" << endl;
        #endif
        
        // Free up memory
        delete[] ipk;
    }

    
    // Store raw measurements in a text file        
    ofstream file;
    file.open("Perfo.txt");    
    file << Perfo;
    file.close();

    // Display the benchmark results
    cout << "\n######################################################################################" << endl;
    cout << "  BENCHMARK RESULTS in seconds (N = " << N << ")" << endl << endl;
        
    cout << "- Issuer.KeyGen:     " << stats(Perfo[0]) << endl;
    cout << "- Holder.Init:       " << stats(Perfo[1]) << endl << endl;
    
    #ifdef USE_ISSUER_SIGNATURE
    cout << "- I.VerCred_Plain:   " << stats(Perfo[3]) << endl;
    cout << "- H.VerCred_Plain:   " << stats(Perfo[4]) << endl << endl;    
    #endif

    #ifdef USE_ISSUER_BLIND_SIGNATURE
    cout << "- Holder.VerCred1:   " << stats(Perfo[2]) << endl;
    cout << "- Issuer.VerCred:    " << stats(Perfo[3]) << endl;
    cout << "- Holder.VerCred2:   " << stats(Perfo[4]) << endl << endl;    
    #endif

    cout << "- Holder.VerPres:    " << stats(Perfo[5]) << endl;
    cout << "- Verifier.Verify:   " << stats(Perfo[6]) << endl << endl;    
    
    cout << "- TOTAL time:        " << stats(Perfo[7]) << endl << endl;    
        
    #ifdef VERBOSE
    cout << "- Prove_Com  trials: " << stats(Perfo[8]) << endl;
    cout << "- Prove_ISIS trials: " << stats(Perfo[9]) << endl << endl;
    #endif

    return 0;
}