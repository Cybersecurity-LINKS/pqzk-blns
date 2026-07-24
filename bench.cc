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
// compare_D - Compare two double values for qsort function.
//=================================================================================
int compare_D(const void *a, const void *b)
{
    double da = *(const double*)a;
    double db = *(const double*)b;

    if (da < db)
    {
        return -1;
    }
    if (da > db)
    {
        return 1;
    }
    return 0;
}


//=================================================================================
// compute_median - Compute the median on an input vector.
//=================================================================================
double compute_median(const vec_D v, const long len)
{
    long i;    
    double sorted[len];

    for (i = 0; i < len; i++) 
    {
        sorted[i] = v[i];
    }

    // Sort the array in ascending order
    qsort(sorted, len, sizeof(double), compare_D);

    if (len % 2 == 0)
    {
        // Even elements: average of the two middle elements
        return (sorted[(len/2) - 1] + sorted[len/2]) / 2.0;
    }
    else
    {
        // Odd elements: the middle element
        return sorted[(len-1)/2];
    }

}


//=================================================================================
// stats - Compute statistics on an input vector, return the results as a string.
//=================================================================================
string stats(const vec_D v)
{
    double  min, max, avg, median, sigma, sum, diff;
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

    // Compute median
    median = compute_median(v, len);
    
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
    sprintf(buf, "min = %-8.5f  max = %-8.5f  avg = %-8.5f  med = %-8.5f  std = %-8.5f", min, max, avg, median, sigma);

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
    uint8_t         seed_crs[SEED_LEN], nonce[NONCE_LEN];
    Vec<string>     attrs;
    mat_zz_p        B_f;
    CRS2_t          crs;
    CRED_t          cred;
    VP_t            VP{};
    long            iter, iter_warm, iter_tot, W, N, valid;
    double          t1, t2, ta, tb;
    string          old_timestamp, new_timestamp;

    uint8_t* ipk = nullptr;
    uint8_t* u = nullptr;
    idx_pub = conv<vec_UL>("[0 5 6 13 14 15]");    // Indexes of disclosed attributes (revealed, i.e. idx)
    idx_hid = Compute_idx_hid(idx_pub);     // Indexes of undisclosed attributes (hidden, i.e. \overline{\idx})
    // NOTE: both are vectors of non-negative integers in ascending order (one could be the empty array)
    // NOTE: in principle, Holder can use different indexes during Issuing and Presentation protocols
 
    W = 100;  // Number of warm-up iterations, to be executed before the actual benchmarking iterations
    N = 1000; // Number of iterations, for benchmarking purposes

    Perfo.SetDims(13, N);

    iter_warm = 0;
    iter      = 0;

    for(iter_tot = 0; iter_tot < (W+N); iter_tot++)
    {       
        #ifdef VERBOSE
        if (iter_tot < W)
        {
            cout << "\n#####################################################################" << endl;
            cout << "  WARM-UP:   " << iter_warm+1 << " of " << W << endl;
            cout << "#####################################################################" << endl;
        }
        else
        {
            cout << "\n#####################################################################" << endl;
            cout << "  ITERATION: " << iter+1 << " of " << N << endl;
            cout << "#####################################################################" << endl;
        }
        
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

            uint8_t* Rho = nullptr;
        
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
            uint8_t*        Rho2 = nullptr;
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
                
        cout << "\n- Verifier.Challenge    (generate a random nonce)" << endl;
        #endif
        ta = GetWallTime();
        GenRandBytes(nonce, NONCE_LEN);
        tb = GetWallTime();
        #ifdef VERBOSE
        cout << "  CPU time: " << (tb - ta) << " s" << endl;    

        cout << "\n- Holder.VerPres        (prove knowledge of signature and attributes)" << endl;
        #endif
        ta = GetWallTime();        
        H_VerPres(VP, cred, nonce, seed_crs, crs, ipk, B_f, attrs, idx_pub);
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
        valid = V_Verify(VP, nonce, seed_crs, crs, B_f, idx_pub);
        tb = GetWallTime();
        delete[] VP.Pi;
        VP.Pi = nullptr;
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
            if (iter_tot < W)
            {
                cout << "  WARM-UP:   " << iter_warm+1 << " of " << W << " - TOT time: " << (t3 - t1) << " s" << endl;
            }
            else
            {
                cout << "  ITERATION: " << iter+1      << " of " << N << " - TOT time: " << (t3 - t1) << " s" << endl;
            }
        #endif


        #ifdef VERBOSE
        cout << "\n=====================================================================\n";
        cout << "  UPDATE CREDENTIAL" << endl;
        cout << "=====================================================================\n";
        #endif

        #ifdef USE_ISSUER_SIGNATURE

            uint8_t*        Rho2 = nullptr;
            STATE_t         state;

            #ifdef VERBOSE
            cout << "\n- Holder.ReqUpd_Plain   (request an updated signature)" << endl;
            #endif
            ta = GetWallTime();
            H_ReqUpd_Plain(&u, old_timestamp, new_timestamp, state, attrs, ipk, cred);
            tb = GetWallTime();
            #ifdef VERBOSE
            cout << "  CPU time: " << (tb - ta) << " s" << endl;
            #endif
            Perfo[10][iter] = tb - ta;
            
        #endif
        // #else
        #ifdef USE_ISSUER_BLIND_SIGNATURE

            #ifdef VERBOSE
            cout << "\n- Holder.ReqUpdate      (request an updated signature)" << endl;
            #endif
            ta = GetWallTime();
            H_ReqUpdate(&u, old_timestamp, new_timestamp, state, attrs, ipk);
            tb = GetWallTime();
            #ifdef VERBOSE
            cout << "  CPU time: " << (tb - ta) << " s" << endl;
            #endif
            Perfo[10][iter] = tb - ta;

        #endif
            
        #ifdef VERBOSE
        cout << "\n- Issuer.UpdateSign     (update signature)" << endl;
        #endif
        ta = GetWallTime();
        I_UpdateSign(&Rho2, B_f, ipk, isk, u, old_timestamp, new_timestamp);
        tb = GetWallTime();
        #ifdef VERBOSE
        cout << "  CPU time: " << (tb - ta) << " s" << endl;
        #endif
        Perfo[11][iter] = tb - ta;
        
        #ifdef VERBOSE
        cout << "\n- Holder.VerCred2       (check signature and store credential)" << endl;
        #endif
        ta = GetWallTime();
        H_VerCred2(cred, ipk, B_f, &Rho2, state);
        tb = GetWallTime();
        #ifdef VERBOSE
        cout << "  CPU time: " << (tb - ta) << " s" << endl;
        #endif
        Perfo[12][iter] = tb - ta;
        assert(cred.valid);


        if (iter_tot < W)
        {
            iter_warm++;
        }
        else
        {
            iter++;
        }

        // Free up memory
        delete[] ipk;
        ipk = nullptr;

    } // end for for loop (iter_tot)

    
    // Store raw measurements in a text file        
    ofstream file;
    file.open("Perfo.txt");    
    file << Perfo;
    file.close();

    // Display the benchmark results
    cout << "\n####################################################################################################" << endl;
    cout << "  BENCHMARK RESULTS in seconds (N = " << N << ")" << endl << endl;
        
    cout << "- Issuer.KeyGen:       " << stats(Perfo[0]) << endl;
    cout << "- Holder.Init:         " << stats(Perfo[1]) << endl << endl;
    
    #ifdef USE_ISSUER_SIGNATURE
    cout << "- I.VerCred_Plain:     " << stats(Perfo[3]) << endl;
    cout << "- H.VerCred_Plain:     " << stats(Perfo[4]) << endl << endl;    
    #endif

    #ifdef USE_ISSUER_BLIND_SIGNATURE
    cout << "- Holder.VerCred1:     " << stats(Perfo[2]) << endl;
    cout << "- Issuer.VerCred:      " << stats(Perfo[3]) << endl;
    cout << "- Holder.VerCred2:     " << stats(Perfo[4]) << endl << endl;    
    #endif

    cout << "- Holder.VerPres:      " << stats(Perfo[5]) << endl;
    cout << "- Verifier.Verify:     " << stats(Perfo[6]) << endl << endl;    
    
    cout << "- TOTAL time:          " << stats(Perfo[7]) << endl;
    cout << "####################################################################################################" << endl << endl;
    
    #ifdef USE_ISSUER_BLIND_SIGNATURE
    cout << "- Prove_Com  trials:   " << stats(Perfo[8]) << endl;
    #endif
    cout << "- Prove_ISIS trials:   " << stats(Perfo[9]) << endl << endl;
    
    #ifdef USE_ISSUER_SIGNATURE
    cout << "- Holder.ReqUpd_Plain: " << stats(Perfo[10]) << endl;
    #endif
    // #else
    #ifdef USE_ISSUER_BLIND_SIGNATURE
    cout << "- Holder.ReqUpdate:    " << stats(Perfo[10]) << endl;
    #endif
    cout << "- Issuer.UpdateSign:   " << stats(Perfo[11]) << endl;
    cout << "- Holder.VerCred2:     " << stats(Perfo[12]) << endl;

    return 0;
}