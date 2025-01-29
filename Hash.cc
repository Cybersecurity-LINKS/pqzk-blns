#include "Hash.h"
#include "Utils.h"


//==============================================================================
// CustomHash - Custom Hash function, implemented using SHAKE256
// 
// Inputs:
// - x:         integer (initial seed)
// - out_len:   length of the output value, in bytes
//
// Output:
// - y:         integer (out_len bytes) 
//==============================================================================
long int CustomHash(const long int x, const size_t out_len) 
{
    EVP_MD_CTX *mdctx;
    const EVP_MD *md;
    char alg[] = "shake256";
    long int y = 0;


    unsigned char* x_arr = new unsigned char[sizeof(long int)];
    unsigned char* y_arr = new unsigned char[out_len];

    memcpy(x_arr, &x, sizeof(long int));

    md = EVP_get_digestbyname(alg);
    if (md == NULL) {
        printf("Unknown message digest %s\n", alg);
        exit(1);
    }

    mdctx = EVP_MD_CTX_new();
    //if(!EVP_DigestInit_ex(mdctx, md, NULL)) {
    if (!EVP_DigestInit_ex2(mdctx, md, NULL)) {          
        // NOTE: EVP_DigestInit_ex2 requires OpenSSL >= 3.x
        printf("Message digest initialization failed.\n");
        EVP_MD_CTX_free(mdctx);
        exit(1);
    }
    if (!EVP_DigestUpdate(mdctx, x_arr, sizeof(x))) {
        printf("Message digest update failed.\n");
        EVP_MD_CTX_free(mdctx);
        exit(1);
    }

    //if (!EVP_DigestFinalXOF(mdctx, y_arr, out_len)) {
    if (!EVP_DigestSqueeze(mdctx, y_arr, out_len)) {
        printf("Message digest squeeze failed.\n");
        EVP_MD_CTX_free(mdctx);
        exit(1);
    }
    // NOTE: EVP_DigestSqueeze (instead of EVP_DigestFinalXOF) requires OpenSSL >= 3.3.0
    
    EVP_MD_CTX_free(mdctx);

    if (out_len > sizeof(y)) {
        memcpy(&y, y_arr, sizeof(long int));
    } else {
        memcpy(&y, y_arr, out_len);
    }

    delete[] y_arr;
    delete[] x_arr;

    return y;
}


//==============================================================================
// Hash_Init -  Initialize the Custom Hash function, implemented using SHAKE256
// 
// Inputs:
// - inputStr:  string containing the input message (initial seed)
//
// Output:
// - mdctx:     digest context in OpenSSL
//==============================================================================
EVP_MD_CTX* Hash_Init(const string& inputStr)
{
    EVP_MD_CTX *mdctx;
    const EVP_MD *md;
    const size_t in_len = inputStr.length(); 
    const char alg[] = "shake256";

    md = EVP_get_digestbyname(alg);
    if (md == NULL) {
        printf("Unknown message digest %s\n", alg);
        exit(1);
    }

    mdctx = EVP_MD_CTX_new();
    if (!EVP_DigestInit_ex2(mdctx, md, NULL)) {   
        printf("Message digest initialization failed.\n");
        EVP_MD_CTX_free(mdctx);
        exit(1);
    }
    if (!EVP_DigestUpdate(mdctx, inputStr.data(), in_len)) {
        printf("Message digest update failed.\n");
        EVP_MD_CTX_free(mdctx);
        exit(1);
    }
    
    return mdctx;
}


//==============================================================================
// Hash_zz_pX - Generate a random polynomial using the Custom Hash function
// 
// Inputs:
// - mdctx:     digest context in OpenSSL
// - n_coeffs:  number of coefficients of the random polynomial (i.e. d_hat)
// - b_coeffs:  number of bytes for each coefficient (i.e. |q_hat|)
//
// Output:
// - out_poly:  random polynomial with n_coeffs coefficients (mod q_hat)
//==============================================================================
void Hash_zz_pX(zz_pX& out_poly, EVP_MD_CTX *mdctx, const long& n_coeffs, const size_t& b_coeffs)
{    
    // NOTE: the current modulus (q_hat or q0) must already be set by the calling function

    int i;
    unsigned char* y_arr = new unsigned char[b_coeffs];   
           
    out_poly.SetLength(n_coeffs);

    for(i=0; i < n_coeffs; i++)
    {
        if (!EVP_DigestSqueeze(mdctx, y_arr, b_coeffs)) {
            printf("Message digest squeeze failed.\n");
            EVP_MD_CTX_free(mdctx);
            exit(1);
        }
        
        out_poly[i] = conv<zz_p>(ZZFromBytes(y_arr, b_coeffs)); 
    }
    
    out_poly.normalize();
    
    delete[] y_arr;   

    // return out_poly;
}


//==============================================================================
// Hash_v_zz_p - Generate a random vector (mod q_hat) using Custom Hash function
// 
// Inputs:
// - mdctx:     digest context in OpenSSL
// - n_elems:   number of elements of the random vector (i.e. 256+n+1)
// - b_num:     number of bytes for each random number (i.e. |q_hat|)
//
// Output:
// - out_vec:   vector of random numbers (modulo q_hat)
//==============================================================================
void Hash_v_zz_p(vec_zz_p& out_vec, EVP_MD_CTX *mdctx, const long& n_elems, const size_t& b_num)
{    
    // NOTE: the current modulus (q_hat or q0) must already be set by the calling function

    long      i;
    unsigned char* y_arr = new unsigned char[b_num];    

    out_vec.SetLength(n_elems);
    
    for(i=0; i < n_elems; i++)
    {
        if (!EVP_DigestSqueeze(mdctx, y_arr, b_num)) {
            printf("Message digest squeeze failed.\n");
            EVP_MD_CTX_free(mdctx);
            exit(1);
        }

        out_vec[i] = conv<zz_p>( ZZFromBytes(y_arr, b_num) );
    }
    
    delete[] y_arr;   

    // return out_vec;
}


//==============================================================================
// Hash_bits -  Generate a random vector of bits using the Custom Hash function
// 
// Inputs:
// - mdctx:     digest context in OpenSSL
// - n_elems:   number of elements of the random vector (i.e. d_hat)
//
// Output:
// - out_bits:  vector of random bits {0, 1}
//==============================================================================
void Hash_bits(vec_ZZ& out_bits, EVP_MD_CTX *mdctx, const long& n_elems)
{    
    long            i, j, k, n_bytes; 
    unsigned char   curr_byte;
    unsigned char*  y_arr;
    
    // Compute the minimum number of bytes needed to fill the vector  
    n_bytes = ceil(n_elems / 8.0);

    y_arr = new unsigned char[n_bytes];

    if (!EVP_DigestSqueeze(mdctx, y_arr, n_bytes)) {
        printf("Message digest squeeze failed.\n");
        EVP_MD_CTX_free(mdctx);
        exit(1);
    }
           
    out_bits.SetLength(n_elems);
    k = 0;

    for(i=0; i < n_bytes; i++)
    {
        curr_byte = y_arr[i];
        
        for(j=0; j < 8; j++)
        {
            if(k < n_elems)
            {
                out_bits[k] = conv<ZZ>( curr_byte & 1 );
                curr_byte = curr_byte >> 1;
            }
            k++;
        }
    } 
    
    delete[] y_arr;

    // return out_bits;
}


//==============================================================================
// Hash_ZZ_xi0 - Generate a random integer modulo (xi0+1) using Custom Hash function
// 
// Inputs:
// - mdctx:      digest context in OpenSSL
// - b_num:      number of bytes of the random integer
//
// Output:
// - out:        random integer modulo (xi0+1), i.e. from 0 to xi0
//==============================================================================

void Hash_ZZ_xi0(ZZ& out, EVP_MD_CTX *mdctx, const size_t& b_num)
{    
    unsigned char* y_arr = new unsigned char[b_num];    

    if (!EVP_DigestSqueeze(mdctx, y_arr, b_num)) {
        printf("Message digest squeeze failed.\n");
        EVP_MD_CTX_free(mdctx);
        exit(1);
    }

    out = (ZZFromBytes(y_arr, b_num)) % (xi0+1);
        
    delete[] y_arr;   

    // return out;
}



//==============================================================================
// Hcrs    -    H_crs, custom Hash function needed in BLNS for crs. 
//              It generates the pair of common random string (crs_ISIS, crs_Com).
// 
// Input:
// - inputStr:  string containing the input message (initial seed)
//
// Output:
// - crs:       structure with the pair (crs_ISIS, crs_Com)
//==============================================================================
void Hcrs(CRS2_t& crs, const string& inputStr)
{
    int                     i, j, n, m1, m2, n256;
    EVP_MD_CTX              *mdctx;
    size_t                  b_coeffs;
       
    mdctx = Hash_Init(inputStr);

    // Create the crs structure  
    crs.SetLength(2); 
    crs[0].SetLength(5); // crs_ISIS
    crs[1].SetLength(9); // crs_Com

    // ###########################  crs_ISIS  #########################################
    {
        // NOTE: elements of all matrices in crs_ISIS are mod q2_hat
        zz_pPush push(q2_hat); 
        // NOTE: backup current modulus q0, temporarily set to q2_hat (i.e., zz_p::init(q2_hat))    

        // Compute the minimum number of bytes to represent each coefficient
        b_coeffs = ceil(log2( conv<double>(q2_hat) ) / 8.0);    
        
        n    = n_ISIS;
        m1   = m1_ISIS;
        m2   = m2_ISIS;
        n256 = 256/d_hat;

        if ( (256 % d_hat) != 0)
        {
            cout << "ERROR! 256 must be divisible by d_hat" << endl;
            assert((256 % d_hat) == 0);
        }

        // Create the crs_ISIS structure, i.e. crs[0]
        // crs[0][0] = A_1;
        // crs[0][1] = A_2;
        // crs[0][2] = B_y;
        // crs[0][3] = B_g;
        // crs[0][4] = b;
            
        // Random generation of A_1 ∈ R^(n x m_1)_q_hat   
        //                      A_2 ∈ R^(n x m_2)_q_hat      
        crs[0][0].SetDims(n, m1);    
        crs[0][1].SetDims(n, m2);

        for(i=0; i<n; i++)
        {
            for(j=0; j<m1; j++)
            {
                Hash_zz_pX(crs[0][0][i][j], mdctx, d_hat, b_coeffs);
            }
            for(j=0; j<m2; j++)
            {
                Hash_zz_pX(crs[0][1][i][j], mdctx, d_hat, b_coeffs);
            }
        }       

        // Random generation of B_y ∈ R^(256/d_hat x m_2)_q_hat  
        crs[0][2].SetDims(n256, m2);

        for(i=0; i<n256; i++)
        {
            for(j=0; j<m2; j++)
            {
                Hash_zz_pX(crs[0][2][i][j], mdctx, d_hat, b_coeffs);
            }
        }

        // Random generation of B_g ∈ R^(tau0^ x m_2)_q_hat
        crs[0][3].SetDims(tau0, m2);

        for(i=0; i<tau0; i++)
        {
            for(j=0; j<m2; j++)
            {
                Hash_zz_pX(crs[0][3][i][j], mdctx, d_hat, b_coeffs);
            }
        }

        // Random generation of b ∈ R^(m_2)_q_hat
        crs[0][4].SetDims(1, m2);
        // NOTE: b is (1 x m_2) matrix, not a vector!

        for(i=0; i<m2; i++)
        {
            Hash_zz_pX(crs[0][4][0][i], mdctx, d_hat, b_coeffs);
        }        
    }

      
    // ###########################  crs_Com  #########################################    
    {
        // NOTE: elements of all matrices in crs_Com are mod q1_hat
        zz_pPush push(q1_hat); 
        // NOTE: backup current modulus q0, temporarily set to q1_hat (i.e., zz_p::init(q1_hat))

        // Compute the minimum number of bytes to represent each coefficient
        b_coeffs = ceil(log2( conv<double>(q1_hat) ) / 8.0);

        n    = n_Com;
        m1   = m1_Com;
        m2   = m2_Com;
        // n256 = 256/d_hat;

        // Create the crs_Com structure, i.e. crs[1]        
        // crs[1][0] = A_1;
        // crs[1][1] = A_2;
        // crs[1][2] = B_y;
        // crs[1][3] = B_g;
        // crs[1][4] = b;
        // crs[1][5] = A_bar_1;
        // crs[1][6] = A_bar_2;
        // crs[1][7] = B_bar_1;
        // crs[1][8] = B_bar_2;
            
        // Random generation of A_1 ∈ R^(n x m_1)_q_hat   
        //                      A_2 ∈ R^(n x m_2)_q_hat      
        crs[1][0].SetDims(n, m1);    
        crs[1][1].SetDims(n, m2); 
        
        for(i=0; i<n; i++)
        {
            for(j=0; j<m1; j++)
            {
                Hash_zz_pX(crs[1][0][i][j], mdctx, d_hat, b_coeffs);
            }
            for(j=0; j<m2; j++)
            {
                Hash_zz_pX(crs[1][1][i][j], mdctx, d_hat, b_coeffs);
            }
        }       

        // Random generation of B_y ∈ R^(256/d_hat x m_2)_q_hat  
        crs[1][2].SetDims(n256, m2);

        for(i=0; i<n256; i++)
        {
            for(j=0; j<m2; j++)
            {
                Hash_zz_pX(crs[1][2][i][j], mdctx, d_hat, b_coeffs);
            }
        }

        // Random generation of B_g ∈ R^(tau0^ x m_2)_q_hat
        crs[1][3].SetDims(tau0, m2);

        for(i=0; i<tau0; i++)
        {
            for(j=0; j<m2; j++)
            {
                Hash_zz_pX(crs[1][3][i][j], mdctx, d_hat, b_coeffs);
            }
        }

        // Random generation of b ∈ R^(m_2)_q_hat
        crs[1][4].SetDims(1, m2);
        // NOTE: b is (1 x m_2) matrix, not a vector!

        for(i=0; i<m2; i++)
        {
            // for(j=0; j<1; j++)
            {
                Hash_zz_pX(crs[1][4][0][i], mdctx, d_hat, b_coeffs);
            }
        }
        
        // Random generation of A_bar_1, B_bar_1 ∈ R^(m_1 x n_1)_q_hat
        crs[1][5].SetDims(m1, n_i);
        crs[1][7].SetDims(m1, n_i);       
        
        for(i=0; i<m1; i++)
        {
            for(j=0; j<n_i; j++)
            {
                Hash_zz_pX(crs[1][5][i][j], mdctx, d_hat, b_coeffs);
                Hash_zz_pX(crs[1][7][i][j], mdctx, d_hat, b_coeffs);
            }
        }

        // Random generation of A_bar_2, B_bar_2 ∈ R^(m_2 x n_2)_q_hat
        crs[1][6].SetDims(m2, n_i);
        crs[1][8].SetDims(m2, n_i);       
        
        for(i=0; i<m2; i++)
        {
            for(j=0; j<n_i; j++)
            {
                Hash_zz_pX(crs[1][6][i][j], mdctx, d_hat, b_coeffs);
                Hash_zz_pX(crs[1][8][i][j], mdctx, d_hat, b_coeffs);
            }
        }
    }

    EVP_MD_CTX_free(mdctx);
    
    // return crs;    
}


//==============================================================================
// HCom1   -    H_Com, custom Hash function needed in BLNS for commitment. 
//              It generates the 1st challenge used in the NIZK proof system.
// 
// Input:
// - inputStr:  string containing the input messages
//
// Output:
// - R_goth:    structure with the pair (R_goth_0, R_goth_1), 3D matrices of {0, 1}
//==============================================================================
void HCom1(R_GOTH_t& R_goth, const string& inputStr)
{
    int         i, j;
    EVP_MD_CTX *mdctx;

    const int  m1 = m1_Com;

    mdctx = Hash_Init(inputStr);  
    
    // Create the R_goth structure  
    R_goth.SetLength(2);   
    R_goth[0].SetDims(256, m1);
    R_goth[1].SetDims(256, m1); 

    // Random generation of R_goth_i ∈ {0, 1}^(256 x m_1 x d_hat)
    for(i=0; i<256; i++)
    {
        for(j=0; j<m1; j++)
        {
            Hash_bits(R_goth[0][i][j], mdctx, d_hat);
            Hash_bits(R_goth[1][i][j], mdctx, d_hat);
        }
    }

    EVP_MD_CTX_free(mdctx);

    // return R_goth;
}


//==============================================================================
// HCom2   -    H_Com, custom Hash function needed in BLNS for commitment. 
//              It generates the 2nd challenge used in the NIZK proof system.
// 
// Input:
// - inputStr:  string containing the input messages
//
// Output:
// - gamma:     matrix of integers modulo q1_hat
//==============================================================================
void HCom2(mat_zz_p& gamma, const string& inputStr)
{
    zz_pPush push(q1_hat); 
    // NOTE: backup current modulus q0, temporarily set to q1_hat (i.e., zz_p::init(q1_hat))

    int         i, n257;
    EVP_MD_CTX *mdctx; 
    
    mdctx = Hash_Init(inputStr); 

    // Compute the minimum number of bytes to represent each coefficient
    const size_t b_coeffs = ceil(log2( conv<double>(q1_hat) ) / 8.0);    
    
    n257 = 256 + d0 + 1; 
    // NOTE: gamma has 256+d+1 columns in Com, while 256+d+3 in ISIS   

    // Random generation of gamma ∈ Z^(tau0 x 256+d0+1)_q_hat
    gamma.SetDims(tau0, n257);

    for(i=0; i<tau0; i++)
    {
        Hash_v_zz_p(gamma[i], mdctx, n257, b_coeffs);
    }

    EVP_MD_CTX_free(mdctx);
    
    // return gamma;
}


//==============================================================================
// HCom3   -    H_Com, custom Hash function needed in BLNS for commitment. 
//              It generates the 3rd challenge used in the NIZK proof system.
// 
// Input:
// - inputStr:  string containing the input messages
//
// Output:
// - mu:        vector with tau0 polynomials with d_hat coefficients modulo q1_hat
//==============================================================================
void HCom3(vec_zz_pX& mu, const string& inputStr)
{
    zz_pPush push(q1_hat); 
    // NOTE: backup current modulus q0, temporarily set to q1_hat (i.e., zz_p::init(q1_hat))
    
    int         i;
    EVP_MD_CTX *mdctx;

    // Compute the minimum number of bytes to represent each coefficient
    const size_t b_coeffs = ceil(log2( conv<double>(q1_hat) ) / 8.0);   

    mdctx = Hash_Init(inputStr);  

    // Random generation of mu ∈ R^(tau0)_q_hat
    mu.SetLength(tau0);

    for(i=0; i<tau0; i++)
    {        
        Hash_zz_pX(mu[i], mdctx, d_hat, b_coeffs);
    }
    
    EVP_MD_CTX_free(mdctx);
        
    // return mu;
}


//==============================================================================
// HCom4   -    H_Com, custom Hash function needed in BLNS for commitment. 
//              It generates the 4th challenge used in the NIZK proof system.
// 
// Input:
// - inputStr:  string containing the input messages
//
// Output:
// - c:         polynomial with d_hat coefficients, c ∈ C ⊂ R^
// NOTE: c without modulo (q1_hat)
//==============================================================================
void HCom4(ZZX& c, const string& inputStr)
{
    zz_pPush push(q1_hat); 
    // NOTE: backup current modulus q0, temporarily set to q1_hat (i.e., zz_p::init(q1_hat))
       
    int         i;
    EVP_MD_CTX *mdctx;    
    ZZ          norm1_c, c_i;    
    ZZX         c_2k;
        
    // Compute the minimum number of bytes to represent each coefficient
    const size_t b_coeffs = ceil(log2(xi0+1) / 8.0);
    
    // Compute (nu0)^(2*k0)
    const ZZ    nu0_2k = power(conv<ZZ>(nu0), 2*k0);
      
    // Initialize the variable norm1_c = ||c^(2k)||_1
    norm1_c = 2*nu0_2k;
    
    mdctx = Hash_Init(inputStr); 

    c.SetLength(d_hat);
    c_2k.SetLength(d_hat*2*k0);

    // Loop to ensure that (2k)√(||c^(2k)||_1 ≤ nu0,  
    // i.e.  ||c^(2k)||_1 ≤ (nu0)^(2k)
    while(norm1_c > nu0_2k)
    {
        // Random generation of c ∈ R^_(xi0+1)
        Hash_ZZ_xi0(c_i, mdctx, b_coeffs);
        // NOTE: generate each coefficient c[i] ∈ [0, xi0], to ensure ||c||∞ ≤ ξ
        
        // c[0] = c_i;
        SetCoeff(c, 0, c_i);       
                
        for(i=1; i<(d_hat/2); i++)
        {
            Hash_ZZ_xi0(c_i, mdctx, b_coeffs);
            
            // c[i] = c_i;
            SetCoeff(c, i, c_i);

            // c[d_hat-i] = -c[i];
            SetCoeff(c, (d_hat-i), -c_i);
            // NOTE: this ensures that σ(c) = c
        }
        c.normalize();
        
        // NOTE: avoid (rare) cases with c == 0
        if (IsZero(conv<zz_pX>(c)))
        {
            continue;
        }
        
        // c_2k = power(c, (2*k0));
        c_2k = c;

        for(i=0; i<(2*k0 - 1); i++)
        {
            // c_2k *= c;
            // c_2k = (c_2k * c) % phi_hat; 
            c_2k = ModPhi_hat(c_2k * c);
        }

        // Compute ||c^(2k)||_1
        norm1_c = 0;

        for(i=0; i<=deg(c_2k); i++)
        {
            // norm1_c = norm1_c + c_2k[i];
            norm1_c += abs(coeff(c_2k, i)); 
        }
    }
    
    EVP_MD_CTX_free(mdctx);
         
    // return c;
}


//==============================================================================
// HISIS1   -   H_ISIS, custom Hash function needed in BLNS for ISIS. 
//              It generates the 1st challenge used in the NIZK proof system.
// 
// Input:
// - inputStr:  string containing the input messages
//
// Output:
// - R_goth:    structure with the pair (R_goth_0, R_goth_1), 3D matrices of {0, 1}
//==============================================================================
// NOTE: HISIS1 is identical to HCom1, apart m1
void HISIS1(R_GOTH_t& R_goth, const string& inputStr)
{
    int                 i, j;
    EVP_MD_CTX *mdctx;

    const int  m1 = m1_ISIS;

    mdctx = Hash_Init(inputStr);  
    
    // Create the R_goth structure  
    R_goth.SetLength(2);   
    R_goth[0].SetDims(256, m1);
    R_goth[1].SetDims(256, m1); 

    // Random generation of R_goth_i ∈ {0, 1}^(256 x m_1 x d_hat) 
    for(i=0; i<256; i++)
    {
        for(j=0; j<m1; j++)
        {
            Hash_bits(R_goth[0][i][j], mdctx, d_hat);
            Hash_bits(R_goth[1][i][j], mdctx, d_hat);            
        }
    }

    EVP_MD_CTX_free(mdctx);

    // return R_goth;
}


//==============================================================================
// HISIS2   -   H_ISIS, custom Hash function needed in BLNS for ISIS. 
//              It generates the 2nd challenge used in the NIZK proof system.
// 
// Input:
// - inputStr:  string containing the input messages
//
// Output:
// - gamma:     matrix of integers modulo q2_hat
//==============================================================================
void HISIS2(mat_zz_p& gamma, const string& inputStr)
{
    zz_pPush push(q2_hat); 
    // NOTE: backup current modulus q0, temporarily set to q2_hat (i.e., zz_p::init(q2_hat))

    int         i, n259;       
    EVP_MD_CTX *mdctx; 
    
    mdctx = Hash_Init(inputStr); 

    // Compute the minimum number of bytes to represent each coefficient
    const size_t b_coeffs = ceil(log2( conv<double>(q2_hat) ) / 8.0);    
    
    n259 = 256 + d0 + 3;
    // NOTE: gamma has 256+d+3 columns in ISIS, while 256+d+1 in Com 

    // Random generation of gamma ∈ R^(tau0 x 256+d+3)_q_hat
    gamma.SetDims(tau0, n259);

    for(i=0; i<tau0; i++)
    {
        Hash_v_zz_p(gamma[i], mdctx, n259, b_coeffs);
    }
    
    EVP_MD_CTX_free(mdctx);

    // return gamma;
}


//==============================================================================
// HISIS3   -   H_ISIS, custom Hash function needed in BLNS for ISIS. 
//              It generates the 3rd challenge used in the NIZK proof system.
// 
// Input:
// - inputStr:  string containing the input messages
//
// Output:
// - mu:        vector with tau0 polynomials with d_hat coefficients modulo q2_hat
//==============================================================================
void HISIS3(vec_zz_pX& mu, const string& inputStr)
// NOTE: HISIS3 is identical to HCom3, apart the modulo  
{     
    zz_pPush push(q2_hat); 
    // NOTE: backup current modulus q0, temporarily set to q2_hat (i.e., zz_p::init(q2_hat))
    
    int         i;
    EVP_MD_CTX *mdctx;

    // Compute the minimum number of bytes to represent each coefficient
    const size_t b_coeffs = ceil(log2( conv<double>(q2_hat) ) / 8.0);   

    mdctx = Hash_Init(inputStr);  

    // Random generation of mu ∈ R^(tau0)_q_hat
    mu.SetLength(tau0);

    for(i=0; i<tau0; i++)
    {        
        Hash_zz_pX(mu[i], mdctx, d_hat, b_coeffs);
    }
    
    EVP_MD_CTX_free(mdctx);
        
    // return mu;
}


//==============================================================================
// HISIS4   -   H_ISIS, custom Hash function needed in BLNS for ISIS. 
//              It generates the 4th challenge used in the NIZK proof system.
// 
// Input:
// - inputStr:  string containing the input messages
//
// Output:
// - c:         polynomial with d_hat coefficients 
// NOTE: c without modulo (q2_hat)
//==============================================================================
void HISIS4(ZZX& c, const string& inputStr)
// NOTE: HISIS4 is identical to HCom4, apart the modulo
{
    zz_pPush push(q2_hat); 
    // NOTE: backup current modulus q0, temporarily set to q2_hat (i.e., zz_p::init(q2_hat))
       
    int         i;
    EVP_MD_CTX *mdctx;    
    ZZ          norm1_c, c_i;    
    ZZX         c_2k;
        
    // Compute the minimum number of bytes to represent each coefficient
    const size_t b_coeffs = ceil(log2(xi0+1) / 8.0);
    
    // Compute (nu0)^(2*k0)
    const ZZ    nu0_2k = power(conv<ZZ>(nu0), 2*k0);
      
    // Initialize the variable norm1_c = ||c^(2k)||_1
    norm1_c = 2*nu0_2k;
    
    mdctx = Hash_Init(inputStr); 

    c.SetLength(d_hat);
    c_2k.SetLength(d_hat*2*k0);

    // Loop to ensure that (2k)√(||c^(2k)||_1 ≤ nu0,  
    // i.e.  ||c^(2k)||_1 ≤ (nu0)^(2k)
    while(norm1_c > nu0_2k)
    {
        // Random generation of c ∈ R^_(xi0+1)
        Hash_ZZ_xi0(c_i, mdctx, b_coeffs);
        // NOTE: generate each coefficient c[i] ∈ [0, xi0], to ensure ||c||∞ ≤ ξ
        
        // c[0] = c_i;
        SetCoeff(c, 0, c_i);       
                
        for(i=1; i<(d_hat/2); i++)
        {
            Hash_ZZ_xi0(c_i, mdctx, b_coeffs);
            
            // c[i] = c_i;
            SetCoeff(c, i, c_i);

            // c[d_hat-i] = -c[i];
            SetCoeff(c, (d_hat-i), -c_i);
            // NOTE: this ensures that σ(c) = c
        }
        c.normalize();
        
        // NOTE: avoid (rare) cases with c == 0
        if (IsZero(conv<zz_pX>(c)))
        {
            continue;
        }
        
        // c_2k = power(c, (2*k0));
        c_2k = c;

        for(i=0; i<(2*k0 - 1); i++)
        {
            // c_2k *= c;
            // c_2k = (c_2k * c) % phi_hat; 
            c_2k = ModPhi_hat(c_2k * c);
        }

        // Compute ||c^(2k)||_1
        norm1_c = 0;

        for(i=0; i<=deg(c_2k); i++)
        {
            // norm1_c = norm1_c + c_2k[i];
            norm1_c += abs(coeff(c_2k, i)); 
        }
    }
    
    EVP_MD_CTX_free(mdctx);
         
    // return c;
}


//==============================================================================
// HM     -     H_M, custom Hash function needed in BLNS for hashing attributes. 
//              It hashes an attribute a_i into a vector of length h0, 
//              with coefficients in the range (−ψ, ψ).
// 
// Input:
// - a_i:       attribute, string of bits of arbitrary length a_i ∈ {0, 1}∗
//
// Output:
// - m_i:       vector with h0 coefficients in the range (−psi0, psi0)
//==============================================================================
void HM(vec_ZZ& m_i, const string& a_i)
{
    long        k, range;
    EVP_MD_CTX *mdctx;
    vec_zz_p    tmp;
    
    mdctx = Hash_Init(a_i);

    // Compute the numerical range of each coefficient
    range = 2*psi0 + 1;

    zz_pPush push(range);
    // NOTE: backup current modulus q0, temporarily set to range (i.e., zz_p::init(range))
    
    // Compute the minimum number of bytes to represent each coefficient
    const size_t b_coeffs = ceil(log2(range) / 8.0);

    // Random generation of m_i (modulo range)
    Hash_v_zz_p(tmp, mdctx, h0, b_coeffs);
    m_i = conv<vec_ZZ>( tmp );

    for(k=0; k<h0; k++)
    {
        m_i[k] = m_i[k] - psi0;
        // NOTE: now each coefficient is in the range (−psi0, psi0)
    }
    
    EVP_MD_CTX_free(mdctx);

    // return m_i;
}
