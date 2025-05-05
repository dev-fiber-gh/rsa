#include <openssl/bn.h>
#include <openssl/err.h>
#include <openssl/rand.h>

#include <time.h>
#include <memory.h>

#include <memory>
#include <vector>
#include <string>
#include <random>
#include <cstdint>
#include <iostream>

#include "primes.h"



std::random_device g_Rd;
std::mt19937 g_Gen(g_Rd());
std::uniform_int_distribution<uint64_t> g_GenRandom(1, 0xffffffffffffffff);



#define MAKESMARTPTR(Type, Alloc, Delloc)\
class Smart##Type {\
    Type* num = nullptr;\
public:\
    Smart##Type()  { num = Alloc(); }\
    ~Smart##Type() { if(num) { Delloc(num); num = nullptr; } }\
    operator const Type*() { return  num; }\
    operator Type*()       { return  num; }\
    operator Type**()      { return &num; }\
};


MAKESMARTPTR(BIGNUM, BN_new, BN_free)
MAKESMARTPTR(BN_CTX, BN_CTX_new, BN_CTX_free)


#define BN_SIG_COMPARE(a, o, b) BN_cmp(a,b) o 0
#define BN_ABS_COMPARE(a, o, b) BN_ucmp(a,b) o 0


void setLog(char* primeLog, const uint32_t primeLogSize, const uint32_t pos, const uint16_t prime) {
    for (uint32_t i = pos; i < primeLogSize ; i += prime) primeLog[i] = 0;
}
int getLog(char* primeLog, const uint32_t primeLogSize, const uint32_t pos) {
    if (pos >= primeLogSize) return true;
    return primeLog[pos];
}

bool SimpleTest(const BIGNUM* v, char* primeLog, const uint32_t primeLogSize, const uint32_t pos) {  
    for (uint32_t i = 0; i < g_PrimesNumber; i++) {
        if(g_PrimeTest[i]) {
            if (!BN_mod_word(v, (uint32_t)g_Primes[i])) {
                g_PrimeTest[i] = 0;
                setLog(primeLog, primeLogSize, pos, g_Primes[i]);
                return false;
            } 
        }
    }
    
    return true;
}
bool Fermat(BIGNUM* v, BIGNUM* vm1, BIGNUM* two, BIGNUM* res, BN_CTX* ctx) {
    BN_copy(vm1, v);
    BN_sub_word(vm1, 1);
    BN_mod_exp(res, two, vm1, v, ctx);

    if (BN_is_one(res)) return true;
    else return false;
}
bool MillerRabin(BIGNUM* v, BIGNUM* vm1, BIGNUM* a,    BIGNUM* b,
                 BIGNUM* r, BIGNUM* two, BIGNUM* p256, BIGNUM* res, 
                 BN_CTX* ctx, const uint32_t trials) {
    BN_copy(vm1, v);
    BN_sub_word(vm1, 1);
    BN_copy(a, vm1);
    BN_copy(b, vm1);

    int k = 0;
    while (!BN_is_zero(b)) {
        if (BN_is_bit_set(b, 0)) BN_lshift1(b, b);
        else break;
        k++;
    }

    BN_lshift(a, a, k);

    uint64_t str[8];

    for (uint32_t x = 0; x < trials; x++) {
        std::cout << x << " ";

        do {
            for (int i = 0; i < 8; i++) str[i] = g_GenRandom(g_Gen);
            BN_bin2bn((unsigned char*)str, 64, r);
        } while ((BN_SIG_COMPARE(r, >, vm1)) || (BN_SIG_COMPARE(r, <, p256)));    

        BN_copy(b, v);
        BN_sub_word(b, 4);
        BN_div(NULL, b, r, b, ctx);
        BN_add_word(b, 2);
        BN_mod_exp(res, b, a, v, ctx);
        
        if (BN_is_one(res) || (BN_SIG_COMPARE(res, ==, vm1))) continue;

        for (int i = 1; i < k; i++) {
            BN_mod_exp(res, res, two, v, ctx);
            if (BN_is_one(res)) return false;
            if (BN_SIG_COMPARE(res, ==, vm1)) goto mainloop;
        }

        return false;

        mainloop:;
    }

    return true;
}
bool isProbablyPrime(BIGNUM* v, BIGNUM* vm1, BIGNUM* a,    BIGNUM* b,
                     BIGNUM* r, BIGNUM* two, BIGNUM* p256, BIGNUM* res,
                     BN_CTX* ctx, char* primeLog, const uint32_t primeLogSize, const uint32_t i) {
    std::cout << "r.";
    if (!getLog(primeLog, primeLogSize, i)) return false;

    std::cout << "s.";
    if (!SimpleTest(v, primeLog, primeLogSize, i)) return false;

    std::cout << "f.";
    if (!Fermat(v, vm1, two, res, ctx)) return false;

    std::cout << "mr.";
    if (!MillerRabin(v, vm1, a, b, r, two, p256, res, ctx, 100)) return false;

    return true;
}

void correct_e(BIGNUM* phi, BIGNUM* e, BIGNUM* res, BN_CTX* ctx) {
    for(;;) {
        BN_gcd(res, phi, e, ctx);
        if (BN_is_one(res)) return;
        else BN_add_word(e, 2);
    }
}

int main() {
    clock_t start = clock();

    const uint32_t BitPerByte = 8;
    const uint32_t BitLength = 6144;
    const uint32_t NumberOfRoots = 2;
    const uint32_t primeLogSize = 50000;
    
    
    
    SmartBIGNUM two;
    SmartBIGNUM p256;
    BN_set_bit(two, 1);
    BN_set_bit(p256, 256);


    char* numberStr;
    char  primeLog[primeLogSize];
    
    

    SmartBIGNUM e;
    SmartBIGNUM magn;
    SmartBIGNUM test_pq;  
    BN_dec2bn((BIGNUM**)magn, std::to_string(BitLength).c_str()); 
    BN_dec2bn((BIGNUM**)e, "167772161"); //Prime with 3 "one" bits
    BN_set_bit(test_pq, (BitLength / NumberOfRoots) - 8); 


    std::vector<unsigned char> rnd;


    SmartBIGNUM p;    
    rnd.resize(BitLength / (NumberOfRoots * BitPerByte));
    do {
        for (uint32_t i = 0; i < rnd.size(); i++) {
            if(RAND_bytes(&rnd[0], (int)rnd.size()) != 1) { std::cerr << "Random generation failed.\n"; return -1; }
        }

        BN_bin2bn((const unsigned char*)&rnd[0], (int)rnd.size(), p);
        rnd.clear();
    } while (BN_SIG_COMPARE(p, <=, test_pq)); 
        
    
    SmartBN_CTX ctx;
    SmartBIGNUM p_magn;
    BN_div(p_magn, NULL, p, magn, ctx);
    
    
    SmartBIGNUM q; 
    SmartBIGNUM pmq;
    rnd.resize(BitLength / (NumberOfRoots * BitPerByte));
    do {
        for (uint32_t i = 0; i < rnd.size(); i++) {
            if(RAND_bytes(&rnd[0], (int)rnd.size()) != 1) { std::cerr << "Random generation failed.\n"; return -1; }
        }
        
        BN_bin2bn((const unsigned char*)&rnd[0], (int)rnd.size(), q);
        rnd.clear();

        BN_sub(pmq, p, q);
        if (BN_ABS_COMPARE(pmq, >, p_magn)) break;
    } while (BN_SIG_COMPARE(q, <=, test_pq));

    
    //set value to odd
    BN_set_bit(p, 0);
    BN_set_bit(q, 0);

    
    
    SmartBIGNUM a;
    SmartBIGNUM b;
    SmartBIGNUM r;
    SmartBIGNUM vm1;
    SmartBIGNUM res;
    
    memset(primeLog, 1, primeLogSize);
    memset(g_PrimeTest, 1, sizeof(g_PrimeTest));
    for (uint32_t i = 0;; i += 2) {
        std::cout << i;
        if (isProbablyPrime(p, vm1, a, b, r, two, p256, res, ctx, primeLog, primeLogSize, i)) {
            std::cout << "\n";
            std::cout << "\n";

            numberStr = BN_bn2dec(p);

            std::cout << numberStr;
            OPENSSL_free(numberStr);

            std::cout << " is a pseudoprime" << "\n";

            break;
        }
        std::cout << " ";
        BN_add_word(p, 2); //odd + even(2) = odd
    }
    std::cout << "\n";

    
    
    memset(primeLog, 1, primeLogSize);
    memset(g_PrimeTest, 1, sizeof(g_PrimeTest)); 
    for (uint32_t i = 0;; i += 2) {
        std::cout << i;
        if (isProbablyPrime(q, vm1, a, b, r, two, p256, res, ctx, primeLog, primeLogSize, i)) {
            std::cout << "\n";
            std::cout << "\n";

            numberStr = BN_bn2dec(q);

            std::cout << numberStr;
            OPENSSL_free(numberStr);

            std::cout << " is a pseudoprime" << "\n";

            break;
        }
        std::cout << " ";
        BN_add_word(q, 2); //odd + even(2) = odd
    }
    std::cout << "\n";


    SmartBIGNUM n;
    BN_mul(n, p, q, ctx);    
    std::cout << "n" << "\n";    
    numberStr = BN_bn2dec(n);
    std::cout << numberStr << "\n" << "\n";
    OPENSSL_free(numberStr);


    ///////////////////////////////////////////////////////////////////
    //Some standards require d < lambda(n), "lambda(n) = lcm(p-1, q-1)"
    SmartBIGNUM pm1; 
    SmartBIGNUM qm1; 
    SmartBIGNUM phi;

    //phi = (p-1)*(q-1)
    BN_copy(pm1, p);
    BN_sub_word(pm1, 1);
    BN_copy(qm1, q);
    BN_sub_word(qm1, 1);
    BN_mul(phi, pm1, qm1, ctx);
    //////////////////////////////////////////////////////////////////

    
    correct_e(phi, e, res, ctx);

    
    SmartBIGNUM d;
    SmartBIGNUM n_magn;
    BN_div(n_magn, NULL, n, magn, ctx); 
    while (true) {
        if (BN_mod_inverse(d, e, phi, ctx) == NULL) { correct_e(phi, e, res, ctx); continue; }
        else if (BN_SIG_COMPARE(d, <=, n_magn))     { correct_e(phi, e, res, ctx); continue; }
        else break;
    }


    std::cout << "e" << "\n";
    numberStr = BN_bn2dec(e);
    std::cout << numberStr << "\n" << "\n";
    OPENSSL_free(numberStr);

    std::cout << "d" << "\n";
    numberStr = BN_bn2dec(d);
    std::cout << numberStr << "\n" << "\n";
    OPENSSL_free(numberStr);



    //Message part
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    std::string message = "Tell me, O Muse, of that ingenious hero who travelled far and wide after he had sacked the famous town of Troy. Many cities did he visit, and many were the nations with whose manners and customs he was acquainted; moreover he suffered much by sea while trying to save his own life and bring his men safely home; but do what he might he could not save his men, for they perished through their own sheer folly in eating the cattle of the Sun-god Hyperion; so the god prevented them from ever reaching home. Tell me, too, about all these things, oh daughter of Jove, from whatsoever source you may know them.";

    SmartBIGNUM m;

    BN_bin2bn((const unsigned char*)&message[0], (int)message.size(), m);
    if(BN_SIG_COMPARE(m, >=, n)) { std::cerr << "Message is too long.\n"; return -1; }
    message.clear();
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////


    //Encrypting
    BN_mod_exp(res, m, e, n, ctx);
    BN_clear(m);

    
    unsigned char randomSizeNumber;
    if(RAND_bytes(&randomSizeNumber, 1) != 1) { std::cerr << "Random generation failed.\n"; return -1; }
    randomSizeNumber = (randomSizeNumber & 31) + 16;
    
    rnd.resize(randomSizeNumber);    
    if(RAND_bytes(&rnd[0], (int)rnd.size()) != 1) { std::cerr << "Random generation failed.\n"; return -1; }    
    BN_bin2bn((const unsigned char*)&rnd[0], (int)rnd.size(), r);
    
    
    //RSA blinding and decrypting
    //////////////////////////////////////////////////////////////////////////////////////////
    
    ///////////blinding///////////////////////
    SmartBIGNUM mask;
    //res = res * r^e mod n    
    BN_mod_exp(mask, r, e, n, ctx);
    BN_mul(res, res, mask, ctx);
    //////////////////////////////////////////


    
    #define CRT

#ifdef CRT
    ///////////////////With CRT///////////////////////
    SmartBIGNUM dq;
    SmartBIGNUM dp;
    SmartBIGNUM qInv;
    SmartBIGNUM m1;
    SmartBIGNUM m2;
    SmartBIGNUM t;

    //dp = e^-1 mod(p-1) = d mod (p-1)
    BN_div(NULL, dp, d, pm1, ctx);

    //dq = e^-1 mod(q-1) = d mod (q-1)
    BN_div(NULL, dq, d, qm1, ctx);
        
    //qInv = q^-1 mod p
    BN_mod_inverse(qInv, q, p, ctx);

    //m1 = res ^dp mod p
    BN_mod_exp(m1, res, dp, p, ctx);

    //m1 = res ^dq mod q
    BN_mod_exp(m2, res, dq, q, ctx);

    //t = m1 - m2
    BN_sub(t, m1, m2);

    //in case if t is negative
    while (BN_is_negative(t))
        BN_add(t, t, p);

    // m = m2 + (qInv * t mod p) * q
    BN_mod_mul(m, t, qInv, p, ctx);
    BN_mul(t, m, q, ctx);
    BN_add(m, t, m2);

#else
    ///////////////Regular(Without CRT)////////////////
    BN_mod_exp(m, res, d, n, ctx); //m = res^d mod n
    ///////////////////////////////////////////////////
#endif


    ////////////////////////////////////////////////////////
    //that block calculating m = m * r^-1 mod n
    if (BN_mod_inverse(mask, r, n, ctx) == NULL) { std::cerr << "Unpredicted Error.\n"; return -1; }
    BN_mod_mul(m, m, mask, n, ctx);
    ///////////////////////////////////////////////////////

    message.resize(BitLength / BitPerByte);
    
    const int len = BN_bn2bin(m, (unsigned char*)&message[0]);
    message.resize((unsigned)len);
    std::cout << message << "\n";
    std::cout << double(clock() - start) / CLOCKS_PER_SEC << " s\n";

    return 0;
}


