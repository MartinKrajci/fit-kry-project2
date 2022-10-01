#include <gmp.h>
#include <stdio.h>
#include <unistd.h>
#include <string>
#include <time.h>
#include <tgmath.h>

#include <iostream>

#define ITERATIONS 40

using namespace std;

/*
* Implementation of Miler-Rabin algorithm used as a test for prime numbers
*/
bool mrTest(mpz_t d, mpz_t prime, gmp_randstate_t rstate, int primeBits)
{

    mpz_t a, prime_sub_4, x, prime_sub_1;
    mpz_inits(a, prime_sub_4, x, prime_sub_1, NULL);

    mpz_sub_ui(prime_sub_4, prime, 4);

    mpz_urandomb(a, rstate, primeBits);
    mpz_mod(a, a, prime_sub_4);
    mpz_add_ui(a, a, 2);

    mpz_powm(x, a, d, prime);
    mpz_sub_ui(prime_sub_1, prime, 1);

    if(mpz_cmp_ui(x, 1) == 0)
    {
        return true;
    }

    if(mpz_cmp(x, prime_sub_1) == 0)
    {
        return true;
    }

    while(mpz_cmp(d, prime_sub_1) == 0)
    {
        mpz_mul(x, x, x);
        mpz_mod(x, x, prime);

        if(mpz_cmp_ui(x, 1) == 0)
        {
            return false;
        }
        if(mpz_cmp(x, prime_sub_1) == 0)
        {
            return true;
        }
    }

    mpz_clears(a, prime_sub_4, x, prime_sub_1, NULL);

    return false;
}

/*
* Generate prime number of given size
*/
bool generatePrime(mpz_t prime, gmp_randstate_t rstate, int primeBits) 
{
    mpz_t prime_mod_res, d, d_mod_res, highest_bit;
    mpz_inits(prime_mod_res, d, d_mod_res, highest_bit, NULL);

    mpz_urandomb(prime, rstate, primeBits);
    mpz_setbit(prime, primeBits-1);

    mpz_mod_ui(prime_mod_res, prime, 2);
    // prime <= 1
    if(!(mpz_cmp_ui(prime, 1) > 0))
    {

        return false;
    }
    // prime is even
    if(mpz_cmp_ui(prime_mod_res, 0) == 0)
    {

        return false;
    }
    if(!(mpz_cmp_ui(prime, 3) > 0))
    {
        return true;
    }

    mpz_sub_ui(d, prime, 1);
    mpz_mod_ui(d_mod_res, d, 2);
    while (mpz_cmp_ui(d_mod_res, 0) == 0)
    {
        mpz_div_ui(d, d, 2);
        mpz_mod_ui(d_mod_res, d, 2);
    }
    

    for(int i = 0; i < ITERATIONS; i++)
    {
        if (!mrTest(d, prime, rstate, primeBits))
        {
            return false;
        }
    }
    return true;
    
}

/*
* Extended euclidean algorithm which helps to find secret exponent
*/
bool modifiedGcd(mpz_t e, mpz_t phi, mpz_t x, mpz_t y)
{
    bool isOne;

    mpz_t phi_temp, e_temp, x_tmp, y_tmp, inv_tmp;
    mpz_inits(phi_temp, e_temp, x_tmp, y_tmp, inv_tmp, NULL);

    if (mpz_cmp_ui(e, 0) == 0)
    {
        mpz_set_ui(x, 0);
        mpz_set_ui(y, 1);
        mpz_clears(phi_temp, e_temp, x_tmp, y_tmp, inv_tmp, NULL);
        return mpz_cmp_ui(phi, 1) == 0;
    }
    mpz_mod(phi_temp, phi, e);
    mpz_set(e_temp, e);
    isOne = modifiedGcd(phi_temp, e_temp, x_tmp ,y_tmp);

    mpz_div(inv_tmp, phi, e);
    mpz_mul(inv_tmp, inv_tmp, x_tmp);
    mpz_sub(x, y_tmp, inv_tmp);
    mpz_set(y, x_tmp);
    return isOne;
    

    mpz_clears(phi_temp, e_temp, x_tmp, y_tmp, inv_tmp, NULL);
}

/*
* Generate exponents (RSA keys)
*/
void generateKeys(mpz_t phi, gmp_randstate_t rstate, mpz_t e)
{
    bool isOne  = false;

    mpz_t e_tmp, phi_temp, x, y, d;
    mpz_inits(e, e_tmp, phi_temp, x, y, d, NULL);

    while (!isOne)
    {
        mpz_urandomm(e, rstate, phi);
        mpz_set(e_tmp, e);
        mpz_set(phi_temp, phi);
        isOne = modifiedGcd(e_tmp, phi_temp, x, y);
    }

    mpz_mod(d, x, phi);

    gmp_printf("%#Zx %#Zx\n", e, d);
    
    mpz_clears(e_tmp, phi_temp, x, y, d, NULL);
}

/*
* Generate modulus, exponents and factors
*/
char* generate(int bits)
{
    int p_bits, q_bits;
    float half_bits;

    mpz_t p, q, n, bit_test, bit_test_res, two, p_sub_one, q_sub_one, phi, e;
    mpz_inits(p, q, n, bit_test, bit_test_res, two, p_sub_one, q_sub_one, phi, e, NULL);

    gmp_randstate_t rstate;
    gmp_randinit_mt(rstate);
    srand(time(0));
    gmp_randseed_ui(rstate, rand());

    mpz_set_ui(two, 2);
    mpz_root(bit_test, two, bits - 1);

    half_bits = bits/2.0;
    if((static_cast<int>(half_bits)) == half_bits)
    {
        p_bits = half_bits;
        q_bits = half_bits;
    }
    else
    {
        p_bits = static_cast<int>(half_bits);
        q_bits = (static_cast<int>(half_bits)) + 1;
    }

    do
    {
        if(!generatePrime(p, rstate, p_bits))
        {
            continue;
        }
        if(!generatePrime(q, rstate, q_bits))
        {
            continue;
        }
        mpz_mul(n, p, q);
        mpz_and(bit_test_res, n, bit_test);
    }
    while(mpz_cmp_ui(bit_test_res, 0) == 0);
    gmp_printf("%#Zx %#Zx %#Zx ", p, q, n);
    mpz_sub_ui(p_sub_one, p, 1);
    mpz_sub_ui(q_sub_one, q, 1);

    mpz_mul(phi, p_sub_one, q_sub_one);
    generateKeys(phi, rstate, e);

    gmp_randclear(rstate);  

    mpz_clears(p, q, n, bit_test, bit_test_res, two, p_sub_one, q_sub_one, phi, e, NULL);
    
    return NULL;
}

/*
* Encode given message
*/
void encode(mpz_t pubExp, mpz_t modulus, mpz_t message)
{
    mpz_t cyphertext;
    mpz_init(cyphertext);
    mpz_powm(cyphertext, message, pubExp, modulus);
    gmp_printf("%#Zx\n", cyphertext);
    mpz_clear(cyphertext);
}

/*
* Decode given cyphertext
*/
void decode(mpz_t privExp, mpz_t modulus, mpz_t cyphertext)
{
    mpz_t message;
    mpz_init(message);
    mpz_powm_sec(message, cyphertext, privExp, modulus);
    gmp_printf("%#Zx\n", message);
    mpz_clear(message);
}

int main(int argc, char *argv[])
{
    int opt;
    while((opt = getopt(argc, argv, ":g:e:d:b:")) != -1) 
    { 
        switch(opt) 
        { 
            case 'g':
                {
                    int bits = stoi(optarg);
                    generate(bits);
                    break;
                }
            case 'e':
                {
                    mpz_t pubExp, modulus, message;
                    mpz_inits(pubExp, modulus, message, NULL);
                    string pubExpStr = string(optarg).erase(0,2);
                    string modulusStr = string(argv[3]).erase(0,2);
                    string messageStr = string(argv[4]).erase(0,2);
                    mpz_set_str(pubExp, pubExpStr.c_str(), 16);
                    mpz_set_str(modulus, modulusStr.c_str(), 16);
                    mpz_set_str(message, messageStr.c_str(), 16);
                    encode(pubExp, modulus, message);
                    mpz_clears(pubExp, modulus, message, NULL);
                    break;
                }
            case 'd':
                {
                    mpz_t privExp, modulus, cyphertext;
                    mpz_inits(privExp, modulus, cyphertext, NULL);
                    string privExpStr = string(optarg).erase(0,2);
                    string modulusStr = string(argv[3]).erase(0,2);
                    string cyphertextStr = string(argv[4]).erase(0,2);
                    mpz_set_str(privExp, privExpStr.c_str(), 16);
                    mpz_set_str(modulus, modulusStr.c_str(), 16);
                    mpz_set_str(cyphertext, cyphertextStr.c_str(), 16);
                    decode(privExp, modulus, cyphertext);
                    mpz_clears(privExp, modulus, cyphertext, NULL);
                    break;
                }
            case 'b':
                return 0; 
            case ':': 
                printf("Missing parameter of argument\n");
                return 1; 
            case '?': 
                printf("Unknown option: %c\n", optopt);
                return 1; 
        } 
    } 

    return 0;
}