# Cryptography
## Implementation of RSA algorithm

## Build
Use `make` command to build this project.

## Run
Tool can be run with 3 different options. All numbers are in hexadecimal format except for size of modulus.

### Generate
Generate 2 prime numbers (used in algorithm), modulus, public exponent (key) and private exponent (key).
```bash
$ ./kry -g size
```
Parameter of `-g` argument is size of modulus in bits.

### Encrypt
Encrypt message with RSA algorithm.
```bash
$ ./kry -e exp mod message
```
Parameters of `-e` are public exponent, modulus and message to encrypt (in form of number).

### Decrypt
Decrypt message with RSA algorithm.
```bash
$ ./kry -d exp mod enc_message
```
Parameters of `-d` are private exponent, modulus and message to decrypt (in form of number).