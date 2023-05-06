Core code for the paper [Exploring Decryption Failures of BIKE: New Class of Weak Keys and Key Recovery Attacks](https://submit.iacr.org/crypto2023/paper/128?cap=hcav128AwcHVhHEVEHnYjcwLymGdRsc)

This program constructs different kinds of weak keys and ciphertexts with gathering properties and simulates the DFR of BIKE. Code for calculating frequency of keys and proof in Appendix is trivial according to our lemmas and theories in the paper.

Installation:
=====
Use make to run the program, g++ is necessary and openmp is optional.

Usage:
=====
The parameters are fixed before compilation.

Users can change the parameter m,out_of (the parameter $\epsilon$ in paper) in main.cpp. 

In main.cpp, if only_key=1, key is special and error is rand; if tot_rand=1, key is random and error is special; if rand_err=1, error is random. 

If eps=1, we output "good" which means the vector is not in overlapping area, while "bad" means the vector is in overlapping area.

The outputs of poly and error are the indexes of nonzero elements, we have tried to test those vectors in [standard BIKE](https://github.com/awslabs/bike-kem) and could obtain the same decryption failure.

To make experiments under small parameters, users can also change the r in util.h.
