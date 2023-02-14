# A weak key in BIKE
use make to run the program, g++ is necessary and openmp is optional \
you can change the parameter m,out_of (represent epsilon in paper) and only_key (1 represents only weak key and 0 represents weak key plus weak ciphertexts) in main.cpp \
If eps=1, we the output "good" means the vector is not in overlapping area, "bad" means the vector is in overlapping area \
The output of poly and error is the indexes of nonzero elements, we have tried to move this input ot standard BIKE and could obtain the same decryption failure 
