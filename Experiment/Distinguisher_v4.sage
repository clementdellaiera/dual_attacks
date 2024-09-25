import numpy as np
import cupy as cp
import time 

start = time.time()
load('Parameters.sage')
modulus , short_vector_list        = load('Dual_short_vectors/dual_short_vectors')
target_list_LWE , target_list_unif = load('Target_samples/target_samples')

print("Loading time : ", time.time() - start , ' sec')

W_len , T_len = len(short_vector_list) , min(len(target_list_LWE) , len(target_list_unif))
modulus = int(modulus)

print("Convertion en np.array...")
start = time.time()
int_W_transpose = np.array([np.array([int(short_vector_list[i][j]) for j in range(rank)]) for i in range(len(short_vector_list))]).transpose()
int_target_list_LWE = np.array([np.array([int(target_list_LWE[i][j][0]) for j in range(rank)]) for i in range(len(target_list_LWE))])
int_target_list_unif = np.array([np.array([int(target_list_unif[i][j][0]) for j in range(rank)]) for i in range(len(target_list_unif))])
print("Conversion en np.array (en s): ", time.time() - start)

print("Shapes check :")
print("Short dual vectors W | target vectors LWE | target vectors unif : ", int_W_transpose.shape ,  int_target_list_LWE.shape , int_target_list_unif.shape)

print("Convertion en cp.array...")
start = time.time()
cp_W_transpose = cp.asarray(int_W_transpose) 
cp_LWE = cp.asarray(int_target_list_LWE)
cp_unif = cp.asarray(int_target_list_unif)
print("Conversion en cp.array (en s): ", time.time() - start)

print("Distingueur GPU...")
start = time.time()
data_LWE = []
start_mul = time.time()
data_LWE = cp.matmul(cp_LWE, cp_W_transpose) % modulus
print("Calcul multiplication #1 GPU (en s): ", time.time() - start_mul)
start_mul = time.time()
data_unif = cp.matmul(cp_unif, cp_W_transpose) % modulus
print("Calcul multiplication #2 GPU (en s): ", time.time() - start_mul)

cp.cuda.Stream.null.synchronize()
start_cos = time.time()
cos_data_LWE , cos_data_unif = cp.cos(2/modulus * cp.pi * data_LWE) , cp.cos(2/modulus * cp.pi * data_unif)
print("Calcul cos GPU (en s): ", time.time() - start_cos)
print('Distinguisher values LWE | uniform : ',cp.mean(cos_data_LWE) , ' | ', cp.mean(cos_data_unif ))
print("Distinguisher running time : ", time.time() - start, ' sec')
