import numpy as np
from time import time 					###

counter = time() 
load('Parameters.sage')
modulus , short_vector_list = load('Dual_short_vectors/dual_short_vectors')
target_list_LWE , target_list_unif = load('Target_samples/target_samples')

print('Load : ' , time() - counter) 		###

W_len , T_len = len(short_vector_list) , min(len(target_list_LWE) , len(target_list_unif))

# Conversion
counter = time() 				###
W_array , T_array_LWE , T_array_unif = np.array(short_vector_list) , np.transpose(np.array(target_list_LWE)[:,:,0]) , np.transpose(np.array(target_list_unif)[:,:,0]) 
print('Conversion : ' , time() - counter) 	###

# Multiplication (computation of inner products)

counter = time()				###	
data_LWE , data_unif = np.matmul(W_array , T_array_LWE) , np.matmul(W_array , T_array_unif )

print('Inner products : ' , time() - counter) 	###
# Computation of cosine

counter = time() 				###
cos_data_LWE , cos_data_unif = np.cos( (2/ modulus) * np.pi * data_LWE) , np.cos((2/ modulus) * np.pi * data_unif)
print('Cosine : ' , time() - counter) 		###

print('Distinguisher values LWE | uniform : ', np.mean(cos_data_LWE) ,' | ', np.mean(cos_data_unif )  )
	
# save((data_LWE , data_unif) , 'Data/distinguisher_values') 
