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
W_tuple = ()
for w in short_vector_list :
	W_tuple += w

W_matrix = matrix(ZZ , W_len , rank , W_tuple)

T_matrix_LWE , T_matrix_unif = target_list_LWE[0] , target_list_unif[0]
for k in range(T_len - 1) :
	T_matrix_LWE , T_matrix_unif = block_matrix([[ T_matrix_LWE , target_list_LWE[k+1] ]]) , block_matrix([[ T_matrix_unif , target_list_unif[k+1] ]])

print('Conversion : ' , time() - counter) 	###

# Multiplication (computation of inner products)

counter = time()				###	
data_LWE , data_unif = W_matrix * matrix(ZZ , T_matrix_LWE) , W_matrix * matrix(ZZ , T_matrix_unif) 

print('Inner products : ' , time() - counter) 	###
# Computation of cosine

counter = time() 				###
cos_data_LWE , cos_data_unif = np.cos( (2/ modulus) * np.pi * data_LWE) , np.cos((2/ modulus) * np.pi * data_unif)
print('Cosine : ' , time() - counter) 		###

test_LWE , test_unif = np.array(data_LWE) , np.array(data_unif)
counter = time() 				###
cos_np_LWE , cos_np_unif = np.cos( (2/ modulus) * np.pi * test_LWE) , np.cos((2/ modulus) * np.pi * test_unif)
print('Cosine  with numpy: ', time() - counter)	###


print('Distinguisher values LWE | uniform : ', np.mean(cos_data_LWE) ,' | ', np.mean(cos_data_unif )  )
print('Distinguisher values LWE | uniform : ', np.mean(cos_np_LWE) ,' | ', np.mean(cos_np_unif )  )

'''
for target_LWE  in target_list_LWE :
	D_LWE = 0
	for w in short_vector_list : 
		ww , tt = vector( ZZ, w ) , vector(target_LWE)
		D_LWE  += np.cos(2 * np.pi * RR(ww.inner_product(tt)))
	data_LWE.append((1/W) * D_LWE)
		
for target_unif in target_list_unif :
	D_unif = 0
	for w in short_vector_list : 
		ww , tt = vector( ZZ, w ) , vector(target_unif)
		D_unif += np.cos(2 * np.pi * RR(ww.inner_product(tt)))		
	data_unif.append((1/W) * D_unif)
'''
	
# save((data_LWE , data_unif) , 'Data/distinguisher_values') 
