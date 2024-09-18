import numpy as np

modulus , short_vector_list = load('Dual_short_vectors/dual_short_vectors')
target_list_LWE , target_list_unif = load('Target_samples/target_samples')

W = len(short_vector_list)

data_LWE , data_unif = [] , [] 

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
	
save((data_LWE , data_unif) , 'Data/distinguisher_values') 
