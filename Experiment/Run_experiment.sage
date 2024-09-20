from time import time

print('*************************************')
print('**   Experiment for Dual Attacks   **')
print('*************************************')
print('')
counter = time()

load('Parameters.sage')
print('Lattice of rank ' ,rank) 
if is_LWE :
	print('Lattice is qary with modulus ' , modulus , ' and index ' , log_covolume )

##
counter = time() - counter
print('Parameters : ' , counter , ' sec')
##
load('LWE_instance_generator.sage')

counter = time() - counter
print('LWE instance : ' , counter, ' sec')
##
load('Dual_short_vectors.sage')

counter = time() - counter
print('Dual short vectors : ' , counter, ' sec')
##
load('Target_samples.sage')

counter = time() - counter 
print('Target samples : ' , counter, ' sec')
print('Size of target samples : ', T)
##
load('Distinguisher_v2.sage')

counter = time() - counter
print('Distinguisher : ' , counter, ' sec')
##

