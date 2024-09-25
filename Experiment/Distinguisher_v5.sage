import numpy as np
import cupy as cp
from time import time 

from matplotlib import pyplot as plt
# Set matplotlib and seaborn plotting style : might need to install seaborn in sage shell by running pip install seaborn
import seaborn as sns
sns.set_style('darkgrid')

############################
##   Loading parameters   ##
############################

start = time()
load('Parameters.sage')
modulus , short_vector_list        = load('Dual_short_vectors/dual_short_vectors')
target_list_LWE , target_list_unif = load('Target_samples/target_samples')

print("Loading time : ", time() - start , ' sec')

W_len , T_len = len(short_vector_list) , min(len(target_list_LWE) , len(target_list_unif))
modulus = int(modulus)

##################################
##   Conversion to cupy array   ##
##################################

number_bound = 10000

counter = time()
W = cp.array(short_vector_list[0:number_bound]) 	
target_LWE , target_unif = cp.array([ target.lift() for target in target_list_LWE])[:,:,0].transpose() , cp.array([ target.lift() for target in target_list_unif])[:,:,0].transpose()

print('Conversion to cupy array : ', time() - counter ,  ' sec.')
counter = time() 	
 	
# cp.array of shape ( W_len , rank )
print(W.shape)# , target_LWE.shape )
print(target_LWE.shape)

counter = time()
data_LWE , data_unif  = cp.matmul(W , target_LWE) % modulus , cp.matmul(W , target_unif) % modulus
print('Computation of inner products : ', time() - counter , ' sec.')
 	
cp.cuda.Stream.null.synchronize()

counter = time()
cos_data_LWE , cos_data_unif = cp.cos(2/modulus * cp.pi * data_LWE) , cp.cos(2/modulus * cp.pi * data_unif)
print('Computation of cosine : ', time() - counter , ' sec.')

counter = time()

AR_distinguisher_LWE , AR_distinguisher_unif =  cp.mean(cos_data_LWE,axis = 0) , cp.mean(cos_data_unif, axis=0 )

print('Mean of distinguisher values LWE | uniform : ', cp.mean(AR_distinguisher_LWE) , ' | ', cp.mean(AR_distinguisher_unif) )
print('Computation of distinguisher values : ', time() - counter , ' sec.')
print('Total execution time : ', time() -start, ' sec.')

################
##   Graphs   ##
################

N_LWE , N_unif = len(AR_distinguisher_LWE) , len(AR_distinguisher_unif)

yy_LWE , yy_unif = np.linspace(0,1,N_LWE) , np.linspace(0,1,N_unif)
#AR_distinguisher_LWE.sort() , AR_distinguisher_unif.sort()

# figure, axis = plt.subplots(2) 

plt.plot( AR_distinguisher_LWE.get()  , yy_LWE  , color = "blue", label="LWE" )
plt.plot( AR_distinguisher_unif.get() , yy_unif , color = "cyan", label="Uniform" )
plt.title("Empirical c.d.f. of distinguisher values")
plt.savefig('Figures/cdf_distinguisher_LWE_rank50_20240925.png', bbox_inches='tight')


























