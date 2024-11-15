##################
##   PACKAGES   ##
##################

from fpylll import FPLLL, IntegerMatrix

import numpy as np
import cupy as cp
from time import time 

from sage.stats.distributions.discrete_gaussian_integer import DiscreteGaussianDistributionIntegerSampler \
    as DiscreteGaussian

from matplotlib import pyplot as plt
# Set matplotlib and seaborn plotting style : might need to install seaborn in sage shell by running pip install seaborn
import seaborn as sns
sns.set_style('darkgrid')

####################
##   Parameters   ##
####################
start = time()

stdev, rank , modulus  = 0.001 , 62 , 931 	# rank can be 20, 30, 40, 62, 64, 66, 68
log_covolume = rank / 2	
m  = rank - log_covolume			# m = n - k

R  = Zmod(modulus)	 			# ZZ mod q

seed_FPLLL = 1337
FPLLL.set_random_seed(seed_FPLLL)

D = DiscreteGaussian(stdev)	

#######################
##   q-ary lattice   ##
#######################

B = IntegerMatrix.random( rank , 'qary' , k = log_covolume , q = modulus)
BB = matrix(ZZ,B)				# basis of L = L_q(A)
A  = BB[0:log_covolume, m:rank]			# A

'''
B_perp = block_matrix( [ [modulus * matrix.identity(ZZ,m) , matrix.zero(ZZ ,m , log_covolume )] , [ - A.transpose() , matrix.identity(ZZ , log_covolume)] ] )
B_perp = IntegerMatrix.from_matrix(B_perp) 	# L perp (rescaled dual)
'''

########################
##   Target samples   ##
########################
T = 1000

target_list_unif , target_list_LWE = [] , []
for i in range(T):
	# s , e = random_matrix(R , 1 , rank ) , random_matrix(R , rank , 1 )
	s = random_matrix(R , 1 , rank )
	e = matrix(ZZ, [D() for _ in range(rank)]).transpose()
	target_list_LWE.append( (s * BB).transpose() + e ) , target_list_unif.append( random_matrix(R , rank , 1 ) ) # BB * s + e )


#######################################
##   Dual short vectors from sieve   ##
#######################################

with open('dual_short_vectors_'+str(rank)+'_931.txt', 'r' , encoding="utf-8") as file :
	data = file.readlines()
	file.close()


short_vector_list = [cp.array(eval(a)) for a in data]
	
#######################
##   Distinguisher   ##
#######################

number_bound = 10000

W = cp.array(short_vector_list[0:number_bound]) 
target_LWE , target_unif = cp.array([ target.lift() for target in target_list_LWE])[:,:,0].transpose() , cp.array([ target.lift() for target in target_list_unif])[:,:,0].transpose()

modulus = int(modulus)
data_LWE , data_unif  = cp.matmul(W , target_LWE) % modulus , cp.matmul(W , target_unif) % modulus


cp.cuda.Stream.null.synchronize()


cos_data_LWE , cos_data_unif = cp.cos(2/modulus * cp.pi * data_LWE) , cp.cos(2/modulus * cp.pi * data_unif)
AR_distinguisher_LWE , AR_distinguisher_unif =  cp.mean(cos_data_LWE,axis = 0) , cp.mean(cos_data_unif, axis=0 )

print('Mean of distinguisher values LWE | uniform : ', cp.mean(AR_distinguisher_LWE) , ' | ', cp.mean(AR_distinguisher_unif) )

################
##   Graphs   ##
################

N_LWE , N_unif = len(AR_distinguisher_LWE) , len(AR_distinguisher_unif)

yy_LWE , yy_unif = np.linspace(0,1,N_LWE) , np.linspace(0,1,N_unif)
AR_distinguisher_LWE.sort() , AR_distinguisher_unif.sort()

# figure, axis = plt.subplots(2) 

plt.plot( AR_distinguisher_LWE.get()  , yy_LWE  , color = "blue", label="LWE" )
plt.plot( AR_distinguisher_unif.get() , yy_unif , color = "cyan", label="Uniform" )
plt.title("Empirical c.d.f. of distinguisher values")
plt.savefig('Figures/cdf_distinguisher_LWE_rank_'+str(rank)+'_modulus'+str(modulus)+'_20241115.png', bbox_inches='tight')
#plt.show()
print('Execution time : ', time() - start , ' sec.')


