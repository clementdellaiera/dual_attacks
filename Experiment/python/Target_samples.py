##   Load lattice basis 
from Parameters import *
from LWE_instance_generator import *
import numpy as np

from sage.stats.distributions.discrete_gaussian_integer import DiscreteGaussianDistributionIntegerSampler \
    as DiscreteGaussian

D = DiscreteGaussian(stdev)
e = np.array([D() % modulus for _ in range(rank)]).transpose()


if is_LWE :
	target_list_unif , target_list_LWE = [] , []
	# R  = Zmod(modulus) 				# ZZ mod q
	m  = rank - log_covolume			# m = n - k
	A = B.submatrix(range(0,int(log_covolume)), range(int(m), rank))  # A in ZZ^{ m x k } such that B = qary Lattice of LWE (A)
	# print(A)
	for i in range(1):
		# s , e = random_matrix(R , 1 , rank ) , random_matrix(R , rank , 1 )
		s = np.array([np.random.randint(modulus-1) for _ in range(rank)])
		_BB = [[0 for _ in range(rank)] for _ in range(rank)]
		_ = A.to_matrix(_BB)
		BB = np.array(_BB)
		v = np.matmul(s, BB) % modulus

		target_list_LWE.append( v.transpose() + e ) 
		target_list_unif.append( np.array([np.random.randint(modulus-1) for _ in range(rank)]) ) # BB * s + e )
		
	# print(target_list_LWE)
		# A_big = [ [Id ] , [A] ]  
		# A_big * s + e = s_0 + e_0 		, A s_1 + e_1
		# B * s + e 	= s_0 + A * s_1 + e_0 	, q s_1 + e_1


# print(len(target_list_LWE[0]))
# print(len(target_list_unif[0]))
# save( ( target_list_LWE , target_list_unif) , 'Target_samples/target_samples' )		
