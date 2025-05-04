###################
##   Libraries   ##
###################

from matplotlib import pyplot as plt
import seaborn as sns
sns.set_style('darkgrid')

import numpy as np

from sage.stats.distributions.discrete_gaussian_integer import DiscreteGaussianDistributionIntegerSampler as DiscreteGaussian

###################
##   Functions   ##
################### 

# Sampler
def discrete_gaussian_vector(s_param , rank , modulo):
	gaussian_sampler = DiscreteGaussian(s_param)
	if modulo == 0 :
		output = vector( ZZ , [ gaussian_sampler() for i in range(rank)] ) 
	else :
		output = vector( Zmod(modulo) , [ gaussian_sampler() for i in range(rank)] ) 
			# vector( ZZ , [ mod_plus(gaussian_sampler() , modulo) for i in range(rank)] )
	return(output)

def uniform_qary_vector(rank, modulo):
	assert modulo > 1 , "line 28 : modulo less than 2"
	output = vector(Zmod(modulo),[ randint( - modulo // 2 , (modulo-1) //2 ) for i in range(rank)])
	return(output)

# Integers to Zmod(q)
def nearest_integer(x_input):
	# input  : x float
	# output : n int defined as argmin d(x,n) for n in ZZ, and if x in ZZ +0.5 then x+0.5
	return(sign(x_input) * floor( abs(x_input) + 0.5 ))
	
def mod_plus(x_input , modulo):
	# input  : x int , modulo int 
	# output : unique y int in { -L , L - 1 } such that x = y mod modulo 
	L = floor( (modulo - 1) / 2 )
	x_mod = x_input % modulo
	if x_mod >= L :
		return(x_mod - modulo)
	else : 	
		return(x_mod)

# Main
def modulo_switch( x_input , modulo_in , modulo_out ):
	# input  : x int , modulo in and out int
	# output : nearest integer of rescaled x  
	x_scaled = ZZ(x_input) * modulo_out / modulo_in # should it be mod_plus(x_input , modulo_in) * modulo_out / modulo_in  
	x_rounded = nearest_integer(x_scaled)
	return(x_rounded)

#####################
##   Simulations   ## 
#####################

def gaussian_distribution_before_and_after_modulo_switching(sample_size):
	# Experiment : 
	# 	1. sample U(Z mod q) and discrete_gaussian(Z,stdev) mod q
	# 	2. Apply modulo swithcing function
	# 	3. Compare and graph empirical cdf
	data_unif =  [ randint(- q // 2 , (q-1) //2) for i in range(sample_size) ]
	data_gaussian = [ D() for i in range(sample_size) ]

	data_unif_modulo_switched , data_gaussian_modulo_switched = [modulo_switch(x , q , N) for x in data_unif] , [modulo_switch(x , q , N) for x in data_gaussian]
 
	data_unif.sort() , data_gaussian.sort()
	data_unif_modulo_switched.sort() , data_gaussian_modulo_switched.sort() 

	fig , axs = plt.subplots(1,2)
	axs[0].plot(data_unif , [i for i in range(sample_size)])
	axs[0].plot(data_gaussian , [i for i in range(sample_size)])
	axs[1].plot(data_unif_modulo_switched , [i for i in range(sample_size)])
	axs[1].plot(data_gaussian_modulo_switched , [i for i in range(sample_size)])

	axs[0].set(xlabel= "modulo q = "+str(q) )
	axs[1].set(xlabel= "modulo N = 2**"+str(L-4)+" = "+str(N) )

	fig.suptitle('Modulo switching for uniform and gaussian distribution')
	plt.show()
	#plt.savefig("Figures/modulo_switching_q"+str(q)+"N"+str(N)+".png")
	return(data_unif , data_gaussian , data_unif_modulo_switched , data_gaussian_modulo_switched)

def score_distribution(A , stdev , stdev_dual):
	target_LWE , target_unif  = [ A * discrete_gaussian_vector(stdev , n , q) +  discrete_gaussian_vector(stdev , m , q)  for i in range(sample_size)] , [ uniform_qary_vector(m, q) for i in range(sample_size)]
	# A * discrete_gaussian_vector(stdev , n , q) +  discrete_gaussian_vector(stdev , m , q)
	A_q = matrix(Zmod(q),A)
	A_q_dual = A_q.transpose().inverse()
	short_dual_vectors_number = 100
	short_dual_vectors = [ A_q_dual * discrete_gaussian_vector(stdev_dual , n , q) for i in range(short_dual_vectors_number)] 

	"""
	cos_LWE , cos_unif			= [ np.cos(2 * np.pi * int(w.inner_product(t))/ q) for w in short_dual_vectors  for t in target_LWE ] , [ np.cos(2 * np.pi * int(w.inner_product(t))/ q) for w in short_dual_vectors for t in target_unif ]
	cos_LWE.sort() , cos_unif.sort()
	"""
	score_LWE , score_unif = [ np.mean([ np.cos(2 * np.pi * int(w.inner_product(t))/ q) for w in short_dual_vectors ])  for t in target_LWE ] , [ np.mean([np.cos(2 * np.pi * int(w.inner_product(t))/ q) for w in short_dual_vectors ]) for t in target_unif ]
	score_LWE.sort() , score_unif.sort()

	plt.clf()
	'''
	plt.plot(cos_LWE  , [i for i in range(sample_size*short_dual_vectors_number)] , color = "cyan" , label="LWE target")
	plt.plot(cos_unif , [i for i in range(sample_size*short_dual_vectors_number)] , color = "gray", label = "Uniform target")
	'''
	plt.plot(score_LWE  , [i for i in range(sample_size)] , color = "cyan" , label="LWE target")
	plt.plot(score_unif , [i for i in range(sample_size)] , color = "gray" , label = "Uniform target")

	plt.title("Distinguisher")
	plt.show()

####################
##   Parameters   ##
####################

q = 937
L = floor(np.log2(q))
N = pow(2,L-4) 

m , n = 10 , 10 

stdev = q // 10
D = DiscreteGaussian(stdev)

sample_size = 10000

####################
##   Experiment   ##
####################


stdev , stdev_dual = 0.05 , 10
short_dual_vectors_number = 100

# Sampling LWE matrix
condition = True	
while condition :
	assert m == n , "line 136 : A not a square matrix" 
	A = matrix( ZZ , m , n , [ randint( - q // 2 , (q-1) //2 ) for i in range(m*n)])
	condition = not(matrix(Zmod(q),A).is_invertible())

A_N = matrix(ZZ, m,n, [modulo_switch(A[i,j],q,N) for i in range(m) for j in range(n)] ) 
A_N_dual = A_N.transpose().inverse()
A_N_perp = A_N.adjugate().transpose() # transpose of the comatrix : det(A_N) * A_N_dual == A_N_perp

target_LWE , target_unif  = [ A * discrete_gaussian_vector(stdev , n , q) +  discrete_gaussian_vector(stdev , m , q)  for i in range(sample_size)] , [ uniform_qary_vector(m, q) for i in range(sample_size)]
short_dual_vectors = [ A_N_perp * discrete_gaussian_vector(stdev_dual , n , N) for i in range(short_dual_vectors_number)] 

target_LWE_modulo_switched , target_unif_modulo_switched = [ [ modulo_switch(t[i] , q , N ) for i in range(n) ] for t in target_LWE ] , [ [ modulo_switch(t[i] , q , N ) for i in range(n) ] for t in target_unif ]
'''
score_LWE , score_unif = [ np.mean([ np.cos(2 * np.pi * int(w.inner_product(t))/ q) for w in short_dual_vectors ])  for t in target_LWE_modulo_switched ] , [ np.mean([np.cos(2 * np.pi * int(w.inner_product(t))/ q) for w in short_dual_vectors ]) for t in target_unif_modulo_switched ]
score_LWE.sort() , score_unif.sort()

# Graphs
plt.clf()	
plt.plot(score_LWE  , [i for i in range(sample_size)] , color = "cyan" , label="LWE target")
plt.plot(score_unif , [i for i in range(sample_size)] , color = "gray" , label = "Uniform target")
plt.title("Distinguisher")
plt.show()

# score_distribution(A , stdev , stdev_dual)

'''























