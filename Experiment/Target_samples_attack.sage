################
##   IMPORT   ##
################

from sage.stats.distributions.discrete_gaussian_integer import DiscreteGaussianDistributionIntegerSampler \
    as DiscreteGaussian
## import numpy as np
import cupy as cp
import numpy as np
from time import time 
from matplotlib import pyplot as plt

# Set matplotlib and seaborn plotting style 
# Might need to install seaborn in sage shell by running pip install seaborn
import seaborn as sns
sns.set_style('darkgrid')

##   Load lattice basis   ##
load('Parameters.sage')
[ is_LWE, stdev , rank , modulus , log_covolume ] , B = load('Lattice_basis/basis_LWE.sobj')

##   Paremeters   ##
dim_guess = 4
stdev = 0.01
D = DiscreteGaussian(stdev)

#########################
##   TARGET SAMPLING   ## 
#########################
# TO DO 
# Vérifier le nombre de zeros sur target et target list
# Prendre en compte la dimension dim_guess dans la partie sur le calcul du distingueur.
# Ecrire un script independant pour le sieve afin de le faire tourner sur Vegeta

if is_LWE :
	BB , m = matrix(ZZ,B) , rank - log_covolume 	# m = n - k 
	A  = BB[0:log_covolume, m:rank]			# A in ZZ^{ m x k } such that B = qary Lattice of LWE (A)
	R  = Zmod(modulus) 				# ZZ mod q
		
	A_guess = random_matrix(R, m, dim_guess )
	B_total = block_matrix( ZZ , [[ matrix.identity(ZZ,m) , A , A_guess ] , [ matrix.zero(ZZ, log_covolume , m) 	, modulus * matrix.identity(ZZ,log_covolume),  matrix.zero(ZZ, log_covolume , dim_guess) ] , [ matrix.zero(ZZ, dim_guess , m), matrix.zero(ZZ, dim_guess , log_covolume )	, modulus*matrix.identity(ZZ, dim_guess ) ]])
	
	s , e = random_matrix(R ,  rank + dim_guess ,1) , matrix(ZZ,rank + dim_guess ,1, [D() for _ in range( rank + dim_guess )]) 
	s_guess = s[ rank:]
	target = B_total * s + e
	
	target_list = []			
	for i in range(T): 
		s_try = random_matrix( R , dim_guess , 1)
		t_try = block_matrix([ [A_guess * s_try ] , [ matrix.zero(R , log_covolume , 1) ] , [ modulus * s_guess ] ] )
		target_list.append( target - t_try) 
			
# Save			
file_path = 'Target_samples/target_samples_attack_isLWE_' + str(is_LWE)+'_stdev_' +str(stdev )+'_rank_' +str(rank )+'_modulus_' +str(modulus )+'_logcovolume' +str(log_covolume) 			
save( ( s_guess , target , target_list) , file_path )		

#####################
### DISTINGUISHER ###
#####################

############################
##   Loading parameters   ##
############################

start = time()
modulus , short_vector_list        = load('Dual_short_vectors/dual_short_vectors.sobj')

#target_list_LWE , target_list_unif = load('Target_samples/target_samples')

print("Loading time : ", time() - start , ' sec')

W_len , T_len = len(short_vector_list) , len(target_list)
modulus = int(modulus)

##################################
##   Conversion to cupy array   ##
##################################

number_bound = 10000 		# Number of samples of dual short vectors used

counter = time()
W = cp.array(short_vector_list[0:number_bound]) 
WW = cp.array(short_vector_list[0:7*number_bound]) 	

print(target)
target_matrix = cp.array([ target_vector.lift() for target_vector in target_list])[0,:,:] #.transpose()


print('Conversion to cupy array : ', time() - counter ,  ' sec.')
counter = time() 	
 	
# cp.array of shape ( W_len , rank )
print('W shape : ' , W.shape)# , target_LWE.shape )
print('target_matrix shape : ', target_matrix.shape)

#print(WW.shape)
#print('target shape : ', target.shape)

counter = time()
print('A DEBUGGER : target est cropped' ,target_matrix[0:rank,:].shape) ## A ENLEVER : les tailles ne sont pas bonnes !
data = cp.matmul(W , target_matrix[0:rank,:]) % modulus 
print('Computation of inner products : ', time() - counter , ' sec.')
 	
cp.cuda.Stream.null.synchronize()

counter = time()
cos_data = cp.cos(2/modulus * cp.pi * data) 
print('Computation of cosine : ', time() - counter , ' sec.')

counter = time()

AR_distinguisher =  cp.mean(cos_data,axis = 0) 

print('Mean of distinguisher values : ', cp.mean(AR_distinguisher) )
print('Computation of distinguisher values : ', time() - counter , ' sec.')
print('Total execution time : ', time() -start, ' sec.')

################
##   Graphs   ##
################

N = len(AR_distinguisher) 
yy = np.linspace(0,1,N) 
AR_distinguisher.sort() 

# figure, axis = plt.subplots(2) 

# Save
plt.plot( AR_distinguisher.get()  , yy  , color = "blue", label="Distinguisher values" )
plt.title("Empirical c.d.f. of distinguisher values")
file_path_figure = 'Figures/cdf_distinguisher_attack_isLWE_' + str(is_LWE)+'_stdev' +str(stdev )+'_rank' +str(rank )+'_modulus' +str(modulus )+'_logcovolume' +str(log_covolume)+'.png'
plt.savefig(file_path_figure , bbox_inches='tight')
plt.show()



























