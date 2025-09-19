####################
##   Librairies   ##
####################

import numpy as np
from sage.stats.distributions.discrete_gaussian_integer import DiscreteGaussianDistributionIntegerSampler #, DiscreteGaussianDistributionLatticeSampler
from FFT_gemini import *        
from time import time

###################
##   Functions   ##
###################

def nearest_integer( float_number ) :
    # input  : float number x
    # output : integer z s.t. |x-z| is minimum, with convention nearest_integer( n + 0.5 ) = n+1 
    return( int( round( float_number ) ) )

def gaussian_sampler( dimension , stdev) :
    D = DiscreteGaussianDistributionIntegerSampler(sigma=stdev)
    return( vector(ZZ, [D() for k in range(dimension)]))

def exp_mod( float_number , q ) :
    a = (2 * I * np.pi ) / q 
    y = np.exp( a * int(float_number)) 
    return( complex(y) )

def short_vectors_get( B_column_basis , stdev , vector_number):
    B_LLL = B_column_basis.LLL()
    # DiscreteGaussianDistributionLatticeSampler takes as input a row basis
    sampler = distributions.DiscreteGaussianDistributionLatticeSampler(B_LLL.transpose() , stdev)
    return( [sampler() for k in range(vector_number)] )


def initialize_table( s_test_enum , short_vector_list , target , A_fft , A_enum ):
    TABLE = { }
    for w_dual_vector in short_vector_list :
        x_j = w_dual_vector[: m ]
        y_j_fft , y_j_enum  = x_j.transpose() * A_fft , x_j.transpose() * A_enum 
        index = nearest_integer( ( p / q ) * y_j_fft)
        if index in TABLE.keys() :
            TABLE[index] +=  exp_mod(x_j.transpose() * target - y_j_enum.transpose() * s_test_enum ,q  )
        else :
            TABLE[index] =  exp_mod(x_j.transpose() * target - y_j_enum.transpose() * s_test_enum , q)
        return(TABLE)

'''
def BitReverse( n_int , bitlen ) : 
    # input  : int n in [0,..,pow(2,bitlen)-1]
    # output : int in [0,..,pow(2,bitlen)-1] such that binary expansion is reversed of binary expansion of n
    bin_rep = bin(n_int)[2:]
    r = bitlen - len(bin_rep)
    if r > 0 :
        bin_rep = r *'0' + bin_rep
    else :
        bin_rep = bin_rep[:bitlen]
    bit_rev = bin_rep[::-1]
    return(int(bit_rev , 2))
'''

def ScoreEval(A , q , N_short_vector , L , target , s_enum_test):
	
	m , n = A.nrows() , A.ncols()
	k_enum , k_lat , k_fft
	assert (k_enum + k_lat + k_fft == n)

	p = pow(2,L) # modulo switching
	
	print('## Evaluating score')
	print('Modulo switching : reducing ', q , ' to ' , p)
	print('FFT')
	print('Dimension k_fft : ', k_fft ,' | size pow(p , k_fft ) : ' , pow(p,k_fft)) 

	A_enum , A_lat , A_fft = A[ : ,  : k_enum ] , A[ : , k_enum : k_enum + k_lat] , A[ : , k_enum + k_lat : ]
	B_lat_dual = block_matrix( [ [ identity_matrix(ZZ , m) , zero_matrix(ZZ,m,k_lat) ] , [ A_lat.lift().transpose() , q * identity_matrix(ZZ,k_lat)] ] )
	
	short_vector_list = short_vectors_get( B_lat_dual , sigma , N_short_vector)
	W = matrix(ZZ, N_short_vector , m + k_lat ,short_vector_list)    
	x_transpose = W[:,:m] 

	y_fft , y_enum = x_transpose * A_fft , x_transpose * A_enum
	c  = x_transpose * target
	u = y_enum * s_enum_test
	
	TABLE = {  }
	labels = [ vector( ZZ , [ nearest_integer( (p/q)*y_fft.lift()[i][j] ) for j in range(k_fft) ])  for i in range(N_short_vector)  ]

	for k in range(N_short_vector) :
		g = tuple(labels[k])
		if g in TABLE.keys() :
			TABLE[g] += exp_mod(c[k] - u[k] , q )
		else :
			TABLE[g] = exp_mod(c[k] - u[k] , q)

	
        
	TABLE_fft = fft_multidim_sparse(TABLE  , p , k_fft)

	score = { v : real(a) for v , a in TABLE_fft.items() }
	
	argmax_score = max(score , key=score.get)
	
	return(argmax_score , score[argmax_score])

####################
##   Experiment   ##
####################        

# LWE parameters
q = 937
m , n = 20 , 10
sigma = 3.0

# Dual sampling and modulo switching parameter
N_short_vector = 1000
L = 3

k_enum , k_lat , k_fft = 1 , 3 , 6
assert (k_enum + k_lat + k_fft == n)


print('LWE instance ')
print(' m , n , q : ' , m , n, q )

# LWE instance
A = random_matrix(  Zmod(q) , m , n)
secret , error = gaussian_sampler(n , sigma) , gaussian_sampler(m , sigma)
target = A * secret + error

print('Secret : ' , secret)

s_enum_test = gaussian_sampler( k_enum , sigma)  
start = time()
print("Score : " , ScoreEval(A , q , N_short_vector , L , target, s_enum_test))
print("Execution time : ", time() - start, ' s')

s0 = secret[:k_enum]
start = time()
print("Score : " , ScoreEval(A , q , N_short_vector , L , target, s_enum_test))
print("Execution time : ", time() - start, ' s')


# Evaluation lazy
# test = [ np.cos( 2*np.pi * I *target.lift().inner_product(w[:m] )/ q) for w in short_vector_list]              
#sum(test)
