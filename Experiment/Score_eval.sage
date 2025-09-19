import numpy as np
from sage.stats.distributions.discrete_gaussian_integer import DiscreteGaussianDistributionIntegerSampler #, DiscreteGaussianDistributionLatticeSampler

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

#############
##   FFT   ##
#############

def addition_tuple(v1, v2, modulus):
    """
    Addition de deux vecteurs (tuples) modulo q.
    """
    resultat = tuple((x + y) % modulus for x, y in zip(v1, v2))
    return(resultat)

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



# ML-DSA
'''
q = pow(2,23) - pow(2,13) + 1   # Solinas prime
rank = 256 

z0 = power_mod(1753 , 2 * BitReverse(randint(0,rank // 2)) + 1 , q)                        # 512-th root of unity in Zmod(q)
zeta = { i : power_mod(z0,BitReverse(i), q) for i in range(1, rank)}

print('q , rank , root of unity : ', q , ' | ' , rank, ' | ' , z0 ) 
'''
# def FFT(T , index ,modulo ):

def FFT(T , modulo):
    assert len(T) == rank , "input of FFT : wrong size in FFT.sage/FFT"
    W = T.copy()
    m , L = 0 , modulo // 2
    while L > 0 :
        start = 0 
        while start < rank : 
            m += 1
            z = zeta[m] # power_mod( z0 , BitReverse(m) , q ) #   
            for j in range(start , start + L ):
                t = z * W[ (j + L) % rank ] % q
                W[ (j + L) % rank] , W[ j % rank ] = (W[j % rank] - t ) % q , (W[j % rank] + t ) % q
            start += 2* L        
        L = floor(L / 2)
    return(W)

def FFT_inverse(W):
    assert len(W) == rank , "input of NTT_inverse : wrong size in NTT.sage/NTT"
    T = W.copy()
    m , L = 256 , 1
    while L < 256 :
        start = 0
        while start < 256 :
            m -= 1
            z = - zeta[m]
            for j in range(start, start + L) :
                t = T[j % rank]
                T[ j % rank] , T[ (j + L) % rank ] = (t + T[ (j + L) % rank ] ) % q , z * (t - T[ (j + L) % rank ] ) % q
                # T[ (j + L) % rank ] = ( z * T[ (j + L) % rank ] )  % q
            start += 2 * L
        L *= 2  
    f = 8347681             # inverse(256) mod q
    return({ j : (f * T[j]) % q for j in range(rank)})

def ScoreEval(A , q , N_short_vector , L , target , s_enum_test):
	
	m , n = A.nrows() , A.ncols()
	k_enum , k_lat , k_fft
	assert (k_enum + k_lat + k_fft == n)

	p = pow(2,L) # modulo switching

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

sigma = 3.0
q = 937
m , n = 20 , 10
N_short_vector = 1000

L = 3

k_enum , k_lat , k_fft = 1 , 3 , 6
assert (k_enum + k_lat + k_fft == n)

A = random_matrix(  Zmod(q) , m , n)
secret , error = gaussian_sampler(n , sigma) , gaussian_sampler(m , sigma)
target = A * secret + error


s_enum_test = gaussian_sampler( k_enum , sigma)  

        
from FFT_gemini import *        

from time import time
start = time()

# TABLE = initialize_table( s_test_enum , short_vector_list , target , A_fft , A_enum )

'''
print('LWE pair : ' , A ,' | ' , target)
print(B_lat_dual)
'''
#print(TABLE)
start = time()
print("Test function" , ScoreEval(A , q , N_short_vector , L , target, s_enum_test))
print("Execution time : ", time() - start, ' s')
# Evaluation lazy
# test = [ np.cos( 2*np.pi * I *target.lift().inner_product(w[:m] )/ q) for w in short_vector_list]              
#sum(test)
