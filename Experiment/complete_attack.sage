###################
##   Libraries   ##
###################

from fpylll import IntegerMatrix
from fpylll.util import gaussian_heuristic
from g6k import Siever 

###################
##   Functions   ##
################### 

def siever(BASIS , RANK):
	##     Dual short vectors of dual/perp lattice via sieve     ##
	##   Modification of g6k example file all_short_vectors.py   ##
	print("Running sieve via g6k")

	short_vector_list = []	
	
	g6k = Siever(BASIS)
	g6k.lll(0, RANK)

	g6k.initialize_local(0, RANK /2 , RANK )
	while g6k.l > 0:
    	# Extend the lift context to the left
    		g6k.extend_left(1)
    	# Sieve
    		g6k()

	with g6k.temp_params(saturation_ratio=.95, saturation_radius=1.7, db_size_base=sqrt(1.7), db_size_factor=5):
    		g6k()

	# Convert all data_base vectors from basis A to cannonical basis and print them 
	# out if they are indeed shorter than 1.7 * gh^2

	gh = gaussian_heuristic([g6k.M.get_r(i, i) for i in range(RANK)])

	data_base = list(g6k.itervalues())
	found = 0

	for x in data_base:
    		v = BASIS.multiply_left(x)
    		l = sum(v_**2 for v_ in v)
    		if l < 1.7 * gh:
        	# print(l/gh, v)
        		found += 1
        		short_vector_list.append(v)
	
	print("Found %d vectors of squared length than 1.7*gh. (expected %f)"%(found, .5 * 1.7**(RANK /2.)))
	
	return(short_vector_list)

def Distinguisher(short_vector_list , target , modulus):
	import numpy as np
	import cupy as cp
	
	W = cp.array(short_vector_list)
	target = cp.array(target) #.transpose() 

	print(W.shape)
	print(target.shape)

	data  = cp.mod(cp.matmul(W , target) , int(modulus))

	cp.cuda.Stream.null.synchronize()

	cos_data = cp.cos(2/modulus * cp.pi * data)

	score = cp.mean(cos_data ,axis = 0)

	return(score)


####################
##   Parameters   ##
#################### 

rank , log_covolume , modulus = 20 , 10 , 937
rank_guess = 3  

####################
##   Experiment   ##
####################

print("Lattice :")
print("q-ary")
print("Rank, index , modulus : ",rank , log_covolume , modulus )

B = IntegerMatrix.random( rank , 'qary' , k = log_covolume , q = modulus)
m, m_guess = rank - log_covolume , rank_guess 


A = matrix(ZZ,B)[0:log_covolume, m:rank]
A_guess = A[: rank_guess]

B_perp = block_matrix( [ [modulus * matrix.identity(ZZ,m) , matrix.zero(ZZ ,m , log_covolume )] , [ - A.transpose() , matrix.identity(ZZ , log_covolume)] ] )
B_perp = IntegerMatrix.from_matrix(B_perp)

B_guess_perp = block_matrix( [ [modulus * matrix.identity(ZZ,m_guess) , matrix.zero(ZZ ,m_guess , log_covolume )] , [ - A_guess.transpose() , matrix.identity(ZZ , log_covolume)] ] )
B_guess_perp = IntegerMatrix.from_matrix(B_guess_perp)

dual_short_vectors = siever(B_guess_perp, log_covolume + rank_guess)

target = [ randint(0,modulus) for i in range(log_covolume+rank_guess) ]

print('test : ',Distinguisher(dual_short_vectors , target , modulus))