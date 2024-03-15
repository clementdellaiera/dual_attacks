##################
##   PACKAGES   ##
##################

from sage.all_cmdline import * 
import time

import numpy as np
from numpy.random import normal, uniform, binomial

from matplotlib import pyplot as plt

from fpylll import *

from itertools import product

from sage.stats.distributions.discrete_gaussian_lattice import DiscreteGaussianDistributionLatticeSampler 
# Documentation at https://doc.sagemath.org/html/en/reference/stats/sage/stats/distributions/discrete_gaussian_lattice.html

###################
##   FUNCTIONS   ##
###################

def Dual(B): # given an integer matrix, returns det(B) , Bv , Lv : determinant, integer matrix such that (1/det(B)) * Bv is a basis of dual lattice, dual lattice
	D = det(B)
	BB = D * B.inverse().transpose()
	Bv = BB.change_ring(ZZ)
	return(D , Bv , (1/D) * Bv.image() )

def hypercube( a , d ): 
	# Returns the list of all integer vectors inside the hypercube 
	# of dimension d and of side 2a+1 (i.e. the max-norm centered ball of Z*d of radius a) 
	L = []
	if (d < 1) or (a < 1):
		return(False)
	elif (d == 1 ):
		for k in range(-a , a+1):
	 		L.append([k])
	else:  
		for v in hypercube(a , d-1):
			for k in range(-a,a+1):
				L.append(v+[k])
	return(L)
		
def Sample_Dual(B,k,T, option):
	# Returns about T**k vectors sampled in the dual
	# Random_with_replacement and Random_Walk : random sample (no control on the size)
	# Sieve option : [| -T/2 , T/2 |] * Bv_k where Bv_k is the basis matrix 
	# of the first k vectors of the BKZ reduction of a dual basis 
	# (obtained with the function Dual)
	D , Bv, Lv = Dual(B)
	Bv_BKZ = Bv.BKZ()
	Bv_k = Bv_BKZ[ 0:k ]
	sample = []
	if (option == "Random_with_replacement"):
		for i in range(T**k):
			sample.append( Bv_k.image().random_element() ) 
			# Basis are maybe in row forms, ie if B = V.basis(), xB in V
	elif (option == "Sieve"):
		for v in hypercube(T//2,k) :
			sample.append( vector(ZZ ,v ) * Bv_k )
	elif (option == "Random_walk"):
		v = 0
		for i in range(T**k):
			u = randint(0,k-1)
			v += Bv_k[u] 
			sample.append(v)	
	return(sample)

def Distinguisher_Sieve(B , t, k, N): 
	# B basis of lattice, X sample of target vector , k , N size of dual vector sample
	# Returns the value of the AR-distinguisher evaluated wrt a sample of short dual vectors
	# obtained via teh sieve method in sample dual
	D , Bv, Lv = Dual(B)
	Bv_BKZ = Bv.BKZ()
	Bv_k = Bv_BKZ[ 0:k ]
	T = 0
	C = 0
	for v in hypercube(N//2,k) :
			w = vector(ZZ ,v ) * Bv_k
			T += np.cos(2 * np.pi * w.inner_product(t))
			C += 1
	return((1/C)*T)   

def Distinguisher_Sieve_preprocessing(B , W , t,k): 
	# Same as Distinguisher_sieve, but allows for external reduction of the basis 
	# and sampling of the dual vectors
	# ie returns sum of cos( 2 * pi * < w, t> )
	# where w ranges W * B_k
	B_k = B[ 0:k ]
	T = 0
	C = 0
	for w in W:
			T += np.cos(2 * np.pi * w.inner_product(t))
			C += 1
	return((1/C)*T)  
	
#################
##   SAMPLER   ##
#################

def LWE(n): 
	# Returns one sample of a LWE distribution (A, As+e) where A in ZZ**(n,n)
	A = random_matrix(ZZ,n,n)
	s = random_vector(ZZ , n  )
	e = random_vector(ZZ , n  )
	sk , pk = s , (A,A*s+e)
	return(sk,pk)
	
def LWE_sample( B, var_1, var_2): 
	# Returns one sample of a target sB +e
	# with s sample via ranom_vector on Z**m and e a N(0,var)
	m , n = B.nrows(), B.ncols()
	D = DiscreteGaussianDistributionLatticeSampler(B, var_1) #random_vector(ZZ , m  )
	s = D()
	e = normal(0 , var_2 , n )
	return(vector(RR, s + e)) #s * B + e))
	
def random_target_sample(B):
	# Returns a real vector sampled uniformly on the fundamental domain [0,1]**m * B 
	# for the lattice with basis B
	m , n = B.nrows(), B.ncols()
	U = vector(RR , uniform(0,1, m))
	return(vector(RR, U * B))	
	
def Sample_Dual_hypercube(B,k,T):
	# Returns about T**k vectors sampled as v * B_k where 
	# v in hypercube [| -T/2 , T/2 |] * k  
	# B_k is the matrix of the first k vectors of B 
	B_k = B[ 0:k ]
	sample = []
	for v in hypercube(T//2,k) :
			sample.append( vector(ZZ ,v ) * B_k )
	return(sample)

def Sample_Dual_random_walk(B,k,T,R):
	m , n = B.nrows() , B.ncols()
	w = vector(ZZ, np.zeros(n))
	sample = [w]
	i , j = 0 , 0
	while (i < T) and ( j < 1000 * T ):
		j += 1
		t , b = randint(0, k-1) , ZZ( 2*binomial(1,0.5)- 1)
		w_succ = w + b * vector(ZZ, B[t])
		if (w_succ.norm() < R**2):
			sample.append( w_succ )
			w = w_succ 
			i += 1
	return(sample)

def Sample_AR_Distinguisher(target , Bv , W , k ):
	# Inputs : 
	# target list of target vectors in RR**n , Bv in RR**(n,n) , W list of dual vectors , 
	# k int number of vectors in Bv, T int number of experiments  
	# Returns : sample of size T of the AR-statistic sum cos(2 * pi * < w , t>))
	# with w sampled as x * Bv_k where Bv_k is the k first vectors of the dual basis Bv
	# and x ranges in hypercube of dimension k and size a/2 i.e.
	# [ (1 / W ) * sum_{w in W} cos( 2 * pi * < w , t > ) ]_{t in target} 
	Bv_k = Bv[0:k]
	X = []
	for t in target:
		cos = []
		for w in W:
			cos.append( np.cos( 2 * np.pi * w.inner_product(vector(RR, t)) ) )
		X.append(np.mean(cos))
	return(X)
	

