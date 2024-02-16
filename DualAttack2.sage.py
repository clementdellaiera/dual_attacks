

# This file was *autogenerated* from the file DualAttack2.sage
from sage.all_cmdline import *   # import sage library

_sage_const_10 = Integer(10); _sage_const_37 = Integer(37); _sage_const_1 = Integer(1); _sage_const_0 = Integer(0); _sage_const_2 = Integer(2); _sage_const_20 = Integer(20); _sage_const_100 = Integer(100); _sage_const_0p1 = RealNumber('0.1'); _sage_const_4 = Integer(4); _sage_const_200 = Integer(200); _sage_const_0p00000001 = RealNumber('0.00000001')
from sage.all_cmdline import * 
import time
import numpy as np
from numpy.random import normal
from fpylll import *

from itertools import product

################
##   BASICS   ##
################

# Generating q-ary lattice of dimension 100 and determinant q**50, where q is a 30 bit prime
n=_sage_const_10 
#AA = IntegerMatrix.random(n , "qary", k = 50 , bits = 30)
## Generating SIS instance
q = _sage_const_37 
AA = random_matrix(Zmod(q),n)

###################
##   FUNCTIONS   ##
###################

def Dual(B): # given an integer matrix, returns det(B) , Bv , Lv : determinant, integer matrix such that (1/det(B)) * Bv is a basis of dual lattice, dual lattice
	D = det(B)
	BB = D * B.inverse().transpose()
	Bv = BB.change_ring(ZZ)
	return(D , Bv , (_sage_const_1 /D) * Bv.image() )

def hypercube( a , d ): 
	# Returns the list of all integer vectors inside the hypercube 
	# of dimension d and of side 2a+1 (i.e. the max-norm centered ball of Z*d of radius a) 
	L = []
	if (d < _sage_const_1 ) or (a < _sage_const_1 ):
		return(False)
	elif (d == _sage_const_1  ):
		for k in range(-a , a+_sage_const_1 ):
	 		L.append([k])
	else:  
		for v in hypercube(a , d-_sage_const_1 ):
			for k in range(-a,a+_sage_const_1 ):
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
	Bv_k = Bv_BKZ[ _sage_const_0 :k ]
	sample = []
	if (option == "Random_with_replacement"):
		for i in range(T**k):
			sample.append( Bv_k.image().random_element() ) 
			# Basis are maybe in row forms, ie if B = V.basis(), xB in V
	elif (option == "Sieve"):
		for v in hypercube(T//_sage_const_2 ,k) :
			sample.append( vector(ZZ ,v ) * Bv_k )
	elif (option == "Random_walk"):
		v = _sage_const_0 
		for i in range(T**k):
			u = randint(_sage_const_0 ,k-_sage_const_1 )
			v += Bv_k[u] 
			sample.append(v)	
	return(sample)

def Distinguisher_Sieve(B , t, k, N): 
	# B basis of lattice, X sample of target vector , k , N size of dual vector sample
	# Returns the value of the AR-distinguisher evaluated wrt a sample of short dual vectors
	# obtained via teh sieve method in sample dual
	D , Bv, Lv = Dual(B)
	Bv_BKZ = Bv.BKZ()
	Bv_k = Bv_BKZ[ _sage_const_0 :k ]
	T = _sage_const_0 
	C = _sage_const_0 
	for v in hypercube(N//_sage_const_2 ,k) :
			w = vector(ZZ ,v ) * Bv_k
			T += np.cos(_sage_const_2  * np.pi * w.inner_product(t))
			C += _sage_const_1 
	return((_sage_const_1 /C)*T)   

	
def LWE(n): 
	# Returns one sample of a LWE distribution (A, As+e) where A in ZZ**(n,n)
	A = random_matrix(ZZ,n,n)
	s = random_vector(ZZ , n  )
	e = random_vector(ZZ , n  )
	sk , pk = s , (A,A*s+e)
	return(sk,pk)
	
def LWE_sample(B,var): 
	# Returns one sample of a target sB +e
	# with s sample via ranom_vector on Z**m and e a N(0,var)
	m , n = B.nrows(), B.ncols()
	s = random_vector(ZZ , m  )
	e = normal(_sage_const_0  , var , n )
	return(vector(RR, s * B + e))
	
def Test(X, epsilon):
	T = (_sage_const_1 /len(XX))* sum(np.cos(_sage_const_2  * np.pi * X)) # Test statistic
	return(T < epsilon )

###############
##   TESTS   ## 
###############

# Parameters

n=_sage_const_10 

# LATTICES
print("** LATTICES **")

A = random_matrix(ZZ,n,n)
#s, e = random_vector(ZZ , n ), random_vector(ZZ , n ) 
#t = A*s+e

print("")
print("** BKZ **")

print("norm of first vector of A")
print(A[_sage_const_1 ].norm())

print("BKZ reduction of A")
print("Norm of the first vector of reduced basis : "+str(A.BKZ()[_sage_const_1 ].norm()))

print("")
print("**  DUAL  **") 

# IntegralLattice returns the lattice generated by a quadratic form over ZZ

LL = IntegralLattice(A.transpose() * A) # for integal lattice, there is a direct method
LLv = LL.dual_lattice()
#print(LL.parent())
#print(LL.dual_lattice().parent())

L = A.image()
Av = A * (A.transpose()*A).inverse()
Lv = Av.image() # problem : this isconsidered a rationnal vector space, not a lattice
w , v = Lv.random_element() , L.random_element()
print(" (w ,v) in ZZ : "+str( w.dot_product(v) in ZZ))

dd = det(A)
B = dd * A.inverse().transpose()
B = B.change_ring(ZZ)
dLv = B.image()
ww , v = (_sage_const_1 /dd) * dLv.random_element() , L.random_element()
print(" (w ,v) in ZZ : "+str( ww.dot_product(v) in ZZ))
vol , Av, Lv = Dual(A)
print("Dual Test : "+str(Lv == (_sage_const_1 /dd) * dLv))
	
ss, (AA,bb) = LWE(_sage_const_20 )	

N = _sage_const_100 	
XX = np.random.normal(_sage_const_0 ,_sage_const_1 ,N)

print("Distinguisher TEST : "+str(Test(XX , _sage_const_0p1 )))

print("")
mm = _sage_const_10  			# rank of lattice
kk = _sage_const_4  				# k first vectors in BKZ reduction of dual basis
TT = _sage_const_10 				# scale of the dual sample
BB = random_matrix(ZZ,mm) 	# Basis of lattice

################################
##    Test for Dual Sampler   ##
################################

print("SAMPLE Test : sample of "+str(TT)+" vectors over dual lattice of rank " + str(mm))

print("Random independant draw with replacement")
SAMPLE = Sample_Dual(BB, kk , TT, "Random_with_replacement")
success = _sage_const_0 
for w in SAMPLE:
	v = BB.image().random_element()
	if (w.dot_product(v) in ZZ):
		success += _sage_const_1 
print("Test (w,v) in ZZ : "+str(_sage_const_100  * success / len(SAMPLE)) + " % success for "+str(len(SAMPLE))+" dual vectors")

print("Random Walk")
SAMPLE = Sample_Dual(BB, kk , TT, "Random_walk")
success = _sage_const_0 
for w in SAMPLE:
	v = BB.image().random_element()
	if (w.dot_product(v) in ZZ):
		success += _sage_const_1 
print("Test (w,v) in ZZ : "+str(_sage_const_100  * success / len(SAMPLE)) + " % success for "+str(len(SAMPLE))+" dual vectors")

print("Sieve")
SAMPLE = Sample_Dual(BB, kk , TT, "Sieve")
success = _sage_const_0 
for w in SAMPLE:
	v = BB.image().random_element()
	if (w.dot_product(v) in ZZ):
		success += _sage_const_1 
print("Test (w,v) in ZZ : "+str(_sage_const_100  * success / len(SAMPLE)) + " % success for "+str(len(SAMPLE))+" dual vectors")
print("")

###############################
##   Test for distinguiser   ##
###############################
print("")
print("DISTINGUISHER test")
tt = LWE_sample(BB,_sage_const_0p1 )
print(Distinguisher_Sieve(BB ,tt, kk,TT))

tt_uniform = vector(RR, normal(_sage_const_0 ,_sage_const_200 ,mm))
print(Distinguisher_Sieve(BB ,tt_uniform, kk,TT))
##################################
##   Hadamard Walsh Transform   ##
##################################

from sympy import fwht
from numpy import cos, linspace
m=_sage_const_10 
t = linspace(_sage_const_0 ,_sage_const_1 ,_sage_const_2 **m)
f = cos(t)
Ff = fwht(f)
e = _sage_const_0p00000001 
print("FAST HADAMARD WALSH TRANSFORM")
print("Check that F**2 = id : "+str(sum(_sage_const_1 *(abs(fwht(Ff)-_sage_const_2 **m * f) < e)) == _sage_const_2 **m ))


