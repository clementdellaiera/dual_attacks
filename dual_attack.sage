##################
##   PACKAGES   ##
##################

from sage.all_cmdline import * 
import time

import numpy as np
from numpy.random import normal, uniform

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
			sample.append( vector(ZZ ,v ) * Bv_k )
	return(sample)

def Sample_AR_Distinguisher(target , Bv , W , k ):
	# target list of target vectors in RR**n , Bv in RR**(n,n) , W list of dual vectors , 
	# k int number of vectors in Bv, T int number of experiments  
	# Returns a sample of size T of the AR-statistic sum cos(2 * pi * < w , t>))
	# with w sampled as x * Bv_k where Bv_k is the k first vectors of the dual basis Bv
	# and x ranges in hypercube of dimension k and size a/2 
	Bv_k = Bv[0:k]
	X = []
	for t in target:
		cos = []
		for w in W:
			cos.append( np.cos( 2 * np.pi * w.inner_product(vector(RR, t)) ) )
		X.append(np.mean(cos))
	return(X)
			
##############	
##   TEST   ##
##############

start_time = time.time()
total_time = 0

mm = 10			# rank of lattice
kk = 5				# k first vectors in BKZ reduction of dual basis
aa = 5				# Size of hypercube
BB = random_matrix(ZZ,mm) 	# Basis of lattice

m , n = BB.nrows(), BB.ncols()

total_time = time.time() - start_time
start_time = time.time()

D , Bv, Lv = Dual(BB)
Bv_BKZ = Bv.BKZ()
Bv_k = Bv_BKZ[ 0:kk ]

BKZ_time =  time.time() - start_time
total_time += BKZ_time
start_time = time.time()

WW = Sample_Dual_hypercube(Bv_BKZ,kk,aa)

WW_time =  time.time() - start_time
total_time += WW_time
start_time = time.time()

D = DiscreteGaussianDistributionLatticeSampler(BB, 0.01) # Sampler discrete gaussian on lattice
sigma = 0.00000000001  					 # Variance of real gaussian used for error in LWE distribution  

print("** START **")
print(" ")
print("Random lattice of rank "+str(mm))
print("BKZ reduction of dual basis, keeping "+str(kk)+" short vectors")
print("Generating "+str(len(WW))+" short vectors in dual lattice in W")
print("W : image of [-"+str(aa//2)+" , "+str(aa//2)+"]^"+str(kk)+ " by Bv[1:"+str(kk) +"]")

# Profile of dual basis
size = []
for i in range(mm):
	size.append(Bv_BKZ[i].norm())
plt.plot(np.linspace(1,mm,mm) ,size)
plt.title("Profile of reduced dual basis")
plt.ylabel("Norm")
plt.xlabel("Ordered set of vectors")
plt.show()	

# EXPERIMENT 1   
# On a sample of T iid targets, compute the CDF of the AR-statistics (W-mean of of cos(2 * pi * <w,t>) ) and draws the CDF's. One sample is draw according to the LWE distribution, the second one according to the image of the uniform on hypercube by B (= basis)  
# For each step, we sample t_lwe = s + e and t_unif = u * B 
# where s, e , u  = gaussian(B, 10000), Normal(0,sigma) , uniform( [0,1]**m * B)
# and evaluate the AR statistics on the set W of short dual vectors
# given by the image of the k-th hypercube by the BKZ reduction of B

T = 50
target_lwe = []
target_unif = []

total_time += time.time() - start_time
start_time = time.time()

for i in range(T):
	ee = vector(RR, normal(0,sigma,m)) * BB #vector(RR, uniform(0,,m))
	target_lwe.append( D()+ee ) 
	target_unif.append(vector(RR, uniform(0,1,mm)) * BB) 
		
data_lwe , data_unif = Sample_AR_Distinguisher(target_lwe , Bv_BKZ , WW , kk) , Sample_AR_Distinguisher(target_unif , Bv_BKZ , WW , kk)

EXP1_time =  time.time() - start_time
total_time += EXP1_time
start_time = time.time()


data_lwe.sort() , data_unif.sort()
L1 , L2 = len(data_lwe) , len(data_unif)
plt.plot(data_lwe, (1/L1) * np.linspace(0,1,L1) , color = "blue", label="LWE",drawstyle="steps")
plt.plot(data_unif, (1/L2) * np.linspace(0,1,L2) , color = "cyan",label="Uniform", drawstyle="steps")
plt.title("Empirical cumulative distribution for AR-distinguisher")
plt.xlabel("LWE vs uniform")
plt.legend()

plt.show()

total_time += time.time() - start_time
start_time = time.time()
	
# EXPERIMENT 2 : 
# CDF of the sample (cos(2 * pi * <w,t> ))_{w in W} for different distribution of t
# t_lwe = s + e , where s is sampled on lattice, e is a centered gaussian of variance sigma
# t_unif = uniform on fundamental domain [0,1]**m * B
ss = D() 				# gaussian vector in lattice 
ee = vector(RR, uniform(0,sigma,n)) 	# vector(RR, normal(0,sigma,n)) * BB# gaussian in RR**n
uu = vector(RR, uniform(0,1,mm)) * BB	# uniform random vector in fundamental domain

success_lwe , success_unif = 0 , 0
cos_lwe, cos_unif = [] , []

for ww in WW :
	if (ww.inner_product(ss) in ZZ):
		success_lwe += 1
	if (ww.inner_product(uu) in ZZ):
		success_unif += 1 
	cos_lwe.append( np.cos( 2 * np.pi * ww.inner_product(ss+ee))) # + ee)))
	cos_unif.append( np.cos( 2 * np.pi * ww.inner_product(uu) ) )

print("Test <w,s> in ZZ for w in W for s LWE / uniform : "+str(round(100 * success_lwe / len(WW),2) ) + " % / "+str(round(100 * success_unif / len(WW) , 2))  +" %")
print("Mean  of cos(2 * pi * < w , s > ) when s is LWE / uniform : "+str(np.mean(cos_lwe)) + " / " + str(np.mean(cos_unif)) )

EXP2_time =  time.time() - start_time
total_time += EXP2_time
start_time = time.time()

cos_lwe.sort() , cos_unif.sort()
plt.plot(cos_lwe, (1/len(WW)) * np.linspace(0,1,len(WW)) , color = "blue",label="LWE")
plt.plot(cos_unif, (1/len(WW)) * np.linspace(0,1,len(WW)) , color = "cyan",label="Uniform")
plt.title("Empirical cumulative distribution function of cos(2 * pi * < w , s >) for w in W")
plt.legend()
plt.xlabel("LWE vs uniform")
plt.show()
		
#print("DISTINGUISHER test")

# Experiment 3 : 
# For each step, we sample t_lwe = s + e and t_unif = u * B 
# where s, e , u  = gaussian(B, 10000), Normal(0,sigma) , uniform( [0,1]**m * B)
# and evaluate the AR statistics on the set W of short dual vectors
# given by the image of the k-th hypercube by the BKZ reduction of B

#T = 100
#X_lwe  , X_unif = [] , []
#for i in range(T):
#	#tt, tt_uniform = LWE_sample(BB,10000,0.1) , random_target_sample(BB)
#	tt, tt_uniform = vector(RR, ss + normal(0 , sigma , n )) , random_target_sample(BB)
#	X_lwe.append(Distinguisher_Sieve_preprocessing(Bv_k , WW , tt, kk))
#	X_unif.append(Distinguisher_Sieve_preprocessing(Bv_k , WW , tt_uniform, kk)) 
#X_lwe.sort()
#X_unif.sort()

#y = (1/T)*np.linspace(0,1,T)

#plt.plot(X_lwe,y,color="blue",label="LWE")
#plt.plot(X_unif,y,color="black",label="Random target")

#plt.title("Distinguisher for LWE and random sampling")
#plt.xlabel("Empirical cumulative distribution functions")
#plt.legend()
#plt.show()
total_time += time.time() - start_time
print(" ")
print("** EXECUTION TIME **")
print(" ")
print("BKZ reduction :              "+str(BKZ_time)+" s ")
print("Short dual vector sampling : "+str(WW_time)+" s ")
print("Short dual vector sampling : "+str(BKZ_time)+" s ")
print("Experiment 1 :               "+str(EXP1_time)+" s ")
print("Experiment 2 :               "+str(EXP2_time)+" s ")
print("TOTAL :                      "+str(total_time)+" s ")
print(" ")
print("** END **")
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
