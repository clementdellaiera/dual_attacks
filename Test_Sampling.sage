load("dual_attack.sage")

##############	
##   TEST   ##
##############
print("** START **")
print(" ")

start_time = time.time()
total_time = 0

mm = 10			# rank of lattice
kk = 5				# k first vectors in BKZ reduction of dual basis
aa = 5				# Size of hypercube
BB = random_matrix(ZZ,mm) 	# Basis of lattice

m , n = BB.nrows(), BB.ncols()

print("Random lattice of rank "+str(mm))
print("BKZ reduction of dual basis, keeping "+str(kk)+" short vectors")

total_time = time.time() - start_time
start_time = time.time()

D , Bv, Lv = Dual(BB)
Bv_BKZ = Bv.BKZ()
Bv_k = Bv_BKZ[ 0:kk ]

BKZ_time =  time.time() - start_time
total_time += BKZ_time
start_time = time.time()

# Sampling of short dual vectors in the dual
print(" ") 
print("**  Generating short vectors in dual lattice in W  **")

WW = Sample_Dual_hypercube(Bv_BKZ,kk,aa)
print("Hypercube - image of [-"+str(aa//2)+" , "+str(aa//2)+"]^"+str(kk)+ " by Bv[1:"+str(kk) +"] : "+str(len(WW))+" vectors")

WW_time =  time.time() - start_time
total_time += WW_time
start_time = time.time()

TT = 100
R = 100 * aa 
WW_random_walk = Sample_Dual_random_walk(Bv_BKZ,kk, TT , R)

print("Random walk - uniform step amongst "+str(kk)+" first vectors of Bv.BKZ() :       "+str(len(WW_random_walk))+" vectors")

WW_rw_time =  time.time() - start_time
total_time += WW_rw_time
start_time = time.time()

D = DiscreteGaussianDistributionLatticeSampler(BB, 0.01) # Sampler discrete gaussian on lattice
sigma = 0.00000000001  					 # Variance of real gaussian used for error in LWE distribution  

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

data_lwe_rw , data_unif_rw = Sample_AR_Distinguisher(target_lwe , Bv_BKZ , WW_random_walk , kk) , Sample_AR_Distinguisher(target_unif , Bv_BKZ ,  WW_random_walk , kk)

EXP1_time =  time.time() - start_time
total_time += EXP1_time
start_time = time.time()


data_lwe.sort() , data_unif.sort()
data_lwe_rw.sort() , data_unif_rw.sort()

L1 , L2 = len(data_lwe) , len(data_unif)
plt.plot(data_lwe, (1/L1) * np.linspace(0,1,L1) , color = "blue", label="LWE",drawstyle="steps")
plt.plot(data_unif, (1/L2) * np.linspace(0,1,L2) , color = "cyan",label="Uniform", drawstyle="steps")
plt.title("Empirical cumulative distribution for AR-distinguisher")
plt.xlabel("LWE vs uniform")
plt.legend()

plt.show()

L1 , L2 = len(data_lwe_rw) , len(data_unif_rw)
plt.plot(data_lwe_rw, (1/L1) * np.linspace(0,1,L1) , color = "blue", label="LWE",drawstyle="steps")
plt.plot(data_unif_rw, (1/L2) * np.linspace(0,1,L2) , color = "cyan",label="Uniform", drawstyle="steps")
plt.title("Empirical cumulative distribution for AR-distinguisher - Random walk")
plt.xlabel("LWE vs uniform - W sampled as a random walk")
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

print("")
print("W sampled from hypercube")
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

cos_lwe.sort() , cos_unif.sort()
plt.plot(cos_lwe, (1/len(WW)) * np.linspace(0,1,len(WW)) , color = "blue",label="LWE")
plt.plot(cos_unif, (1/len(WW)) * np.linspace(0,1,len(WW)) , color = "cyan",label="Uniform")
plt.title("Empirical cumulative distribution function of cos(2 * pi * < w , s >) for w in W")
plt.legend()
plt.xlabel("LWE vs uniform - W hypercube")
plt.show()

print("")
print("W sampled as random walk")
success_lwe , success_unif = 0 , 0
cos_lwe, cos_unif = [] , []

for ww in WW_random_walk :
	if (ww.inner_product(ss) in ZZ):
		success_lwe += 1
	if (ww.inner_product(uu) in ZZ):
		success_unif += 1 
	cos_lwe.append( np.cos( 2 * np.pi * ww.inner_product(ss+ee))) # + ee)))
	cos_unif.append( np.cos( 2 * np.pi * ww.inner_product(uu) ) )

print("Test <w,s> in ZZ for w in W for s LWE / uniform : "+str(round(100 * success_lwe / len(WW_random_walk),2) ) + " % / "+str(round(100 * success_unif / len(WW_random_walk) , 2))  +" %")
print("Mean  of cos(2 * pi * < w , s > ) when s is LWE / uniform : "+str(np.mean(cos_lwe)) + " / " + str(np.mean(cos_unif)) )

EXP2_time =  time.time() - start_time
total_time += EXP2_time
start_time = time.time()

cos_lwe.sort() , cos_unif.sort()
plt.plot(cos_lwe, (1/len(WW_random_walk)) * np.linspace(0,1,len(WW_random_walk)) , color = "blue",label="LWE")
plt.plot(cos_unif, (1/len(WW_random_walk)) * np.linspace(0,1,len(WW_random_walk)) , color = "cyan",label="Uniform")
plt.title("Empirical cumulative distribution function of cos(2 * pi * < w , s >) for w in W")
plt.legend()
plt.xlabel("LWE vs uniform - W random walk")
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
print("Short dual vector sampling (cube): "+str(WW_time)+" s ")
print("Short dual vector sampling (random walk): "+str(WW_rw_time)+" s ")
print("Experiment 1 :               "+str(EXP1_time)+" s ")
print("Experiment 2 :               "+str(EXP2_time)+" s ")
print("TOTAL :                      "+str(total_time)+" s ")
print(" ")
print("** END **")

