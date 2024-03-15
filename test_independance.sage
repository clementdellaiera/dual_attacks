load("dual_attack.sage")
import pandas as pd
# To do : time package _ time.perf() , seminaire , bernstein

####################
##   PARAMETERS   ##
####################

mm = 10				# rank of lattice
kk = 5				# k first vectors in BKZ reduction of dual basis
aa = 5				# Size of hypercube
BB = random_matrix(ZZ,mm) 	# Basis of lattice

D , Bv, Lv = Dual(BB)
Bv_BKZ = Bv.BKZ()

WW = Sample_Dual_hypercube(Bv_BKZ,kk,aa)

P = DiscreteGaussianDistributionLatticeSampler(BB, 0.01) # Sampler discrete gaussian on lattice
sigma = 0.00000000001  					 # Variance of real gaussian used for error in LWE distribution

def Summary(X):
	print("Min , Max      : "+str(min(X)) +" , "+str(max(X)))
	print("Mean, Variance : "+str(np.mean(X)) + " , "+str(np.var(X)))
	print("Quantiles : "+str(np.quantile(X,0.25))+" , "+ str(np.median(X))+" , "+str(np.quantile(X,0.75)))
	
###############################
##   DISTRIBUTION OF <w,t>   ##
###############################

print("**   START   **")

start_time_perf , start_time_proc = time.perf_counter() , time.process_time()

ee = vector(RR, normal( 0 , sigma , mm )) * BB
target_lwe , target_unif = P() + ee , vector(RR, uniform(0,1,mm)) * BB 

data_lwe , data_unif = [] , []
data_cos_lwe , data_cos_unif = [] , []

for w in WW:
	data_lwe.append(  w.inner_product(vector( RR, target_lwe  )))
	data_unif.append( w.inner_product(vector( RR, target_unif )))
	data_cos_lwe.append(  np.cos( 2 * np.pi * w.inner_product(vector( RR, target_lwe  ))))
	data_cos_unif.append( np.cos( 2 * np.pi * w.inner_product(vector( RR, target_unif ))))
	
################
##   GRAPHS   ##
################

data_lwe.sort() , data_unif.sort()
yy = (1/len(WW)) * np.linspace(0,1,len(WW))

figure, axis = plt.subplots(2, 2) 

axis[0,0].plot(data_lwe , yy , color = "blue", label="LWE" )
axis[0,0].plot(data_unif, yy , color = "cyan", label="Uniform" )
axis[0,1].plot(data_lwe , yy , color = "blue", label="LWE" )
axis[0,0].set_title("Empirical c.d.f. of < w , t > for w in W")
axis[1,0].scatter(np.cos(data_lwe)  , np.sin(data_lwe) )
axis[1,1].scatter(np.cos(data_unif) , np.sin(data_unif) )
#plt.legend()
#plt.xlabel("LWE vs uniform")
plt.show()

print("Time : "+str(time.perf_counter()-start_time_perf) + " sec")
print("Processing  : "+str(time.perf_counter()-start_time_proc) + " ms")
print("")
print("LWE")
Summary(data_lwe)
print("")
print("Unif")
Summary(data_unif)
print("**   END   **")

######################
##   INDEPENDANCE   ##
######################

# Calcul de la statistique max_{u,v} | F_n_(X,Y)(u,v) - F_n_X(u) * F_n_Y(v) | pour construire un test d'indépendance non paramétrique 
# Sur une tirage de cos(2 * pi * <w,t> ) où w est choisi sur deux valeurs (ici les deux premiers vecteurs de la base du dual reduite par BKZ
# Tentative d'afficher les fdr empiriques à 2 variables
 
T = 50
data_lwe_0 , data_lwe_1 , data_unif_0 , data_unif_1 = [] , [] , [] , []

w0 , w1 = vector(ZZ , Bv_BKZ[1]) , vector( ZZ , Bv_BKZ[2] )

for i in range(T):
	ee = vector(RR, normal( 0 , sigma , mm )) * BB
	target_lwe , target_unif = P() + ee  , vector(RR, uniform(0,1,mm)) * BB 
	data_lwe_0.append(np.cos( 2 * np.pi * w0.inner_product(vector( RR, target_lwe  )))) , data_unif_0.append(np.cos( 2 * np.pi * w0.inner_product(vector( RR, target_lwe  ))))
	data_lwe_1.append(np.cos( 2 * np.pi * w1.inner_product(vector( RR, target_lwe  )) )) , data_unif_1.append(np.cos( 2 * np.pi * w1.inner_product(vector( RR, target_lwe  ))))

T_max , h = 0 , 100
x_max , y_max = max(data_lwe_0) , max(data_lwe_1)
for k in range(h):
	for l in range(h):
		u , v = x_max * k / h , y_max * l / h
		F_X , F_Y = np.mean( 1*(np.array(data_lwe_0) < u ) ) , np.mean( 1*(np.array(data_lwe_1) < v ) )	
		F_X_Y = np.mean( 1 * (np.array(data_lwe_0) < u ) * ( np.array(data_lwe_1) < v ) ) 
		T = abs(F_X_Y - F_X * F_Y )
		if ( T > T_max ) :
			T_max = T 

print("Kolmogorov independance statistics : "+str(T_max)) 

# Graph of Empirical CDF

z = (1/T) * np.linspace(0,1,T)

xy_plane , cdf_X_Y , cdf_X_times_cdf_Y = matrix.zero(T) , matrix.zero(T) , matrix.zero(T)

for k in range(T):
	for l in range(T):
		#xy_plane[k,l] = (data_lwe_0[k] , data_lwe_1[l]) 
		cdf_X_Y[k,l] = min( k , l ) #  * np.linspace(0,1,T)
		cdf_X_times_cdf_Y[k,l] = k * l 
cdf_X_Y , cdf_X_times_cdf_Y = (1/T) * cdf_X_Y , (1/T**2) * cdf_X_times_cdf_Y



























