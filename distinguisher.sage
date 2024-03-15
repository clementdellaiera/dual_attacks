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
	
def character(w,t):
	# Given w in dual lattice and t in vector space 
	# Returns exp(2 * i * pi * < w , t >) where the circle is identified as (-0.5 , 0.5)
	# via exp( 2* i * pi * x ) -> floor( x- 0.5 ) + 0.5
	# Lvee x V -> (-0.5 , 0.5) = RR / ZZ
	inner_product = w.inner_product( vector(RR, t) )
	return(  ((inner_product + 0.5 ) % 1) + 0.5 )
	
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
	data_lwe.append( character(w, target_lwe   ))
	data_unif.append(character(w, target_unif ))
	
################
##   GRAPHS   ##
################

# Set matplotlib and seaborn plotting style : might need to install seaborn in sage shell by running pip install seaborn
import seaborn as sns
sns.set_style('darkgrid')

data_lwe.sort() , data_unif.sort()
yy = (1/len(WW)) * np.linspace(0,1,len(WW))

figure, axis = plt.subplots(2) 

axis[0].plot(data_lwe , yy , color = "blue", label="LWE" )
axis[0].plot(data_unif, yy , color = "cyan", label="Uniform" )
axis[0].set_title("Empirical c.d.f. of < w , t > for w in W")
hh = len(data_lwe) //100
#sns.lineplot(x = np.array(data_lwe) , y=np.array(yy) , ax = axis[0])
sns.displot(np.array(data_lwe), bins = hh , ax = axis[1] , color="blue")
sns.displot(np.array(data_unif), bins = hh , ax = axis[1] , color = "cyan")

#axis[1].scatter(0*data_lwe , data_lwe , s = 1, marker='o', color = "blue")
#axis[1].scatter(0*data_unif, data_unif , s = 1, marker='o', color = "cyan")
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

#f, ax = plt.subplots(2,1, figsize=(12,6), sharex=True, sharey=True) 
#sns.lineplot(data_lwe, yy)#, hue="region", style="event")#, data=fmri)

#k1 = sns.distplot(np.array(data_lwe) , bins=20 )#, ax=ax[0][0] ) 
#k2 = sns.distplot(np.array(data_unif), bins=20 )#, ax=ax[1][0] )
