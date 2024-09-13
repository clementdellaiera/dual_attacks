

# This file was *autogenerated* from the file distinguisher.sage
from sage.all_cmdline import *   # import sage library

_sage_const_10 = Integer(10); _sage_const_5 = Integer(5); _sage_const_0p01 = RealNumber('0.01'); _sage_const_0p00000000001 = RealNumber('0.00000000001'); _sage_const_0p25 = RealNumber('0.25'); _sage_const_0p75 = RealNumber('0.75'); _sage_const_0p5 = RealNumber('0.5'); _sage_const_1 = Integer(1); _sage_const_0 = Integer(0); _sage_const_2 = Integer(2); _sage_const_100 = Integer(100)
load("dual_attack.sage")
import pandas as pd

# To do : time package _ time.perf() , seminaire , bernstein

####################
##   PARAMETERS   ##
####################

mm = _sage_const_10 				# rank of lattice
kk = _sage_const_5 				# k first vectors in BKZ reduction of dual basis
aa = _sage_const_5 				# Size of hypercube
BB = random_matrix(ZZ,mm) 	# Basis of lattice

D , Bv, Lv = Dual(BB)
Bv_BKZ = Bv.BKZ()

WW = Sample_Dual_hypercube(Bv_BKZ,kk,aa)

P = DiscreteGaussianDistributionLatticeSampler(BB, _sage_const_0p01 ) # Sampler discrete gaussian on lattice
sigma = _sage_const_0p00000000001   					 # Variance of real gaussian used for error in LWE distribution

def Summary(X):
	print("Min , Max      : "+str(min(X)) +" , "+str(max(X)))
	print("Mean, Variance : "+str(np.mean(X)) + " , "+str(np.var(X)))
	print("Quantiles : "+str(np.quantile(X,_sage_const_0p25 ))+" , "+ str(np.median(X))+" , "+str(np.quantile(X,_sage_const_0p75 )))
	
def character(w,t):
	# Given w in dual lattice and t in vector space 
	# Returns exp(2 * i * pi * < w , t >) where the circle is identified as (-0.5 , 0.5)
	# via exp( 2* i * pi * x ) -> floor( x- 0.5 ) + 0.5
	# Lvee x V -> (-0.5 , 0.5) = RR / ZZ
	inner_product = w.inner_product( vector(RR, t) )
	return(  ((inner_product + _sage_const_0p5  ) % _sage_const_1 ) + _sage_const_0p5  )
	
###############################
##   DISTRIBUTION OF <w,t>   ##
###############################

print("**   START   **")

start_time_perf , start_time_proc = time.perf_counter() , time.process_time()

ee = vector(RR, normal( _sage_const_0  , sigma , mm )) * BB
target_lwe , target_unif = P() + ee , vector(RR, uniform(_sage_const_0 ,_sage_const_1 ,mm)) * BB 

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
yy = (_sage_const_1 /len(WW)) * np.linspace(_sage_const_0 ,_sage_const_1 ,len(WW))

figure, axis = plt.subplots(_sage_const_2 ) 

axis[_sage_const_0 ].plot(data_lwe , yy , color = "blue", label="LWE" )
axis[_sage_const_0 ].plot(data_unif, yy , color = "cyan", label="Uniform" )
axis[_sage_const_0 ].set_title("Empirical c.d.f. of < w , t > for w in W")

#sns.lineplot(x = np.array(data_lwe) , y=np.array(yy) , ax = axis[0])
histogrammes , axis_2 = plt.subplots(_sage_const_2 ) 
hh = len(data_lwe) //_sage_const_100 
sns.displot(np.array(data_lwe), bins = hh , ax = axis_2[_sage_const_0 ] , color="blue")
sns.displot(np.array(data_unif), bins = hh , ax = axis_2[_sage_const_1 ] , color = "cyan")

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
