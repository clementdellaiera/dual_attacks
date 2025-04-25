# Libraries

from matplotlib import pyplot as plt
import seaborn as sns
sns.set_style('darkgrid')

import numpy as np

from sage.stats.distributions.discrete_gaussian_integer import DiscreteGaussianDistributionIntegerSampler as DiscreteGaussian

# Parameters 

q = 937
L = floor(np.log2(q))
N = pow(2,L-4) 

stdev = q // 10
D = DiscreteGaussian(stdev)

def nearest_integer(x_input):
	return(sign(x_input) * floor( abs(x_input) + 0.5 ))
	
def mod_plus(x_input , modulo):
	L = floor( (modulo - 1) / 2 )
	x_mod = x_input % modulo
	if x_mod >= L :
		return(x_mod - modulo)
	else : 	
		return(x_mod)

def modulo_switch( x_input , modulo_in , modulo_out ):
	x_scaled = x_input * modulo_out / modulo_in
	x_rounded = nearest_integer(x_scaled)
	return(x_rounded)

# Simulations 

sample_size = 10000

data_unif =  [ randint(- q // 2 , q //2) for i in range(sample_size) ]
data_gaussian = [ D() for i in range(sample_size) ]

data_unif_modulo_switched , data_gaussian_modulo_switched = [modulo_switch(x , q , N) for x in data_unif] , [modulo_switch(x , q , N) for x in data_gaussian]
 
data_unif.sort() , data_gaussian.sort()
data_unif_modulo_switched.sort() , data_gaussian_modulo_switched.sort() 

# Graphs
fig , axs = plt.subplots(1,2)
axs[0].plot(data_unif , [i for i in range(sample_size)])
axs[0].plot(data_gaussian , [i for i in range(sample_size)])
axs[1].plot(data_unif_modulo_switched , [i for i in range(sample_size)])
axs[1].plot(data_gaussian_modulo_switched , [i for i in range(sample_size)])

axs[0].set(xlabel= "modulo q = "+str(q) )
axs[1].set(xlabel= "modulo N = 2**"+str(L-4)+" = "+str(N) )

fig.suptitle('Modulo switching for uniform and gaussian distribution')

plt.savefig("Experiment/Figures/modulo_switching_q"+str(q)+"N"+str(N)+".png")





























