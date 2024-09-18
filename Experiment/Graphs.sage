import numpy as np
from matplotlib import pyplot as plt
# Set matplotlib and seaborn plotting style : might need to install seaborn in sage shell by running pip install seaborn
import seaborn as sns
sns.set_style('darkgrid')

################
##   GRAPHS   ##
################

data_LWE , data_unif = load('Data/distinguisher_values')
N_LWE , N_unif = len(data_LWE) , len(data_unif)

data_LWE.sort() , data_unif.sort()
yy_LWE , yy_unif = (1/N_LWE) * np.linspace(0,1,N_LWE) , (1/N_unif) * np.linspace(0,1,N_unif)

# figure, axis = plt.subplots(2) 

plt.plot(data_LWE , yy_LWE , color = "blue", label="LWE" )
plt.plot(data_unif, yy_unif , color = "cyan", label="Uniform" )
plt.title("Empirical c.d.f. of distinguisher values")
plt.savefig('Figures/cdf_distinguisher_LWE.png', bbox_inches='tight')

