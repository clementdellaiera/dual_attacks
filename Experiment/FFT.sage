import numpy as np
from matplotlib import pyplot as plt

modulus = 937

points = [k for k in range(modulus)]
data = [np.cos( 2*np.pi * k / modulus ) for k in range(modulus)]
data_fft = np.fft.fft(data)

fix , axs = plt.subplots(2)
axs[0].plot(points , data , 'ro', markersize = 0.1 )
axs[1].plot(points , data_fft , 'ro' , markersize = 0.1 ) 

plt.show()

data = [ 1 for k in range(modulus)]
data_fft = np.fft.fft(data)

def score(t):
	values = [np.cos(2*np.pi * w * t / modulus) for w in range(modulus)]
	return(sum(values))
target = [score(k) for k in range(modulus)]

fix , axs = plt.subplots(2)
axs[0].plot(points , target , 'ro', markersize = 0.1 )
axs[1].plot(points , data_fft , 'ro' , markersize = 0.1 ) 

plt.show()

"""
In dimension 1
Say q = 2**n
a not invertible mod q , a = 2**k
L 	= { x in ZZ | x = a y mod q }
Lperp 	= { w in ZZ | a * w = 0 mod q } = 2**(n-k)ZZ

Say q prime
a random mod q
L 	= { ( x_0 , x_1 ) in ZZ^2 | x_0 = y_0 + a * y_1 mod q , x_1 = y_1 mod q }
Lperp 	=  {( w_0 , w_1 ) in ZZ^2 | w_0  ,  - a* w_0 + w_1 = 0 mod q }
"""
