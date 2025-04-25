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
