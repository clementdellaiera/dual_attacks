""" Modification of the script all_short_vectors.py from g6k github repository available at https://github.com/fplll/g6k"""
### In order to use 	: conda activate g6k
### Once done 		: conda deactivate

##################
##   Packages   ##
##################

from math import sqrt
from fpylll import IntegerMatrix
from fpylll.util import gaussian_heuristic


try:
    from g6k import Siever
except ImportError:
    raise ImportError("g6k not installed. Please run './setup.py install' from ../")

#from simple_pump import pump # Pourquoi cela ne fonctionne pas si décommenté ? Pourquoi le script s'éxécute quand même lorsque commenté ?

# Set up the instance

n = 30
A = IntegerMatrix.random(n, "qary", k=n/2, bits=8)
g6k = Siever(A)
g6k.lll(0, n)

g6k.initialize_local(0, n/2, n)
while g6k.l > 0:
    # Extend the lift context to the left
    g6k.extend_left(1)
    # Sieve
    g6k()

with g6k.temp_params(saturation_ratio=.95, saturation_radius=1.7, 
                     db_size_base=sqrt(1.7), db_size_factor=5):
    g6k()

# Convert all db vectors from basis A to cannonical basis and print them 
# out if they are indeed shorter than 1.7 * gh^2

gh = gaussian_heuristic([g6k.M.get_r(i, i) for i in range(n)])

db = list(g6k.itervalues())
found = 0

for x in db:
    v = A.multiply_left(x)
    l = sum(v_**2 for v_ in v)
    if l < 1.7 * gh:
        print(l/gh, v)
        found += 1

print("found %d vectors of squared length than 1.7*gh. (expected %f)"%(found, .5 * 1.7**(n/2.)))
