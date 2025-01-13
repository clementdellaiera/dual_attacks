from fpylll import IntegerMatrix, FPLLL
from fpylll.util import gaussian_heuristic
from g6k import Siever 
from time import time
import sys

from tqdm import tqdm
from multiprocessing import Pool

FPLLL.set_random_seed(1337)		# Seed set for reproducibility of experiments


# Parameters

# stdev, rank , modulus  = 1, 20 , 931  # stdev parameters for discrete gaussian error in LWE, dimension of lattice , modulus of LWE instances
rank , modulus  = Integer(sys.argv[1]), Integer(sys.argv[2])
log_covolume = rank // 2			# log_covolume = k = n-m | IntegerMatrix takes for qary lattices k in det L = q**k as parameters | basis of type [ [I_{n-k} , A ] , [0 , q I_k]] with size(A) = (n-k , k)
# T = 1000 				# Size of target samples


# Instance generation
B = IntegerMatrix.random( rank , 'qary' , k = log_covolume , q = modulus)
print("B : [OK]")

start = time()
## Basis of dual/perpendicular lattice
m = rank - log_covolume
A = matrix(ZZ,B)[0:log_covolume, m:rank]
"""
If B = [[I_m | A ] , [0 | q I_k]] , then
    B is a basis of Lq(A) 
    Bv = [[I_m | 0 ],[ -(1/q) * A | (1/q) * I_k ]] 	is a basis of dual(Lq(A))
    Bperp = [[q * I_m | 0 ] , [ - A | I_k ]]	is a basis of perp(Lq(A)) 
"""
B_perp = block_matrix( [ [modulus * matrix.identity(ZZ,m) , matrix.zero(ZZ ,m , log_covolume )] , [ - A.transpose() , matrix.identity(ZZ , log_covolume)] ] )
B_perp = IntegerMatrix.from_matrix(B_perp)

print("B_perp : [OK]")

###############################################################
##     Dual short vectors of dual/perp lattice via sieve     ##
##   Modification of g6k example file all_short_vectors.py   ##
###############################################################
short_vector_list = []	
	
g6k = Siever(B_perp)
g6k.lll(0,rank)
print("LLL(B_perp) : [OK]", flush=True)
print(g6k.params.saturation_radius,g6k.params.saturation_ratio,g6k.params.db_size_base, g6k.params.db_size_factor, flush=True)

# sat_radius = float(sys.argv[3])
# sat_ratio = .5
# with g6k.temp_params(saturation_ratio=sat_ratio, saturation_radius=sat_radius, db_size_base=sqrt(2), db_size_factor=5):
g6k.initialize_local(0, rank /2, rank)
while g6k.l > 0:
    _start = time()
    # Extend the lift context to the left
    g6k.extend_left(1)
    # Sieve
    # g6k()
    # g6k(alg="gauss" if g6k.n<45 else "hk3")
    g6k(alg="gauss")
    # g6k(alg="hk3")
    _duration = time() - _start
    print("Progressive sieving : dim={} in {:.4f} sec".format(g6k.n, _duration), flush=True)


sat_radius = float(sys.argv[3])
sat_ratio = .8
with g6k.temp_params(saturation_ratio=sat_ratio, saturation_radius=sat_radius, 
                     db_size_base=sqrt(2), db_size_factor=5):
    _start = time()
    # g6k(alg="gauss" if g6k.n<45 else "hk3")
    g6k(alg="gauss")
    _duration = time() - _start
    print("Final sieving : dim={} in {:.4f} sec".format(g6k.n, _duration))

# Convert all data_base vectors from basis A to cannonical basis and print them 
# out if they are indeed shorter than 1.7 * gh^2
print("Vecteurs mis sous forme de liste...", flush=True)
gh = gaussian_heuristic([g6k.M.get_r(i, i) for i in range(rank)])


from itertools import islice
# Define a function to yield chunks of a specified size from an iterable
def chunks(iterable, size):
    # Create an iterator from the input iterable
    iterator = iter(iterable)   
    # Loop over the iterator, taking the first element in each iteration
    for first in iterator:
        # Yield a list consisting of the first element and the next 'size-1' elements from the iterator
        yield [first] + list(islice(iterator, size - 1))

# # Exemple sur la liste suivante
# my_list = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
# # Convert the generator returned by the chunks function into a list of chunks
# chunked_list = list(chunks(my_list, 10))
# print(chunked_list)
## on trouve : [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10], [11, 12, 13, 14, 15]]

chunks_size = 500000
data_base = list(chunks(g6k.itervalues(), chunks_size))
print('#database = {}'.format(chunks_size*(len(data_base)-1) + len(data_base[-1])), flush=True)

def sauvegarde(x):
    res = []
    v = B_perp.multiply_left(x)
    l = sum(v_**2 for v_ in v)
    if l < sat_radius * gh:
        res.append(v)
    return res


print("Sauvegarde des vecteurs retenus...", flush=True)
_start = time()
# chunked_db = list(chunks(data_base, chunks_size))
nb_cpu = int(sys.argv[4])
res = []

for db in tqdm(data_base):
    pool = Pool(nb_cpu)
    short_vector_list = pool.map(sauvegarde, db)
    short_vector_list = [ent for sublist in short_vector_list for ent in sublist]
    pool.close()
    res += short_vector_list


file_name = '/local/opt/tuonguye/Dual_short_vectors/dual_short_vectors_{}_{}.txt'.format(sys.argv[1], sys.argv[2])
with open(file_name, "w") as f:
    f.write(str(res))
_duration = time() - _start
print("Ecriture dans le fichier de sauvegarde : [OK] en {:.4f} sec".format(_duration), flush=True)


print("Found %d vectors of squared length less than %.2f*gh. (expected %f)"%(len(res), sat_radius, sat_ratio * .5 * sat_radius**(rank /2.)))
print('Execution time : {:.4f} sec'.format(time() - start))





# save((modulus ,short_vector_list) , 'Dual_short_vectors/dual_short_vectors_{}_{}'.format(sys.argv[1], sys.argv[2]))
