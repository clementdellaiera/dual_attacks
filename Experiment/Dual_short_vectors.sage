from fpylll import IntegerMatrix
from fpylll.util import gaussian_heuristic
from g6k import Siever 
from time import time

start = time()
##   Load lattice basis 
[ is_LWE , stdev , rank , modulus , log_covolume ] , B = load('Lattice_basis/basis_LWE.sobj')

## Basis of dual/perpendicular lattice
if is_LWE :
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
else :
	assert False , "Code not done yet"
	det_B , inverse_matrix = det(B) ,  
	Bv = B.inverse().transpose()

###############################################################
##     Dual short vectors of dual/perp lattice via sieve     ##
##   Modification of g6k example file all_short_vectors.py   ##
###############################################################
short_vector_list = []	
	
g6k = Siever(B_perp)
g6k.lll(0,rank)

g6k.initialize_local(0, rank /2, rank)
while g6k.l > 0:
    # Extend the lift context to the left
    g6k.extend_left(1)
    # Sieve
    g6k()

with g6k.temp_params(saturation_ratio=.95, saturation_radius=1.7, 
                     db_size_base=sqrt(1.7), db_size_factor=5):
    g6k()

# Convert all data_base vectors from basis A to cannonical basis and print them 
# out if they are indeed shorter than 1.7 * gh^2

gh = gaussian_heuristic([g6k.M.get_r(i, i) for i in range(rank)])

data_base = list(g6k.itervalues())
found = 0

for x in data_base:
    v = B_perp.multiply_left(x)
    l = sum(v_**2 for v_ in v)
    if l < 1.7 * gh:
        # print(l/gh, v)
        found += 1
        short_vector_list.append(v)

print("Found %d vectors of squared length than 1.7*gh. (expected %f)"%(found, .5 * 1.7**(rank /2.)))
print('Execution time : ', time() - start , ' sec.')

save((modulus ,short_vector_list) , 'Dual_short_vectors/dual_short_vectors')
