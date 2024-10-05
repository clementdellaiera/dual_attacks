##   Load lattice basis 
load('Parameters.sage')
[ is_LWE, stdev , rank , modulus , log_covolume ] , B = load('Lattice_basis/basis_LWE.sobj')

from sage.stats.distributions.discrete_gaussian_integer import DiscreteGaussianDistributionIntegerSampler \
    as DiscreteGaussian

dim_guess = 4
stdev = 0.01
D = DiscreteGaussian(stdev)

if is_LWE :
	BB , m = matrix(ZZ,B) , rank - log_covolume 	# m = n - k 
	A  = BB[0:log_covolume, m:rank]			# A in ZZ^{ m x k } such that B = qary Lattice of LWE (A)
	R  = Zmod(modulus) 				# ZZ mod q
		
	A_guess = random_matrix(R, m, dim_guess )
	B_total = block_matrix( ZZ , [[ matrix.identity(ZZ,m) , A , A_guess ] , [ matrix.zero(ZZ, log_covolume , m) 	, modulus * matrix.identity(ZZ,log_covolume),  matrix.zero(ZZ, log_covolume , dim_guess) ] , [ matrix.zero(ZZ, dim_guess , m), matrix.zero(ZZ, dim_guess , log_covolume )	, modulus*matrix.identity(ZZ, dim_guess ) ]])
	
	s , e = random_matrix(R ,  rank + dim_guess ,1) , matrix(ZZ,rank + dim_guess ,1, [D() for _ in range( rank + dim_guess )]) 
	s_guess = s[ rank:]
	target = B_total * s + e
	
	target_list = []			
	for i in range(T): 
		s_try = random_matrix( R , dim_guess , 1)
		t_try = block_matrix([ [A_guess * s_try ] , [ matrix.zero(R , log_covolume , 1) ] , [ modulus * s_guess ] ] )
		target_list.append( target - t_try) 
			
save( ( s_guess , target , target_list) , 'Target_samples/target_samples_attack' )		
