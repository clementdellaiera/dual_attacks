##   Load lattice basis 
load('Parameters.sage')
[ is_LWE, stdev , rank , modulus , log_covolume ] , B = load("Lattice_basis/basis_isLWE_True_stdev_0.500000000000000_rank_50_modulus_931_logcovolume25.sobj")
# B = load('Lattice_basis/basis_LWE_50.sobj')

from sage.stats.distributions.discrete_gaussian_integer import DiscreteGaussianDistributionIntegerSampler \
    as DiscreteGaussian

stdev = 0.01

D = DiscreteGaussian(stdev)
# e = matrix(ZZ, [D() for _ in range(rank)]).transpose()


if is_LWE :
	target_list_unif , target_list_LWE = [] , []
	R  = Zmod(modulus) 				# ZZ mod q
	m  = rank - log_covolume			# m = n - k
	BB = matrix(ZZ,B) 
	A  = BB[0:log_covolume, m:rank]			# A in ZZ^{ m x k } such that B = qary Lattice of LWE (A)
	for i in range(T):
		# s , e = random_matrix(R , 1 , rank ) , random_matrix(R , rank , 1 )
		s = random_matrix(R , 1 , rank )
		e = matrix(ZZ, [D() for _ in range(rank)]).transpose()
		target_list_LWE.append( (s * BB).transpose() + e ) , target_list_unif.append( random_matrix(R , rank , 1 ) ) # BB * s + e )
		
		# A_big = [ [Id ] , [A] ]  
		# A_big * s + e = s_0 + e_0 		, A s_1 + e_1
		# B * s + e 	= s_0 + A * s_1 + e_0 	, q s_1 + e_1
		
save( ( target_list_LWE , target_list_unif) , 'Target_samples/target_samples' )		



size_db_random_s = 100
if is_LWE :
	target_list = []
	R  = Zmod(modulus) 				# ZZ mod q
	m , k = rank - log_covolume , log_covolume			# m = n - k
	BB = matrix(ZZ,B) 
	A  = BB[0:log_covolume, m:rank]			# A in ZZ^{ m x k } such that B = qary Lattice of LWE (A)
	A0 = BB[0:k, m:rank]
	A1 = BB[k:log_covolume, m:rank]
	s0 = random_matrix(R , 1 , k )
	s1 = random_matrix(R , 1 , rank - k )
	e = matrix(ZZ, [D() for _ in range(rank)]).transpose()
	target_list_LWE.append( (s0 * BB).transpose() + e )
	for i in range(size_db_random_s):
		# s , e = random_matrix(R , 1 , rank ) , random_matrix(R , rank , 1 )
		s_tilde = random_matrix(R , 1 , rank )
		target_list.append( (s * BB).transpose() + e )  # BB * s + e )
		
		# A_big = [ [Id ] , [A] ]  
		# A_big * s + e = s_0 + e_0 		, A s_1 + e_1
		# B * s + e 	= s_0 + A * s_1 + e_0 	, q s_1 + e_1
