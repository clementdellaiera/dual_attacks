##   Load lattice basis 
load('Parameters.sage')
[ is_LWE , rank , modulus , log_covolume ] , B = load('Lattice_basis/basis_LWE.sobj')


if is_LWE :
	target_list_unif , target_list_LWE = [] , []
	R  = Zmod(modulus) 				# ZZ mod q
	m  = rank - log_covolume			# m = n - k
	BB = matrix(ZZ,B) 
	A  = BB[0:log_covolume, m:rank]			# A in ZZ^{ m x k } such that B = qary Lattice of LWE (A)
	for i in range(T):
		s , e = random_matrix(R , 1 , rank ) , random_matrix(R , rank , 1 )
		target_list_LWE.append( (s * BB).transpose()  + e) , target_list_unif.append( random_matrix(R , rank , 1 ) ) # BB * s + e )
		
		# A_big = [ [Id ] , [A] ]  
		# A_big * s + e = s_0 + e_0 		, A s_1 + e_1
		# B * s + e 	= s_0 + A * s_1 + e_0 	, q s_1 + e_1
		
save( ( target_list_LWE , target_list_unif) , 'Target_samples/target_samples' )		
