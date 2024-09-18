# Scripts

Parameters.sage
	No input
	Declaration of variables
		rank , modulus 	: int , int 	| dimension of lattice , modulus of LWE instances
		log_covolume 	: int		| in the case of a LWE lattice, this is k = n-m 
		is_LWE 		: bool		| True if lattice is qary

LWE_instance_generator.sage
	Input 	| Variables from Parameters.sage
	Output	| [ is_LWE , rank , modulus , log_covolume ] , B : parameter_set ( bool , int , int , int ) , basis of lattice  
			to 'Lattice_basis/basis_LWE.sobj'
			
Dual_short_vectors.sage
	Input   | [ is_LWE , rank , modulus , log_covolume ] , B : parameter_set ( bool , int , int , int ) , basis of lattice  
			from Lattice_basis/basis_LWE.sobj'
	Output  | modulus , short_vector_list : int , list of vectors in dual lattice
			to 'Dual_short_vectors/dual_short_vectors.sobj'

Target_samples.sage
	Input	| Variables from Parameters.sage 
		| [ is_LWE , rank , modulus , log_covolume ] , B : parameter_set ( bool , int , int , int ) , basis of lattice  
			from Lattice_basis/basis_LWE.sobj'	
	Output	| target_list_LWE , target_list_unif 
			to 'Target_samples/target_samples'

Distinguisher.sage
	Input 	| modulus , short_vector_list : int , list of vectors in dual lattice
			from 'Dual_short_vectors/dual_short_vectors.sobj'
		| target_list_LWE , target_list_unif 
			from 'Target_samples/target_samples'
	Output	| data_LWE , data_unif : list of floats , list of floats  
			to 'Data/distinguisher_values'
Graphs.sage
	Input 	| data_LWE , data_unif 
			from 'Data/distinguisher_values')
	Output	| 
	
# Folders

Lattice_basis
	[ is_LWE , rank , modulus , log_covolume ] , B
Dual_short_vectors
	modulus , short_vector_list 
Target_samples
	target_list_LWE , target_list_unif
