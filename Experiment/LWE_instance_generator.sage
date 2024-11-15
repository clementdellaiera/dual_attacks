load('Parameters.sage')
from fpylll import IntegerMatrix

# Instance generation
if is_LWE :
	B = IntegerMatrix.random( rank , 'qary' , k = log_covolume , q = modulus)
else :
	B = random_matrix(ZZ, rank)
	
Param = [ is_LWE , stdev , rank , modulus , log_covolume ]	

# Save to folder
file_path = 'Lattice_basis/basis_isLWE_' + str(is_LWE)+'_stdev_' +str(stdev )+'_rank_' +str(rank )+'_modulus_' +str(modulus )+'_logcovolume' +str(log_covolume) 			
save((Param , B) , file_path)
