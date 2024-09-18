load('Parameters.sage')
from fpylll import IntegerMatrix

# Instance generation
if is_LWE :
	B = IntegerMatrix.random( rank , 'qary' , k = log_covolume , q = modulus)
else :
	B = random_matrix(ZZ, rank)
	
Param = [ is_LWE , rank , modulus , log_covolume ]	
# Save to folder
save((Param , B) , 'Lattice_basis/basis_LWE')
