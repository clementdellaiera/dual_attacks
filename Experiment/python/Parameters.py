# Parameters

stdev, rank , modulus  = 0.2, 30 , 37  # stdev parameters for discrete gaussian error in LWE, dimension of lattice , modulus of LWE instances
log_covolume = rank / 2			# log_covolume = k = n-m | IntegerMatrix takes for qary lattices k in det L = q**k as parameters | basis of type [ [I_{n-k} , A ] , [0 , q I_k]] with size(A) = (n-k , k)

is_LWE = True

T = 1000 				# Size of target samples

