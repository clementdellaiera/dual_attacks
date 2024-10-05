################################
##   INFORMATION ON DATASET   ##
################################

# Dual vectors
def Dual_vectors_info():
	W = load('Dual_short_vectors/dual_short_vectors')
	number = len(W[1])
	length_dual_vectors = [matrix(ZZ,w).norm() for w in W[1]]
	min_W , max_W , mean_W = min(length_dual_vectors) , max(length_dual_vectors) , mean(length_dual_vectors)
