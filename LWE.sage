load('Distinguisher/dual_attack.sage')

def character(w,t):
	# Given w in dual lattice and t in vector space 
	# Returns exp(2 * i * pi * < w , t >) where the circle is identified as (-0.5 , 0.5)
	# via exp( 2* i * pi * x ) -> floor( x- 0.5 ) + 0.5
	# Lvee x V -> (-0.5 , 0.5) = RR / ZZ
	inner_product = w.inner_product( vector(RR, t) )
	return(  ((inner_product + 0.5 ) % 1) + 0.5 )

q , n , s = 37 , 50 , 3  # 32768 , 600 , 3
A = random_matrix( Zmod(q) , n)

# Lattices basis
B = block_matrix([ [ matrix.identity(ZZ,n) , A.lift()], [matrix.zero(ZZ,n) , q * matrix.identity(ZZ , n)]]) 

Gauss_L_A = DiscreteGaussianDistributionLatticeSampler(B.BKZ(), s)

det_B , B_v , L_A_dual  = Dual(B)
B_v_basis = (1/det_B) * B_v.BKZ()
Gauss_L_A_dual = DiscreteGaussianDistributionLatticeSampler(B_v_basis, s)
# Sampling of short vectors

W , W_size = [] , 1000	
ee = DiscreteGaussianDistributionLatticeSampler( matrix.identity(ZZ , 2*n) , s)()

target_LWE , target_UNIF = Gauss_L_A() + ee , vector(RR, uniform(0,1,2*n)) * B 
data_inner_products_LWE , data_inner_products_UNIF = [] , []

for i in range(W_size) :
	w = Gauss_L_A_dual()
	W.append( w ) , data_inner_products_LWE.append( character(w,target_LWE) )  , data_inner_products_UNIF.append( character(w,target_UNIF) ) 

################
##   GRAPHS   ##
################
	
import seaborn as sns
sns.set_style('darkgrid')

data_inner_products_LWE.sort() , data_inner_products_UNIF.sort()
yy = (1/W_size) * np.linspace(0,1,W_size)

figure, axis = plt.subplots(2) 

axis[0].plot(data_inner_products_LWE , yy , color = "blue", label="LWE" )
axis[0].plot(data_inner_products_UNIF, yy , color = "cyan", label="Uniform" )
axis[0].set_title("Empirical c.d.f. of < w , t > for w in W")

plt.show()
