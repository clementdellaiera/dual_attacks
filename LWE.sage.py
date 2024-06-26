

# This file was *autogenerated* from the file LWE.sage
from sage.all_cmdline import *   # import sage library

_sage_const_0p5 = RealNumber('0.5'); _sage_const_1 = Integer(1); _sage_const_37 = Integer(37); _sage_const_50 = Integer(50); _sage_const_3 = Integer(3); _sage_const_2 = Integer(2); _sage_const_0 = Integer(0); _sage_const_1000 = Integer(1000)
load('dual_attack.sage')

def character(w,t):
	# Given w in dual lattice and t in vector space 
	# Returns exp(2 * i * pi * < w , t >) where the circle is identified as (-0.5 , 0.5)
	# via exp( 2* i * pi * x ) -> floor( x- 0.5 ) + 0.5
	# Lvee x V -> (-0.5 , 0.5) = RR / ZZ
	inner_product = w.inner_product( vector(RR, t) )
	return(  ((inner_product + _sage_const_0p5  ) % _sage_const_1 ) + _sage_const_0p5  )

q , n , s = _sage_const_37  , _sage_const_50  , _sage_const_3   # 32768 , 600 , 3
A = random_matrix( Zmod(q) , n)

# Lattices basis

B = matrix(ZZ , _sage_const_2 *n , _sage_const_2 *n )
B[_sage_const_0 :n , _sage_const_0 :n] , B[ n:_sage_const_2 *n , _sage_const_0 :n] , B[n:_sage_const_2 *n , n:_sage_const_2 *n]= matrix.identity(ZZ , n) , A.lift().transpose() , q * matrix.identity(ZZ , n)
"""B =  [ I_n  ,    A  ]
	[  0   , q I_n ] 
	basis of L_A """

Gauss_L_A = DiscreteGaussianDistributionLatticeSampler(B.BKZ(), s)

det_B , B_v , L_A_dual  = Dual(B)
B_v_basis = (_sage_const_1 /det_B) * B_v.BKZ()
Gauss_L_A_dual = DiscreteGaussianDistributionLatticeSampler(B_v_basis, s)
# Sampling of short vectors

W , W_size = [] , _sage_const_1000 	
ee = DiscreteGaussianDistributionLatticeSampler( matrix.identity(ZZ , _sage_const_2 *n) , s)()

target_LWE , target_UNIF = Gauss_L_A() + ee , vector(RR, uniform(_sage_const_0 ,_sage_const_1 ,_sage_const_2 *n)) * B 
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
yy = (_sage_const_1 /W_size) * np.linspace(_sage_const_0 ,_sage_const_1 ,W_size)

figure, axis = plt.subplots(_sage_const_2 ) 

axis[_sage_const_0 ].plot(data_inner_products_LWE , yy , color = "blue", label="LWE" )
axis[_sage_const_0 ].plot(data_inner_products_UNIF, yy , color = "cyan", label="Uniform" )
axis[_sage_const_0 ].set_title("Empirical c.d.f. of < w , t > for w in W")

plt.show()

