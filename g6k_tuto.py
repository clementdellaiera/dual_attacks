import g6k
from fpylll import IntegerMatrix

n , m , b = 10 , 5 , 5 
A = IntegerMatrix.random( n , "qary" , m , bits = b )
print(A)
