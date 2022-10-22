
from numpy import array, zeros, dot
from numpy.linalg import inv, norm

#----- MODULO MATEMATICAS -----#


# Jacobiano
def Jacobiano(F, U):

    dx = 1e-3
    dim = len(U)
    Jab = zeros( (dim,dim) )
    x = zeros( dim )

    for i in range(dim):

        x[i] = dx
        Jab[:,i] = ( F(U + x) - F(U - x) )/( 2*dx )

    return Jab
# ------------

# Metodo de Newton
def Newton(F, U0):

    dim = len(U0)
    dx = array(zeros(dim))
    b  = array(zeros(dim))
    
    eps = 1
    it = 0
    it_max = 10000

    while ( eps > 1e-8 ) and ( it <= it_max ):
   
        it = it + 1
        Jab = Jacobiano(F,U0)
        b = F(U0)
        dx = dot( inv(Jab),b )
        U0  = U0 - dx
        eps = norm( dx )
        
    return U0
# ------------
