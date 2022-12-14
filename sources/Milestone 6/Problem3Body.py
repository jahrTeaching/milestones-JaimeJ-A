
from numpy import zeros, size, linspace, array, log10, absolute, sqrt
from numpy.linalg import eig

from Matematicas import Newton, Jacobiano

#----- MODULO PROBLEMA 3 CUERPOS -----#


# 3 Body Problem
def Problem3Body(U, t, mu):
    
    x, y = U[0], U[1]       
    dxdt, dydt = U[2], U[3]     

    r1 = sqrt( (x + mu)**2 + y**2 )       
    r2 = sqrt( (x - 1 + mu)**2 + y**2 )    

    dvxdt = 2*dydt + x - ((1 - mu)*(x + mu))/(r1**3) - mu*(x + mu - 1)/(r2**3)
    dvydt = -2*dxdt + y -((1 - mu)/(r1**3) + mu/(r2**3))*y

    return array([dxdt, dydt, dvxdt, dvydt]) 
# -----------------

# Lagrange
def Lagrange_Points(U0, NL, mu):

    L_P = zeros([5,2])

    def F(y):
        x = zeros(4)
        x[0:2] = y
        x[2:4] = 0
        F_x = Problem3Body(x, 0, mu)
        return F_x[2:4]
        
    for i in range(NL):
        L_P[i,:] = Newton(F, U0[i,0:2])

    return L_P
# ----------------

# Lagrange estabilidad
def Lagrange_Points_Stability(U0, mu):

    def F(y):
        return Problem3Body(y, 0 , mu)

    A = Jacobiano(F, U0)
    values, vectors = eig(A)

    return values
# ---------------




