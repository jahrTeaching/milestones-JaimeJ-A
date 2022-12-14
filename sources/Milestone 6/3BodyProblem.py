
from numpy import zeros, size, linspace, array, log10, absolute, sqrt
from numpy.linalg import eig

#----- MODULO PROBLEMA 3 CUERPOS -----#


# 3 Body Problem
def Problem3Body(U, t, mu):
    
    x, y = U[0], U[1]       # The first two components are the position coordinates
    vx, vy = U[2], U[3]     # The last two components velocity 

    r1 = sqrt((x + mu)**2 + y**2)        # Distance from the third body to the mass m1
    r2 = sqrt((x - 1 + mu)**2 + y**2)    # Distance from the third body to the mass m2

    dxdt = vx
    dydt = vy

    dvxdt = 2*vy + x - ((1 - mu)*(x + mu))/(r1**3) - mu*(x + mu - 1)/(r2**3)
    dvydt = -2*vx + y -((1 - mu)/(r1**3) + mu/(r2**3))*y


    return array([dxdt, dydt, dvxdt, dvydt])
# -----------------

# Lagrange
def Lagrange_Points_Calculation(U0, NL, mu):

    LP = zeros([5,2])

    def F(Y):
        
        X = zeros(4)
        X[0:2] = Y
        X[2:4] = 0
        FX = CR3BP(X, 0, mu)
        return FX[2:4]
        
    for i in range(NL):
        LP[i,:] = newton(F, U0[i,0:2])

    return LP

def Lagrange_Points_Stability(U0, mu):

    def F(Y):
        return CR3BP(Y, 0 , mu)

    A = jacobiano(F, U0)
    values, vectors = eig(A)

    return values
