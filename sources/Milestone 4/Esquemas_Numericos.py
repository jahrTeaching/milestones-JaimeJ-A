
import sys
sys.path.append('.')
from scipy.optimize import newton, fsolve

#----- MODULO ESQUEMAS NUMERICOS -----#


# Euler        
def Euler(U, dt:float, t, F):
    Euler.__name__ = "EULER"

    return U + dt*F(U,t)
# ------------

# Runge-Kutta orden 4
def RK4(U, dt, t, F):
    RK4.__name__ = "RUNGE-KUTTA 4"

    k1 = F( U, t )
    k2 = F( U + (k1 * dt)/2, t + dt/2 )
    k3 = F( U + (k2 * dt)/2, t + dt/2 )
    k4 = F( U + (k3 * dt)  , t + dt   )

    return U + ( dt/6 )*( k1 + 2*k2 + 2*k3 + k4 )
# ------------

# Crank-Nicolson
def Crank_Nicolson(U, dt, t, F):
    Crank_Nicolson.__name__ = "CRANK-NICOLSON"
    
    def CN(x):
        return x - U - dt/2 * ( F(U, t) + F(x,t) )

    return newton( CN, U )
# ------------

# Euler Inverso
def Inverse_Euler(U, dt, t, F):
    Inverse_Euler.__name__ = "EULER INVERSO"

    def IE(x):
        return x - U - dt*F( x, t )

    return newton( IE, U )
# ------------

# Euler Inverso
def Leap_Frog (U1, U2, dt, t, F): 
    Leap_Frog.__name__ = "LEAP FROG"

    return U1 + 2*dt*F(U2,t)
# ------------
















