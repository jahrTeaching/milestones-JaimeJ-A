
from Matematicas import Newton

#----- MODULO ESQUEMAS NUMERICOS -----#


# Euler        
def Euler(U, dt, t, F):

    return U + dt*F(U,t)
# ------------

# Runge-Kutta orden 4
def RK4(U, dt, t, F):

    k1 = F( U, t )
    k2 = F( U + (k1 * dt)/2, t + dt/2 )
    k3 = F( U + (k2 * dt)/2, t + dt/2 )
    k4 = F( U + (k3 * dt)  , t + dt   )

    return U + ( dt/6 )*( k1 + 2*k2 + 2*k3 + k4 )
# ------------

# Crank-Nicolson
def Crank_Nicolson(U, dt, t, F):
    
    def CN(x):
        return x - U - dt/2*( F(U, t) + F(x,t) )

    return Newton( CN, U )
# ------------

# Euler Inverso
def Inverse_Euler(U, dt, t, F):

    def IE(x):
        return x - U - dt*F( x, t )

    return Newton( IE, U )
# ------------