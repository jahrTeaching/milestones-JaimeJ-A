
from numpy import array

#----- MODULO OSCILADOR -----#


# Oscilador Lineal
def Oscilador(U, t):

    x = U[0]; dx_dt = U[1]

    return array( [ dx_dt, -x] )
# ------------































