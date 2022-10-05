
from numpy import array

#----- MODULO KEPLER -----#


# Kepler(orbita)
def Kepler(U, t):

    x = U[0]; y = U[1]; dx_dt = U[2]; dy_dt = U[3]
    d = (x**2 + y**2)**1.5

    return array( [ dx_dt, dy_dt, -x/d, -y/d] )
# ------------