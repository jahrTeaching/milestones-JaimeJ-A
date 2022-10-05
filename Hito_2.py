
from numpy import array, zeros, linspace
import matplotlib.pyplot as plt

from Esquemas_Numericos import Euler, RK4, Crank_Nicolson, Inverse_Euler
from EDO import Cauchy
from Kepler import Kepler

#----- MODULO PRINCIPAL HITO 2 -----#


# Variables
T = 20          # Duracion
dt = 0.005      # Paso integracion


n = int(T/dt)                 
t  = linspace(0,T,n)            
U0 = array( [1,0,0,1] )
U = U0

E_Temporal = [ Euler, RK4, Crank_Nicolson, Inverse_Euler ]
E_Temporal_Plot = ['EULER','RUNGE-KUTTA 4','CRANK-NICOLSON','EULER INVERSO']

# Graficas
for i in range (4):

    U = Cauchy( Kepler, t, U0, E_Temporal[i] )
 
    print( U[:, len(t)-1] )
    plt.title(f'{E_Temporal_Plot[i]}')
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.plot( U[0,:], U[1,:] )

    plt.savefig('Plots/' + E_Temporal_Plot[i]+ ' ' + str(dt)+'.png')
    plt.show()
