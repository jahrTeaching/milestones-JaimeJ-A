
from numpy import array, zeros, linspace, size
import matplotlib.pyplot as plt

from Esquemas_Numericos import Euler, RK4, Crank_Nicolson, Inverse_Euler
from EDO import Cauchy, Error_Richardson, Error_Module, Convergence_rate
from Kepler import Kepler

#----- MODULO PRINCIPAL HITO 3 -----#


# Variables
T = 15          # Duracion
dt = 0.01      # Paso integracion

N = int(T/dt)                 
t  = linspace(0,T,N)            
U0 = array( [1,0,0,1] )


# Errores
E_Euler          = Error_Richardson(Kepler, t, U0, Euler, 1)
E_RK4            = Error_Richardson(Kepler, t, U0, RK4, 4)
E_Crank_Nicolson = Error_Richardson(Kepler, t, U0, Crank_Nicolson, 2)
E_Inverse_Euler  = Error_Richardson(Kepler, t, U0, Inverse_Euler, 1)


# Graficas
E_Temporal       = [ Euler, RK4, Crank_Nicolson, Inverse_Euler ]
E_Temporal_Error = [ E_Euler, E_RK4, E_Crank_Nicolson, E_Inverse_Euler ]
E_Temporal_Plot  = ['EULER','RUNGE-KUTTA 4','CRANK-NICOLSON','EULER INVERSO']

for i in range(4): 

    E_Temporal_Error_Module = Error_Module( E_Temporal_Error[i] )

    plt.title(f'ERROR {E_Temporal_Plot[i]}, dt = {dt} s')
    plt.plot(t, E_Temporal_Error[i][:,0],"r", label = "x")
    plt.plot(t, E_Temporal_Error[i][:,1],"b", label = "y")
    plt.plot(t, E_Temporal_Error_Module, 'c', label = 'Module')
    plt.xlabel("t (s)")
    plt.ylabel("Error")
    plt.legend(loc = "upper left")
    plt.grid()

    plt.savefig('Plots/' + E_Temporal_Plot[i]+ ' ' + str(dt)+'.png')
    plt.show()

    [log_E, log_N, log_E_lineal, log_N_lineal, order] = Convergence_rate(Kepler, t, U0, E_Temporal[i])
    
    plt.plot(log_N, log_E, "b", label = E_Temporal_Plot[i])
    plt.plot(log_N_lineal, log_E_lineal, "r", label = 'Linear regression')
    plt.legend(loc ='lower left')
    plt.xlabel("log(N)")
    plt.ylabel("log(U2-U1)")
    plt.title(f'{E_Temporal_Plot[i]}, order = {order}')
    plt.plot()
    plt.grid()
    plt.savefig('Plots/'+ E_Temporal_Plot[i]+ ' ' + str(dt)+'.png')
    plt.show() 
