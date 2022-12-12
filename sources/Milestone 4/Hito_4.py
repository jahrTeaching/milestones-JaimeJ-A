
from numpy import array, zeros, linspace, size, transpose
import matplotlib.pyplot as plt

from Esquemas_Numericos import Euler, RK4, Crank_Nicolson, Inverse_Euler, Leap_Frog
from EDO import Cauchy, Region_Estabilidad
from Oscilador import Oscilador

#----- MODULO PRINCIPAL HITO 4 -----#


# Variables
T = 30              # Duracion
dt = 0.001          # Paso integracion

N = int(T/dt)                 
t  = linspace(0,T,N)

# Condiciones iniciales
U0 = array( [0,1] )


# Graficas
E_Temporal = [ Euler, RK4, Crank_Nicolson, Inverse_Euler, Leap_Frog ]

for j in range(5): 

        # Oscilador lineal

        Osc = Cauchy(Oscilador, t, U0, E_Temporal[j])

        plt.title(f'Oscilador con {E_Temporal[j].__name__}; dt = {dt}')
        plt.xlabel("Tiempo (s)")
        plt.ylabel("x")
        plt.grid()

        plt.plot(t, Osc[:,0], label = 'dt =' + str(dt) + ' s')
        plt.legend(loc ='lower left')    
        
        plt.savefig('Plots/Oscilador/' + E_Temporal[j].__name__+ ' ' + str(dt)+'.png')
        plt.close()
    

        # Región absoluta de estabilidad
    
        Reg, R, I = Region_Estabilidad(E_Temporal[j])

        zer = zeros(100)
        fig = plt.figure()
        ax = fig.add_subplot()
        plt.title(f'Región de estabilidad absoluta del {E_Temporal[j].__name__} ')
        plt.plot(R,zer,'k-')                         ## X axis
        plt.plot(zer,I,'k-')                         ## Y axis
        plt.grid()
        plt.contour(R, I, transpose(Reg), levels = [0, 1], colors = ['r'], linewidth = 2 )
        plt.contourf(R, I, transpose(Reg), levels = [0, 1], colors =['grey'])            ## Stability region

        if E_Temporal[j] != Crank_Nicolson:
            ax.set_aspect('equal', adjustable = 'box')

        plt.xlabel("Re")
        plt.ylabel("Im")

        plt.savefig('Plots/Regiones/' + E_Temporal[j].__name__ + ' ' + str(dt)+ '.png')
        plt.close()