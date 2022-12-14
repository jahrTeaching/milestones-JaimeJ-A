from numpy import array, zeros, linspace, size, transpose, reshape, around
from numpy.random import random
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d

from Esquemas_Numericos import Euler, RK4, Crank_Nicolson, Inverse_Euler, Leap_Frog, RKE
from EDO import Cauchy
from Problem3Body import Problem3Body, Lagrange_Points, Lagrange_Points_Stability

from random import random


#----- MODULO PRINCIPAL HITO 6 -----#


# VARIABLES
T = 1000         # Duracion
N = int(1e6)                 
t  = linspace(0,T,N)


## CR3BP SOLUCION ##

mu = 3.0039e-7  # Earth - Sun
#mu = 1.2151e-2 # Earth - Moon

def F(U, t):
    return Problem3Body(U, t, mu)


## PUNTOS LAGRANGE ##

NLP = 5 # Numbero Puntos Lagrange (5)

# Condiciones Iniciales
U0 = zeros([NLP,4])
U0[0,:] = array([0.8, 0.6, 0, 0])
U0[1,:] = array([0.8, -0.6, 0, 0])
U0[2,:] = array([-0.1, 0, 0, 0])
U0[3,:] = array([0.1, 0, 0, 0])
U0[4,:] = array([1.01, 0, 0, 0])


LagrangePoints = Lagrange_Points(U0, NLP, mu)
U0LP = zeros(4)
U0SLP = zeros(4)
eps = 1e-3*random()
Lagrange_Points_List = array([1,2,3,4,5])

Temp_Schemes = [Euler, RK4, Crank_Nicolson, Inverse_Euler, Leap_Frog, RKE]
T_S_list = ['Euler', 'RK4', 'CrankNicolson', 'InverseEuler', 'LeapFrog', 'Embedded_RK']

for k in range(6):

    for i in range(NLP):

        sel_LP = i + 1

        if sel_LP == 5:
            label = 'L2'
        elif sel_LP == 4:
            label = 'L1'
        elif sel_LP == 3:
            label = 'L3'
        elif sel_LP == 2:
            label = 'L5'
        elif sel_LP == 1:
            label = 'L4'

        U0LP[0:2] = LagrangePoints[sel_LP-1,:] + eps
        U0LP[2:4] = eps
        U0SLP[0:2] = LagrangePoints[sel_LP-1,:] 
        U0SLP[2:4] = 0

        ## ESTABILIDAD ##

        TS = Temp_Schemes[k]
        T_S = T_S_list[k]

        U_LP = Cauchy(F, t, U0LP, TS)
        eig = Lagrange_Points_Stability(U0SLP, mu)
        print(around(eig.real,8))

        ## GRAFICAS ##
        fig, (ax1, ax2) = plt.subplots(1, 2)
        ax1.plot(U_LP[:,0], U_LP[:,1],'-',color = "r")
        ax1.plot(-mu, 0, 'o', color = "g")
        ax1.plot(1-mu, 0, 'o', color = "b")

        for i in range(NLP):
            ax1.plot(LagrangePoints[i,0], LagrangePoints[i,1] , 'o', color = "k")

        ax1.set_xlim(-9,9)
        ax1.set_ylim(-9,9)
        ax1.set_title("Orbital system view")

        ax2.plot(U_LP[:,0], U_LP[:,1],'-',color = "r")
        ax2.plot(LagrangePoints[sel_LP-1,0], LagrangePoints[sel_LP-1,1] , 'o', color = "k")


        ax2.set_title("Lagrange point view")
        ax2.set_xlim(LagrangePoints[sel_LP-1,0]-0.25,LagrangePoints[sel_LP-1,0]+0.25)
        ax2.set_ylim(LagrangePoints[sel_LP-1,1]-0.25,LagrangePoints[sel_LP-1,1]+0.25)
        fig.suptitle(f"Earth-Sun - CR3BP ({T_S} ) - Orbit around the {label} point with t = " + str(t[N-1])+'s')

        for ax in fig.get_axes():
            ax.set(xlabel='x', ylabel='y')
            ax.grid()

        plt.savefig('Plots/' + label +' '+ T_S +'.png')
        #plt.show()

