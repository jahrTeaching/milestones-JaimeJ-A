from numpy import array, zeros, linspace, size, transpose, reshape
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d

from Esquemas_Numericos import Euler, RK4, Crank_Nicolson, Inverse_Euler, Leap_Frog
from EDO import Cauchy
from N_Bodies import N_Bodies

#----- MODULO PRINCIPAL HITO 4 -----#


# VARIABLES
T = 100          # Duracion
dt = 0.1        # Paso integracion

N = int(T/dt)                 
t  = linspace(0,T,N)


# CONDICIONES INICIALES
Nb = 4 # Bodies
Nc = 3 # Coord

U0 = zeros(Nb*Nc*2)
U_0 = reshape(U0, (Nb, Nc, 2) )
r0 = reshape(U_0[:, :, 0], (Nb, Nc))
v0 = reshape(U_0[:, :, 1], (Nb, Nc))

# body 1 
r0[0,:] = [ 2, 2, 0]
v0[0,:] = [ -0.4, 0, 0]
# body 2 
r0[1,:] = [ -2, 2, 0] 
v0[1,:] = [ 0, -0.4, 0]
# body 3 
r0[2, :] = [ -2, -2, 0 ] 
v0[2, :] = [ 0.4, 0., 0. ] 
# body 4 
r0[3, :] = [ 2, -2, 0 ] 
v0[3, :] = [ 0., 0.4, 0. ]

# SOLUCIÓN
U = Cauchy(N_Bodies, t, U0, RK4)

# GRAFICAS
colors = ['b','r','g','m']

U_s  = reshape( U, (N, Nb, Nc, 2) ) 
r   = reshape( U_s[:, :, :, 0], (N, Nb, Nc) )

    # 2D
for i in range(Nb):
    plt.figure(1)
    plt.plot(r[:, i, 0], r[:, i, 1], colors[i])

plt.title(f'N = {Nb} body problem, xy')
plt.xlabel("x")
plt.ylabel("y")
plt.grid()
plt.savefig('Plots/'+ '{Nb}' + '2D' +'.png')
plt.show()

    # 3D 
fig = plt.figure(2)
ax = fig.add_subplot(projection='3d')

for i in range(Nb):

    ax.plot(r[:, i, 0], r[:, i, 1], r[:, i, 2], colors[i])

plt.title(f'N = {Nb} body problem, 3D')
plt.xlabel("x")
plt.ylabel("y")
plt.ylabel("z")
plt.savefig('Plots/'+ '{Nb}' + '3D' +'.png')
plt.show()