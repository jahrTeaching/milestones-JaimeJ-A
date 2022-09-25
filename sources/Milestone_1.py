"""
Milestone 1
"""

from numpy import array, zeros
import matplotlib.pyplot as plt
from scipy.optimize import fsolve


"Integración"
tf = 20 
N = 2000
dt = tf/N

"Condiciones Iniciales"
U0 = array( [ 1, 0, 0, 1] ) 
U = zeros((len(U0) , N))
U[:,0] = U0



"EULER"

for i in range(1, N):
    F = array( [U[2,i-1], U[3,i-1], -U[0,i-1]/(U[0,i-1]**2+U[1,i-1]**2)**(1.5), -U[1,i-1]/(U[0,i-1]**2+U[1,i-1]**2)**(1.5)] )
    U[:,i] = U[:,i-1] + dt*F

plt.figure(1)
plt.title("Euler")
plt.plot( U[0,:], U[1,:] )
plt.show()



"RUNGE-KUTTA 4"

k1 = zeros(len(U0)); k2 = zeros(len(U0)); k3 = zeros(len(U0)); k4 = zeros(len(U0))
U_2 = zeros(len(U0)); U_3 = zeros(len(U0)); U_4 = zeros(len(U0))

for i in range(1, N):
    k1 = array( [U[2,i-1], U[3,i-1], -U[0,i-1]/(U[0,i-1]**2 + U[1,i-1]**2)**(1.5), -U[1,i-1]/(U[0,i-1]**2 + U[1,i-1]**2)**(1.5)] )
    U_2 = U[:,i-1] + dt*k1[:]/2
    k2 = array( [U_2[2], U_2[3], -U_2[0]/(U_2[0]**2 + U_2[1]**2)**(1.5), -U_2[1]/(U_2[0]**2 + U_2[1]**2)**(1.5)] )
    U_3 = U[:,i-1] + dt*k2[:]/2
    k3 = array( [U_3[2], U_3[3], -U_3[0]/(U_3[0]**2 + U_3[1]**2)**(1.5), -U_3[1]/(U_3[0]**2 + U_3[1]**2)**(1.5)] )
    U_4 = U[:,i-1] + dt*k3[:]
    k4 = array( [U_4[2], U_4[3], -U_4[0]/(U_4[0]**2 + U_4[1]**2)**(1.5), -U_4[1]/(U_4[0]**2 + U_4[1]**2)**(1.5)] )
    U[:,i] = U[:,i-1] + dt*(k1+2*k2+2*k3+k4)/6

plt.figure(2)
plt.title("Runge-Kutta 4")
plt.plot( U[0,:], U[1,:] )
plt.show()



"CRANK-NICHOLSON"

for i in range(1,N):
    F = array( [U[2,i-1], U[3,i-1], -U[0,i-1]/(U[0,i-1]**2 + U[1,i-1]**2)**(1.5), -U[1,i-1]/(U[0,i-1]**2 + U[1,i-1]**2)**(1.5)] )
    def CN(x):
        return [x[0] - U[0,i-1] - (x[2] + F[0])*dt/2,
                x[1] - U[1,i-1] - (x[3] + F[1])*dt/2,
                x[2] - U[2,i-1] - (-x[0]/((x[0]**2+x[1]**2)**(1.5)) + F[2])*dt/2,
                x[3] - U[3,i-1] - (-x[1]/((x[0]**2+x[1]**2)**(1.5)) + F[3])*dt/2 ]
    U[:,i] = fsolve(CN, [U[:,i-1]])
        

plt.figure(3)
plt.title("Crank-Nicolson")
plt.plot( U[0,:], U[1,:] )
plt.show()
