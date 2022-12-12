
from Esquemas_Numericos import Leap_Frog

from numpy import zeros, float64, size, linspace, array, log10, round_, absolute
from numpy.linalg import norm
from sklearn.linear_model import LinearRegression

#----- MODULO EDO(Ecuaciones diferenciales ordinarias)-----#


# Problema  Cauchy
def Cauchy(F, t, U0, E_Temporal):

    n, nv = len(t)-1, len(U0)
 
    U = zeros((n+1,nv), dtype = float64)
    dt = t[1] - t[0]
    U[0,:] = U0

    if E_Temporal == Leap_Frog:

        U[1,:] = U[0,:] + dt*F(U[0,:],t[0])

        for i in range(1,n):

            U1 = U[i-1, :]
            U2 = U[i, :]
            U[i+1, :] = E_Temporal(U1, U2, dt, t[i], F)

    else:
        for i in range(n):

            U[i+1,:] = E_Temporal(U[i,:], t[i+1] - t[i],t[i],F)
            
    return U
# ------------

# Error de Richardson
def Error_Richardson(F, t, U0, E_Temporal, q):

    N  = size(t)         
    Error = zeros((N,size(U0)))

    T  = t[N-1]             
    t1 = t                
    t2 = linspace(0,T,2*N)
    
    U1 = Cauchy(F, t1, U0, E_Temporal)
    U2 = Cauchy(F, t2, U0, E_Temporal)

    for i in range(N):
        Error[i,:] = ( U2[2*i,:] - U1[i,:] )/( 1 - (1/(2**q)) )

    return Error 

def Error_Module(Error_E_Temporal):

    Module_Error = zeros( size(Error_E_Temporal[:,0]) )

    for i in range(0,size(Error_E_Temporal[:,0])):
     Module_Error[i] = (Error_E_Temporal[i,0]**2 + Error_E_Temporal[i,1]**2)**(1/2)

    return Module_Error
# ------------

# Convergencia
def Convergence_rate(F, t, U0, E_Temporal):  

    n = size(t)
    T = t[n-1]
    t1 = t
 
    U1 = Cauchy(F, t1, U0, E_Temporal)
   
    m = 2 
    log_E = zeros(m)
    log_N = zeros(m) 
    Error = zeros(m)
    n = 2*n
    
    for i in range (0,m):

        t2 = linspace(0, T, (2**i)*n)
        U2 = Cauchy(F, t2, U0, E_Temporal)

        Error[i] = norm(U2[int((2**i)*n-1),:] - U1[int((2**i)*n/2-1),:])
        log_E[i] = log10(Error[i])
        log_N[i] = log10((2**i)*n)

        U1 = U2
        print(i)

        for j in range(0,m):

         if (abs(log_E[j]) > 12):

             break

    j = min(j, m)
    
    # Regresión lineal y calculo de la pendiente 
    
    reg = LinearRegression().fit(log_N[0:j+1].reshape((-1, 1)),log_E[0:j+1]) 
    order = round_(abs(reg.coef_),1)

    log_N_lineal = log_N[0:j+1]
    log_E_lineal = reg.predict(log_N[0:j+1].reshape((-1, 1)))

    return [log_E, log_N, log_E_lineal, log_N_lineal, order]
# ------------

# Region de estabilidad
def Region_Estabilidad(E_Temporal): 

    N = 100
    x, y = linspace(-5, 5, 100), linspace(-5, 5, 100)
    rho =  zeros( (N, N),  dtype=float64)

    for i in range(N): 
      for j in range(N):

          w = complex(x[i], y[j])
         
          if E_Temporal.__name__ == "LEAP FROG":
            rho[i, j] = abs(E_Temporal( 1., 1., 1., 0., lambda u, t: w*u )) 

          else:
            rho[i, j] = abs(E_Temporal( 1., 1., 0., lambda u, t: w*u )) 

    return rho, x, y 
# ------------




















