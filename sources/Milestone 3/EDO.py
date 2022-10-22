
from numpy import zeros, float64, size, linspace, array, log10, round_
from numpy.linalg import norm
from sklearn.linear_model import LinearRegression

#----- MODULO EDO(Ecuaciones diferenciales ordinarias)-----#


# Problema  Cauchy
def Cauchy(F, t, U0, E_Temporal):

    n = len(t)-1
    nv = len(U0)

    U = array(zeros([n+1,nv]))

    U[0,:] = U0

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
   
    m = 7 
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
    
    # Regresión lineal y cálculo de la pendiente 
    
    reg = LinearRegression().fit(log_N[0:j+1].reshape((-1, 1)),log_E[0:j+1]) 
    order = round_(abs(reg.coef_),1)

    log_N_lineal = log_N[0:j+1]
    log_E_lineal = reg.predict(log_N[0:j+1].reshape((-1, 1)))

    return [log_E, log_N, log_E_lineal, log_N_lineal, order]
