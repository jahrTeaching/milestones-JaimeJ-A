
from numpy import zeros, float64

#----- MODULO EDO(Ecuaciones diferenciales ordinarias)-----#



def Cauchy(F,t,U0,E_Temporal):

    n, nv = len(t)-1,len(U0)

    U = zeros((nv,n+1), dtype=float64)

    U[:,0] = U0

    for i in range(n):

        U[:,i+1] = E_Temporal(U[:,i],t[i+1] - t[i],t[i],F)

    return U





















