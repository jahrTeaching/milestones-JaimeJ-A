
from numpy import array, reshape, zeros
from numpy.linalg import norm

#----- MODULO N-CUERPOS -----#


# N-Cuerpos
def N_Bodies(U,t):

    #(Nb, Nc) = (2,3)                # 2 cuerpos  
    #(Nb, Nc) = (2,3)                # 3 cuerpos
    (Nb, Nc) = (4,3)                 # 4 cuerpos

    U_s = reshape( U, (Nb,Nc,2)) 
    F = zeros(len(U))
    F_s = reshape( F, (Nb,Nc,2)) 
    
    r = reshape(U_s[:,:,0], (Nb,Nc))    # Posiciones
    v = reshape(U_s[:,:,1], (Nb,Nc))    # Velocidades
    
    drdt = reshape(F_s[:,:,0], (Nb,Nc)) # Velocidades
    dvdt = reshape(F_s[:,:,1], (Nb,Nc)) # Aceleraciones
    
    for i in range(Nb):
        drdt[i,:] = v[i,:]
        for j in range(Nb):
            if j != i:
                d = r[j,:] - r[i,:]
                dvdt[i,:] = dvdt[i,:] + d[:]/( norm(d)**3)
            
    return F
# ------------































