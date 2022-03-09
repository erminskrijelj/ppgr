import numpy as np
from numpy import linalg as LA

#print(np.__version__)

if __name__ =='__main__':
    M=np.array([[-3,-1,1],[3,-1,1],[1,1,1]],dtype=np.float64)
    #print("M={}".format(M))
    M1=M.transpose()
    #print("M1={}".format(M1))
    
    d=np.array([[-1],[1],[1]])
    x=LA.inv(M1).dot(d)
    #print(x)
    [alfa]=x[0]
    #alfa=round(alfa,5)
    [beta]=x[1]
    #beta=round(beta,5)
    [gama]=x[2]
    #gama=round(gama,5)
    #print(P)
    P=([alfa,beta,gama],[alfa,beta,gama],[alfa,beta,gama])
    g=M1
    g=g*P
    
    #DODATO RADI ISPISA
    #g=g*3

    g=g.round(4)
    print(g)
    #print(alfa,beta,gama)

    N=np.array([[-2,-1,1],[2,-1,1],[2,1,1]],dtype=np.float64)

    N1=N.transpose()
    #print("M1={}".format(M1))
    
    dprim=np.array([[-2],[1],[1]])
    xprim=LA.inv(N1).dot(dprim)

    [ni]=xprim[0]
    [psi]=xprim[1]
    [ro]=xprim[2]
    #print(ni,psi,ro)
    print()
    P2=([ni,psi,ro],[ni,psi,ro],[ni,psi,ro])
    h=N1
    h=h*P2
    h=h.round(4)
    
    print(f"H={h}\n")
    ginv=LA.inv(g)
    #print(ginv)    
    #DODATO RADI ISPISA
    #ginv=ginv*6
    
    f=h.dot(ginv)
    f=f.round(4)
    print(f"F:{f}\n")

    #######################################################################


