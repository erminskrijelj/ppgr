import numpy as np

def ParametriKamere(T):

    T0=T[:,:-1]
    #prvi korak iz algoritma resenja problema
    if np.linalg.det(T0)<0:
        T=-T
    #TC=0 odatle trazimo c1 c2 c3 c4
    c1=np.linalg.det([[T[0][1],T[0][2],T[0][3]],
                      [T[1][1],T[1][2],T[1][3]],
                      [T[2][1],T[2][2],T[2][3]]])
    
    c2=-np.linalg.det([[T[0][0],T[0][2],T[0][3]],
                       [T[1][0],T[1][2],T[1][3]],
                       [T[2][0],T[2][2],T[2][3]]])
    

    c3=np.linalg.det([[T[0][0],T[0][1],T[0][3]],
                       [T[1][0],T[1][1],T[1][3]],
                       [T[2][0],T[2][1],T[2][3]]])
    
    c4=-np.linalg.det([[T[0][0],T[0][1],T[0][2]],
                       [T[1][0],T[1][1],T[1][2]],
                       [T[2][0],T[2][1],T[2][2]]])
    
    c1=np.round(c1/c4)
    c2=np.round(c2/c4)
    c3=np.round(c3/c4)

    C=np.array([c1,c2,c3])

    Q,R=np.linalg.qr(np.linalg.inv(T0))
    #dobijamo gornjw trougaonuu i Q koja je ortogonalna
    
    #2. korak iz algoritma resenja problema
    if R[0][0]<0:
        R[0]= - R[0]
        Q[:,0]=-Q[:,0]
    
    if R[1][1]<0:
        R[1]= - R[1]
        Q[:,1]=-Q[:,1]
    #T0=KA
    #K=T0 Q --matrica kalibracije
    if R[2][2]<0:
        R[2]= - R[2]
        Q[:,2]=-Q[:,2]
    #A je transponovano Q ili Q^-1 --matrica rotacije
    A=Q.T
    
    #K je R^-1 
    K=np.linalg.inv(R)
    #3.korak algoritma
    if K[2][2]!=1:
        K=K/K[2][2]
    # K=K.round()
    # A=A.round()

    return K,A,C

def CameraDLP(o,s):
    n=len(o)
    A=np.zeros((2*n,12))

    for i in range(0,n):
        x=o[i][0]
        y=o[i][1]
        z=o[i][2]
        t=o[i][3]

        u=s[i][0]
        v=s[i][1]
        w=s[i][2]
    
        tmp1=np.array([0,0,0,0,-w*x,-w*y,-w*z,-w*t,v*x,v*y,v*z,v*t])
        tmp2=np.array([w*x,w*y,w*z,w*t,0,0,0,0,-u*x,-u*y,-u*z,-u*t])

        A[2 * i]=tmp1
        A[2*i+1]=tmp2

    U,D,Vt=np.linalg.svd(A)

    T=Vt[-1].reshape(3,4)
    T=T.round(4)
    #matrica 3 x 4
    return T

if __name__=='__main__':

    n=4 #moj indeks je 194/2018
    M1=np.array([460,280,250,1])
    M2=np.array([50,380,350,1])
    M3=np.array([470,500,100,1])
    M4=np.array([380,630,50*n,1])
    M5=np.array([30*n,290,0,1])
    M6=np.array([580,0,130,1])

    M1p=np.array([288,251,1])
    M2p=np.array([79,510,1])
    M3p=np.array([470,440,1])
    M4p=np.array([520,590,1])
    M5p=np.array([365,388,1])
    M6p=np.array([365,20,1])

    org=np.array([M1,M2,M3,M4,M5,M6])
    slike=np.array([M1p,M2p,M3p,M4p,M5p,M6p])

    T=np.array([[5,-1-2*n,3,18-3*n],
                [0,  -1 , 5 , 21],
                [0, -1, 0 ,1]])


    # T=np.array([[1,  5 , 7 , -2],
    #             [3,  0 , 2 , 3  ],
    #             [4 , -1, 0 , 1  ]])
    print(f"n= {n}")
    print("T:\n",T)
    K,A,C=ParametriKamere(T)
    print("K :",K)
    print("A :",A)
    print("C :",C)
    print("--------------------------------------")
   
    T1=CameraDLP(org,slike)
    T1=T1/T1[0][0]
    print("T: ",T1.round(4))