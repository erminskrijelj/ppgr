import numpy as np

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
    #matrica 3 x 4
    return T

def afinizuj(x):
    x[0]=x[0]/x[2]
    x[1]=x[1]/x[2]
    x[2]=1
    
def zaokruzi(x):
    x[0]=round(x[0])
    x[1]=round(x[1])
    x[2]=round(x[2])

def cross(a, b):
    c = [a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2],a[0]*b[1] - a[1]*b[0]]
    return c

def nevidljivo(p1,p2,p3,p5,p6,p7,p8):
    p1.append(1)
    p2.append(1)
    p3.append(1)
    p5.append(1)
    p6.append(1)
    p7.append(1)
    p8.append(1)
#
    xb1=cross((cross(p2,p6)),(cross(p1,p5)))
    xb2=cross((cross(p2,p6)),(cross(p3,p7)))
    xb3=cross((cross(p1,p5)),(cross(p3,p7)))
    xb=[]
    xb.append(round((xb1[0]+xb2[0]+xb3[0])/3.0))
    xb.append(round((xb1[1]+xb2[1]+xb3[1])/3.0))
    xb.append(round((xb1[2]+xb2[2]+xb3[2])/3.0))

    yb=cross((cross(p5,p6)),(cross(p7,p8)))
    p4=cross((cross(p8,xb)),(cross(p3,yb)))
    #print(p4)
    afinizuj(p4)
#    print(p4)
    zaokruzi(p4)
# Ako zelimo ispis u afinim koordinatama  
    p4.remove(1)
    return p4



if __name__=='__main__':

    #originali su one koordinate koje su merene na sceni 
    M1=np.array([5,160,215,1])
    M2=np.array([15,310,0,1])
    M3=np.array([200,190,0,1])
    M4=np.array([140,270,145,1])
    M5=np.array([90,0,25,1])
    M6=np.array([290,0,25,1])

    M1p=np.array([768,497,1])
    M2p=np.array([1055,1022,1])
    M3p=np.array([639,1154,1])
    M4p=np.array([914,761,1])
    M5p=np.array([410,908,1])
    M6p=np.array([70,1107,1])

    org=np.array([M1,M2,M3,M4,M5,M6])
    slike=np.array([M1p,M2p,M3p,M4p,M5p,M6p])

    T1=CameraDLP(org,slike)
    T1=T1/T1[0][0]
    print("T: \n",T1)
