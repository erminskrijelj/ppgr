import numpy as np
import math

def Euler2A(fi,teta,psi):
    
    Rx=np.array([[1,0,0],
                [0,math.cos(fi),-math.sin(fi)],
                [0,math.sin(fi), math.cos(fi)]])

    Ry=np.array([[math.cos(teta),0,math.sin(teta)],
                [0, 1, 0],
                [-math.sin(teta),0,math.cos(teta)]])

    Rz=np.array([[math.cos(psi),-math.sin(psi),0],
                [math.sin(psi),math.cos(psi),0],
                [0, 0, 1]])

    A=(Rz.dot(Ry)).dot(Rx)
    return A


def crossProduct(r,q):
    p=np.zeros((1,3))

    p[0][0]=r[1]*q[2]-r[2]*q[1]
    p[0][1]=r[2]*q[0]-r[0]*q[2]
    p[0][2]=r[0]*q[1]-r[1]*q[0]

    return p

def dotProduct(r,q):
    return r[0]*q[0] +r[1]*q[1]+r[2]*q[2]

###nastavak drugog
def axisAngle(A):
    #Proveravamo za matricu A da li je matrica rotacije
    u1=A.dot(A.T).round(5)==np.eye(3)
    
    u2=np.linalg.det(A)==1.0

    if not u1.all() or not u2:
        print("Matica A nije matrica kretanja")
        return -1,-1

    X=A-np.eye(3)

    r=np.array([X[0][0],X[0][1],X[0][2]])
    q=np.array([X[1][0],X[1][1],X[1][2]])
    
    p=crossProduct(r,q)
    #slucaj kada su r i q lin. zavsini
    if p[0][0]==0 and p[0][1]==0 and p[0][2]==0:
        q=np.array([X[2][0],X[2][1],X[2][2]])
    
    #p=crossProduct(r,q)

    normap=math.sqrt(p[0][0]**2+p[0][1]**2+p[0][2]**2)
    #dobijamo jedinicni vektor
    p=p*(1/normap)

    normar=math.sqrt(r[0]**2+r[1]**2+r[2]**2)
    
    u=r*(1/normar)

    uprim=A.dot(u)

    cosfi=dotProduct(u,uprim)
    fi=math.acos(cosfi)

    return p,fi

#3 RODRIGEZ

def rodrigez (p,fi):
    
    if (p[0][0]**2+p[0][1]**2+p[0][2]**2).round()!=1.0:
        print("VEKTOR P NIJE JEDINICNI!!!", )
        return -1
    #R=ppt+N+M
    R=p.T.dot(p) #ppt
    Q=np.eye(3)-R 
    R=R+math.cos(fi)*Q
    #N=math.cos(fi)*Q
    px=np.array([[0,-p[0][2],p[0][1]],
                 [p[0][2],0,-p[0][0]],
                 [-p[0][1],p[0][0],0]])

    R=R+math.sin(fi)*px
        #M
    return R


def A2Euler(A):
    a31=A[2][0]
    if ((A.dot(A.T).round(5))!=np.eye(3)).all():
        print("Matrica A nije ortogonalna!!!")
        return -1,-1,-1
    if np.linalg.det(A)!=1.0:
        print("Matrica A nije ortogonalna!!!")
        return -1,-1,-1

    if a31<1:
        if a31>-1:
            psi=math.atan2(A[1][1],A[0][0])
            teta=math.asin(-a31)
            fi=math.atan2(A[2][1],A[2][2])
        else:
            #Ox3=-Oz
            psi=math.atan2(-A[0][1],A[1][1])
            teta=math.pi/2
            fi=0
    else:
        psi=math.atan2(-A[0][1],A[1][1])     
        teta=-math.pi/2
        fi=0

    return fi,teta,psi


def axisAngle2Q(p,fi):

    if (p[0][0]**2+p[0][1]**2+p[0][2]**2)!=1.0:
        print("vektorr p nije jedinicni!!!")
        norma=math.sqrt(p[0][0]**2+p[0][1]**2+p[0][2]**2)
        #normalizujemo
        p=p*(1/norma)

    w=math.cos(fi/2)
    #|w|=1 --identitet tj za fi=0 --->q=+/- 1
    r=math.sin(fi/2)*p

    #print("-------------")
    #print(f"r: {r} , w: {w}")

    #q=np.array([r,w])
    
    #q=np.array([q[0][0][0],q[0][0][1],q[0][0][2],q[1]])
    q=np.array([r[0][0],r[0][1],r[0][2],w])

    return q


def Q2AngleAxis(q):

    
    if ((q[0]**2 +q[1]**2+q[2]**2+q[3]**2).round(3))!=1.0:
        print("Kvaternion nije jedinicni!!!")
        norma=math.sqrt(q[0]**2 +q[1]**2+q[2]**2+q[3]**2)
        #normalizujemo
        q=q*(1/norma)

    w=q[3]

    if w<0:
        q=-q
    
    fi=2*math.acos(w)

    if abs(w)==1:
        p=np.array([1,0,0])
    else:
        v=np.array([q[0],q[1],q[2]])
        normav=math.sqrt(v[0]**2+v[1]**2+v[2]**2)
        p=v*(1/normav)
    
    return p,fi

if __name__=='__main__':
    #1:
    print("Euler2A:")
    #A=Euler2A(math.acos(1/4),-math.asin(1/3),math.atan(2))
    #A = Euler2A(math.acos(2/8), math.atan(2/5), math.asin(1/4))
    

    #A=Euler2A(-math.acos(1/4),-math.asin(8/9),math.atan(4))
    
    #RODRIGEZ
    #p=[1.0,0.0,0.0]
    #fi=math.pi/2  -->dobicemo rotaciju oko X
    #primer
    print(f"PRIMER: {math.acos(2/7)} , {math.atan(1/4)} , {math.asin(1/5)}",)
    A = Euler2A(math.acos(2/7), math.atan(1/4), math.asin(1/5))
    

    print(A)
    print()
    #2:
    print("AxisAngle:")

    p,fi1=axisAngle(A)
    print("prava:  ",p)
    print("ugao: ",round(fi1,7))
    #izlaz 1,-1 oznacava da A nije matrica rotacije
    print()

    #3:RODRIGEZ
    print("Rodrigez:")
    R=rodrigez(p,fi1)
    print(R)
    print()

    #4:
    print("A2Euler:")
    fi,teta,psi=A2Euler(A)
    print("ugao fi: ",fi)
    print("ugao teta: ",teta)
    print("ugao psi: ",psi)
    print()

    #5:
    print("AxisAngle2Q:")
    q=axisAngle2Q(p,fi1)
    print("q: ",q)
    print()

    #6
    print("Q2AngleAxis:")
    p,fi=Q2AngleAxis(q)
    print(f"p: {p}  fi: {fi}")     
    #<----2