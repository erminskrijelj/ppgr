import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
%matplotlib notebook
#from mpl_toolkits.mplot3d import Axes3D
#tacke leveslike
L1=[816,111,1]
L2=[950,164,1]
L3=[987,125,1]
L4=[854,78,1]
L5=[790,303,1]
L6=[911,358,1]
L7=[952,317,1]
L8=[0,0,1] # nevidljiva tacka
L9=[320,345,1]
L10=[452,369,1]
L11=[511,271,1]
L12=[387,249,1]
L13=[364,559,1]
L14=[477,582,1]
L15=[526,487,1]
L16=[0,0,1]
L17=[137,549,1]
L18=[435,761,1]
L19=[815,382,1]
L20=[547,252,1]
L21=[175,656,1]
L22=[449,860,1]
L23=[806,487,1]
L24=[0,0,1]


R1=[910,444,1]
R2=[807,562,1]
R3=[917,613,1]
R4=[1012,492,1]
R5=[0,0,1]
R6=[773,771,1]
R7=[864,824,1]
R8=[956,700,1]
R9=[298,71,1]
R10=[251,119,1]
R11=[370,139,1]
R12=[413,89,1]
R13=[0,0,1]
R14=[288,326,1]
R15=[395,343,1]
R16=[431,287,1]
R17=[0,0,1]
R18=[137,319,1]
R19=[525,531,1]
R20=[742,347,1]
R21=[0,0,1]
R22=[160,425,1]
R23=[529,645,1]
R24=[735,453,1]

#rekonstrukcija skrivenih tacaka

L8=np.cross(np.cross(np.cross(np.cross(L5,L6),np.cross(L1,L2)),L7),
            np.cross(np.cross(np.cross(L2,L6),np.cross(L3,L7)),L4))
        
L8=L8/L8[2]
L8=[np.round(L8[0]),np.round(L8[1]),L8[2]]

L16=np.cross(np.cross(np.cross(np.cross(L13,L14),np.cross(L9,L10)),L15),
            np.cross(np.cross(np.cross(L9,L13),np.cross(L10,L14)),L12))
      
L16=L16/L16[2]
L16=[np.round(L16[0]),np.round(L16[1]),L16[2]]


L24=np.cross(np.cross(np.cross(np.cross(L21,L22),np.cross(L19,L20)),L23),
            np.cross(np.cross(np.cross(L17,L21),np.cross(L18,L22)),L20))
        
L24=L24/L24[2]
L24=[np.round(L24[0]),np.round(L24[1]),L24[2]]

R5=np.cross(np.cross(np.cross(np.cross(R1,R2),np.cross(R3,R4)),R6),
            np.cross(np.cross(np.cross(R2,R6),np.cross(R3,R7)),R1))
      
R5=R5/R5[2]
R5=[np.round(R5[0]),np.round(R5[1]),R5[2]]

R13=np.cross(np.cross(np.cross(np.cross(R10,R14),np.cross(R11,R15)),R9),
           np.cross(np.cross(np.cross(R9,R10),np.cross(R11,R12)),R14))
       
R13=R13/R13[2]
R13=[np.round(R13[0]),np.round(R13[1]),R13[2]]

##
R17=np.cross(np.cross(np.cross(np.cross(R19,R20),np.cross(R23,R24)),R18),
            np.cross(np.cross(np.cross(R18,R19),np.cross(R22,R23)),R20))      
R17=R17/R17[2]
R17=[np.round(R17[0]),np.round(R17[1]),R17[2]]


R21=np.cross(np.cross(np.cross(np.cross(R19,R20),np.cross(R23,R24)),R22),
             np.cross(np.cross(np.cross(R18,R19),np.cross(R22,R23)),R24))
       
R21=R21/R21[2]
R21=[np.round(R21[0]),np.round(R21[1]),R21[2]]

L=[L1,L2,L3,L4,L5,L6,L7,L8,L9,L10,L11,L12,L13,L14,L15,L16,L17,L18,L19,L20,L21,L22,L23,L24]
R=[R1,R2,R3,R4,R5,R6,R7,R8,R9,R10,R11,R12,R13,R14,R15,R16,R17,R18,R19,R20,R21,R22,R23,R24]

#8 tacaka na osnovu kojih pravimo fundamentalnu matricu FF
LL=np.array([L1,L2,L3,L4,L6,L7,L9,L10])
RR=np.array([R1,R2,R3,R4,R6,R7,R9,R10])


def jednacina(l,d):
    a1=l[0]
    a2=l[0]
    a3=l[0]

    b1=d[0]
    b2=d[0]
    b3=d[0]

    return np.array([a1*b1,a2*b1,a3*b1,
                     a1*b2,a2*b2,a3*b2,
                     a1*b3,a2*b3,a3*b3])

# def MapThread(l,d):
#     matrica=jednacina(l[0],d[0])
#     n=len(l)
#     for i in range(1,n):
#         tmp=jednacina(l[i],d[i])
#         matrica=np.concanate((matrica,tmp),axis=0)

#     return matrica

#svd dekompozicija
###TRAZIMO EPIPOLOVE

def epipolovi(FF):
    U,D,Vt=np.linalg.svd(FF)
    #F*e1=0
    # print("U:\n",U)
    # print("D:\n",D)
    # print("Vt\n:",Vt)
    print("################# EPIPOLOVI ###############")

    print("e1:")
    e1=Vt[:][-1]
    #print(e1)
    e1=(1/e1[2])*e1 #afinee
    print(e1)
    
    #FT *e2=0
    print("e2:")
    e2=U.T[:][-1]
    #treca kolona tj treca vrsta od U.T
    #print(e2)
    e2=(1/e2[2])*e2 #afinee
    print(e2)

    #zarad postizanje uslova det(FF)=0  
    D=np.array([[D[0],0,0],
                [0,D[1],0],
                [0,0,0]])
    #np.diag([1,1,0])*D

    print("D:\n",D)

    FF1=U.dot(D).dot(Vt)
    # print(FF)
    print("FF1:\n",FF1)

    #FF i FF1 imaju iste matrice U i Vt samim tim i iste epipolove
    print("det FF:\n",np.linalg.det(FF))
    print("det FF1:\n",np.linalg.det(FF1))
    return e2,FF1

#Triangulacija
#kanonska matrica kamere
def triangulacija(e2,FF1):
    print("#################TRIANGULACIJA###############")
    T1=np.array([[1,0,0,0],
                 [0,1,0,0],
                 [0,0,1,0]])

    print("T1:") #kanonska matrica kamere
    print(T1)
    #matrica vektorskog mnozenja
    E2=np.array([[0,-e2[2],e2[1]],
                 [e2[2],0,-e2[0]],
                 [-e2[1],e2[0],0]])
    print("E2\n",E2)

    #matrica kamere T2
    tmp=E2.dot(FF1)
    T2=np.array([[tmp[0][0],tmp[0][1],tmp[0][2],e2[0]],
                 [tmp[1][0],tmp[1][1],tmp[1][2],e2[1]],
                 [tmp[2][0],tmp[2][1],tmp[2][2],e2[2]]])

    print("T2:\n",T2)
    return T1,T2

#za svaku jednacinu dobijmao sistem od 4 jednacine sa 4 homogene nepoznate
def jednacine(l,r,T1,T2):

    return np.array([l[1]*T1[2]-l[2]*T1[1],
                    -l[0]*T1[2]+l[2]*T1[0],
                     r[1]*T2[2]-r[2]*T2[1],
                    -r[0]*T2[2]+r[2]*T2[0]])

def tridkoordinate(l,r,T1,T2):
    U,D,Vt=np.linalg.svd(jednacine(l,r,T1,T2))

    P=Vt[-1] #poslednja vrsta svd 
    P=P/P[3] #iz homogenih pravimo afine
    return P[:-1] 

def mapThread2(l,d,T1,T2):
    rekonstruisane=[]
    for i in range(len(l)):
        rekonstruisane.append(tridkoordinate(l[i],d[i],T1,T2))
        print(f"{i}. rekonstrukcija:\n{(tridkoordinate(l[i],d[i],T1,T2))}")
    return rekonstruisane

if __name__=='__main__':

    jednacine8=[]
    for i in range(8):
        jednacine8.append(jednacina(LL[i],RR[i]))

    U,D,Vt=np.linalg.svd(jednacine8)
    FF=np.zeros((3,3))
    #FF=(Vt[:][-1]) # poslednja kolona 
    Vt=Vt[:][-1]
    for i in range(3):
        for j in range(3):
            FF[i][j]=Vt[3*i+j]
    #koeficijenti matrice FF su poslednja kolona matrice Vt
    print("Fundamentalna matrica:")
    print(FF)

    print("Determinta Fundamentalne matrica:")
    print(np.linalg.det(FF))

    #C-check proveravamo da li R.T F L=0
    C=np.ones((8,1))
    for i in range(8):
        C[i][0]=(RR[i].T.dot(FF).dot(LL[i]))
    print("Provera uslova yT F x = 0:\n",C.T)

    #print(np.linalg.det(FF))
    # treba da bude blizu nulee

    e2,FF1=epipolovi(FF)
    
    T1,T2=triangulacija(e2,FF1)
    
    ljedn=jednacine(L2,R2,T1,T2)
    mjedn=np.matrix(ljedn)
    print("MAtrica iz jednacina:\n",mjedn)
    #za L2 i R2
    
    #tridkoordinatte nam vraca 3D za date tacke
    M=tridkoordinate(L2,R2,T1,T2)
    print("Poslednja vrsta Vt matrice iz SVD-a:\n",M)

    rekonstruisane=mapThread2(L,R,T1,T2)

    #Mnozimo z koordinate sa nekoliko stotina zato sto nismo radili normalizaciju
    tmp=np.eye(3)
    tmp[2][2]=400    
    rek400=np.zeros((len(L),3))
   
    for i in range(len(L)):   
        rek400[i]=tmp.dot(rekonstruisane[i])
        print(rek400[i])
        #rek400[i]=tmp.dot(rekonstruisane[i][:,0])
    #print(rek400)

# ####ISCRTAVANJE 
    keks=np.array([[1,2],[2,3],[3,4],[4,1],
                   [5,6],[6,7],[7,8],[8,5],
                   [1,5],[2,6],[3,7],[4,8]])

    caj=np.array([[9,10],[10,11],[11,12],[12,9],
                  [13,14],[14,15],[15,16],[16,13],
                  [9,13],[10,14],[11,15],[12,16]])

    modem=np.array([[17,18],[18,19],[19,20],[20,17],
                    [21,22],[22,23],[23,24],[24,21],
                    [17,21],[18,22],[19,23],[20,24]])

    
    fig = plt.figure()
    ax = plt.axes(projection='3d')

    for ivica in keks:
        ax.plot3D([rek400[ivica[0] - 1][0],rek400[ivica[1] - 1][0]], [rek400[ivica[0] - 1][1], rek400[ivica[1] - 1][1]], [rek400[ivica[0] - 1][2], rek400[ivica[1] - 1][2]], 'red')

    for ivica in caj:
        ax.plot3D([rek400[ivica[0] - 1][0],rek400[ivica[1] - 1][0]], [rek400[ivica[0] - 1][1], rek400[ivica[1] - 1][1]], [rek400[ivica[0] - 1][2], rek400[ivica[1] - 1][2]], 'blue')

    for ivica in modem:
        ax.plot3D([rek400[ivica[0] - 1][0],rek400[ivica[1] - 1][0]], [rek400[ivica[0] - 1][1], rek400[ivica[1] - 1][1]], [rek400[ivica[0] - 1][2], rek400[ivica[1] - 1][2]], 'green')
