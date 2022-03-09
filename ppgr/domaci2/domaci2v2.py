import numpy as np
from numpy import linalg as LA
from scipy.linalg import svd
import matplotlib.pyplot as plt
import math

#TEST PRIMERI

originali=[[-3,-1,1],[3,-1,1],[1,1,1],[-1,1,1]]
slike=[[-2,-1,1],[2,-1,1],[2,1,1],[-2,1,1]]

originali2=[[-3,-1,1],[3,-1,1],[1,1,1],[-1,1,1],[1,2,3],[-8,-2,1]]
slike2=[[-2,-1,1],[2,-1,1],[2,1,1],[-2,1,1],[2,1,4],[-16,-5,4]]

originalitest1=[[1,1,1],[5,2,1],[6,4,1],[-1,7,1]]
sliketest1=[[0,0,1],[10,0,1],[10,5,1],[0,5,1]]

originalitest2=[[1,1,1],[5,2,1],[6,4,1],[-1,7,1],[3,1,1]]
sliketest2=[[0,0,1],[10,0,1],[10,5,1],[0,5,1],[3,-1,1]]

test2=[[-1,-4,1],[-4,1,1],[-8,-5,1],[-1,0,1],[1,-2,1]]
test2p=[[3,-3,1],[5,1,1],[6,0,1],[5,2,1],[3,-13,1]]


def naivni(originali,slike):
    
    M=np.array([[originali[0][0],originali[1][0],originali[2][0]],
                [originali[0][1],originali[1][1],originali[2][1]],
                [originali[0][2],originali[1][2],originali[2][2]]])
    
    
    D=np.array([originali[3][0],originali[3][1],originali[3][2]])
    x=np.linalg.solve(M,D)
    #print(x)
    alfa=x[0]
    beta=x[1]
    gama=x[2]

    tmp1=np.array([alfa*originali[0][0],alfa*originali[0][1],alfa*originali[0][2]])
    tmp2=np.array([beta*originali[1][0],beta*originali[1][1],beta*originali[1][2]])
    tmp3=np.array([gama*originali[2][0],gama*originali[2][1],gama*originali[2][2]])
    g=np.column_stack([tmp1,tmp2,tmp3])

    N=np.array([[slike[0][0],slike[1][0],slike[2][0]],
                [slike[0][1],slike[1][1],slike[2][1]],
                [slike[0][2],slike[1][2],slike[2][2]]])
    #N1=N.transpose()
    Dprim=np.array([slike[3][0],slike[3][1],slike[3][2]])
    xp=np.linalg.solve(N,Dprim)
    
    ni=xp[0]
    psi=xp[1]
    teta=xp[2]

    tmp1p=np.array([ni*slike[0][0],ni*slike[0][1],ni*slike[0][2]])
    tmp2p=np.array([psi*slike[1][0],psi*slike[1][1],psi*slike[1][2]])
    tmp3p=np.array([teta*slike[2][0],teta*slike[2][1],teta*slike[2][2]])

    h=np.column_stack([tmp1p,tmp2p,tmp3p])
    #print(f"H={h}\n")
    f=np.dot(h,np.linalg.inv(g))
    return f

def dlt(originali,slike):
    n=len(originali)

    A=np.zeros((2*n,9))
    for i in range (0,n):
        x1=originali[i][0]
        x2=originali[i][1]
        x3=originali[i][2]

        x1p=slike[i][0]
        x2p=slike[i][1]
        x3p=slike[i][2]
        
        tmp1=np.array([0, 0, 0,-x3p*x1,-x3p*x2,-x3p*x3, x2p*x1, x2p*x2, x2p*x3 ])
        tmp2=np.array([ x3p*x1,x3p*x2,x3p*x3,0,0,0,-x1p*x1,-x1p*x2,-x1p*x3])
        #print(tmp1)
        #print(tmp2)
        A[2*i]=tmp1
        A[2*i+1]=tmp2
    
    #print(A)
    U,D,Vt=LA.svd(A)
    P=Vt[-1].reshape(3,3)
    return P

def normalizuj(originali):
    #racunacmo teziste C sa koordinatama x i y
    x=0
    y=0
    n=len(originali)
    for i in range(n):
        x=x+(originali[i][0]/originali[i][2])
        y=y+(originali[i][1]/originali[i][2])
    #prosecna vrednost-za teziste
    x=x/n
    y=y/n

    #srednje rastojanje
    d=0.0
    #matrica transliranja
    for i in range(n):
        tr1=float(originali[i][0]/originali[i][2])-x
        tr2=float(originali[i][1]/originali[i][2])-y

        d=d+math.sqrt(tr1**2+tr2**2)

    d=d/float(n)
    #racunamo lambdu
    l=float(math.sqrt(2))/d
    return np.array([[l,0,-l*x],
                     [0,l,-l*y],
                     [0,0,1]])

def dlt_mod(originali,slike):
    T=normalizuj(originali)
    Tp=normalizuj(slike)

    M=T.dot(np.transpose(originali))
    M=np.transpose(M)
    Mp=Tp.dot(np.transpose(slike))
    Mp=np.transpose(Mp)

    Ptmp=dlt(M,Mp)
    P=(LA.inv(Tp)).dot(Ptmp).dot(T)
    return P

def uporedi_dlt_naivni(originali,slike):
    print('Poredjenje naivnog i dlt algoritma')

    Pna=naivni(originali,slike)
    #Pna=Pna.round(5)
    print(f'Naivni:\n{Pna.round(5)}')

    Pdlt=dlt(originali,slike)
    #Pdlt=Pdlt.round(5)
    print(f'DLT:\n{Pdlt.round(5)}')

    Pdlt2=(Pdlt/Pdlt[0][0])*Pna[0][0]
    #Pdlt2=Pdlt2.round(5)
    print(f'Matrica DLT-a nakon sredjivanja:\n{Pdlt2.round(5)}')
    #return Pdlt.round(4)==Pna.round(4)
    
    #dobijamo isto preslikavanje

def uporedi_dlt_dlt_mod(originali,slike):
    print('Poredjenje dlt-a i modifikovanog dlt-a(sa normalizacijom)')
    Pdlt=dlt(originali,slike)
    #Pdlt=Pdlt.round(5)
    print(f'Pdlt:\n{Pdlt} ')
    Pdlt_m=dlt_mod(originali,slike)
    #Pdlt_m=Pdlt_m.round(5)
    print(f'Pdlt_modifikovan:\n{Pdlt_m}')
    Pna=naivni(originali,slike)
    #Pna=Pna.round(5)

    Pdlt_m2=(Pdlt_m/Pdlt_m[0,0])*Pdlt[0,0]
    #Pdlt_m2=Pdlt_m2.round(5)
    print(f'Pdlt_mod nakon sredjivanja:\n{Pdlt_m2}')
    #return Pdlt.round(4)==Pdlt_m.round(4)
    #nema promene(razlika)
 
def promeni_koordinate(test,testp):
    C1=np.array([[0,1, 2],
                 [-1,0,3],
                 [ 0,0,1]])

    C2=np.array([[1,-1,5],
                 [1,1,-2],
                 [ 0,0,1]])

    org=[]
    slk=[]
    
    for i in range(0,len(test)):
        org.append(np.dot(C1,test[i]))
        slk.append(np.dot(C2,testp[i]))

    #pravimo numpy matrice od org-a i slk-a
    org=np.array(org)
    slk=np.array(slk)

    #Provera kod DLT-a
    Pdlt=dlt(test,testp)
    #Pdlt=Pdlt.round(5)
    Pdltnovo=dlt(org,slk)

    #C2^-1*nAdlp*C1
    #P1 kao AdlpStari 
    P1=np.dot(np.linalg.inv(C2),Pdltnovo)
    P1=np.dot(P1,C1)
    P1=(Pdlt[0]/P1[0])*P1
    #P1=P1.round(5)

    print("Uporedjivanje DLT-a nakon promene koordinata(primene C1 i C2):")
    print(f"{P1}\n{Pdlt}") #primecujemo da se razlikuju
    
    #Provera kod modifikovanog DLT-a:
    Pdlt_mod=dlt_mod(test,testp)
    Pdlt_mod_novo=dlt_mod(org,slk)
    #Pdlt_mod=Pdlt_mod.round(5)
    
    P2=np.dot(np.linalg.inv(C2),Pdlt_mod_novo)
    P2=np.dot(P2,C1)
    P2=(Pdlt_mod[0]/P2[0])*P2
    #P2=P2.round(5)
    print("Uporedjivanje modifikovanog DLT-a nakon promene koordinata(primene C1 i C2):")
    print(f"{P2}\n{Pdlt_mod}") #nema promene
    

if __name__=='__main__':

    originali=[[-3,-1,1],[3,-1,1],[1,1,1],[-1,1,1]]
    originali2=[[-3,-1,1],[3,-1,1],[1,1,1],[-1,1,1],[1,2,3],[-8,-2,1]]

    slike=[[-2,-1,1],[2,-1,1],[2,1,1],[-2,1,1]]
    slike2=[[-2,-1,1],[2,-1,1],[2,1,1],[-2,1,1],[2,1,4],[-16,-5,4]]

    test1=[[2,0,1],[-2,1,1],[-1,-4,1],[0,2,1]]
    test1p=[[-2,1,1],[2,-1,1],[1,-2,1],[3,-1,1]]

    testprosireni=[[2,0,1],[-2,1,1],[-1,-4,1],[0,2,1],[2,2,1]]
    testprosirenip=[[-2,1,1],[2,-1,1],[1,-2,1],[3,-1,1],[-12,1,1]]
    #koordinate nastale modifikovanjem
    test2=[[-1,-4,1],[-4,1,1],[-8,-5,1],[-1,0,1],[1,-2,1]]
    test2p=[[3,-3,1],[5,1,1],[6,0,1],[5,2,1],[3,-13,1]]

    #uporedi_dlt_naivni(test1,test1p)
    #uporedi_dlt_dlt_mod(testprosireni,testprosirenip)
    #promeni_koordinate(testprosireni,testprosirenip)
    Pn=naivni(test1,test1p)
    Pn=Pn/Pn[0][0]

    Pdlt=dlt(test1,test1p)
    Pdlt=Pdlt/Pdlt[0][0]

    Pdlt_mod=dlt_mod(test1,test1p)
    Pdlt_mod=Pdlt_mod/Pdlt_mod[0][0]

    print("1.) za prve cetiri tacke(y1,y2,y3,y4)) -->(y1p,y2p,y3p,y4p)")
    print(f"Naivni:\n{Pn}")
    print(f"DLT:\n{Pdlt}")
    print(f"DLT-modifikovan:\n{Pdlt_mod}")


    Pdlt=dlt(testprosireni,testprosirenip)
    Pdlt=Pdlt/Pdlt[0][0]

    Pdlt_mod=dlt_mod(testprosireni,testprosirenip)
    Pdlt_mod=Pdlt_mod/Pdlt_mod[0][0]

    print("2.) za prve cetiri tacke(y1,y2,y3,y4,y5)) -->(y1p,y2p,y3p,y4p,y5p)")
    print(f"DLT:\n{Pdlt}")
    print(f"DLT-modifikovan:\n{Pdlt_mod}")
    
    Pdlt_modn=dlt_mod(test2,test2p)
    Pdlt_modn=Pdlt_modn/Pdlt_modn[0][0]

    print("3.)  DLP-mod za (yn1,yn2,yn3,yn4,yn5)) -->(yn1p,yn2p,yn3p,yn4p,yn5p)")
    print(Pdlt_modn)