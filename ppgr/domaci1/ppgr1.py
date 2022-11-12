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

if __name__ == '__main__':

    p1=[501,404]  
    p2=[285,673]
    p3=[110,542]
    p5=[542,193]
    p6=[267,436]
    p7=[44,312]
    p8=[353,130]

    p4=nevidljivo(p1,p2,p3,p5,p6,p7,p8)
    print(p4)  #dobijeno resenje [349,339]


