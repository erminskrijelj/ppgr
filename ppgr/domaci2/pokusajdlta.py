    n=int(input("Unesi n:"))
    print(n)
    A=np.zeros((2*n,9))
    #print(A)
    for i in range (0,n):
        print(f'Unesite koordinate tacke M{i}:')
        Xtmp=[]
        x1=int(input())
        Xtmp.append(x1)
        x2=int(input())
        Xtmp.append(x2)
        x3=int(input())
        Xtmp.append(x3)
        X=np.array(Xtmp)  
        print(f'Unesite koordinate tacke M{i} prim:')
        Xptmp=[]
        x1=int(input())
        Xptmp.append(x1)
        x2=int(input())
        Xptmp.append(x2)
        x3=int(input())
        Xptmp.append(x3)
        Xp=np.array(Xptmp) 

        tmp1=np.array([0, 0, 0,-Xp[2]*X[0],-Xp[2]*X[1],-Xp[2]*X[2], Xp[1]*X[0], Xp[1]*X[1], Xp[1]*X[2] ])
        tmp2=np.array([ Xp[2]*X[0],Xp[2]*X[1],Xp[2]*X[2],0,0,0,-Xp[0]*X[0],-Xp[0]*X[1],-Xp[0]*X[2]])
        A[2*i]=tmp1
        A[2*i+1]=tmp2
        
    #Nula=np.array([[0],[0]])
    
    #Pij=LA.inv(M1).dot(Nula)
    #print(f"({A})")
    U,D,Vt=LA.svd(A)
    #print(U)
    #print(Vt[-1])
    P=Vt[-1].reshape(3,3)