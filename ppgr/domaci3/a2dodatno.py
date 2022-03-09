    if a31==1:
        teta=-math.pi/2
        fi=0
        psi=math.atan2(-A[0][1],A[1][1])
    elif a31==-1:
        teta=math.pi/2
        fi=0
        psi=math.atan2(-A[0][1],A[1][1])
    else:
        teta=math.asin(-a31)
        fi=math.atan2(A[2][1],A[2][2])
        psi=math.atan2(A[1][0],A[0][0])