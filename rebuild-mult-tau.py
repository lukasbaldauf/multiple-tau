import sys
import numpy as np
from numba import njit
# usage
# pytrgon rebuild-mult-tay.py <Temp> <Volume> <Outfile>


@njit
def push_D(D,w,k):
    for j in range(p-1,0,-1):
        D[k,j]=D[k,j-1]
    D[k,0]=w

@njit
def add(D,C,N,A,M,kmax,w,k=0):
    if k == S:
        return
    if k>kmax[0]:
        kmax[0]=k
    push_D(D,w,k)
    A[k]+=w
    M[k]+=1
    if M[k]==m:
        add(D,C,N,A,M,kmax,A[k]/m,k+1)
        A[k]=0
        M[k]=0
    if k==0:
        for j in range(p):
            if D[k,j]>-1e10:
                C[k,j]+=D[k,0]*D[k,j]
                N[k,j]+=1
    else:
        for j in range(p_m,p):
            if D[k,j]>-1e10:
                C[k,j]+=D[k,0]*D[k,j]
                N[k,j]+=1

def rebuild(C,N):
    corr=[]
    time=[]
    
    ci=0
    for t in range(p):
        corr.append(C[0,t]/N[0,t])
        time.append(ci)
        ci+=1
    
    
    for s in range(1,S):
        for i in range(p_m,p):
            if N[s,i]>0:
                corr.append(C[s,i]/N[s,i])
                time.append(i*m**s)
                ci+=1
    return time,corr



m=2
p=16
S=40
p_m=int(p/m)

print(*sys.argv)

T=float(sys.argv[1]) # K
V_w=float(sys.argv[2])
writefile=sys.argv[3]
V=V_w*(1e-9)**3 # nm^3 -> m^3
kb=1.38064852e-23

corrs=[]
for iters in range(1,10):
    if iters not in [1,5,9]: # consider only off diagonal elements
        D=np.load(f"D-mat{iters}.npy")
        C=np.load(f"C-mat{iters}.npy")
        N=np.load(f"N-mat{iters}.npy")
        A=np.load(f"A-mat{iters}.npy")
        M=np.load(f"M-mat{iters}.npy")
        kmax=np.load(f"kmax-mat{iters}.npy")
        t,corr=rebuild(C,N)
        corrs.append(corr)
        print(f"Matrix {iters} rebuilt")
        
        
corr_av=V/(kb*T)*np.mean(corrs,axis=0)
np.savetxt(writefile,np.c_[t,corr_av],header='Temp: '+str(T)+'\tVol: '+str(V_w)+'\n\t t \t\t G(t)',fmt=['%15.3f','%15.4e'])
print(f'G(t) written to {writefile}.')
