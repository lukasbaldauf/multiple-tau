"""
A multiple-tau algorithm for calculating time correlation functions
on the fly, as described in 
J. Chem. Phys. 133, 154103 (2010); http://dx.doi.org/10.1063/1.3491098

The code is rewritten for Python from the source below
https://blogs.upm.es/compsoftmatter/software/multiple-tau-correlator/
"""


import sys
import numpy as np
from numba import njit


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

x = np.loadtxt(sys.argv[1],comments=["#","@"])

for iters in range(1,10):
    try:
        D=np.load(f"D-mat{iters}.npy")
        C=np.load(f"C-mat{iters}.npy")
        N=np.load(f"N-mat{iters}.npy")
        A=np.load(f"A-mat{iters}.npy")
        M=np.load(f"M-mat{iters}.npy")
        kmax=np.load(f"kmax-mat{iters}.npy")
        print("Matrices loaded.")
    
    except FileNotFoundError:
        D=np.ones([S,p])*(-1e10) # data array
        C=np.zeros([S,p]) # correlator array
        N=np.zeros([S,p]) # counter
        A=np.zeros(S) # accumulator
        M=np.zeros(S) # counter
        kmax=np.zeros(1,dtype=int)
        print("Initialized matrices.")
    
    for val in x[:,iters]:
        add(D,C,N,A,M,kmax,val)
    
    
    np.save(f"D-mat{iters}",D)
    np.save(f"C-mat{iters}",C)
    np.save(f"N-mat{iters}",N)
    np.save(f"A-mat{iters}",A)
    np.save(f"M-mat{iters}",M)
    np.save(f"kmax-mat{iters}",kmax)
    print(f"DONE run {iters}. Saved files to X-mat{iters}.npz")

