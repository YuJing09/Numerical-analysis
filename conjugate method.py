# -*- coding: utf-8 -*-
"""
Created on Wed May 13 20:43:12 2020

@author: User
"""
import numpy as np
def cgm(A,x0,b,tol):
    r0=b-A.dot(x0)
    p=r0
    infinitynorm=1
    iteration=1
    r=[r0,r0]
    x=x0
    while infinitynorm>tol:
          r[:]=r[::-1]
          xi=x.copy()
          alpha=(np.dot(r[0].T,r[0]))/(np.dot(p.T,A).dot(p))
          x=xi+alpha*p
          r[1]=r[0]-alpha*np.dot(A,p)
          beta=(np.dot(r[1].T,r[1]))/(np.dot(r[0].T,r[0]))
          p=r[1]+beta*p
          iteration+=1
          infinitynorm=np.linalg.norm(np.abs(x-xi),ord=np.inf)
          
          res=np.linalg.norm(np.dot(A,x)-b,ord=np.inf)
          
    return x,iteration,res
N=160 
A=np.zeros([N,N])
b=np.zeros([N,1])
tol=10**(-6)
x0=np.zeros([N,1])
for i in range(N):
    A[i,i]=2*(i+1)
    if i < N-1:
        A[i,i+1]=0.5*(i+1)
    if i > 0:
        A[i,i-1]=0.5*(i+1)
    b[i]=1.5*(i+1)-6
x,it,res=cgm(A,x0,b,tol)

print('iteration : %3d\nresidual norms: %.8f' % (it,res))

