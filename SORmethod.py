# -*- coding: utf-8 -*-
"""
Created on Thu May 26 12:57:46 2022

@author: 88697
"""

import numpy as np

def sormethod(A,b,x0,w=1,tol=10**-6):
    x=x0.copy()
    infinitynorm=1
    itera=0
    while infinitynorm>tol and 0<=itera<100:
        x1=x.copy()
        itera+=1
        for i in range(len(x0)):
            x[i]=(1-w)*x[i]+(w/A[i][i])*(b[i]-np.sum(np.fromiter((A[i][j]*x[j] for j in range(0,i)),dtype=np.float32))-np.sum(np.fromiter((A[i][j]*x[j] for j in range(i+1,len(x0))),dtype=np.float32)))

            infinitynorm=np.linalg.norm(np.abs(x1-x))
            
            res=np.linalg.norm(np.dot(A,x)-b,ord=np.inf)
    return res,itera,x
if __name__=='__main__':
    N=160
    A=np.zeros([N,N])
    for i in range(N):
        A[i,i]=2*(i+1)
        if(i+2)<=N-1:
            A[i,i+2]=0.5*(i+1)
        if(i-2)>=0:
            A[i,i-2]=0.5*(i+1)
        if (i+4)<=N-1:
            A[i,i+4]=0.25*(i+1)
        if (i-4)>=0:
            A[i,i-4]=0.25*(i+1)
    b=np.array(np.math.pi)
    b=np.tile(b,[N,1])
    x0=np.zeros([N,1])
    tol=10**-6
    w=0.6,0.8,1.,1.2,1.4
    for wi in w:
        res,iteration,x=sormethod(A,b,x0,wi,tol)
        print('omega:%.1f \n numbber of iteration :%d \n residdual norms :%.10f' % (wi,iteration,res))