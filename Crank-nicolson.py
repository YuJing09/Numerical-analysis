# -*- coding: utf-8 -*-
"""
Created on Thu May 26 14:35:42 2022

@author: 88697
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
fig=plt.figure(figsize=(8,8))
ax=Axes3D(fig)
N=30
y=np.arange(51)*0.01
x=np.linspace(0,2,N)
x=x[1:-1]
linx,liny=np.meshgrid(x,y)
A=np.zeros([N-2,N-2])
B=np.zeros([N-2,N-2])
tol=10**-7
h=2/(N-1)
k=0.01
r=k/(2*h**2)
def exasol(t,x):
    return np.reshape(np.exp(-1*t*(np.pi)**2/4)*np.sin(np.pi*x/2),[-1,1])
def ut0(x):
    return np.reshape(np.sin(np.pi*x/2),[-1,1])

for i in range(len(A)):
    A[i,i]=1+2*r
    if i>0:
        A[i,i-1]=-r
    if i<(N-2)-1:
        A[i,i+1]=-r
for i in range(len(B)):
    B[i,i]=1-2*r
    if i>0:
        B[i,i-1]=r
    if i<N-2-1:
        B[i,i+1]=r
w=1
u=ut0(x)
b=B.dot(ut0(x))
infinitynorm=1
x0=np.zeros([N-2,1])
while k<0.51:
    while infinitynorm>tol: 
        xi=x0.copy()
        for i in range(len(x0)):
            x0[i]=(1-w)*x0[i]+w*(b[i])/(A[i,i])
            if i>0:
                x0[i]-=w*A[i,i-1]*x0[i-1]/A[i,i]
            if i<N-2-1:
                x0[i]-=w*A[i,i+1]*x0[i+1]/A[i,i]
            infinitynorm=np.linalg.norm(x0-xi,ord=np.inf,axis=0)[0]
            
            
    u=np.concatenate((u,x0),axis=1)
    k+=0.01
    b=B.dot(x0)
    infinitynorm=1
    x0=np.zeros([N-2,1])
    print(k)
ureal=ut0(x)
for i in range(50):
    ureal=np.concatenate((ureal,exasol(0.01*(i+1),x)),axis=1)
ax.plot_surface(linx,liny,u.swapaxes(1,0),rstride=2,cstride=2,cmap=plt.get_cmap('rainbow'))
plt.show()
fig=plt.figure(figsize=(8,8))
ax=Axes3D(fig)
ax.plot_surface(linx,liny,ureal.swapaxes(1,0),rstride=2,cstride=2,cmap=('rainbow'))
       
            
            