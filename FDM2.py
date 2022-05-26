# -*- coding: utf-8 -*-
"""
Created on Thu May 26 15:32:11 2022

@author: 88697
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from SORmethod import sormethod
import time

def f(x,y,e):
    return (2*e*np.exp(-x/e)-x+1)*(y-y**2)+2*e**2*(x-1)*(np.exp(-x/e)-1)

fig=plt.figure(figsize=(8,8))

ax=Axes3D(fig)
N=100
e=10**-4
x=np.linspace(0,1,N)
y=np.linspace(0,1,N)
h=1/(N-1)
lx,ly=np.meshgrid(x,y)
lxx=lx.reshape([-1,1])
lyy=ly.reshape([-1,1])
allxy=np.concatenate((lxx,lyy),axis=1)
mask0=np.greater(allxy,0)
mask1=np.greater(1,allxy)
mask=np.concatenate([mask0,mask1],axis=1).all(axis=1)
allxymask=allxy[mask]
fxy=f(allxymask[:,0],allxymask[:,1],e)*h**2
x=np.zeros([(N-2)**2,1])
t0=time.time()

w=1
tol=10**-15
itera=0

t1=time.time()
infn=1
while infn>tol:
    itera+=1
    x0=x.copy()
    for i in range(len(x)):
        x[i]=(1-w)*x[i]+w*fxy[i]/(4*e**2+h**2)
        if(i%(N-2))!=0:
            x[i]-=w*(-e**2)*x[i-1]/(4*e**2+h**2)
        if(i-(N-2)>=0):
            x[i]-=w*(-e**2)*x[i-(N-2)]/(4*e**2+h**2)
        if (i+N-2)<=(len(x)-1):
            x[i]-=w*(-e**2)*x[i+N-2]/(4*e**2+h**2)
        if i%(N-2)!=(N-3):
            x[i]-=w*(-e**2)*x[i+1]/(4*e**2+h**2)
    infn=np.linalg.norm(np.abs(x-x0),ord=np.inf)
x=x.reshape([N-2,N-2])
zero=np.zeros([N-2,1])
u=np.concatenate([zero,x,zero],axis=1)
zero=np.zeros([1,N])
u=np.concatenate([zero,u,zero],axis=0)
ax.plot_surface(lx,ly,u,rstride=2,cstride=2,cmap=plt.get_cmap('rainbow'))
ax.contourf(lx,ly,u,zdir='z',offset=-0.1,cmap=plt.get_cmap('rainbow'))










