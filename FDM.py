# -*- coding: utf-8 -*-
"""
Created on Thu May 26 13:23:39 2022

@author: 88697
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
fig=plt.figure(figsize=(8,8))
ax=Axes3D(fig)
n=11
m=11
xi=np.linspace(0,2,m)
yi=np.linspace(0,1,n)
k=1/(n-1)
h=2/(m-1)
g=k/h
def f(x,y,k):
    u=-k**2*(x**2+y**2)*np.exp(x*y)
    return u
def exactsol(x,y):
    return np.exp(x*y)
linx,liny=np.meshgrid(xi,yi)
lx=np.reshape(linx,[-1,1])
ly=np.reshape(liny,[-1,1])
lxy=np.concatenate((lx,ly),axis=-1)
lxyin=np.reshape(lxy,[n,m,2])
lxyin=lxyin[1:-1]
lxyin=lxyin[:,1:-1]
lxyin=lxyin.reshape([(m-2)**2,2])
fxy=f(lxyin[:,0],lxyin[:,1],k)

w=3/(2+((4-(np.cos((np.pi/(m-1)))+np.cos(np.pi/(n-1)))**2))**0.5)
x=np.zeros([(m-2)**2,1])
tol=10**-6
infinitynorm=1
res=1
inte=0

while infinitynorm>tol:
    print(infinitynorm)
    inte+=1
    x0=x.copy()
    fxyt=np.zeros_like(fxy)
    for i in range(len(lxyin)):
        x[i]=(1-w)*x[i]
        if i%(m-2)==0:
            fxyt[i]+=g**2
        else :
            x[i]+=w*(g**2)*x[i-1]/(2*(g**2+1))
        if i-(m-2)<0:
            fxyt[i]+=1
        else: 
            x[i]+=w*x[i-(m-2)]/(2*(g**2+1))
        if i%(m-2)==(m-3):
            fxyt[i]+=np.exp(2*lxyin[:,1][i])*g**2
            
        else:
            x[i]+=w*(g**2)*x[i+1]/(2*(g**2+1))
        if i+(m-2)>len(lxyin)-1:
            fxyt[i]+=np.exp(lxyin[:,0][i])
          
        else:
            x[i]+=w*x[i+(m-2)]/(2*(g**2+1))
        fxyt[i]+=fxy[i]
        
        x[i]+=w*fxyt[i]/(2*(g**2+1))
    infinitynorm=np.linalg.norm(np.abs(x-x0),ord=np.inf) 
x1=x.reshape([m-2,n-2])
ones=np.ones([1,m-2])
x1=np.concatenate((ones,x1),axis=0)
ones=np.ones([m-1,1])
x1=np.concatenate((ones,x1),axis=1)
y2=np.exp((xi))
x1=np.concatenate((x1,y2[None,:-1]),axis=0)
x2=np.exp((2*yi))
x1=np.concatenate((x1,x2[:,None]),axis=1)
sol=exactsol(lx,ly)
ax.plot_surface(ly.reshape([m,n]),lx.reshape([m,n]),x1,rstride=2,cstride=2,cmap=plt.get_cmap('rainbow'))
plt.show()
for j in range(len(x)):
    print("x:%f y:%f approximation:%f solution:%f"% (lx[j],ly[j],x1.reshape([-1,1])[j],sol[j]))
      
