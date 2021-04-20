# -*- coding: utf-8 -*-
"""
Created on Thu Apr 15 10:27:28 2021

@author: henri and julien
"""
""" The goal of this code is to compute the velocity around an airfoil with the linear vortex sheet
panel method"""

import numpy as np
from matplotlib import pyplot as plt


def y_e(x,c):
    return 0.278536*(np.sqrt(x/c)) - 0.148567*(x/c) + 0.006397*(x/c)**2 - 0.22098*(x/c)**3 + 0.081084*(x/c)**4
def y_i(x,c):
    return - 0.190361*(np.sqrt(x/c)) + 0.161628*(x/c) - 0.341176*(x/c)**2 + 0.897631*(x/c)**3 - 0.531252*(x/c)**4

c = 1 
#x = np.linspace(0,1,1000)
#y = y_e(x,c)
#yy = y_i(x,c)
#plt.figure(figsize=(15,4))
#plt.plot(x,yy)
#plt.plot(x,y)
#plt.show()



"""Discretization of the airfoil"""
class Panel:
    def __init__(self, xi, xf, yi, yf):
        self.xi=xi
        self.yi=yi
        self.xf=xf
        self.yf=yf
        
        self.length = np.sqrt( (xf-xi)**2 + (yf-yi)**2)
        
        self.xc = (xi+xf)/2
        self.yc = (yi+yf)/2
        
        if xf-xi > 0: #upper side
            self.angle = np.arccos(-(yf-yi)/self.length) #angle defined between the normal to the panel and the x-axis
        else: #lower side
            self.angle = -np.arccos(-(yf-yi)/self.length)
            
        self.nx = np.cos(self.angle)
        self.ny = np.sin(self.angle)
        
        
def discretization(n):
    #n*2 -1 panels
    #the panels are constructed clockwise
    x1=np.linspace(1,0,n,endpoint=False)
    y1=y_i(x1,c)
    x2=np.linspace(0,1,n)
    y2=y_e(x2,c)
    x=np.append(x1,x2)
    y=np.append(y1,y2) 
    
    plt.plot(x,y,'-o') 
    
    panels = np.empty(n*2-1, dtype=object)
    for i in range(0, np.size(x)-1):
        panels[i] = Panel(x[i],x[i+1],y[i],y[i+1])
    N=2*n-1
    return N,x,y,panels
    
N,x,y,panels = discretization(30)

gamma=np.zeros(len(x)) #gamma[0] must be equal to -gamma[-1]


###print for debugging###
for i in range(0,np.size(panels)): 
#    print(i)
#    print("xi",panels[i].xi)
#    print(panels[i].length)
#    print(np.rad2deg(panels[i].angle))
#    print(i)
#    print("nx",panels[i].nx)
#    print("ny",panels[i].ny)
        
        
    
    
            
        
        
        
        
        
        
    
