# -*- coding: utf-8 -*-
"""
Created on Thu Apr 15 10:27:28 2021

@author: henri and julien
"""
""" The goal of this code is to compute the velocity around an airfoil with the linear vortex sheet
panel method"""

import numpy as np
import scipy as sc
from matplotlib import pyplot as plt


def y_e(x,c):
    return 0.278536*(np.sqrt(x/c)) - 0.148567*(x/c) + 0.006397*(x/c)**2 - 0.22098*(x/c)**3 + 0.081084*(x/c)**4
def y_i(x,c):
    return - 0.190361*(np.sqrt(x/c)) + 0.161628*(x/c) - 0.341176*(x/c)**2 + 0.897631*(x/c)**3 - 0.531252*(x/c)**4

c = 1 

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
        
#        if xf-xi > 0: #upper side
#            self.angle = np.arccos(-(yf-yi)/self.length) #angle defined between the normal to the panel and the x-axis
#        else: #lower side
#            self.angle = -np.arccos(-(yf-yi)/self.length)
#            
#        self.nx = np.cos(self.angle)
#        self.ny = np.sin(self.angle)
        
#        self.angle = (np.arctan(np.abs(yf-yi)/np.abs(xf-xi)))
#        self.nx = np.sin(self.angle)*np.sign(xf-xi)
#        self.ny = np.cos(self.angle)*np.sign(yf-yi)
        
        
        
#        self.normal = np.array([self.nx,self.ny]) #y-axis of the panel
        
        self.longitudinal = np.array([self.xf-self.xi,self.yf-self.yi]) #x-axis of the panel
        
        self.normal = np.array([-self.longitudinal[1],self.longitudinal[0]])/self.length #y-axis of the panel
        
def discretization(n):
    #(n-1)*2 panels
    theta = np.linspace(0,0.95*np.pi,n)
    x1= (1-np.cos(theta))/(1-np.cos(theta[-1]) )
    x1=np.flip(x1)
    x2=np.flip(x1)
    x2=np.delete(x2,0)
    
    y1=y_i(x1,c)
    y2=y_e(x2,c)
    
    x=np.append(x1,x2)
    y=np.append(y1,y2)
    
#    print(x)
#    print(y)
    
    print("size x",np.size(x))
    
    print("size y",np.size(y))
    
    print("n",n)         
    
    plt.plot(x,y,'-o') 
    
    panels = np.empty((2*n)-2, dtype=object)
    for i in range(0, (2*n)-2):
#        print("yi",y[i],i)
#        print("yi+1",y[i+1],i+1)
        panels[i] = Panel(x[i],x[i+1],y[i],y[i+1])
    N=(2*n)-2
    return N,x,y,panels
    
N,x,y,panels = discretization(20)

print("print")
print(len(x))
print(len(y))
print(N)
"""Creation of left and side of the equation"""


gamma=np.zeros(N+1) #gamma[0] must be equal to -gamma[-1]
A=np.zeros((N+1,N+1))


def u(x,y,a,b):
    return -1/(2*np.pi) *( ( (a*y)/2 * np.log((x-b)**2 +y**2) - (a*x+b)*np.arctan((x-b)/y) ) - ( (a*y)/2 * np.log((x+b)**2 +y**2) - (a*x+b)*np.arctan((x+b)/y) ) )

def v(x,y,a,b):
    return -1/(2*np.pi) *( (a*(-(x-b)+y*np.arctan((x-b)/y)) + (1/2)*(a*x+b)*np.log((x-b)**2+y**2)) - (a*(-(x+b)+y*np.arctan((x+b)/y)) + (1/2)*(a*x+b)*np.log((x+b)**2+y**2)) )

def contribution(x,y,b):
    A=(-1/(2*np.pi))*(y/2)*np.log(((x-b)**2 + y**2)/((x+b)**2 + y**2))
    B=x*( (np.arctan((x-b)/y)) - (np.arctan((x+b)/y)) )*(-1/(2*np.pi))
    C=( (np.arctan((x-b)/y)) - (np.arctan((x+b)/y)) )*(-1/(2*np.pi))
    contribu = np.array([(A-B-C)/(2*b) , (B-A-C)/(2*b)])
    
    F= (-1/(2*np.pi))*(2*b + y*( np.arctan((x-b)/y) - np.arctan((x+b)/y) ) )
    G= (-1/(2*np.pi))*np.log(((x-b)**2 + y**2)/((x+b)**2 + y**2))
    contribv = np.array([(F+(G/2)*(x+1) )/(2*b), ((G/2)*(1-x) - F)/(2*b)])
    
    return np.flip(contribu), np.flip(contribv) #(gamma left, gamma right)


for i in range(0,len(panels)):
    panel_receiver=panels[i]
    for j in range(0,len(panels)):
        panel_sender=panels[j]
        if (i!=j):
            a=np.array([panel_receiver.xc-panel_sender.xc,panel_receiver.yc-panel_sender.yc])
            
            #cos = np.dot(a,b)/ ( np.linalg.norm(a)*np.linalg.norm(b) )
            cos = np.dot(a,panel_sender.normal)/ ( np.linalg.norm(a)*np.linalg.norm(panel_sender.normal))
            angle = np.arccos(cos)
            #print(angle-np.pi/2)
            if(angle >= np.pi/2):
                angle=np.pi-angle
            xp = np.sin(angle) *np.linalg.norm(a)* np.sign(np.dot(a,panel_sender.longitudinal))
            yp = np.cos(angle) *np.linalg.norm(a)* np.sign(np.dot(a,panel_sender.normal))
#            if(i==8):
#                print(np.rad2deg(angle))
#                print("panel (i,j)","(",i,",",j,")","xp",xp,"yp",yp,"dot x",np.sign(np.dot(a,panel_sender.longitudinal)),"dot y",np.sign(np.dot(a,panel_sender.normal)))
#               taila=[panel_sender.xc,panel_sender.yc]
#                plt.quiver(*taila,a[0],a[1],scale=1,scale_units='xy',angles = 'xy',color=['g', 'r', 'k'])
#                tailb=[panel_sender.xi,panel_sender.yi]
#                plt.quiver(*tailb,panel_sender.longitudinal[0],panel_sender.longitudinal[1],scale=1,scale_units='xy',angles = 'xy',color=['g', 'r', 'k'])
#                tailc=[panel_sender.xc,panel_sender.yc]
#                plt.quiver(*tailc,panel_sender.normal[0],panel_sender.normal[1],scale=1,scale_units='xy',angles = 'xy',color=['g', 'r', 'k'])
#                plt.xlim([-0.5,1.5])
#                plt.ylim([-1,1])               
            
            
            
            
#            a=np.array([panel_receiver.xc-panel_sender.xc,panel_receiver.yc-panel_sender.yc])
#            print("a",a,"j",j)
#            b=np.array([panel_sender.xf-panel_sender.xi,panel_sender.yf-panel_sender.yi])
#            print("a norm",np.linalg.norm(a),"j:",j)
#            print("panel_receiver.xc",panel_receiver.xc)
#            print("panel_receiver.yc",panel_receiver.yc)
#            print("panel_sender.xc",panel_sender.xc)
#            print("panel_receiver.yc",panel_sender.yc)
#            cos = np.dot(a,b)/ ( np.linalg.norm(a)*np.linalg.norm(b) )
#            sin = np.sqrt(1-cos**2) 
#            xp= np.linalg.norm(a)*cos
#            yp= np.linalg.norm(a)*sin 
            
            contribu, contribv = contribution(xp,yp,panel_sender.length/2)
            contribu = contribu * np.dot(panel_receiver.normal, panel_sender.longitudinal/panel_sender.length)
            contribv = contribv * np.dot(panel_receiver.normal,panel_sender.normal)
               
            #print(contribu[0]+contribv[0])
            #print(contribu[1]+contribv[1])
            A[i,j]+= contribu[0]+contribv[0]
            A[i,j+1]+= contribu[1]+contribv[1] 
            #print("a",a,"b",b)
 
#Equation to have gamma0=-gamma(N+1)       
A[-1,0]=1 
A[-1,-1]= 1



"""Creation of right hand side of the equation"""
Uinf=20*np.array([1,0])
bright=np.zeros(N+1)

for k in range(0,len(panels)):
    panel=panels[k]
    bright[k]=-np.dot(Uinf,panel.normal)


"""Solving the system"""


gamma_solve = np.linalg.solve(A,bright)
#gamma_solve, residuals, rank, s = np.linalg.lstsq(A,bright)


print(gamma_solve/20)

plt.figure()
plt.plot(np.arange(0,len(gamma_solve)),gamma_solve/20,"-o")    





"""Print for debugging"""

####print for debugging###
#plt.figure()
#for i in range(0,np.size(panels)):
#    panel = panels[i]
#    print("Panel",i+1)
#    print("xi",panels[i].xi)
#    print("xf",panels[i].xf)
#    print("yi",panels[i].yi)
#    print("yf",panels[i].yf)
#    print("xc",panels[i].xc)
#    print("yc",panels[i].yc)
#    plt.scatter(panels[i].xc,panels[i].yc)
#    print(panels[i].length)
#    print(np.rad2deg(panels[i].angle))
#    print("nx",panels[i].nx)
#    print("ny",panels[i].ny)
#    taild=[panel.xi,panel.yi]
#    plt.quiver(*taild,panel.longitudinal[0],panel.longitudinal[1],scale=1,scale_units='xy',angles = 'xy',color=['g', 'r', 'k'])
#    tailc=[panel.xc,panel.yc]
#    plt.quiver(*tailc,-panel.longitudinal[1]/panel.length,panel.longitudinal[0]/panel.length,scale=1,scale_units='xy',angles = 'xy',color=['g', 'r', 'k'])
#    plt.xlim([-1,2])
#    plt.ylim([-1,1])



##influence of one panel test on a single panel
#
#xc1=1
#xc2=0
#yc1=1
#yc2=2
#xi=0
#xf=2
#yi=1
#yf=2
#nx=-1
#ny=-1
#nxp=0
#nyp=1
#normalp=np.array([nxp,nyp]) #normal du panel qui influence
#normal=np.array([nx,ny]) #normal du panel sur lequel on calcule l'influence
#a=np.array([xc2-xc1,yc2-yc1])
#b=np.array([xf-xi,yf-yi])
#cos = np.dot(a,b)/ ( np.linalg.norm(a)*np.linalg.norm(b) )
#sin = np.sqrt(1-cos**2) 
#x= np.linalg.norm(a)*cos
#y= np.linalg.norm(a)*sin  
#
#print(x,y)  
#
##compute u(x,y)
#u=10
#v=20
#
#
#contribution = u*np.dot(a,normal) + v*np.dot(normalp,normal) 




        
        
        
        
        
        
    
