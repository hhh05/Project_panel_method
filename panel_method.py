# -*- coding: utf-8 -*-
"""
Created on Thu Apr 15 10:27:28 2021

@author: henri and julien
"""
""" The goal of this code is to compute the velocity around an airfoil with the linear vortex sheet
panel method"""

"""Comment: 200 or 300 hundreds panels necessary"""

import numpy as np
from matplotlib import pyplot as plt


def y_e(x,c):
    return 0.278536*(np.sqrt(x/c)) - 0.148567*(x/c) + 0.006397*(x/c)**2 - 0.22098*(x/c)**3 + 0.081084*(x/c)**4
def y_i(x,c):
    return - 0.190361*(np.sqrt(x/c)) + 0.161628*(x/c) - 0.341176*(x/c)**2 + 0.897631*(x/c)**3 - 0.531252*(x/c)**4

c = 1 
x = np.linspace(0,1,200)
y = y_e(x,c)
yy = y_i(x,c)
plt.figure(figsize=(15,4))
plt.plot(x,yy)
plt.plot(x,y)
plt.show()