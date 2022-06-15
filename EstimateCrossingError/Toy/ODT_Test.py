#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 10 14:13:34 2022

@author: tquah
"""
import numpy as np
import matplotlib.pyplot as plt
from ErrorPropagationFunctions import *
plt.close('all')
#define other ones as needed
xmin = 1
xmax = 100
npts = 5



x = np.linspace(1,100,5)

 
y1 = -np.exp(-0.1*x) #this is lam
y2 = np.zeros_like(x) # imaggine this is dis

plt.scatter(x,y1)
plt.scatter(x,y2)