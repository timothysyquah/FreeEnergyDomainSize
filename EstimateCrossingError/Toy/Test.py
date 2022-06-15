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



data = np.loadtxt('data.dat')

plt.figure()
plt.plot(data[:,0],data[:,1],'ok')
plt.plot(data[:,0],data[:,2],'^g')


param1, param2, error1, error2, xcross, crosserror = EstimateIntersection(data[:,0],data[:,1],data[:,0],data[:,2])
x = np.linspace(np.min(data[:,0]),np.max(data[:,0]))
plt.plot(x,param1[0]+x*param1[1],'-.',color = 'k')
plt.plot(x,param2[0]+x*param2[1],'--',color = 'g')
plt.errorbar(xcross,param2[0]+xcross*param2[1],xerr = crosserror,marker = 's',color = 'r',capsize = 5.0)






# plt.axvline(x0)