#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 10:43:34 2022

@author: tquah
"""
import numpy as np
import matplotlib.pyplot as plt
import os
plt.close('all')
path = '/home/tquah/DataStore'
file1 = 'urgCL_DISPhase.dat'
file2 = 'eCL_DISPhase.dat'
def return_sorted_data(data,loc):

    sort = np.argsort(data[loc,1])
    tdatax = data[loc,1]
    tdatay = data[loc,6]
    tdatayerr = data[loc,7]
    tdatax = tdatax[sort]
    tdatay = tdatay[sort]
    tdatayerr = tdatayerr[sort]
    return tdatax,tdatay,tdatayerr



data1= np.loadtxt(os.path.join(path,file1))
data2 = np.loadtxt(os.path.join(path,file2))
N = 100
Cunique = np.unique(np.array(list(data1[:,0])+list(data2[:,0])))

for i in range(0,len(Cunique)):
    loc1 = np.where(data1[:,0] == Cunique[i])[0]
    loc2 = np.where(data2[:,0] == Cunique[i])[0]
    
    plt.figure()
        
        
    plt.title(f'C = {Cunique[i]*np.sqrt(N)}')
    tdatax,tdatay,tdatayerr =  return_sorted_data(data1,loc1)
    plt.errorbar(tdatax,tdatay,yerr = tdatayerr,marker = '^',color = 'b')
    tdatax,tdatay,tdatayerr =  return_sorted_data(data2,loc2)
    plt.errorbar(tdatax,tdatay,yerr = tdatayerr,marker = 's',color = 'r')