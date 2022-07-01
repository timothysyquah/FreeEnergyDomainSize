#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 17 12:08:41 2022

@author: tquah
"""
import matplotlib.pyplot as plt
import os 
import numpy as np
data = np.loadtxt('/home/tquah/DataStore/dtPOCL.dat')
plt.close('all')
plt.figure(figsize = (6,6))
plt.errorbar(data[:,0],data[:,1],yerr = data[:,2],label = r'$\beta \mu $',marker = '^')
plt.errorbar(data[:,0],data[:,3],yerr = data[:,4],label = r'$\beta P/c$',marker = '+')
plt.errorbar(data[:,0],data[:,5],yerr = data[:,6],label = r'$\beta F/n$',marker = 's')
plt.xscale('log')
plt.legend()
plt.xlabel(r'$dt$')
plt.ylabel(r'$\left< x \right>$')
plt.savefig('/home/tquah/Figures/dtFPo.pdf',dpi = 300)
plt.figure(figsize = (6,6))
# plt.errorbar(data[:,0],data[:,1],yerr = data[:,2],label = r'$\beta \mu $',marker = '^')
# plt.errorbar(data[:,0],data[:,3],yerr = data[:,4],label = r'$\beta P/c$',marker = '+')
plt.errorbar(data[:,0],data[:,5],yerr = data[:,6],label = r'$\beta F/n$',marker = 's')
plt.xscale('log')
plt.ylim(3.20,3.24)
plt.legend()
plt.xlabel(r'$dt$')
plt.ylabel(r'$\left< x \right>$')
plt.savefig('/home/tquah/Figures/dtzoomFPo.pdf',dpi = 300)