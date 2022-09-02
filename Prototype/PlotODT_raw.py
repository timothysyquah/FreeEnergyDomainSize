#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 26 15:00:42 2022

@author: tquah
"""
import numpy as np
import matplotlib.pyplot as plt

plt.close('all')
odt = np.loadtxt('ODT.dat')
plt.figure(figsize = (6,6))
plt.scatter(odt[2:,0],odt[2:,1],marker = 's',color = 'r',label = 'CL-FFT')
plt.scatter(odt[0:2,0],odt[0:2,1],marker = '*',color = 'r',label = 'CL-FFT (2 Period)')
plt.scatter(odt[:,0],odt[:,2],marker = '^',color = 'b',label = 'L-PFT')
plt.axhline(np.array([12.0382]),color = 'k',linestyle = '--',label = 'MFT')

plt.legend()
plt.xscale('log')
plt.xlabel('$C$')
plt.ylabel(r'$\chi N_{\text{ODT}}$')
plt.tight_layout()
plt.savefig('/home/tquah/Figures/ODT.pdf',dpi = 300)