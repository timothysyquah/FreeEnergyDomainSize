#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 20:15:48 2022

@author: tquah
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from scipy.interpolate import interp1d
from scipy.optimize import fsolve
# Extend ScalarFormatter
class MyScalarFormatter(ScalarFormatter):
    # Override '_set_format' with your own
    def _set_format(self):
        self.format = '%.2f'  # Show 2 decimals


def sortarray(array,col):
    sort= np.argsort(array[:,col])
    array = array[sort,:]
    return array

partialfts = np.loadtxt('PartialFTS_F.dat')
cl = np.loadtxt('CL_F.dat')
# cl = np.loadtxt('PartialFTS_t1_F.dat')

C = np.unique(list(partialfts[:,0])+list(cl[:,0]))
# C = np.array([50.])
plt.close('all')
N = 100
for i in range(0,len(C)):
    
    fig, axs = plt.subplots(2,2,figsize = (10,10))
    fig.suptitle(f'C = {C[i]*np.sqrt(N)}')

    locp = np.where(partialfts[:,0]==C[i])[0]
    locc = np.where(cl[:,0]==C[i])[0]
    
    temp_pfts = partialfts[locp,:]
    temp_pfts = sortarray(temp_pfts,1)

    temp_cl = cl[locc,:]
    temp_cl = sortarray(temp_cl,1)
    
    #plot free eenergy
    axs[0,0].set_title('Free Energy')

    axs[0,0].errorbar(temp_pfts[:,1]*N,temp_pfts[:,2],yerr =temp_pfts[:,3],marker = 'o',color = 'r',label = 'PartialFTS')
    axs[0,0].axhline(y=0.0, color='b', linestyle='--')
    axs[0,0].errorbar(temp_cl[:,1]*N,temp_cl[:,2],yerr = temp_cl[:,3],marker = '^',color = 'k',label = 'CL')
    
    
    
    def stupidodt(x,y):
        dx = x[1:]-x[0:-1]
        
        dy= y[1:]-y[0:-1]
        m = dy/dx
        b = y[0:-1]-m*x[0:-1]
        x0 = -b/m
        logic = np.ma.greater(x0,x[0:-1]).data
        logic *= np.ma.greater(x[1:],x0).data
        loc = np.where(logic)[0]
        return x0[loc]
        
        
    pfts_odt = stupidodt(temp_pfts[:,1]*N,temp_pfts[:,2])
    
    cl_odt = stupidodt(temp_cl[:,1]*N,temp_cl[:,2])
    
    for j in range(0,len(pfts_odt)):
        axs[0,0].axvline(x=pfts_odt[j], color='y', linestyle='--')
        axs[1,0].axvline(x=pfts_odt[j], color='y', linestyle='--')
        axs[1,1].axvline(x=pfts_odt[j], color='y', linestyle='--')
        axs[0,1].axvline(x=pfts_odt[j], color='y', linestyle='--')
    for j in range(0,len(cl_odt)):

        axs[0,0].axvline(x=cl_odt[j], color='m', linestyle='--')
        axs[1,0].axvline(x=cl_odt[j], color='m', linestyle='--')
        axs[1,1].axvline(x=cl_odt[j], color='m', linestyle='--')
        axs[0,1].axvline(x=cl_odt[j], color='m', linestyle='--')

    
    
    axs[0,0].legend()
    custom_formatter = MyScalarFormatter(useMathText=True)
    axs[0,0].yaxis.set_major_formatter(custom_formatter)

    fig.text(0.5, 0.04, r'$\chi N$', ha='center')
    axs[0,1].set_title('SCFT "Stress" Operator')
    axs[0,1].errorbar(temp_pfts[:,1]*N,temp_pfts[:,4]-temp_pfts[0,4],yerr =temp_pfts[:,5],marker = 's',color = 'r',label = 'PFTS-DIS')
    axs[0,1].errorbar(temp_pfts[:,1]*N,temp_pfts[:,6]-temp_pfts[0,6],yerr =temp_pfts[:,7],marker = '^',color = 'r',label = 'PFTS-LAM')
    axs[0,1].axhline(y=0.0, color='b', linestyle='--')
    axs[0,1].errorbar(temp_cl[:,1]*N,temp_cl[:,4]-temp_cl[0,4],yerr =temp_cl[:,5],marker = 's',color = 'k',label = 'CL-DIS')
    axs[0,1].errorbar(temp_cl[:,1]*N,temp_cl[:,6]-temp_cl[0,6],yerr =temp_cl[:,7],marker = '^',color = 'k',label = 'CL-LAM')
    axs[0,1].legend()
    axs[0,1].yaxis.set_major_formatter(custom_formatter)

    #structure factor
    axs[1,0].set_title('Order Parameter')
    axs[1,0].errorbar(temp_pfts[:,1]*N,temp_pfts[:,10],yerr =temp_pfts[:,11],marker = 's',color = 'r',label = 'PFTS-LAM',linestyle = ':')
    axs[1,0].errorbar(temp_pfts[:,1]*N,temp_pfts[:,12],yerr =temp_pfts[:,13],marker = '^',color = 'r',label = 'PFTS-DIS',linestyle = '--')
    # axs[2].axhline(y=0.0, color='b', linestyle='--')
    axs[1,0].errorbar(temp_cl[:,1]*N,temp_cl[:,10],yerr =temp_cl[:,11],marker = 's',color = 'k',label = 'CL-LAM',linestyle = ':')
    axs[1,0].errorbar(temp_cl[:,1]*N,temp_cl[:,12],yerr =temp_cl[:,13],marker = '^',color = 'k',label = 'CL-DIS',linestyle = '--')
    axs[1,0].legend()
    

    
    
    #structure factor
    axs[1,1].set_title('$\Delta$ Order Parameter')
    axs[1,1].errorbar(temp_pfts[:,1]*N,np.abs(temp_pfts[:,8]),yerr =temp_pfts[:,9],marker = 's',color = 'r',label = 'PFTS',linestyle = ':')
    # axs[2].errorbar(temp_pfts[:,1]*N,temp_pfts[:,12],yerr =temp_pfts[:,13],marker = '^',color = 'r',label = 'PFTS-LAM',linestyle = '--')
    # axs[2].axhline(y=0.0, color='b', linestyle='--')
    axs[1,1].errorbar(temp_cl[:,1]*N,np.abs(temp_cl[:,8]),yerr =temp_cl[:,9],marker = 's',color = 'k',label = 'CL',linestyle = '--')
    # axs[2].errorbar(temp_cl[:,1]*N,temp_cl[:,12],yerr =temp_cl[:,13],marker = '^',color = 'k',label = 'CL-LAM',linestyle = '--')
    axs[1,1].legend()

