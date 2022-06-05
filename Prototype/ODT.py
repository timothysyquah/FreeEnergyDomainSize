#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  1 19:26:09 2022

@author: tquah
"""
from Functions import *
import argparse
import matplotlib.pyplot as plt
import numpy as np
import os 
plt.close('all')
os.chdir('/home/tquah/DataStore')
if __name__ == "__main__":
    IDIR = os.getcwd()
    parser = argparse.ArgumentParser(description='Calculate ODT')
    #only process one method at a time!
    
    
    
    parser.add_argument('-i','--input',action = 'store',nargs='+',default = [''], help = 'input files', type = str)
    parser.add_argument('-o','--output',action = 'store',default = '', help = 'output file', type = str)
    parser.add_argument('-c','--column',action = 'store',default = [],nargs='+', help = 'Interprets this as C, chi, Free Energy and Free Energy Error', type = int)

    parser.add_argument('-t','--tol',action = 'store',default = 1e-6, help = 'tolerance', type = float)
    parser.add_argument('-p','--plot',action = 'store',default = False, help = 'Plot?', type = bool)
    args = parser.parse_args()
    data_dictionary = dict()
    
    args.input = ["CL_DISPhase.dat","CL_LAMPhase.dat"]
    # args.input = ["PFTS_DISPhase.dat","PFTS_LAMPhase.dat"]
    markerlist = ['^','s']
    args.column = [0,1,-2,-1]

    for file in args.input:
        if "DIS" in file:
            data_dictionary['DIS'] = np.loadtxt(file)
    
        elif "LAM" in file:
            data_dictionary['LAM'] = np.loadtxt(file)
            
        else:
            print('Only Supports LAM and DIS ODT')
    
    header = list(data_dictionary)
    C=[]
    chi = []
    # F = []
    # Ferror = []
    for i in range(0,len(header)):
        C += list(np.unique(data_dictionary[header[i]][:,args.column[0]]))
    C = np.unique(np.array(C))

    for i in range(0,len(C)):
        plt.figure()
        plt.title(f"C={10*C[i]}")
        for j in range(0,len(header)):

            loc = np.where(C[i]==data_dictionary[header[j]][:,0])[0]
            plt.errorbar(100*data_dictionary[header[j]][loc,args.column[1]],\
                         data_dictionary[header[j]][loc,args.column[2]],\
                         yerr = data_dictionary[header[j]][loc,args.column[3]],label = header[j],marker = markerlist[j])
        plt.legend()
        plt.xlabel(r'$\chi N$')
        plt.ylabel(r'$\beta A_{ex}/n$')
        plt.tight_layout()
        plt.savefig(f'/home/tquah/Figures/C{10*C[i]}_PFTS.pdf',dpi = 300)

    
    if args.plot:
        plt.errorbar(data_dictionary['DIS'][1,:],data_dictionary['DIS'][-2,:] )
        plt.errorbar(data_dictionary['LAM'][1,:],data_dictionary['LAM'][-2,:] )
