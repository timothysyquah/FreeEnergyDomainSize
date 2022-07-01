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
    parser.add_argument('-p','--plot',action = 'store',default = True, help = 'Plot?', type = bool)
    args = parser.parse_args()
    data_dictionary = dict()
    
    args.input = ["CL_LAMPhase.dat","CL_DISPhase.dat","sCL_DISPhase.dat"]
    # args.input = ["PFTS_DISPhase.dat","PFTS_LAMPhase.dat"]
    
    Fdat =  np.loadtxt('FSCFT.dat')
    sort = np.argsort(Fdat[:,0])
    Fdat = Fdat[sort,:]
    
    markerlist = ['^','s','p']
    args.column = [0,1,-2,-1]

    for file in args.input:
        print(file)
        if "sCL" in file:
            data_dictionary['sCL'] = np.loadtxt(file)
    
        elif "CL" in file:
            if "sCL" not in file and "LAM" not in file:
                
                data_dictionary['CL'] = np.loadtxt(file)
            if "LAM" in file:
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
    name = ['LAM',r'$L_i = 60 \: l$',r'$L_i = 30 \: l$']
    for i in range(0,len(C)):
        plt.figure()
        plt.title(f"C={10*C[i]}")
        for j in range(0,len(header)):
            H = []
            


            loc = np.where(C[i]==data_dictionary[header[j]][:,0])[0]
            
            for k in range(len(data_dictionary[header[j]][loc,args.column[1]])):
                loc1 = np.where(data_dictionary[header[j]][loc,args.column[1]][k]==Fdat[:,0])[0]
                H.append(Fdat[loc1])
            H = np.vstack(H)
            H[:,3] = 0
            plt.errorbar(100*data_dictionary[header[j]][loc,args.column[1]],\
                         data_dictionary[header[j]][loc,args.column[2]]-H[:,3]*100,\
                         yerr = data_dictionary[header[j]][loc,args.column[3]],label = name[j],marker = markerlist[j])

                
        plt.legend()
        plt.xlabel(r'$\chi N$')
        plt.ylabel(r'$(\beta A_{ex})/n$')
        plt.tight_layout()
        plt.savefig(f'/home/tquah/Figures/C{10*C[i]}_FTS.pdf',dpi = 300)

    
    # if args.plot:
    #     plt.figure()
    #     plt.errorbar(data_dictionary['DIS'][1,:],data_dictionary['DIS'][-2,:] )
    #     plt.errorbar(data_dictionary['LAM'][1,:],data_dictionary['LAM'][-2,:] )

