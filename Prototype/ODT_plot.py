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
period = 2
os.chdir(f'/home/tquah/DataStore/{period}Period')
if __name__ == "__main__":
    IDIR = os.getcwd()
    parser = argparse.ArgumentParser(description='Calculate ODT')
    #only process one method at a time!
    
    
    logi = 0
    parser.add_argument('-i','--input',action = 'store',nargs='+',default = [''], help = 'input files', type = str)
    parser.add_argument('-o','--output',action = 'store',default = '', help = 'output file', type = str)
    parser.add_argument('-c','--column',action = 'store',default = [],nargs='+', help = 'Interprets this as C, chi, Free Energy and Free Energy Error', type = int)

    parser.add_argument('-t','--tol',action = 'store',default = 1e-6, help = 'tolerance', type = float)
    # parser.add_argument('-p','--plot',action = 'store',default = True, help = 'Plot?', type = bool)
    args = parser.parse_args()
    data_dictionary = dict()
    stress = False
    OP = True
    Job = 'FTS'
    args.input = ["CL_DISPhase.dat","CL_LAMPhase.dat"]
    # args.input = ["PartialFTS_DISPhase.dat","PartialFTS_LAMPhase.dat"]
    # 
    # Fdat =  np.loadtxt('FSCFT.dat')
    # sort = np.argsort(Fdat[:,0])
    # Fdat = Fdat[sort,:]
    
    markerlist = ['^','s']
    markerlist1 = ['>','D']

    color = ['r','b']
    color1 = ['m','c']  

    args.column = [0,1,6,7]
    N = 100
    sqrtN = np.sqrt(N)
    for file in args.input:
        if "DIS" in file:
            data_dictionary['DIS'] = np.loadtxt(file)
            # data_dictionary['DIS'][:,0] *= 10
            # data_dictionary['DIS'][:,1] *= 100

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
        fig, ax1 = plt.subplots(figsize = (9,6))
        ax2 = ax1.twinx()

        plt.title(f"C={sqrtN*C[i]}")
        alist = []
        for j in range(0,len(header)):
            H = []
            


            loc = np.where(C[i]==data_dictionary[header[j]][:,0])[0]
            
            # for k in range(len(data_dictionary[header[j]][loc,args.column[1]])):
                # loc1 = np.where(data_dictionary[header[j]][loc,args.column[1]][k]==Fdat[:,0])[0]
                # H.append(Fdat[loc1])
            # H = np.vstack(H)
            # H *=0
            norm = data_dictionary['DIS'][loc,args.column[2]]
            norm*=1
            ax1.errorbar(N*data_dictionary[header[j]][loc,args.column[1]],\
                         data_dictionary[header[j]][loc,args.column[2]]-norm,\
                         yerr = data_dictionary[header[j]][loc,args.column[3]],label = header[j],marker = markerlist[j],color = color[j],zorder = 2)
                
            asym = np.sqrt(np.square(data_dictionary[header[j]][loc,8]-data_dictionary[header[j]][loc,10])\
                    +np.square(data_dictionary[header[j]][loc,8]-data_dictionary[header[j]][loc,12])\
                    +np.square(data_dictionary[header[j]][loc,10]-data_dictionary[header[j]][loc,12])\
                    +np.square(data_dictionary[header[j]][loc,14])+np.square(data_dictionary[header[j]][loc,16])+np.square(data_dictionary[header[j]][loc,18]))

            alist.append(asym)
    
                # print(asym/np.max(asym))
                
            if OP:
                ax2.errorbar(N*data_dictionary[header[j]][loc,args.column[1]],\
                             data_dictionary[header[j]][loc,-2],\
                             yerr = data_dictionary[header[j]][loc,-1],marker = markerlist1[j],linestyle = '-.',zorder = 2,color =color1[j])
        
        
        if stress:
            maxn = np.max(np.array(alist))
            for j in range(0,len(header)):
                norm = data_dictionary['DIS'][loc,args.column[2]]
                norm*=1

                for l in range(0,len(alist[j])):
                    ax1.scatter(N*data_dictionary[header[j]][loc,args.column[1]][l],\
                                 data_dictionary[header[j]][loc,args.column[2]][l]-norm[l],s = 500.0\
                                 ,zorder = 1,marker = 'o',color = 'k',alpha = alist[j][l]/maxn)
    
        ax1.axvline(np.array([12.0382]),color = 'k',linestyle = '--')
        
        ax1.set_xlabel(r'$\chi N$')
        ax1.set_ylabel(r'$\beta F/n -\beta F_{DIS}/n$')
        ax1.set_ylabel(r'$\beta F/n -\beta F_{DIS}/n$')
        ax2.set_ylabel(r'$\left< \Psi \right>$')
        
        if stress:
            plt.scatter([],[],s = 500.0\
                          ,zorder = 1,marker = 'o',color = 'k',alpha =1.0,label = f'Isotropic Measure-{maxn:0.0e}')
        plt.errorbar([],[]\
                      ,zorder = 1,marker = 's',color = 'b',label = 'LAM')
        plt.errorbar([],[],\
                      zorder = 1,marker = '^',color = 'r',label = 'DIS')
        plt.plot([],[],\
                      zorder = 1,color = 'k',linestyle = '--', label = 'SCFT-ODT')

        plt.legend(loc='center left', bbox_to_anchor=(1.15, 0.5))
        plt.tight_layout()
        plt.savefig(f'/home/tquah/Figures/C{sqrtN*C[i]}_{Job}_N100_period_{period}.pdf',dpi = 300)

    
    # if args.plot:
    #     plt.figure()
    #     plt.errorbar(data_dictionary['DIS'][1,:],data_dictionary['DIS'][-2,:] )
    #     plt.errorbar(data_dictionary['LAM'][1,:],data_dictionary['LAM'][-2,:] )

