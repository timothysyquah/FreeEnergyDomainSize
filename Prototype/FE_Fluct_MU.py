#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  1 14:07:29 2022

@author: tquah
"""
import glob
import os 
from Functions import *
import sys
sys.path.append('/home/tquah/PolyFTS_ALL/PolyFTS/tools') #gives us access to Kris' tools
from stats import *
import argparse
import numpy as np

longdes = \
    "FECOMBINE-Combines Stress Components/Mu then uses stats.py to determine cutoff in free energy\
      FELIMIT-Uses mu warmup as cutoff"

"""
Here we will outline the code before writing it
1) use argparse to specify
* Directory Structure
* Specify Keywords 
* Tolerance (free energy)
* Activate Filters-remove bad CL....
* Keywords for averaging
* N-and definitions of C
2) Load data and characterize it based on directory structure and keywords
* may need to do cleaning step here
3) Average the operators 
* Gather volscl
* Gather mu
* Find the series that is shorter
* Compute vectorized Free Energy
* Push through the tools in stats
** Get back averaged values with standard error of the mean
"""
def F(_diag,_mu,_N,_C):
    _P = np.sum(_diag,axis = 1)/3/_C*np.power(_N,3/2)
    _P = (_diag[:,0]+_diag[:,1]+_diag[:,2])/3/_C*np.power(_N,3/2)
    return _P+_mu,_P

# don't include the operator track in operators
# os.chdir("/media/tquah/TIMMY/Projects/PodProjects/Projects/BenchmarkComparePFTS/ODT_main_workflow/CL_PartialFTS_ODT/CL_PartialFTS")
if __name__ == "__main__":
    
    IDIR = os.getcwd()
    parser = argparse.ArgumentParser(description='Tool to Calculate Free Energy')
    #do not mix simulation methods needs "" or else does something strange
    parser.add_argument('-d', '--dirs', action='store', default="",help='list of directories that contain each phase point',type=str)
    parser.add_argument('-f','--file',action = 'store',default = '', help = 'File to read with averages and error', type = str)
    parser.add_argument('-kl','--keyloc',action = 'store',nargs='+',default = [], help = 'keywords location', type = int)
    parser.add_argument('-kt','--keytxt',action = 'store',nargs='+',default = [''], help = 'keywords label', type = str)
    parser.add_argument('-n','--N',action = 'store',default = 1, help = 'Degree of Polymerization', type = float)
    parser.add_argument('-c','--Cref',action = 'store',default = None, help = 'Cref', type = float)
    parser.add_argument('-ed','--exportpath',action = 'store',default = '', help = 'path to export', type = str)
    parser.add_argument('-e','--exportname',action = 'store',default = '', help = 'name to export', type = str)
    parser.add_argument('-eft','--exporttraj',action = 'store',default = False, help = 'export free energy trajectory', type = bool)
    parser.add_argument('-v','--verbose',action = 'store',default = True, help = 'Output More Details', type = bool)
    # parser.add_argument('-at','--averagingtype',action = 'store',default = 'FECOMBINE', help = longdes, type = bool)

    args = parser.parse_args()
    
    # args.dirs = "C*/CL/chi*"
    # args.file = "model1_operators.dat"
    # args.keyloc = [0,2]
    # args.keytxt = ['C','chi']
    # args.exportpath = './'
    # args.exportname = 'operator_averages.dat'
    path = os.path.join(args.dirs,args.file)
    filelist = glob.glob(path)
    epath = os.path.join(args.exportpath,args.exportname)
    print(filelist)
    paramlist = GetParameters(filelist,args.keyloc)    
    # we only allow a single parameter per line
    param_array = np.vstack(paramlist)
    # print(param_array)
    # We know the chemical potential is the hardest to average
    operatorlist = ['ChemicalPotential.Real','StressXX_VolScl.Real','StressYY_VolScl.Real','StressZZ_VolScl.Real']
    Clock = True
    if 'C' not in args.keytxt  and args.Cref==None:
        print("Error C is not specified")
        sys.exit(1)
    elif 'C' not in args.keytxt:
        C = args.Cref*np.sqrt(args.N)
        Clock = False
    else:
        Cindex = args.keytxt.index('C')
    data_array = []
    for i in range(len(filelist)):
        if Clock:
            C = param_array[i,Cindex]*np.sqrt(args.N)
        #get list of stuff

        data,warm =  GetOperators(filelist[i],operatorlist)
        
        
        
        Farray,P = F(data[:,1:],data[:,0],args.N,C)
        # print(param_array)
        # print(args.keytxt)
        # print(Cindex)
        warmupdata3, mueq, idx = autoWarmupMSER_(data[:,0])

        warmupdata1, Feq, idx = autoWarmupMSER_(Farray)
        warmupdata2, Peq, idx = autoWarmupMSER_(P)
        
        nsamples,bounds_F,mean_F,semcc_F,kappa_F,unbiasedvar_F,autocor_F = doStats(np.array([]), Farray[idx:])
        # print(mean_F)
        nsamples,bounds_P,mean_P,semcc_P,kappa_P,unbiasedvar_P,autocor_P = doStats(np.array([]), P[idx:])
        nsamples,bounds_Sx,mean_Sx,semcc_Sx,kappa_Sx,unbiasedvar_Sx,autocor_Sx = doStats(np.array([]), data[:,1])
        nsamples,bounds_Sy,mean_Sy,semcc_Sy,kappa_Sy,unbiasedvar_Sy,autocor_Sy = doStats(np.array([]), data[:,2])
        nsamples,bounds_Sz,mean_Sz,semcc_Sz,kappa_Sz,unbiasedvar_Sz,autocor_Sz = doStats(np.array([]), data[:,3])

        
        
        
        nsamples,bounds_mu,mean_mu,semcc_mu,kappa_mu,unbiasedvar_mu,autocor_mu = doStats(np.array([]), mueq)
        data_array.append(np.array([mean_mu,semcc_mu,mean_P,semcc_P,mean_F,semcc_F]))
        
        if args.exporttraj:
            outputlist = []
            outputlist.append(data[:,0])
            outputlist.append(P)
            outputlist.append(Farray)
            array = np.vstack(outputlist).transpose()
            bpath = PrunePath(filelist[i],-1)
            Fpath = os.path.join(bpath,'F.dat')
            print(Fpath)
            np.savetxt(Fpath,array,header = 'mu P F')
        if args.verbose:
            print('------------------------------------------------------------------------------')
            print(f"C = {C} N = {args.N}")
            print("Warmup")
            print(f'F: {len(warmupdata1)} P: {len(warmupdata2)} MU: {len(warmupdata3)}')
            print("Averages and Errors")
            print(f'F: {mean_F}+-{semcc_F} P: {mean_P}+-{semcc_P} Mu: {mean_mu}+-{semcc_mu}')
            print("Averages and Errors")
            print(f'Sx: {mean_Sx}+-{semcc_Sx} Sy: {mean_Sy}+-{semcc_Sy} Sz: {mean_Sz}+-{semcc_Sz}')

            print('------------------------------------------------------------------------------')

    stackarray = np.hstack((param_array,np.vstack(data_array)))
    args.keytxt.append('ChemicalPotential')
    args.keytxt.append('ChemicalPotential.Error')
    args.keytxt.append('P')
    args.keytxt.append('P.Error')
    args.keytxt.append('F')
    args.keytxt.append('F.Error')
    header = Make_Header(args.keytxt)
    np.savetxt(epath,stackarray,header = header)
    
"""
**will be in another script
4) Estimate ODT
* Using Mean and Standard Error of mean
Simplest Algo is to use Mean/SEM(+-)
These three lines capture the most probable transition
"""



"""
5) Output Results
* Write Free Energy File given parameters output Chemical Potential, Pressure, Free Energy
* Write ODT file compiles all things scanned
"""

