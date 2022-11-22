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

def FEALL(_data,_N,_C,_idxdummy):
    _Farray,_P = F(_data[:,1:],_data[:,0],_N,_C)
    # print(param_array)
    # print(args.keytxt)
    # print(Cindex)
    _warmupdata1, _Feq, _idx_F = autoWarmupMSER_(_Farray)
    _warmupdata2, _Peq, _idx = autoWarmupMSER_(_P)
    _warmupdata3, _mueq, _idx = autoWarmupMSER_(_data[:,0])
    
    _nsamples,_bounds_F,_mean_F,_semcc_F,_kappa_F,_unbiasedvar_F,_autocor_F = doStats(np.array([]), _Feq)
    _nsamples,_bounds_mu,_mean_mu,_semcc_mu,_kappa_mu,_unbiasedvar_mu,_autocor_mu = doStats(np.array([]), _mueq)

    # print(mean_F)
    _nsamples,_bounds_P,_mean_P,_semcc_P,_kappa_P,_unbiasedvar_P,_autocor_P = doStats(np.array([]), _Peq)
    _nsamples,_bounds_Sx,_mean_Sx,_semcc_Sx,_kappa_Sx,_unbiasedvar_Sx,_autocor_Sx = doStats(np.array([]), _data[_idx_F:,1])
    _nsamples,_bounds_Sy,_mean_Sy,_semcc_Sy,_kappa_Sy,_unbiasedvar_Sy,_autocor_Sy = doStats(np.array([]), _data[_idx_F:,2])
    _nsamples,_bounds_Sz,_mean_Sz,_semcc_Sz,_kappa_Sz,_unbiasedvar_Sz,_autocor_Sz = doStats(np.array([]), _data[_idx_F:,3])
    return _mean_mu,_semcc_mu,_mean_P,_semcc_P,_mean_F,_semcc_F,_idx_F,\
        _warmupdata1,_warmupdata2,_warmupdata3,\
            _mean_Sx,_semcc_Sx,_mean_Sy,_semcc_Sy,_mean_Sz,_semcc_Sz

def FEMU(_data,_N,_C,_idxdummy):
    _warmupdata3, _mueq, _idx = autoWarmupMSER_(_data[:,0])
    _Farray,_P = F(_data[_idx:,1:],_mueq,_N,_C)

    _nsamples,_bounds_F,_mean_F,_semcc_F,_kappa_F,_unbiasedvar_F,_autocor_F = doStats(np.array([]), _Farray)
    _nsamples,_bounds_mu,_mean_mu,_semcc_mu,_kappa_mu,_unbiasedvar_mu,_autocor_mu = doStats(np.array([]), _mueq)

    # print(mean_F)
    _nsamples,_bounds_P,_mean_P,_semcc_P,_kappa_P,_unbiasedvar_P,_autocor_P = doStats(np.array([]), _P)
    _nsamples,_bounds_Sx,_mean_Sx,_semcc_Sx,_kappa_Sx,_unbiasedvar_Sx,_autocor_Sx = doStats(np.array([]), _data[_idx:,1])
    _nsamples,_bounds_Sy,_mean_Sy,_semcc_Sy,_kappa_Sy,_unbiasedvar_Sy,_autocor_Sy = doStats(np.array([]), _data[_idx:,2])
    _nsamples,_bounds_Sz,_mean_Sz,_semcc_Sz,_kappa_Sz,_unbiasedvar_Sz,_autocor_Sz = doStats(np.array([]), _data[_idx:,3])
    return _mean_mu,_semcc_mu,_mean_P,_semcc_P,_mean_F,_semcc_F,_idx,\
        _warmupdata3,_warmupdata3,_warmupdata3,\
            _mean_Sx,_semcc_Sx,_mean_Sy,_semcc_Sy,_mean_Sz,_semcc_Sz

def FEChoose(_data,_N,_C,_idx):
    _Farray,_P = F(_data[_idx:,1:],data[_idx:,0],_N,_C)

    _nsamples,_bounds_F,_mean_F,_semcc_F,_kappa_F,_unbiasedvar_F,_autocor_F = doStats(np.array([]), _Farray)
    _nsamples,_bounds_mu,_mean_mu,_semcc_mu,_kappa_mu,_unbiasedvar_mu,_autocor_mu = doStats(np.array([]), data[_idx:,0])

    # print(mean_F)
    _nsamples,_bounds_P,_mean_P,_semcc_P,_kappa_P,_unbiasedvar_P,_autocor_P = doStats(np.array([]), _P)
    _nsamples,_bounds_Sx,_mean_Sx,_semcc_Sx,_kappa_Sx,_unbiasedvar_Sx,_autocor_Sx = doStats(np.array([]), _data[_idx:,1])
    _nsamples,_bounds_Sy,_mean_Sy,_semcc_Sy,_kappa_Sy,_unbiasedvar_Sy,_autocor_Sy = doStats(np.array([]), _data[_idx:,2])
    _nsamples,_bounds_Sz,_mean_Sz,_semcc_Sz,_kappa_Sz,_unbiasedvar_Sz,_autocor_Sz = doStats(np.array([]), _data[_idx:,3])
    return _mean_mu,_semcc_mu,_mean_P,_semcc_P,_mean_F,_semcc_F,_idx,\
        _data[:_idx,0],_data[:_idx,0],_data[:_idx,0],\
            _mean_Sx,_semcc_Sx,_mean_Sy,_semcc_Sy,_mean_Sz,_semcc_Sz

# don't include the operator track in operators
# os.chdir("/media/tquah/TIMMY/Projects/PodProjects/Projects/BenchmarkComparePFTS/ODT_main_workflow/CL_PartialFTS_ODT/CL_PartialFTS")

FUNCTION_MAP = {'FEALL' : FEALL,
                'FEMU' : FEMU,
                "FEChoose":FEChoose}



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
    parser.add_argument('-efn','--exporttrajname',action = 'store',default = "F.dat", help = 'export free energy trajectory', type = str)
    parser.add_argument('-opf','--OPFile',action = 'store',default = "", help = 'Order Parameter File', type = str)

    parser.add_argument('-at','--averagingtype',default = FEALL, help = "Different Averaging Methods",choices=FUNCTION_MAP.keys())
    parser.add_argument('-w','--warmup',default = 0, help = "Warmup",type = int)
    parser.add_argument('--ignorestatus', action='store', default=[], nargs='+',help='status to ignore')




    parser.add_argument('-ep','--export_all',action = 'store_true', help = 'ExportAll')
    parser.add_argument('-eft','--exporttraj',action = 'store_true', help = 'export free energy trajectory')
    parser.add_argument('-v','--verbose',action = 'store_true', help = 'Output More Details')
    parser.add_argument('-T','--sta',action = 'store_true', help = 'Stress Tensor Analysis')
    parser.add_argument('-O','--OP',action = 'store_true',help = 'Order Parameter Analysis')
    # parser.add_argument('-t', '--plottype', action='store', default='simplecolors',help='type of plot to generate')
    # parser.add_argument('-a','--autowarmup',default=False,dest='autowarmup',action='store_true',help='Use MSER-5 method to automate warmup detection')


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
    param_array_used = []
    for i in range(len(filelist)): 
        bpath = PrunePath(filelist[i],-1)

        try:
            statuspath = os.path.join(bpath,'STATUS')
            status = np.loadtxt(statuspath)
            if status in args.ignorestatus:
                print(f'{bpath} Status {status} IGNORE')
                continue
        except:
            print('No Status-No Stuff')
            continue
        
        
        
        
        if Clock:
            C = param_array[i,Cindex]*np.sqrt(args.N)
        #get list of stuff
        operatorlist = ['ChemicalPotential.Real','StressXX_VolScl.Real','StressYY_VolScl.Real','StressZZ_VolScl.Real']
        try: 
            data,warm =  GetOperators(filelist[i],operatorlist)
            func = FUNCTION_MAP[args.averagingtype]
    
            
            mean_mu,semcc_mu,mean_P,semcc_P,mean_F,semcc_F, idx_F,\
                warmupdata1,warmupdata2,warmupdata3,\
                    mean_Sx,semcc_Sx,mean_Sy,semcc_Sy,mean_Sz,semcc_Sz= func(data,args.N,C,args.warmup)
        except:
            print("Some Issue")
            continue
        datalist = []
        
        
        datalist+=[mean_mu,semcc_mu,mean_P,semcc_P,mean_F,semcc_F]




        if args.export_all:
            if args.exporttraj:
                outputlist = []
                outputlist.append(data[:,0])
                Farray,P = F(data[:,1:],data[:,0],args.N,C)

                outputlist.append(P)
                outputlist.append(Farray)
                    
                array = np.vstack(outputlist).transpose()
                Fpath = os.path.join(bpath,args.exporttrajname)
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

            


        
        

        if args.sta:
            def StressTensor(Sxx,Syy,Szz,Sxy,Syz,Sxz):
                A = np.zeros((3,3))
                A[0,0] = Sxx
                A[1,1] = Syy
                A[2,2] = Szz
                A[0,1] = Sxy
                A[1,0] = Sxy
                A[1,2] = Syz
                A[2,1] = Syz
                A[0,2] = Sxz
                A[2,0] = Sxz
                return A
    
                ChemicalPotential
                
                
            operatorlist = ['StressXX.Real','StressYY.Real','StressZZ.Real','StressXY.Real','StressYZ.Real','StressXZ.Real']
    
            data,warm =  GetOperators(filelist[i],operatorlist)
            warmupdata1, Sxx, idx = autoWarmupMSER_(data[:,0])
            warmupdata2, Syy, idx = autoWarmupMSER_(data[:,1])
            warmupdata3, Szz, idx = autoWarmupMSER_(data[:,2])
            warmupdata4, Sxy, idx = autoWarmupMSER_(data[:,3])
            warmupdata5, Syz, idx = autoWarmupMSER_(data[:,4])
            warmupdata6, Sxz, idx = autoWarmupMSER_(data[:,5])
            
            nsamples,bounds_Sxx,mean_Sxx,semcc_Sxx,kappa_Sxx,unbiasedvar_Sxx,autocor_Sxx = doStats(np.array([]), Sxx)
            nsamples,bounds_Syy,mean_Syy,semcc_Syy,kappa_Syy,unbiasedvar_Syy,autocor_Syy = doStats(np.array([]), Syy)
            nsamples,bounds_Szz,mean_Szz,semcc_Szz,kappa_Szz,unbiasedvar_Szz,autocor_Szz = doStats(np.array([]), Szz)
            nsamples,bounds_Sxy,mean_Sxy,semcc_Sxy,kappa_Sxy,unbiasedvar_Sxy,autocor_Sxy = doStats(np.array([]), Sxy)
            nsamples,bounds_Syz,mean_Syz,semcc_Syz,kappa_Syz,unbiasedvar_Syz,autocor_Syz = doStats(np.array([]), Syz)
            nsamples,bounds_Sxz,mean_Sxz,semcc_Sxz,kappa_Sxz,unbiasedvar_Sxz,autocor_Sxz = doStats(np.array([]), Sxz)

            Tensor = StressTensor(mean_Sxx,mean_Syy,mean_Szz,mean_Sxy,mean_Syz,mean_Sxz)
            Eig =np.diag(np.real(np.linalg.eig(Tensor)[0]))
            MSE = np.sum(np.square(Tensor-Eig))

            datalist+= [mean_Sxx,semcc_Sxx,\
                         mean_Syy,semcc_Syy,\
                         mean_Szz,semcc_Szz,\
                         mean_Sxy,semcc_Sxy,\
                         mean_Syz,semcc_Syz,\
                         mean_Sxz,semcc_Sxz,MSE]
            if args.verbose:
                print("Averages and Errors")
                print(f'Sx: {mean_Sxx}+-{semcc_Sxx} Sy: {mean_Syy}+-{semcc_Syy} Sz: {mean_Szz}+-{semcc_Szz}')
                print(f'Sxy: {mean_Sxy}+-{semcc_Sxy} Syz: {mean_Syz}+-{semcc_Syz} Sxz: {mean_Sxz}+-{semcc_Sxz}')



        if args.OP:
            fullOPpath = os.path.join(bpath,args.OPFile)
            data_load = np.loadtxt(fullOPpath)
            nsamples,bounds_OP,mean_OP,semcc_OP,kappa_OP,unbiasedvar_OP,autocor_OP = doStats(np.array([]), data_load[idx_F:,1])
            datalist+=[mean_OP,semcc_OP]
            if args.verbose:
                print("Averages and Errors")
                print(f'OP: {mean_OP}+-{semcc_OP}')


        data_array.append(np.array(datalist))
        param_array_used.append(np.array(param_array[i,:]))
    if args.export_all:
        args.keytxt.append('ChemicalPotential')
        args.keytxt.append('ChemicalPotential.Error')
        args.keytxt.append('P')
        args.keytxt.append('P.Error')
        args.keytxt.append('F')
        args.keytxt.append('F.Error')

        if args.sta: 
            operatorlist = ['StressXX.Real','StressYY.Real','StressZZ.Real','StressXY.Real','StressYZ.Real','StressXZ.Real']
            for operator in operatorlist:
                args.keytxt.append(operator.split('.')[0])
                args.keytxt.append(operator.split('.')[0]+'.Error')
            args.keytxt.append('MSE')
        
        if args.OP:
            args.keytxt.append('OP')
            args.keytxt.append('OP.Error')
        print('------------------------------------------------------------------------------')

        stackarray = np.hstack((np.vstack(param_array_used),np.vstack(data_array)))
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

