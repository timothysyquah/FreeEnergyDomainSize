#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  1 14:07:29 2022

@author: tquah
"""
import os
import re
import glob
import numpy as np
import sys
sys.path.append('/home/tquah/PolyFTS_ALL/PolyFTS/tools') #gives us access to Kris' tools

from stats import *
from scipy import ndimage ,signal,stats
from Functions_Analysis import *
from scipy.stats import mode
import matplotlib.pyplot as plt

filelist = glob.glob('C*/**/chi*/')



uniqueC  = []
uniquechi = []
datadict = dict()
datadict['PartialFTS'] = []
datadict['PartialFTS_t1'] = []

datadict['CL'] = []
figure = True
OP = True
Stress = True

LAM_Files = ['LAM.dat','model2_operators.dat','model2_OrientationPersistenceOP.dat']
DIS_Files = ['DIS.dat','model1_operators.dat','model1_OrientationPersistenceOP.dat']




for file in filelist:
    print(file)
    
    #handle status
    try:
        path = os.path.join(file,'STATUS')
        if np.loadtxt(path) in [1,3]:
            continue
    except:
        continue

    tmplist = []
    
    ele = file.split('/')
    Ctemp = float(extract_value(ele[0])[0])
    chitemp = float(extract_value(ele[2])[0])
    tmplist+=[Ctemp,chitemp]
    try:
        tmplist+=FDiff(file,LAM_Files[0],DIS_Files[0],2,split=False,verbose = False) #dire,file1,file2,col,verbose = False
    except:
        print('F Fail')
        continue
    # print(tmplist)
    try:
        tmplist+=SDiff(file,LAM_Files[1],DIS_Files[1])
    except:
        print('Stress Fail')
        continue
    # print(tmplist)

    try:
        tmplist+=FDiff(file,LAM_Files[2],DIS_Files[2],1,split = True,verbose = False) #dire,file1,file2,col,verbose = False
    except:
        print('OP Fail')
        continue
    # print(tmplist)

    datadict[ele[1]].append(np.array(tmplist))


    
    if np.sum(np.isnan(np.array(tmplist)))>0:
        break

    
    
    
for key in list(datadict):
    datadict[key] = np.vstack(datadict[key])
    np.savetxt(f'{key}_F.dat',datadict[key])