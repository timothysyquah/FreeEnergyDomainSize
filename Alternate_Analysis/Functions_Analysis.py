#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 31 11:13:53 2022

@author: tquah
"""

from scipy import ndimage ,signal,stats
import numpy as np
import re

import sys
sys.path.append('/home/tquah/PolyFTS_ALL/PolyFTS/tools') #gives us access to Kris' tools
from stats import *

from scipy.stats import mode
import os


def autoWarmupMSER_(data, debug=False):
    # Algorithm:
    # - Load all data
    # - Compute SEM with d samples truncated from beginning. The samples are supposed to be batch-averaged. This will be the case in most simulation data anyway.
    # - Find the value of d that minimizes SEM. This is the warmup length.
    loc = np.where(np.isnan(data)==False)
    data = data[loc]

    # MSER-5 - pre-block-average the data in chunks of 5 steps, and truncate integer counts of those steps
    m = 5 # Block size
    N = len(data)
    # Data will be padded in such a way that block average is not modified.
    # We are short by m-N%m points to make N/m an integer.
    # Then do the block averaging
    if N%m != 0:
        #dataint = np.pad(data, (0,m-N%m), mode='mean', stat_length=N%m)
        #databa = np.mean(dataint.reshape(int(N/m), m), axis=1)
        databa = np.mean(np.pad(data, (0,m-N%m), mode='mean', stat_length=N%m).reshape(int(N/m)+1, m), axis=1)
        N = N + m - N%m
    else:
        databa = np.mean(data.reshape(int(N/m), m), axis=1)

    if debug:
        fileout = open("debugMSER.dat","w")
    SEMlist=[]
    fullsize = len(databa)
    while True:
      # It is slow to do a full correlation-corrected statistical analysis in the inner loop
      # Rather than computing everything, just use correlation biased SEM
#      (n,(min,max),mean,semcc,kappa,unbiasedvar,autocor)=doStats(dummy,databa,False,False)
      sem = stats.sem(databa)
#      var = unbiasedvar * (n-1)/n # The MSER method uses the biased sample variance as a measure of homogeneity of the data
#      VARlist.append(var)
#      SEMlist.append(semcc)
      SEMlist.append(sem)
      if debug:
          fileout.write("{0} {1}\n".format(len(databa),sem))
      # Ensure we don't go below 25% of data remaining or 5 samples or that we don't iterate more than 1000 times.
      Nrem = len(databa)
      if Nrem <= 0.25*fullsize or Nrem < 5:
          break
      # TO DO: replace with masked array to avoid repeated copy
      databa = databa[1:] # Delete the first m elements

    if debug:
        fileout.close()

    # Find the index that minimizes the variance
    idx = np.argmin(SEMlist) # idx is therefore the number of warmup blocks to be removed; m*idx is # samples
    idx = idx*m

    if idx > 0:
        warmupdata = data[:idx-1]
    else:
        warmupdata = np.empty([0])
    proddata   = data[idx:]

    return warmupdata, proddata, idx

def GaussianFilter_DyDx_Peak(_y):
    _smeary =  \
        ndimage.gaussian_filter(_y,(np.max(_y)-np.min(_y)))
    #dx not necessary
    _dydx = (_smeary[1:]-_smeary[0:-1])
    _peaks = signal.find_peaks(np.abs(_dydx))
    return _peaks,_dydx,_smeary

def RemoveCorruptedOperatorRows(_path,newfile = "operators_clean.dat"):
    
    _op = open(_path,"r")
    _rows = _op.readlines()
    _neles = []
    for _row in _rows:
        _neles.append(len(_row.split(' ')))
    _op.close()
    _results = mode(np.array(_neles))
    _loc = np.where(np.max(_results[1])==_results[1])[0]
    _mele = _results[0][_loc][0]
    
    _lst = []
    _lst.append(_rows[0])
    for _row in _rows:
        if len(_row.split(' '))==_mele:
            _lst.append(_row)
    _op = open(newfile,'w')
    _op.writelines(_lst)
    _op.close()

def Generate_Test_Case(_minx,_maxx,_dx,_meany,_widthy,_numjumps):
    _jump = np.sort(np.random.randint(_minx,_maxx,_numjumps))
    _x = np.arange(_minx,_maxx+1e-6,_dx)
    _pts = len(_x)
    _y = np.zeros(_pts)
    _R = np.random.normal(_meany,_widthy,_pts)
    _state = 0

    for _j in _jump:
        if _state == 0:
            _state = 1
        else:
            _state = 0
        
        
        _loc = np.where(_j==_x)[0][0]
        _y[_loc:] = _state

    _y+=_R
    return _x,_y

def PickPeak_Histogram(_peaks,_dydx,_pts,_nstds=3.):
    _nobs,_minmax,_mean,_var,_skew,_kurt = stats.describe(np.abs(_dydx[_peaks[0]]))
    _loc = np.where((np.abs(_dydx[_peaks[0]])>=_mean+_nstds*np.sqrt(_var)))[0].tolist()
    _loc+= np.where((np.abs(_dydx[_peaks[0]])<=_mean-_nstds*np.sqrt(_var)))[0].tolist()
    _loc = np.unique(np.array(_loc))
    #need to switch to actual locations
    _aloc = []
    _aloc.append(0)
    for _i in range(len(_loc)):
        _j = np.where(_dydx[_peaks[0][_loc[_i]]]==_dydx)[0][0]
        _aloc.append(_j)
    _aloc.append(_pts-1)
    return _aloc,_loc,_nobs,_minmax,_mean,_var,_skew,_kurt

def AutomateCutoff(_y,_desiredvalue,_loc,_nstd = 3.,_tol = 1e-6):
    _n = len(_loc)
    
    _mean = np.zeros(_n-1)
    _sem = np.zeros(_n-1)
    for i in range(0,len(_loc)-1):
        _mean[i] = np.mean(_y[_loc[i]:_loc[i+1]])
        _sem[i] = stats.sem(_y[_loc[i]:_loc[i+1]])
    _minmean = _mean - (_nstd * _sem) - _tol
    _maxmean = _mean + (_nstd * _sem) + _tol
    _bounds = np.where((_desiredvalue>=_minmean) & (_desiredvalue<=_maxmean))[0]
    _ylist = []
    _rloc = []
    for _b in _bounds:
        _ylist+=_y[_loc[_b]:_loc[_b+1]].tolist()
        _rloc.append((_loc[_b],_loc[_b+1]))

    _ycorrect = np.array(_ylist)
    return _ycorrect,_rloc,_bounds,_mean,_sem,_minmean,_maxmean
def GetRows(_M,_rloc):
    _A = []
    for _pair in _rloc:
        _A.append(_M[_pair[0]:_pair[1],:])
    return np.vstack(_A)



def extract_value(string):
    return re.findall(r"[-+]?\d*\.\d+|\d+", string) #im lazy so i made a function 

def PrunePath(path,skip):
    splitpath = path.split('/')
    newpath = ''
    i = 0
    maxcount = len(splitpath)
    for ele in splitpath:
        newpath+=ele
        newpath+='/'
        i+=1
        if i==maxcount-skip:
            break
        elif skip==-1 and i==maxcount-1:
            break
    return newpath
        
def ConvertToTensor(array):
    a = np.zeros((3,3))
    a[0,0] = array[0,0]
    a[1,1] = array[1,0]
    a[2,2] = array[2,0]

    a[0,1] = array[3,0]
    a[0,2] = array[5,0]

    a[1,0] = array[3,0]
    a[1,2] = array[4,0]

    a[2,0] = array[5,0]
    a[2,1] = array[4,0]
    return a
def GetLocation(_path,_keylist):
    _list = _path.split('/')
    _vallist = []
    for _key in _keylist:
        _temp = extract_value(_list[_key])
        _templist = []
        if len(_temp)>1:
            for _ele in _temp:
                _templist.append(float(_ele))
            _vallist.append(np.array(_templist))
        else:
            _vallist.append(float(_temp[0]))
    return _vallist
def GetParameters(_filelist,_keyloc):
    _param_list = []
    for _file in _filelist:
        _param_list.append(GetLocation(_file,_keyloc))
    return _param_list

def GetOperators_Fancy(_file,_operators,_operator_track):
    _op = open(_file,'r')
    _col_list = []
    _col_track = decideColumn(_op,_operator_track)
    _warm,_prod,_idx = autoWarmupMSER(_op, _col_track)
    _warm,_data = extractData(_op,_col_track,_idx)
    _datalist = []
    _warmdatalist = []
    _datalist.append(_data)
    _warmdatalist.append(_warm)
    for _operate in _operators:
        _col = decideColumn(_op,_operate)
        _warm,_data = extractData(_op,_col,_idx)
        _datalist.append(_data)
        _warmdatalist.append(_warm)
    _op.close()
    _operators = [_operator_track]+_operators
    return _operators, np.vstack(_datalist).transpose(), np.vstack(_warmdatalist).transpose()

def GetOperators(_file,_operators):
    _op = open(_file,'r')
    _col_list = []
    _datalist = []
    _warmdatalist = []
    for _operate in _operators:
        _col = decideColumn(_op,_operate)
        _warm,_data = extractData(_op,_col,0)
        _datalist.append(_data)
        _warmdatalist.append(_warm)
    _op.close()
    return np.vstack(_datalist).transpose(), np.vstack(_warmdatalist).transpose()


def FDiff(dire,file1,file2,col,split=False,verbose = False):
    #file 1- file2
    templist = []
    path1 = os.path.join(dire,file1)
    path2 = os.path.join(dire,file2)
    t1 = np.loadtxt(path1)
    t2 = np.loadtxt(path2)
    #remove nan
    # loc = np.where(np.isnan(t1)==False)[0]
    # t1 = t1[loc,:]
    # loc = np.where(np.isnan(t2)==False)[0]
    # t2 = t2[loc,:]

    
    
    l1 = len(t1[:,col])  
    l2 = len(t2[:,col])
    diff = l1-l2
    if diff>0:
        warmupdata, proddata, idx = autoWarmupMSER_(t1[diff:,col]-t2[:,col])
    elif diff<0:
        warmupdata, proddata, idx = autoWarmupMSER_(t1[:,col]-t2[np.abs(diff):,col])
    else:
        warmupdata, proddata, idx = autoWarmupMSER_(t1[:,col]-t2[:,col])
    _nsamples,_bounds_,_mean_,_semcc_,_kappa_,_unbiasedvar_,_autocor_ = doStats(np.array([]), proddata)
    templist+= [_mean_,_semcc_]

    if verbose:
        print(f'Warmup Time: {len(warmupdata)} Production Time: {len(proddata)} Mean: {_mean_} Error: {_semcc_}')
    if split:
        
        
        warmupdata, proddata, idx = autoWarmupMSER_(t1[:,col])
        # print(proddata)
        _nsamples,_bounds_,_mean_,_semcc_,_kappa_,_unbiasedvar_,_autocor_ = doStats(np.array([]), proddata)
        templist+= [_mean_,_semcc_]
        warmupdata, proddata, idx = autoWarmupMSER_(t2[:,col])
        _nsamples,_bounds_,_mean_,_semcc_,_kappa_,_unbiasedvar_,_autocor_ = doStats(np.array([]), proddata)
        templist+= [_mean_,_semcc_]


    return templist
def SDiff(dire,file1,file2):
    lst = []

    operatorlist = ['StressXX.Real','StressYY.Real','StressZZ.Real','StressXY.Real','StressYZ.Real','StressXZ.Real']
    stress_path = os.path.join(dire,file1)
    data,warm =  GetOperators(stress_path,operatorlist)
    n = len(data)
    Astack = np.dstack(((data[:,0],data[:,3],data[:,5]),(data[:,3],data[:,1],data[:,4]),(data[:,5],data[:,4],data[:,2]))).reshape((n,3,3))
    value,vector = np.linalg.eig(Astack)
    Pmean = np.real(np.sum(value,axis = 1))
    warmupdata1, Prod1, idx = autoWarmupMSER_(Pmean)
    nsamples,bounds_S,mean_S,semcc_S,kappa_S,unbiasedvar_S,autocor_S = doStats(np.array([]), Prod1)
    lst+=[mean_S,semcc_S]
    stress_path = os.path.join(dire,file2)
    data,warm =  GetOperators(stress_path,operatorlist)
    n = len(data)
    Astack = np.dstack(((data[:,0],data[:,3],data[:,5]),(data[:,3],data[:,1],data[:,4]),(data[:,5],data[:,4],data[:,2]))).reshape((n,3,3))
    value,vector = np.linalg.eig(Astack)
    Pmean = np.real(np.sum(value,axis = 1))
    warmupdata1, Prod2, idx = autoWarmupMSER_(Pmean)
    nsamples,bounds_S,mean_S,semcc_S,kappa_S,unbiasedvar_S,autocor_S = doStats(np.array([]), Prod2)
    lst+=[mean_S,semcc_S]
    
    
    
    return lst




def Make_Header(_list):
    _str = ''
    for _ele in _list:
        _str+=_ele+' '
    return _str
