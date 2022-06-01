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

import re


"""
Here we will outline the code before writing it
1) use argparse to specify
* Directory Structure
* Specify Keywords 
* Tolerance (free energy)
* Activate Filters-remove bad CL....
* Keywords for averaging
* N-and definitions of C
* This tool will not plot!
"""
def extract_value(string):
    return re.findall(r"[-+]?\d*\.\d+|\d+", string) #im lazy so i made a function 

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

"""
2) Load data and characterize it based on directory structure and keywords
* may need to do cleaning step here
"""

if __name__ == "__main__":
    
    IDIR = os.getcwd()
    parser = argparse.ArgumentParser(description='Tool to Calculate Free Energy')
    #do not mix simulation methods
    parser.add_argument('-d', '--dirs', action='store', default="",help='list of directories that contain each phase point',type=str)
    parser.add_argument('-f','--file',action = 'store',default = '', help = 'File to read with averages and error', type = str)
    parser.add_argument('-kl','--keyloc',action = 'store',nargs='+',default = [], help = 'keywords location', type = int)
    parser.add_argument('-kt','--keytxt',action = 'store',nargs='+',default = [''], help = 'keywords label', type = str)
    
    args = parser.parse_args()
    
    path = os.path.join(args.dirs,args.file)
    filelist = glob.glob(path)
    print(filelist)
    paramlist = GetParameters(filelist,args.keyloc)    
    # we only allow a single parameter per line
    param_array = np.vstack(paramlist)
    print(param_array)

"""
3) Average the operators 
* Gather volscl
* Gather mu
* Find the series that is shorter
* Compute vectorized Free Energy
* Push through the tools in stats
** Get back averaged values with standard error of the mean
"""

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

