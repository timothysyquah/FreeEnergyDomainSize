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

if __name__ == "__main__":
    
    IDIR = os.getcwd()
    parser = argparse.ArgumentParser(description='Tool to Calculate Free Energy')
    parser.add_argument('-d', '--dirs', action='store', default="",help='list of directories that contain each phase point')
    parser.add_argument('-f','--file',action = 'store',default = '', help = 'File to read with averages and error', type = str)
    parser.add_argument('-k','--keyword',action = 'store',nargs='+',default = [""], help = 'keywords', type = str)
    
    args = parser.parse_args()

    path = os.path.join(args.dirs,args.file)
    filelist = glob.glob(path)

"""
2) Load data and characterize it based on directory structure and keywords
* may need to do cleaning step here
"""
    

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

