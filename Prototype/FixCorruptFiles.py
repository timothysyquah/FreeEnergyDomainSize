#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  1 17:29:01 2022

@author: tquah
"""

from Functions import *
import glob
import os
import argparse
if __name__ == "__main__":
    IDIR = os.getcwd()
    parser = argparse.ArgumentParser(description='Tool to Clean Corrupted Files')
    parser.add_argument('-d', '--dirs', action='store', default="",help='list of directories that contain each phase point',type=str)
    parser.add_argument('-i','--input',action = 'store',default = '', help = 'input file', type = str)
    parser.add_argument('-o','--output',action = 'store',default = '', help = 'output file', type = str)
    args = parser.parse_args()

    path = os.path.join(args.dirs,args.input)

    file_list = glob.glob(path)
    for file in file_list:
        directory = PrunePath(file,1)
        outpath = os.path.join(directory,args.output)
        RemoveCorruptedOperatorRows(file,outpath)