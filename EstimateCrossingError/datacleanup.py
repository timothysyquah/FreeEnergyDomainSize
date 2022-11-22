#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 19 19:06:57 2022

@author: tquah
"""
import pickle
import matplotlib.pyplot as plt
import numpy as np



with open('/home/tquah/DataStore/ODT_CL.pkl', 'rb') as f:
    cl = pickle.load(f)
with open('/home/tquah/DataStore/ODT_PFTS.pkl', 'rb') as f:
    pfts = pickle.load(f)
    


f = open('/home/tquah/DataStore/ODT_CL_clean.pkl',"wb")
pickle.dump(cl,f)
f.close()

f = open('/home/tquah/DataStore/ODT_PFTS_clean.pkl',"wb")
pickle.dump(pfts,f)
f.close()
