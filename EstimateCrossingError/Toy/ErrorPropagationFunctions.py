#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 11:01:45 2022

@author: tquah
"""
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats.distributions import t
def Linear(x,a,b):
    return a+b*x

#only does linear
def ErrorRegression(x,y,alpha):
    param,pcov = curve_fit(Linear,x,y)
    sigma = np.sqrt(np.diag(pcov))
    # n = len(x)
    # p = len(param)
    # dof = max(0,n-p)
    # tval = t.ppf(1.-alpha/2.,dof)
    # error = sigma*tval
    return param,sigma #error
def Intersection(b1,b2,m1,m2):
    db = (b2-b1)
    dm = (m2-m1)
    return db,dm,-db/dm
def IntersectionError(b1,b2,m1,m2,db1,db2,dm1,dm2):
    db,dm, xcross = Intersection(b1,b2,m1,m2)    
    var_b= (np.square(db1)+np.square(db2))/np.square(db)
    var_m= (np.square(dm1)+np.square(dm2))/np.square(dm)
    var = (var_b+var_m)
    return xcross, xcross*np.sqrt(var)

def EstimateIntersection(x1,y1,x2,y2,alpha = 0.05):
    param1,error1 = ErrorRegression(x1,y1,alpha)
    param2,error2 = ErrorRegression(x2,y2,alpha)
    xcross, crosserror = IntersectionError(param1[0],param2[0],param1[1],param2[1],\
                                           error1[0],error1[0],error2[1],error2[1])    
    return param1, param2, error1, error2, xcross, crosserror
