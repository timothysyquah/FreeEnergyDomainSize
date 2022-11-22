#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 19 19:06:57 2022

@author: tquah
"""
import sys
sys.path.append('/home/tquah/FreeEnergyDomainSize/Alternate_Analysis')
from Functions_Analysis import *
sys.path.append('/home/tquah/PolyFTS_ALL/PolyFTS/tools/')
from stats import * 
import pickle
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
plt.close('all')

def LoadPKL(_file):
    with open(_file, 'rb') as f:
        return pickle.load(f)


def UnpackHeader(_dict):
    lst = list(_dict)
    rlst = []
    for ele in lst:
        rlst.append(np.array(ele))
    return np.vstack(rlst)
    
def convert_array_tuplelist(array):
    #assume rows are each 
    shape = np.shape(array)
    tuple_list = []
    for i in range(0,shape[0]):
        tuple_list.append(tuple(array[i,:].tolist()))
    return tuple_list

def StackArray(dic,tuplist):
    data_array = []

    for i in range(0,len(tuplist),1):
        FLAM = dic[tuplist[i]]['LAM'][:,2]
        FDIS = dic[tuplist[i]]['DIS'][:,2]
        try:
            diff = FLAM-FDIS
        except:
            if len(FLAM)>len(FDIS):
                diff = FLAM[1:]-FDIS
            else:
                diff = FLAM-FDIS[1:]

        # temp = np.zeros((len(diff),3))
        # temp[:,0] = tuplist[i][0]
        # temp[:,1] = tuplist[i][1]
        # temp[:,2] = diff
        warm,prod,idx = autoWarmupMSER_(diff)
        nsamples,(minx,maxx),mean,semcc,kappa,unbiasedvar,autocor = doStats(warm, prod)
        templist = list(tuplist[i])+[mean,semcc]
        data_array.append(np.array(templist))
    return np.vstack(data_array)



def linear(x,a,b):
    return a*x+b


def weighted_least_squares(x,y,sigmay):
    vary = np.square(sigmay)
    x.reshape(1,len(x))
    y.reshape(len(y),1)
    A = np.vstack((x,np.ones_like(x))).T
    w = np.diag(1/vary)
    
    ATWA = np.linalg.multi_dot((A.T,w,A))
    ATWy = np.linalg.multi_dot((A.T,w,y))
    beta = np.dot(np.linalg.inv(ATWA),ATWy)
    error = np.square(np.diag(np.linalg.inv(ATWA)))
    return beta,error
    

def linear_regression(data_array):
    param,pcov = curve_fit(linear,data_array[:,0],data_array[:,1],sigma = data_array[:,2])
    return param, np.sqrt(np.diag(pcov))

#plot structure factor determine C from first tuple element
# should improve it using https://en.wikipedia.org/wiki/Weighted_least_squares
#https://stats.stackexchange.com/questions/15511/estimating-the-intersection-of-two-lines
cl  = LoadPKL( '/home/tquah/DataStore/ODT_CL.pkl')
pfts  = LoadPKL('/home/tquah/DataStore/ODT_PFTS.pkl')

cl_header_array = UnpackHeader(cl)
pfts_header_array = UnpackHeader(pfts)


cl_dataset = np.loadtxt('CL.dat')
cl_dataset[:,0]/=10
cl_dataset[:,1]/=100
cl_dataset[:,1]= np.round(cl_dataset[:,1],4)
cl_tup_list = convert_array_tuplelist(cl_dataset)


cl_array = StackArray(cl,cl_tup_list)
cl_uniqueC= np.unique(cl_array[:,0])

plot_cl = []

N = 2000

for i in range(0,len(cl_uniqueC)):
    loc = np.where(cl_array[:,0]==cl_uniqueC[i])[0]
    param,pcov = linear_regression(cl_array[loc,1:])
    # print(-param[1]/param[0])
    # print(cl_uniqueC[i])
    # plot_cl.append(np.array([cl_uniqueC[i],-param[1]/param[0]]))
    beta, error = weighted_least_squares(cl_array[loc,1],cl_array[loc,2],cl_array[loc,3])
    print(cl_uniqueC[i])
    print(f'direct: {beta[0]} +- {error[0]}')
    print(f'direct: {beta[1]} +- {error[1]}')
    print(f'{param[0]} +- {pcov[0]}')
    print(f'{param[1]} +- {pcov[1]}')


    if error[0] == np.inf:
        x0 = np.nan
        x0_err = np.nan
        plot_cl.append(np.array([cl_uniqueC[i],-param[1]/param[0],x0,x0_err]))
        continue
    else:
        # e1 = np.random.normal(beta[0],error[0],N)
        # e2 = np.random.normal(beta[1],error[1],N)
        
        # intercept = -param[1]/param[0]
        # x0 = np.mean(intercept)
        # x0_err = np.std(intercept)
        plot_cl.append(np.array([cl_uniqueC[i],-beta[1]/beta[0],-beta[1]/beta[0],error[1]/error[0]]))
    
    
    
    

    plt.figure()
    
    plt.title(f'C = {cl_uniqueC[i]}')
    
    plt.plot(cl_array[loc,1]*100,cl_array[loc,2],'ok')

    xmin = np.min(cl_array[loc,1])
    xmax = np.max(cl_array[loc,1])
    x = np.linspace(xmin,xmax,100)
    # plt.plot(x*100,param[0]*x+param[1],'r')
    # # xmat = np.ones((N,100))*x
    # for j in range(0,N,1):
    #     plt.plot(x*100,e1[j]*x+e2[j],'k',alpha = 0.2,zorder = 5)


# plt.plot(plot_cl[:,0],plot_cl[:,1])                  
pfts_dataset = np.loadtxt('PFTS.dat')
pfts_dataset[:,0]/=10
pfts_dataset[:,1]/=100
pfts_dataset[:,1]= np.round(pfts_dataset[:,1],4)
pfts_tup_list = convert_array_tuplelist(pfts_dataset)


pfts_array = StackArray(pfts,pfts_tup_list)

pfts_uniqueC= np.unique(pfts_array[:,0])
plot_pfts = []

for i in range(0,len(pfts_uniqueC)):
    loc = np.where(pfts_array[:,0]==pfts_uniqueC[i])[0]
    # param,pcov = linear_regression(pfts_array[loc,1:])
    beta, error = weighted_least_squares(pfts_array[loc,1],pfts_array[loc,2],pfts_array[loc,3])

    # print(-param[1]/param[0])
    print(pfts_uniqueC[i])
    # plot_pfts.append(np.array([pfts_uniqueC[i],-param[1]/param[0]]))
    # e1 = np.random.normal(param[0],pcov[0],N)
    # e2 = np.random.normal(param[1],pcov[1],N)
    # intercept = -e2/e1
    x0 = -beta[1]/beta[0]
    x0_err = error[1]/error[0] #np.std(intercept)/np.sqrt(N)
    plot_pfts.append(np.array([pfts_uniqueC[i],-beta[1]/beta[0],x0,x0_err]))
    
    
    print(f'{param[0]} +-{pcov[0]}')
    print(f'{param[1]} +-{pcov[1]}')

    print(f'{x0} +- {x0_err}')
    
    
    
    
    
    
    
    
    

plt.figure()
    
plot_cl = np.vstack(plot_cl)     
#plt.errorbar(plot_cl[:,0]*10,plot_cl[:,2]*100,yerr = plot_cl[:,3]*100,color = 'r', marker = '^')                  
plt.errorbar(plot_cl[:,0]*10,plot_cl[:,1]*100,color = 'k', marker = '^',label = 'CL')                  

plot_pfts = np.vstack(plot_pfts)   
#plt.errorbar(plot_pfts[:,0]*10,plot_pfts[:,2]*100,yerr = plot_cl[:,3]*100,color = 'b', marker = '^')                  
plt.errorbar(plot_pfts[:,0]*10,plot_pfts[:,1]*100,color = 'g', marker = '^',label = 'PartialFTS')                       

# plt.plot(plot_pfts[:,0]*10,plot_pfts[:,1],'ok')                  
plt.xscale('log')
plt.xlabel('$C$')
plt.ylabel('$(\chi N)_{ODT}$')
plt.legend()

plt.savefig('/home/tquah/Figures/ODT_cross.pdf')

plt.figure()
    
plot_cl = np.vstack(plot_cl)     
plt.errorbar(plot_cl[:,0]*10,plot_cl[:,2]*100,yerr = plot_cl[:,3]*100,color = 'r', marker = '^',label = 'CL')             
#plt.errorbar(plot_cl[:,0]*10,plot_cl[:,1]*100,color = 'k', marker = '^')                  

plot_pfts = np.vstack(plot_pfts)   
plt.errorbar(plot_pfts[:,0]*10,plot_pfts[:,2]*100,yerr = plot_cl[:,3]*100,color = 'b', marker = '^',label = 'PartialFTS')                
#plt.errorbar(plot_pfts[:,0]*10,plot_pfts[:,1]*100,color = 'g', marker = '^')                  
plt.legend()
# plt.plot(plot_pfts[:,0]*10,plot_pfts[:,1],'ok')                  
plt.xscale('log')
plt.xlabel('$C$')
plt.ylabel('$(\chi N)_{ODT}$')

plt.savefig('/home/tquah/Figures/ODT_MC_cross.pdf')
# cl_C = np.unique(cl_header_array[:,0])
# pfts_C = np.unique(pfts_header_array[:,0])


# plot difference in order parameter

# for i in range(0,len(cl_C)):
#     plt.figure()
#     plt.title(fr'$C = {cl_C[i]}$')
#     loc = np.where(cl_C[i]==cl_header_array[:,0])[0]
#     tuple_list = convert_array_tuplelist(cl_header_array[loc,:])
#     for j in range(0,len(tuple_list)):
#         try:
#             OPLAM = cl[tuple_list[j]]['model2_OrientationPersistenceOP'][:,1]
#             OPDIS = cl[tuple_list[j]]['model1_OrientationPersistenceOP'][:,1]
#             # warmup,data,idx = autoWarmupMSER_(OPLAM[1:])
#             # plt.scatter(tuple_list[j][1]*np.ones_like(data),data)
#             # warmup,data,idx = autoWarmupMSER_(OPDIS[1:])
#             # plt.scatter(tuple_list[j][1]*np.ones_like(data),data)
#             warmup,data,idx = autoWarmupMSER_(OPLAM[1:]-OPDIS[1:])
#             plt.scatter(tuple_list[j][1]*np.o

#         except:
#             continue
#             # deltaOP = 
#     plt.ylabel('OP')

# plot          difference in free energy
# for i in range(0,len(cl_C)):
#     plt.figure()
#     plt.title(fr'$C = {10*cl_C[i]}$')
#     loc = np.where(cl_C[i]==cl_header_array[:,0])[0]
#     tuple_list = convert_array_tuplelist(cl_header_array[loc,:])
#     for j in range(0,len(tuple_list)):
#         try:
            
            
#             if cl[tuple_list[j]]['ST'] ==1:
#                 continue
#             FLAM = cl[tuple_list[j]]['LAM'][:,2]
#             FDIS = cl[tuple_list[j]]['DIS'][:,2]
#             # warmup,data,idx = autoWarmupMSER_(OPLAM[1:])
#             # plt.scatter(tuple_list[j][1]*np.ones_like(data),data)
#             # warmup,data,idx = autoWarmupMSER_(OPDIS[1:])
#             # plt.scatter(tuple_list[j][1]*np.ones_like(data),data)
#             warmup,data,idx = autoWarmupMSER_(FLAM[1:]-FDIS[1:])
#             # plt.scatter(tuple_list[j][1]*np.ones_like(data),data)
#             nsamples,(minx,maxx),mean,semcc,kappa,unbiasedvar,autocor = doStats([[]], data)
#             plt.errorbar(tuple_list[j][1],mean,yerr = semcc,linestyle = 'none',zorder = 1,marker = 's',color = 'k')

#         except:
#             continue
#             # deltaOP = 
#     plt.ylabel('F')



