# -*- coding: utf-8 -*-
"""
Created on Thu Feb  2 15:39:10 2023

@author: DEBOJYOTI PANDIT
"""
def simp11(n,phi,dl):   
 import numpy as np
 import matplotlib.pyplot as plt   
 import math
 #dl=1/(n-1)  # interval length
 x=np.zeros((n,1));
 y=np.zeros((n,1));
 x[0]=0;
 y[0]=0;
 
 x[1]=dl/2*(np.cos(phi[0])+np.cos(phi[1]));
 y[1]=dl/2*(np.sin(phi[0])+np.sin(phi[1]));
 
 for j in range(2,n):
    x[j]=x[j-2]+dl/3*(np.cos(phi[j-2])+4*np.cos(phi[j-1])+np.cos(phi[j]));
    y[j]=y[j-2]+dl/3*(np.sin(phi[j-2])+4*np.sin(phi[j-1])+np.sin(phi[j]));
    

    
 return (x,y)
