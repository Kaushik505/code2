# -*- coding: utf-8 -*-
"""
Created on Thu Feb  2 15:39:10 2023

@author: DEBOJYOTI PANDIT
"""
def simp2(n,phi,dl,inc):   
 import numpy as np
 import matplotlib.pyplot as plt   
 import math
 #dl=1/(n-1)  # interval length
 x=np.zeros((n,1));
 y=np.zeros((n,1));
 x[0]=inc;

 
 x[1]=x[0]+dl/2*((phi[0])+(phi[1]));
 
 
 for j in range(2,n):
    x[j]=x[j-2]+dl/3*((phi[j-2])+4*(phi[j-1])+(phi[j]));
   
    

    
 return (x)
