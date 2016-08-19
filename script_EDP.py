# -*- coding: utf-8 -*-
"""
@author: Ahmed REBAI
"""

import numpy as np
import matplotlib.pyplot as plt


def Explicit_Euler(q0,ti,tf,N,fun):
    """
    This function computes and returns a N+1 dimension array of the quantity q 
    containing the numerical solution at instants ti + i*h, with h is the time 
    step defined by h = ( tf - ti ) / N
    
    For the input arguments:
    fun: function of the quantity q_n and the time t
    ti : the initial time
    tf : the final time
    q0 : the initial condition at ti
    N  : number of time steps between initial time ti and final time tf
    """
    
    h = ( tf - ti ) / N
    q = np.zeros((1,N+1), dtype = float )
    q[0,0] = q0
    
    for i in range(N):
        q[0,i+1] = q[0,i] + h * fun(q[0,i]) 

    return q
    
fun = lambda arg1: -4*arg1;
q0 = 1
ti = 0
tf = 3
N = 24
q = Explicit_Euler(q0,ti,tf,N,fun)
time = np.linspace(0,3,25)
plt.figure(1)
plt.plot(time,q[0],'r-')
plt.xlabel("Time in seconds")
plt.ylabel("q physical quantity")
plt.title("Radioactive decay partial differential equation numerical resolution")
plt.savefig("numerical_solution.png")