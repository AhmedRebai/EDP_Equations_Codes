# -*- coding: utf-8 -*-
"""
Created on Sun Aug 21 22:49:00 2016

@author: Ahmed Rebai
"""
import numpy as np
import matplotlib.pyplot as plt

def Runge_Kutta4(fun,u0,ti,tf,N):
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
        
    h = (tf-ti)/N
    q = np.zeros((1,N+1),dtype=float)
    q[0,0] = q0
    
    for i in range(0,n):
        k1 = fun(q[0,i],ti)
        k2 = fun(q[0,i]+0.5*h*k1,ti+0.5*h)
        k3 = fun(q[0,i]+0.5*h*k2,ti+0.5*h)
        k4 = fun(q[0,i]+h*k3,ti+h)
        q[0,i+1] = q[0,i] + (h/6)*(k1+2*k2+2*k3+k4)
        ti = ti + h
        
    return q

fun = lambda arg1,arg2: -4*arg1;
q0 = 1
ti = 0
tf = 3
N = 100
q = Runge_Kutta4(fun,q0,ti,tf,N)
time = np.linspace(0,3,n+1)
plt.figure(1)
plt.plot(time,q[0],'r-')
plt.xlabel("Time in seconds")
plt.ylabel("q physical quantity")
plt.title("Radioactive decay partial differential equation numerical resolution")
plt.savefig("numerical_solution_Runge_Kutta4.png")