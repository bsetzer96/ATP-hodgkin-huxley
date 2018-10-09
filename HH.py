# -*- coding: utf-8 -*-
"""
Created on Thu Oct  4 13:08:13 2018

@author: bever
"""

import numpy as np
import matplotlib.pyplot as plt

E_Na=60
g_Na=120
E_K=-77
g_K=35
E_L= -54.4
g_L=.3
I_app=15
y_0=np.array([-80,.1,.4,.4])
t1=0
t2=10
h=0.01
T=int((t2/h)+1)
t=np.linspace(t1,t2,T)
y=np.zeros((4,T))
y[:,0]=y_0

p = np.array([E_Na,g_Na, E_K,g_K, E_L,g_L,I_app])

def HH(state,p):
    E_Na=p[0]
    g_Na=p[1]
    E_K=p[2]
    g_K=p[3]
    E_L=p[4]
    g_L=p[5]
    I_app=p[6]
    
    v=state[0]
    m=state[1]
    h=state[2]
    n=state[3]
    
    I_Na=-g_Na*(m**3)*h*(v-E_Na)
    I_K=-g_K*(n**4)*(v-E_K)
    I_L=g_L*(v-E_L)
    
    beta_m=4*np.exp((-v-65)/18)
    alpha_m=(2.5-0.1*(v+65))/(np.exp(2.5-0.1*(v+65))-1)
    beta_n=0.125*np.exp((-v-65)/80)
    alpha_n=(0.1-0.01*(v+65))/((np.exp((1-0.1*(v-65))/80))-1)
    beta_h=1/(np.exp(3-0.1*(v+65))+1)
    alpha_h=0.07*np.exp((-v-65)/20)
    
    dv= I_Na+I_K+I_L+I_app
    dm= alpha_m*(1-m)-beta_m*m
    dh= alpha_n*(1-n)-beta_n*n
    dn= alpha_h*(1-h)-beta_h*h
    
    dx=np.array([dv,dm,dh,dn])
    return dx

for i in np.arange(T-1):
    k1= h*HH(y[:,int(i)],p) #first step
    k2=h*HH(y[:,int(i)]+k1/2,p) #second step
    k3=h*HH(y[:,int(i)]+k2/2,p) #third step
    k4=h*HH(y[:,int(i)]+k3,p) #fourth step
    y[:,int(i+1)]=y[:,int(i)]+(1/6)*(k1+2*k2+2*k3+k4) #RK order 4

plt.plot(t,y[0,:])    
    