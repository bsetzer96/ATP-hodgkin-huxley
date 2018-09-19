# -*- coding: utf-8 -*-
"""
Spyder Editor
neurophysiological model for burst suppression
This is a temporary script file.
"""

import numpy as np
import math
import matplotlib.pyplot as plt

##time
t1=0            #start time
t2=1500           #end time
delta=0.05      #step size
h=delta         #runge kutta step
T= (t2-t1)/delta+1 #length of time vector
t=np.linspace(int(t1),int(t2),int(T)) #time vector


## paramenters
y0_py=np.array([-70, 0, 0, 0, 0, 2])    #initial condition
y0_fs= np.array([-70, 0, 0, 0])
N=1                   #number of neurons
E_Na= 50            #Sodium equilibrium potential
E_k= -100             #potassium equilibrium potential
g_Na= 100             #Sodium conductance
g_K= 80               #potassium conductance
g_K_ATP = 0.15        #potassium ATP conductance
I_app_py= 1.8        #applied current
J_ATP= 2              #production rate of ATP (linked to metabolism)
ATP_max= 0.6          # ?
K_m= 6*(10**(-8))     # govern the NA-ATP pump dynamics
F= 8.8*(10**(-5))     # help govern the NA-ATP pump dynamics

p_py=np.array([E_Na, E_k, g_Na, g_K, g_K_ATP, I_app_py, J_ATP, ATP_max, K_m, F])
p_py=np.transpose(p_py)

I_app_fs=.5
p_fs=np.array([E_Na, E_k, g_Na, g_K, I_app_fs])
p_fs=np.transpose(p_fs)

##preallocate space
y_py=np.zeros((6,int(T)))
y_py[:,0]=y0_py

y_fs=np.zeros((4,int(T)))
y_fs[:,0]=y0_fs

## Hodgkin huxley network of neurons with ATP metabolism
def hodghuxATP(state, N,p):
    #N = number of neurons in network
    v= state[0]
    m= state[1]
    n= state[2]
    h= state[3]
    Na= state[4]
    ATP= state[5]
    
    E_Na= p[0]
    E_K= p[1]
    g_Na= p[2]
    g_K = p[3]
    g_K_ATP = p[4]
    I_app= p[5]
    J_ATP = p[6]
    ATP_max = p[7]
    K_m = p[8]
    F = p[9]
    
    alpha_m = ((0.32)*(v+54))/(1-math.exp(-(v+54)/4))
    beta_m = ((0.28)*(v+27))/(math.exp((v+27)/5)-1)
    alpha_h = (0.128)*math.exp(-(v+50)/18)
    beta_h = 4/(1+math.exp(-(v+27)/5))
    #in the paper it calls n m... but it has a different alpha function
    alpha_n = ((0.032)*(v+52))/(1-math.exp(-(v+52)/5))
    beta_n = 0.5*math.exp(-(v+57)/40)
    z = 1/(1+6*ATP)
    
    mdot = alpha_m*(1-m) - beta_m*m
    hdot= alpha_h*(1-h) - beta_h*h
    ndot = alpha_n*(1-n)-beta_n*n

    
    I_Na = g_Na*(m**3)*h*(v-E_Na)
    I_K = g_K*(n**4)*(v-E_K)
    I_leak = 0.1*(v+61)
    I_K_ATP = g_K_ATP*z*(v-E_K)
    
    Nadot = F*I_Na-3*K_m*(Na**3)*ATP
    ATPdot = J_ATP*(ATP_max - ATP) - K_m*(Na**3)*ATP
    
    vdot = I_app-I_Na-I_K-I_K_ATP-I_leak
    
    dx= np.array([vdot, mdot, ndot, hdot, Nadot, ATPdot])
    
    return dx

def hodghuxFS(state, N, p):
    v= state[0]
    m= state[1]
    n= state[2]
    h= state[3]
    
    E_Na= p[0]
    E_K= p[1]
    g_Na= p[2]
    g_K = p[3]
    I_app= p[4]

    
    alpha_m = ((0.32)*(v+54))/(1-math.exp(-(v+54)/4))
    beta_m = ((0.28)*(v+27))/(math.exp((v+27)/5)-1)
    alpha_h = (0.128)*math.exp(-(v+50)/18)
    beta_h = 4/(1+math.exp(-(v+27)/5))
    alpha_n = ((0.032)*(v+52))/(1-math.exp(-(v+52)/5))
    beta_n = 0.5*math.exp(-(v+57)/40)
    
    mdot = alpha_m*(1-m) - beta_m*m
    hdot= alpha_h*(1-h) - beta_h*h
    ndot = alpha_n*(1-n)-beta_n*n

    
    I_Na = g_Na*(m**3)*h*(v-E_Na)
    I_K = g_K*(n**4)*(v-E_K)
    I_leak = 0.1*(v+61)
 
    vdot = I_app-I_Na-I_K-I_leak
    
    dx= np.array([vdot, mdot, ndot, hdot])
    
    return dx

##Runge kutta order 4 integration of pyramidal cell
for n in np.arange(T-1):
    k1=h*hodghuxATP(y_py[:,int(n)],N,p_py)
    k2=h*hodghuxATP(y_py[:,int(n)]+k1/2,N,p_py)
    k3=h*hodghuxATP(y_py[:,int(n)]+k2/2,N,p_py)
    k4=h*hodghuxATP(y_py[:,int(n)]+k3,N,p_py)
    y_py[:,int(n+1)]=y_py[:,int(n)]+(1/6)*(k1+2*k2+2*k3+k4)


plt.plot(t,y_py[0,:])
plt.show 

##Runge kutta order 4 integration of FS interneuron cell
for n in np.arange(T-1):
    k1=h*hodghuxFS(y_fs[:,int(n)],N,p_fs)
    k2=h*hodghuxFS(y_fs[:,int(n)]+k1/2,N,p_fs)
    k3=h*hodghuxFS(y_fs[:,int(n)]+k2/2,N,p_fs)
    k4=h*hodghuxFS(y_fs[:,int(n)]+k3,N,p_fs)
    y_fs[:,int(n+1)]=y_fs[:,int(n)]+(1/6)*(k1+2*k2+2*k3+k4)
    
plt.plot(t,y_fs[0,:])

## Q: Can we go over the set up of the code.... what different functions should I have? 
#What is done seperately what can be done together?