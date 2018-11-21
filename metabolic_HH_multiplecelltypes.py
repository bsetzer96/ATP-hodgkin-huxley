# -*- coding: utf-8 -*-
"""
Created on Mon Nov 19 20:12:04 2018

@author: bever
"""


import numpy as np
import matplotlib.pyplot as plt


## basic values for integration -------------------------------
t1=0            #start time
t2=8000           #end time
delta=0.05      #step size
h=delta         #runge kutta step
T= (t2-t1)/delta+1 #length of time vector
t=np.linspace(int(t1),int(t2),int(T)) #time vector
n= 10               #number of pyramidal cells
m= 4                 #number of in cells
o= 2                #number of TC cells
p = 2            #number of RE cells
ATP_scale=1 

#Pyramidal cell initial condition
y0_py=np.zeros(7*n)    #initial condition
y0_py[:n]=-70 #first n are voltage 
y0_py[n*4:n*5] = 10 #Na 0
y0_py[n*5:n*6]=2 #2nd to last n are ATP

#IN cell 
y0_in= np.zeros(4*m)
y0_in[:m]=-70 

#TC cell
y0_tc= np.zeros(4*o)
y0_tc[:o]=-70

#RE cell
y0_re=np.zeros(4*p)
y0_re[:p]=-70

#synapses
#y0_x = np.zeros(n+m)

def current_Na(m,h,V):
    ENa = 50
    gNa = 100
    am = 0.32*(V+54)/(1-np.exp(-1*(V+54)/4))
    bm = 0.28*(V+27)/(np.exp((V+27)/5)-1)
    ah = 0.128*np.exp(-1*(V+50)/18)
    bh = 4/(1+np.exp(-1*(V+27)/5))
    mdot= am*(1-m)-bm*m
    hdot = ah*(1-h)-bh*h
    current = gNa*m*m*m*h*(V-ENa)
    return [mdot, hdot, current]

def current_K(m,V):
    EK = -100
    gK = 80
    am = 0.032*(V+52)/(1-np.exp(-1*(V+52)/5))
    bm = 0.5*(np.exp(-1*(V+57)/40))
    mdot = am*(1-m)-bm*m
    current = gK*m*m*m*m*(V-EK)
    return [mdot, current]

def current_leak(gL,Erev,V):
    current = gL*(V-Erev)
    return [current]

def current_M(m,V):
    EM = -100
    gM = 2
    Qs = 3.209
    am = Qs*(0.0001)*(V+30)/(1-np.exp(-1*(V+30)/9))
    bm = -1*Qs*(0.0001)*(V+30)/(1-np.exp((V+30)/9))
    mdot = am*(1-m)-bm*m
    current = gM*m*(V-EM)
    return [mdot, current]


def nrn_IN(m,h,n,w,V):
    Erev=-67
    gL=0.1
    I_Na=current_Na(m,h,V)
    I_K=current_K(n,V)
    I_M=current_M(w,V)
    I_leak=current_leak(gL,Erev,V)
    Vdot= Iapp-I_Na-I_K-I_M-I_leak
    return Vdot

def current_Na_des96(m,h,V):
    ENa = 50
    gNa = 90
    Vtraub = -35
    Vt = V-Vtraub
    am = 0.32 * (13-Vt)/(np.exp((13-Vt)/4) - 1)
    bm = 0.28 * (Vt-40)/(np.exp((Vt-40)/5) - 1)
    ah = 0.128 * np.exp((17-Vt)/18)
    bh = 4/(1 + np.exp((40-Vt)/5))
    mdot= am*(1-m)-bm*m
    hdot = ah*(1-h)-bh*h
    current = gNa*m*m*m*h*(V-ENa);
    return [mdot, hdot, current]

def current_K_Des96: ### not in folderrr
    
def current_h_des96(c1, p0, o1, V, Ca):
    cac	= 0.002
    k2	= 0.0004 #;// (1/ms)		: inverse of time constant
    Pc	= 0.01 #;//			: half-activation of CB protein dependence
    k4	= 0.001 #;//	(1/ms)		: backward binding on Ih
    nca	= 4 #;//			: number of binding sites of ca++
    nexp	= 1 #;//			: number of binding sites on Ih channels
    ginc	= 2 #;//			: augmentation of conductance with Ca++
    taum	= 20 #;//	(ms)		: min value of tau
    shift = 0 #;//	(mV)	
    p1 = 1-p0 #;
    tadj = 1
    Eh = -40
    gh = 0.025
    h_inf = 1 / ( 1 + exp((v+75-shift)/5.5) )
    tau_s = (taum + 1000 / ( exp((v+71.5-shift)/14.2) + exp(-(v+89-shift)/11.6) ) ) / tadj;
    alpha = h_inf / tau_s
    beta  = ( 1 - h_inf ) / tau_s
    k1ca = k2 * pow((cai/cac),nca)
    k3p = k4 * pow((p1/Pc),nexp)
    c1dot= beta*o1-alpha*c1
    p0dot = k2*(1-p0)-k1ca*p0
    o1dot = k4*(1-c1-o1)-k3p*o1
    current = gh*(o1+(2*(1-c1-o1)))*(V-Eh);
    return current

def current_T

def nrn_RE(m,h,n,w,V):
    g_L=0.05
    E_L= -90
    I_Na=current_Na_des96(m,h,V)
    I_K=
    I_TRE=
    I_leak=current_leak(g_L,E_L, V)
    current=Iapp-I_Na-I_K-I_TRE-I_leak
    







