# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 17:40:20 2020

@author: bever

attempt at converting to Brian2
"""
from brian2 import *

eq_FS='''

dV/dt=(Iapp-Isyn-IL-INa-IK-Iran)/Cm : volt
Isyn=IsynFS+ IsynPY : amp * meter ** -2
IsynFS : amp * meter ** -2
IsynPY : amp * meter ** -2

IL=gL*(V-VL) : amp * meter ** -2
Iran=ransc*randn() : amp * meter ** -2 (constant over dt)
INa=gNa*m**3*h*(V-VNa) : amp * meter ** -2
    alpha_m = ((0.32/(mV*ms))*(V+54*mV))/(1-exp(-(V+54*mV)/(4*mV)))   : Hz
    beta_m = ((0.28/(mV*ms))*(V+27*mV))/(exp((V+27*mV)/(5*mV))-1) : Hz
    dm/dt = alpha_m*(1-m) - beta_m*m : 1
    alpha_h = (0.128/ms)*exp(-(V+50*mV)/(18*mV)) : Hz
    beta_h = (4/ms)/(1+exp(-(V+27*mV)/(5*mV))) : Hz
    dh/dt= (alpha_h*(1-h) - beta_h*h) : 1
IK=gK*n**4*(V-VK) : amp * meter ** -2
    alpha_n = ((0.032/(mV*ms))*(V+52*mV))/(1-exp(-(V+52*mV)/(5*mV))) : Hz
    beta_n = (0.5/ms)*exp(-(V+57*mV)/(40*mV)): Hz
    dn/dt = alpha_n*(1-n)-beta_n*n : 1
    


'''


##Constants :

Cm = 1.0* ufarad * cm ** -2
gL=0.1 * msiemens * cm **-2
VL=-61*mV

gNa=100 * msiemens * cm **-2
VNa=50*mV

gK=80 * msiemens * cm **-2
VK=-100*mV

Iapp = 1.5 *amp*(10**-6) * cm ** -2 #*amp*meter**-2 #

ransc=0.5*(amp*(10**-6) * cm ** -2)

if __name__=='__main__' :
    

    FS=NeuronGroup(1,eq_FS,method='rk4')
    FS.V = '-65*mvolt'
    FS.h = '0.9932525096349422'
    FS.m = '0.022083418537191028'
    FS.n = '0.05182107474228339'
    
    V1=StateMonitor(FS,'V',record=True)
    INa=StateMonitor(FS,'INa',record=[0])
    IK=StateMonitor(FS,'IK',record=[0])
    Iran=StateMonitor(FS, 'Iran', record=[0])
    
    all_neurons= FS
    #all_synapses= S_PYPY
    all_monitors= V1, INa, IK, Iran
    
    
    run(0.5*second)
    
    figure()
    plot(V1.t/second,V1.V[0]/mV)
    xlabel('Time (s)')
    ylabel('Membrane potential (mv)')
    title('PY cell')

    figure()
    ina=all_monitors[1]
    plot(ina.t/second,ina.INa[0]/(amp*(10**-6) * cm**-2))
    title('INa')  
    xlabel('Time (s)')
    ylabel('INa current (amp*(10^-6) * cm^-2)')
    
    figure()
    ik=all_monitors[2]
    plot(ik.t/second,ik.IK[0]/(amp*(10**-6) * cm**-2))
    xlabel('Time (s)')
    ylabel('IK current (amp*(10^-6) * cm^-2)')
    title('IK')  
     
    figure()
    iran=all_monitors[3]
    plot(ik.t/second,iran.Iran[0]/(amp*(10**-6) * cm**-2))
    xlabel('Time (s)')
    ylabel('Iran current (amp*(10^-6) * cm^-2)')
    title('Iran')  
    
    clear_cache('cython')
