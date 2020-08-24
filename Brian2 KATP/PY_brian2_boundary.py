# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 18:15:25 2020

@author: bever
"""


from brian2 import *


prefs.codegen.target = 'numpy'
prefs.codegen.target = 'cython'


defaultclock.dt = 0.01*ms

eq_PY='''
dV/dt=(Iapp-Isyn-IL-INa-IK-IKATP-Iran)/Cm : volt
Isyn=IsynFS+ IsynPY : amp * meter ** -2
IsynFS : amp * meter ** -2
IsynPY : amp * meter ** -2
Iapp : amp * meter ** -2

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
    
IKATP=gKATP*z*(V- VK)  : amp*meter **-2
    z = (1*mM)/(1*mM+ATP/0.06) : 1
    dNa/dt= -(F*mM/ms)*INa/(amp*(10**-6) * cm ** -2)-((3/ms)*K_m*ATP*(Na/mM)*(Na/mM)*(Na/mM)) : mole * meter**-3
    dATP/dt =  0*(mM/ms) : mole * meter**-3

'''
#z = (1*mM)/(1*mM+ATP/0.06) : 1
#milisemins* cm^2 * mV
#Iapp micro amps per cm^2

##Constants :
Cm = 1.0* ufarad * cm ** -2
gL=0.1 * msiemens * cm **-2
VL=-66*mV

gNa=100 * msiemens * cm **-2
VNa=50*mV

gK=80 * msiemens * cm **-2
VK=-100*mV

gKATP=0.56 *msiemens*cm**-2 #0*msiemens*cm**-2 #

ransc=0*(amp*(10**-6) * cm ** -2)

F = 8.8*(10**(-5)) #units?
K_m =  6*(10**(-8)) 


#for atp
#for atp in np.arange(3)+1:
n=100
atp=1
#Iapp =  4*amp*(10**-6) * cm ** -2 #*amp*meter**-2 #
#atp = 1.5883100381194408
PY=NeuronGroup(1,eq_PY,method='rk4',threshold='V>20*mV', reset='V=0*mV' )
PY.V = '-65*mvolt'
PY.h = '0.9932525096349422'
PY.m = '0.022083418537191028'
PY.n = '0.05182107474228339'
PY.Na = '((2*0.0004)*(2-atp)/(atp*2*0.00000006))**(1.0/3)*(mM)'
PY.ATP = 'atp*(mM)'
PY.Iapp = 'np.arange(n)*amp*(10**-6) * cm ** -2'

#ATP = 2*2*0.0004/(2*0.0004+2*0.00000006*Na*Na*Na)
   # Na = ((2*0.0004)*(2-ATP)/(ATP*2*0.00000006))**(1/3)

#PY.h = '0+0.05*rand()'
#PY.m = '0+0.05*rand()'
#PY.n = '0+0.05*rand()'
#PY.m = '0+0.05*rand()'

V1=StateMonitor(PY,'V',record=[0])
event_mon = SpikeMonitor(PY) #see if I can record variables
#    
run(0.5*second)

#spike_dict = event_mon.spike_trains()
#spike=0
#neur=0
##find the first spike in the neuron group
#while spike ==0:
#    #which variable actually records the spikes???
#    x=spike_dict[neur]/ms
#    
#    neur=neur+1
    

figure()
plot(V1.t/second,V1.V[0]/volt)
xlabel('Time (s)')
ylabel('Membrane potential (V)')
title('PY cell')
    
    #    figure()
    #    plot(I1.t/second,I1.IL[0],label='L')
    #    plot(I1.t/second,I2.INa[0],label='Na')
    #    plot(I1.t/second,I3.IK[0],label='K')
    #    plot(I1.t/second,I4.IAR[0],label='AR')
    #    title('Synaptic currents')
#    legend()