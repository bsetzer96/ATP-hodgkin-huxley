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
    dATP/dt =  (0.0008/ms)*(2.00*mM-ATP)-(K_m*ATP*(Na/mM)*(Na/mM)*(Na/mM))/ms : mole * meter**-3
'''
#z = (1*mM)/(1*mM+ATP/0.06) : 1
#milisemins* cm^2 * mV
#Iapp micro amps per cm^2
#(mmole/cm**3)

##Constants :
Cm = 1.0* ufarad * cm ** -2
gL=0.1 * msiemens * cm **-2
VL=-66*mV

gNa=100 * msiemens * cm **-2
VNa=50*mV

gK=80 * msiemens * cm **-2
VK=-100*mV

gKATP=0.56 *msiemens*cm**-2 #0*msiemens*cm**-2 #

Iapp = 1.5*(amp*(10**-6) * cm ** -2) #*amp*meter**-2 #

ransc=0.5*(amp*(10**-6) * cm ** -2)

F = 8.8*(10**(-5)) #units?
K_m =  6*(10**(-8)) 



if __name__=='__main__' :
#    start_scope()
#    Vrev_inp=0*mV
#    taurinp=0.1*ms # don't think i need most of this
#    taudinp=0.5*ms
#    tauinp=taudinp
#    Vhigh=0*mV
#    Vlow=-80*mV
#    ginp_IB=0* msiemens * cm **-2
#    ginp=0* msiemens * cm **-2
        
    PY=NeuronGroup(1,eq_PY,method='rk4')
    PY.V = '-65*mvolt'
    PY.h = '0.9932525096349422'
    PY.m = '0.022083418537191028'
    PY.n = '0.05182107474228339'
    PY.Na = '12*(mM)'
    PY.ATP = '1.5883100381194408*(mM)'
    
    
    
    #PY.h = '0+0.05*rand()'
    #PY.m = '0+0.05*rand()'
    #PY.n = '0+0.05*rand()'
    #PY.m = '0+0.05*rand()'
    
    V1=StateMonitor(PY,'V',record=True)
    Na=StateMonitor(PY, 'Na' , record=True)
    ATP=StateMonitor(PY, 'ATP' , record=True)
    #ATP1=StateMonitor(PY, 'ATP', record=True)
   # R1=SpikeMonitor(PY,record=True)
    INa=StateMonitor(PY,'INa',record=[0])
    IK=StateMonitor(PY,'IK',record=[0])
    IKATP=StateMonitor(PY, 'IKATP', record=[0])
    Iran=StateMonitor(PY, 'Iran', record=[0])
    
    all_neurons= PY
    #all_synapses= S_PYPY
    all_monitors= V1, ATP,Na, INa, IK, IKATP,Iran
    
    
    run(0.5*second)
    
    figure()
    plot(V1.t/second,V1.V[0]/mV)
    xlabel('Time (s)')
    ylabel('Membrane potential (mv)')
    title('PY cell')
    figure()
    atp=all_monitors[1]
    plot(atp.t/second,atp.ATP[0]/mM)
    title('ATP')  
    xlabel('Time (s)')
    ylabel('ATP concentration (mM)')
    figure()
    na=all_monitors[2]
    plot(na.t/second,na.Na[0]/mM)
    title('Na')  
    xlabel('Time (s)')
    ylabel('Na concentration (mM)')
    figure()
    ina=all_monitors[3]
    plot(ina.t/second,ina.INa[0]/(amp*(10**-6) * cm**-2))
    title('INa')  
    xlabel('Time (s)')
    ylabel('INa current (amp*(10^-6) * cm^-2)')
    figure()
    ik=all_monitors[4]
    plot(ik.t/second,ik.IK[0]/(amp*(10**-6) * cm**-2))
    xlabel('Time (s)')
    ylabel('IK current (amp*(10^-6) * cm^-2)')
    title('IK')  
    figure()
    ikatp=all_monitors[5]
    plot(ik.t/second,ikatp.IKATP[0]/(amp*(10**-6) * cm**-2))
    xlabel('Time (s)')
    ylabel('IKATP current (amp*(10^-6) * cm^-2)')
    title('IKATP')  
    figure()
    iran=all_monitors[6]
    plot(ik.t/second,iran.Iran[0]/(amp*(10**-6) * cm**-2))
    xlabel('Time (s)')
    ylabel('Iran current (amp*(10^-6) * cm^-2)')
    title('Iran')  
    clear_cache('cython')
    #    figure()
    #    plot(I1.t/second,I1.IL[0],label='L')
    #    plot(I1.t/second,I2.INa[0],label='Na')
    #    plot(I1.t/second,I3.IK[0],label='K')
    #    plot(I1.t/second,I4.IAR[0],label='AR')
    #    title('Synaptic currents')
#    legend()