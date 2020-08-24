# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 15:52:46 2020

@author: bever
"""


from brian2 import *

start_scope()

from scipy import signal
from FS_brian2 import *
from PY_brian2 import *


prefs.codegen.target = 'numpy'
prefs.codegen.target = 'cython'

defaultclock.dt = 0.01*ms

def generate_network(N_PY,N_FS):
    FS=NeuronGroup(N_FS,eq_FS,method='rk4')
    FS.V = '-110*mvolt+10*rand()*mvolt'
    FS.h = '0+0.05*rand()'
    FS.m = '0+0.05*rand()'
    FS.n = '0+0.05*rand()'
    
    PY=NeuronGroup(N_PY,eq_PY,method='rk4')
    PY.V = '-65*mvolt'
    PY.h = '0.9932525096349422'
    PY.m = '0.022083418537191028'
    PY.n = '0.05182107474228339'
    PY.Na = '12*(mole * meter**-3)'
    PY.ATP = '1.5883100381194408*(mole * meter**-3)'
    PY.Iapp = '5 *amp*(10**-6) * cm ** -2'
    
    #'''_post=0 : amp * meter ** -2 (summed)'''   
     ##Synapses
    eq_syn=  '''_post=(g_i/a_i)*s_i*(V_post-E_i) : amp * meter ** -2 (summed)
        ds_i/dt=b_i*(1+tanh((V_pre/4)/volt))*(1-s_i)/ms-(s_i/tau_i) : 1
        g_i : siemens* meter**-2
        V_i : volt
        b_i : 1
        E_i : volt
        tau_i : second
        a_i : 1
    '''
    
    S_PYPY=Synapses(PY,PY,model='IsynPY'+eq_syn,method='exact')
    S_PYPY.connect() 
    S_PYPY.g_i=(0.1/(N_PY-1))* msiemens * cm **-2
    S_PYPY.tau_i=2*ms
    S_PYPY.E_i=0*mV
    S_PYPY.a_i= 1
    S_PYPY.b_i = 5
    
    S_FSPY=Synapses(FS,PY,model='IsynFS'+eq_syn,method='exact')
    S_FSPY.connect()
    S_FSPY.g_i=(0.64/(N_FS-1)) * msiemens * cm **-2
    #S_RSSI.g_i=0* msiemens * cm **-2
    S_FSPY.b_i=2
    S_FSPY.tau_i=5*ms
    S_FSPY.a_i= 1 #4
    S_FSPY.E_i= -80*mV
    
    
    S_PYFS=Synapses(PY,FS,model='IsynPY'+eq_syn,method='exact')
    S_PYFS.connect()
    S_PYFS.g_i=(2/(N_PY-1))* msiemens * cm **-2
#    S_SIRS.g_i=0.1* msiemens * cm **-2
    S_PYFS.tau_i=2*ms
    S_PYFS.a_i= 1 #9
    S_PYFS.b_i=5
    S_PYFS.E_i= 0*mV
    
    S_PYPY=Synapses(PY,PY,model='IsynPY'+eq_syn,method='exact')
    S_PYPY.connect()
    S_PYPY.g_i=(0.1/(N_PY-1))* msiemens * cm **-2
    #S_SISI.g_i=5* msiemens * cm **-2
    S_PYPY.tau_i=2*ms
    S_PYPY.E_i=0*mV
    S_PYPY.a_i= 1
    S_PYPY.b_i = 5
    
    S_FSFS=Synapses(FS,FS,model='IsynFS'+eq_syn,method='exact')
    S_FSFS.connect() #j= 'i' is one-to-one connections, () all to all #p = 0.2
    S_FSFS.g_i=1/(N_FS-1)* msiemens * cm **-2
    #S_SISI.g_i=5* msiemens * cm **-2
    S_FSFS.tau_i=5*ms
    S_FSFS.E_i=-80*mV  
    S_FSFS.a_i= 1 #3
    S_FSFS.b_i = 2
    
    V1=StateMonitor(PY,'V',record=True)
    V2=StateMonitor(FS,'V',record=True)
    
    all_neurons=PY, FS
    all_synapses= S_FSPY, S_PYFS, S_PYPY, S_FSFS
    all_monitors=V1,V2

    return all_neurons,all_synapses,all_monitors    

#j='i' 
#look up how to plot connectivity and save
    
if __name__=='__main__':
    net = Network()
    N_PY, N_FS=10,4
    [all_neurons, all_synapses, all_monitors]= generate_network(N_PY,N_FS)
    net.add(all_neurons)
    net.add(all_synapses)
    net.add(all_monitors)
    net.run(1*second, report='text',report_period=300*second)
    figure()
    V1=all_monitors[0]
    V2=all_monitors[1]
    plot(V1.t/second,V1.V[0]/volt)
    xlabel('Time (s)')
    ylabel('Membrane potential (V)')
    title('PY cell')
    
    
    # How do i control connections?
    # how do i set up random connections
    # we don't know what properties of the randomness are important 
    
    
    
    