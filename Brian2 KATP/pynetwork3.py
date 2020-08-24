# -*- coding: utf-8 -*-
"""
Created on Thu May  7 12:59:15 2020

@author: bever
"""


from brian2 import *

start_scope()

from scipy import signal

#from FS_brian2 import *
from PY_brian2_working import *


prefs.codegen.target = 'numpy'
prefs.codegen.target = 'cython'


defaultclock.dt = 0.01*ms

def generate_network(N_PY):
    
    PY=NeuronGroup(N_PY,eq_PY,method='rk4')
    PY.V = '(-65)*mvolt' #+0.05*rand()
    PY.h = '0.9932525096349422'
    PY.m = '0.022083418537191028'
    PY.n = '0.05182107474228339'
    PY.Na = '(12)*mM'
    PY.ATP = '(1.5883100381194408)*mM'
    #mole*meter**-3
#    PY1 = PY[0]
#    PY2 = PY[1]
#    PY1.Iapp = '5 *amp*(10**-6) * cm ** -2'
#    PY2.Iapp = '5 *amp*(10**-6) * cm ** -2'
    
 #'''_post=0 : amp * meter ** -2 (summed)'''   
     ##Synapses
     
    #xdot_AMPA = 5*(1+np.tanh(v_py/4))*(1-x_AMPA)-x_AMPA/2 #py cells
    #xdot_GABA = 2*(1+np.tanh(v_fs/4))*(1-x_GABA)-x_GABA/tau_GABA  #fs cells
    eq_syn=  '''_post=(g_i/a_i)*s_i*(V_post-E_i) : amp * meter ** -2 (summed)
        ds_i/dt=b_i*(1+tanh((V_pre/4)/mV))*(1-s_i)/ms-(s_i/tau_i) : 1
        g_i : siemens* meter**-2
        V_i : volt
        b_i : 1
        E_i : volt
        tau_i : second
        a_i : 1
    '''
    #is v_pre/4 too small?
    
    S_PYPY=Synapses(PY,PY,model='IsynPY'+eq_syn,method='rk4')
    S_PYPY.connect(i=[0,1,2], j=[1,2,0]) #j='i' #i=0, j=1
    S_PYPY.g_i=(0.3/(N_PY-1))* msiemens * cm **-2
    S_PYPY.tau_i=2*ms
    S_PYPY.E_i=0*mV
    S_PYPY.a_i= 1
    S_PYPY.b_i = 5 #maybe b_i is too high
    S_PYPY.s_i=0.0001
    
    V1=StateMonitor(PY,'V',record=True)
    Na=StateMonitor(PY, 'Na' , record=True)
    ATP=StateMonitor(PY, 'ATP' , record=True)
    #ATP1=StateMonitor(PY, 'ATP', record=True)
    AMPA=StateMonitor(S_PYPY, 'IsynPY', record=True)
    S= StateMonitor(S_PYPY, 's_i', record=True)
   # R1=SpikeMonitor(PY,record=True)
    INa=StateMonitor(PY,'INa',record=[0])
    IK=StateMonitor(PY,'IK',record=[0])
    
    all_neurons= PY
    all_synapses= S_PYPY
    all_monitors= V1, ATP,Na, INa, IK, S, AMPA
    #am=AMPA

    return all_neurons,all_synapses,all_monitors  

#j='i' 
#look up how to plot connectivity and save
    
if __name__=='__main__':
    pynet=Network()
    N_PY=3
    [all_neurons, all_synapses, all_monitors]= generate_network(N_PY)
    pynet.add(all_neurons)
    pynet.add(all_synapses)
    pynet.add(all_monitors)
    pynet.run(20*second, report='text',report_period=60*second)
    
    figure()
    V1=all_monitors[0]
    plot(V1.t/second,V1.V[0]/mV)
    plot(V1.t/second, V1.V[1]/mV)
    plot(V1.t/second, V1.V[2]/mV)
    xlabel('Time (s)')
    ylabel('Membrane potential (mV)')
    title('PY cell')
    figure()
    am=all_monitors[6]
    plot(am.t/second,am.IsynPY[1]/second)
    xlabel('Time (s)')
    ylabel('')
    title('AMPA')
    figure()
    s=all_monitors[5]
    plot(am.t/second,s.s_i[0]/second)
    xlabel('Time (s)')
    ylabel('')
    title('s_i')   
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
    clear_cache('cython')
    
    # How do i control connections?
    # how do i set up random connections
    # we don't know what properties of the randomness are important 
    
    
    #make sure ampa current is working, make sure spike rate going up
    #check post-synaptic ampa current , does it look like an AMPA current
    #set one cell to Iapp=0 to see if it's getting ampa
    #check with Amelie about speed
    #could use timestep of 0.05
    
    
    #Q: does the 'threshold' make it spike, or just help to record spikes
    #why does it break when I add Iapp into equation
    
    #Check units of currents
    #ask about looping through for iapp/katp curve
    
    
    
    