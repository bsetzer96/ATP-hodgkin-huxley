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

#range of atp and iapp values
atpval=(np.arange(200)+1)*0.01#0-2
iappval=np.arange(150)*0.1#0-15

#atp = 1.5883100381194408
duration=0.2*second

atp=0.1
Iapp = 1.2*amp*(10**-6) * cm ** -2 #*amp*meter**-2 #

PY=NeuronGroup(1,eq_PY,method='rk4',threshold='V>20*mV')
PY.V = '-65*mvolt'
PY.h = '0.9932525096349422'
PY.m = '0.022083418537191028'
PY.n = '0.05182107474228339'
PY.Na = '((2*0.0004)*(2-atp)/(atp*2*0.00000006))**(1.0/3)*(mM)'
PY.ATP = 'atp*(mM)'
#PY.Iapp = 'iapp*amp*(10**-6) * cm ** -2'

#ATP = 2*2*0.0004/(2*0.0004+2*0.00000006*Na*Na*Na)
   # Na = ((2*0.0004)*(2-ATP)/(ATP*2*0.00000006))**(1/3)
#PY.h = '0+0.05*rand()'
#PY.m = '0+0.05*rand()'
#PY.n = '0+0.05*rand()'
#PY.m = '0+0.05*rand()'
V1=StateMonitor(PY,'V',record=[0])
event_mon = SpikeMonitor(PY) #see if I can record variables

net = Network(PY, V1, event_mon) # create network object so we don't include objects later
#Net = Network(NG, monitor) #include neuron variables and monitors, stores net variable


net.run(duration) 
net.store()
bndry=np.zeros([np.size(atpval),2])
p=-1
#for each ATP loop through Iapp values
for atp in atpval:
    p=p+1
    spk=0
    k=0
    #while Iapp isn't high enough to spike 
    while spk==0 and k<150 :
        #changing Iapp on each run
        net.restore()
        Iapp =  iappval[k]*amp*(10**-6) * cm ** -2 
        #rerun
        net.run(duration)
        #see if there were any spikes
        spike_dict = event_mon.spike_trains()
#        print(atp)
#        print(Iapp)
#        print(spike_dict[0])
        if np.size(spike_dict[0])>0:
            #if there were spikes, record value and move on
            bndry[p,:]= [atp, Iapp]
            spk=1
        
        k=k+1
        #if it doesnt spike for any Iapp record at high value
        if k==150:
            bndry[p,:]= [atp, 100]
            

figure()
plot(bndry[:,1], bndry[:,0])
xlabel('Iapp')
ylabel('ATP')
    

#figure()
#plot(V1.t/second,V1.V[0]/volt)
#xlabel('Time (s)')
#ylabel('Membrane potential (V)')
#title('PY cell')
    
    
    
#net = Network(P_stn, MP_stn, P_gpe, MP_gpe, P_dbs, MP_dbs) # create network object so we don't include objects later
#net.run(duration) # generate spikes
#Net = Network(NG, monitor) #include neuron variables and monitors, stores net variable
#net.store()
#net.run()
    #can run a for loop seperately
    #for parameter
#
#for n_str in str_freqs:
#    # restore network
#    net.restore()
    #str_freqs = [0, 10, 20, 30, 40, 50, 60] * Hz
     
    #run duration again
