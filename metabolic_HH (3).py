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
t2=150           #end time
delta=0.05      #step size
h=delta         #runge kutta step
T= (t2-t1)/delta+1 #length of time vector
t=np.linspace(int(t1),int(t2),int(T)) #time vector


## paramenters
y0_py=np.array([-70, 0, 0, 0, 0, 2])    #initial condition
y0_fs= np.array([-70, 0, 0, 0])
y0_x = [0,0]
n=1                   #number of pyramidal cells
m=1                 #number of fs cells

E_Na= 50            #Sodium equilibrium potential
E_k= -100            #potassium equilibrium potential
E_AMPA= 0           #AMPA EP
E_GABA= -80       #GABA EP
g_Na= 100             #Sodium conductance
g_K= 80               #potassium conductance
g_K_ATP = 0.15        #potassium ATP conductance
g_AMPA= 2               #AMPA excitatory connection current max
g_GABA= 0.64        #GABA inhibitory connection current max
I_app_py= 1.8       #applied current to pyramidal cells
I_app_fs=.5          #applied current to fs cells
J_ATP= 2              #production rate of ATP (linked to metabolism)
ATP_max= 0.6          # ?
K_m= 6*(10**(-8))     # govern the NA-ATP pump dynamics
F= 8.8*(10**(-5))     # help govern the NA-ATP pump dynamics
tau_GABA= 5

p=np.array([E_Na, E_k, E_AMPA, E_GABA, g_Na, g_K, g_K_ATP, g_AMPA, g_GABA, I_app_py, I_app_fs, J_ATP, ATP_max, K_m, F, tau_GABA])
p=np.transpose(p)


##preallocate space
y=np.zeros((6*n+4*m+2,int(T)))
y[:6,0]=y0_py
y[6:10,0]=y0_fs
y[10:,0]=y0_x



## Hodgkin huxley network of neurons with ATP metabolism
def hodghuxATP(state, N,p):
    #N = number of neurons in network
    state_py=state[:6]
    state_fs=state[6:10]
    state_con=state[10:]
    
    v_py= state_py[0]
    m_py= state_py[1]
    n_py= state_py[2]
    h_py= state_py[3]
    Na= state_py[4]
    ATP= state_py[5]
    
    v_fs= state_fs[0]
    m_fs= state_fs[1]
    n_fs= state_fs[2]
    h_fs= state_fs[3]
    
    x_AMPA = state_con[0]
    x_GABA = state_con[1]
    
    E_Na= p[0]
    E_K= p[1]
    E_AMPA = p[2]
    E_GABA= p[3]
    g_Na= p[4]
    g_K = p[5]
    g_K_ATP = p[6]
    g_AMPA = p[7]
    g_GABA = p[8]
    I_app_py= p[9]
    I_app_fs = p[10]
    J_ATP = p[11]
    ATP_max = p[12]
    K_m = p[13]
    F = p[14]
    tau_GABA = p[15] 
    
    
    #hm
    alpha_m_py = ((0.32)*(v_py+54))/(1-math.exp(-(v_py+54)/4))
    beta_m_py = ((0.28)*(v_py+27))/(math.exp((v_py+27)/5)-1)
    alpha_h_py = (0.128)*math.exp(-(v_py+50)/18)
    beta_h_py = 4/(1+math.exp(-(v_py+27)/5))
    alpha_n_py = ((0.032)*(v_py+52))/(1-math.exp(-(v_py+52)/5))
    beta_n_py = 0.5*math.exp(-(v_py+57)/40)
    z = 1/(1+6*ATP)
    
    alpha_m_fs = ((0.32)*(v_fs+54))/(1-math.exp(-(v_fs+54)/4))
    beta_m_fs = ((0.28)*(v_fs+27))/(math.exp((v_fs+27)/5)-1)
    alpha_h_fs = (0.128)*math.exp(-(v_fs+50)/18)
    beta_h_fs = 4/(1+math.exp(-(v_fs+27)/5))
    alpha_n_fs = ((0.032)*(v_fs+52))/(1-math.exp(-(v_fs+52)/5))
    beta_n_fs = 0.5*math.exp(-(v_fs+57)/40)
    
    #gating fuctions 
    mdot_py = alpha_m_py*(1-m_py) - beta_m_py*m_py
    hdot_py= alpha_h_py*(1-h_py) - beta_h_py*h_py
    ndot_py = alpha_n_py*(1-n_py)-beta_n_py*n_py
    
    mdot_fs = alpha_m_fs*(1-m_fs) - beta_m_fs*m_fs
    hdot_fs= alpha_h_fs*(1-h_fs) - beta_h_fs*h_fs
    ndot_fs = alpha_n_fs*(1-n_fs)-beta_n_fs*n_fs
    
    #connectivity functions
    I_AMPA = g_AMPA*x_AMPA*(v_fs-E_AMPA)
    I_GABA=g_GABA*x_GABA*(v_py-E_GABA)

    #pyramidal cell currents
    I_Na_py = g_Na*(m_py**3)*h_py*(v_py-E_Na)
    I_K_py = g_K*(n_py**4)*(v_py-E_K)
    I_leak_py = 0.1*(v_py+61)
    I_K_ATP = g_K_ATP*z*(v_py-E_K)
    
    #FS cell currents
    I_Na_fs = g_Na*(m_fs**3)*h_fs*(v_fs-E_Na)
    I_K_fs = g_K*(n_fs**4)*(v_fs-E_K)
    I_leak_fs = 0.1*(v_fs+61)
 
    #pyramidal cell ATP dif eqs
    Nadot = F*I_Na_py-3*K_m*(Na**3)*ATP
    ATPdot = J_ATP*(ATP_max - ATP) - K_m*(Na**3)*ATP
    
    #connection gating 
    xdot_AMPA = 5*(1+math.tanh(v_py/4))*(1-x_AMPA)-x_AMPA/2
    xdot_GABA = 2*(1+math.tanh(v_fs/4))*(1-x_GABA)-x_GABA/tau_GABA
    
    #voltage equations for both cells
    vdot_py = I_app_py-I_Na_py-I_K_py-I_K_ATP-I_leak_py-I_GABA
    vdot_fs = I_app_fs-I_Na_fs-I_K_fs-I_leak_fs -I_AMPA
    
    #state variable arrays 
    dx_py= np.array([vdot_py, mdot_py, ndot_py, hdot_py, Nadot, ATPdot])
    dx_fs= np.array([vdot_fs, mdot_fs, ndot_fs, hdot_fs])
    dx_con=np.array([xdot_AMPA, xdot_GABA])
    dx= np.concatenate([dx_py, dx_fs, dx_con])
    return dx

##Runge kutta order 4 integration of pyramidal cell
for n in np.arange(T-1):
    k1= h*hodghuxATP(y[:,int(n)],n,p)
    k2=h*hodghuxATP(y[:,int(n)]+k1,n,p)
    k3=h*hodghuxATP(y[:,int(n)]+k2/2,n,p)
    k4=h*hodghuxATP(y[:,int(n)],n,p)
    y[:,int(n+1)]=y[:,int(n)]+(1/6)*(k1+2*k2+2*k3+k4)

v_py=y[0,:]
v_fs=y[6,:]

plt.plot(t,v_py)
plt.plot(t,v_fs)
plt.legend(['Pyramidal Cell', 'FS Cell'])
plt.show 



## Q: Can we go over the set up of the code.... what different functions should I have? 
#What is done seperately what can be done together?