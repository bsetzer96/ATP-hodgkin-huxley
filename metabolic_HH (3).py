# -*- coding: utf-8 -*-
"""
Spyder Editor
neurophysiological model for burst suppression
This is a temporary script file.
"""

import numpy as np
import matplotlib.pyplot as plt

## basic values for integration -------------------------------
t1=0            #start time
t2=150           #end time
delta=0.05      #step size
h=delta         #runge kutta step
T= (t2-t1)/delta+1 #length of time vector
t=np.linspace(int(t1),int(t2),int(T)) #time vector
n= 10               #number of pyramidal cells
m=4                 #number of fs cells

## parameters   ------------------------------------------------------
y0_py=np.zeros(6*n)    #initial condition
y0_py[:n]=-70 #first n are voltage 
y0_py[n*5:]=2 #last n are ATP
y0_fs= np.zeros(4*m)
y0_fs[:m]=-70 
y0_x = np.zeros(n+m)


E_Na= 50            #Sodium equilibrium potential
E_k= -100            #potassium equilibrium potential
E_AMPA= 0           #AMPA EP
E_GABA= -80       #GABA EP
g_Na= 100             #Sodium conductance
g_K= 80               #potassium conductance
g_K_ATP = 0.15        #potassium ATP conductance
g_AMPA_py= 0.1/(n-1)               #AMPA excitatory connection current max
g_AMPA_fs=2/(n-1)
g_GABA_py= 0.64/(m-1)       #GABA inhibitory connection current max
g_GABA_fs=1/(m-1)
I_app_py= 1.8       #applied current to pyramidal cells
I_app_fs=.5          #applied current to fs cells
J_ATP= 2              #production rate of ATP (linked to metabolism)
ATP_max= 0.6          # ?
K_m= 6*(10**(-8))     # govern the NA-ATP pump dynamics
F= 8.8*(10**(-5))     # help govern the NA-ATP pump dynamics
tau_GABA= 5

#Adding noise to parameters
p=np.array([E_Na, E_k, E_AMPA, E_GABA, g_Na, g_K, g_K_ATP, g_AMPA_py, g_AMPA_fs, \
            g_GABA_py, g_GABA_fs, I_app_py, I_app_fs, J_ATP, ATP_max, K_m, F, tau_GABA])
p=np.transpose(p)


##preallocate space
y=np.zeros((7*n+5*m,int(T)))
y[:6*n,0]=y0_py
y[6*n:(6*n+4*m),0]=y0_fs
y[(6*n+4*m):,0]=y0_x



## Hodgkin huxley network of neurons with ATP metabolism---------------------
def hodghuxATP(state, n, m,p):
    #state = vector of all changing variables!
    #       [v_py--m_py--n_py--h_py--Na--ATP | v_fs m_fs n_fs h_fs | x_AMPA x_GABA]
    #n = number of pyramidal neurons
    #m = number of fast spiking interneurons
    #p = vector of all parameters 
    #       [E_Na, E_k, E_AMPA, E_GABA, g_Na, g_K, g_K_ATP, g_AMPA_py, g_AMPA_fs, g_GABA_py, g_GABA_fs, I_app_py, I_app_fs, J_ATP, ATP_max, K_m, F, tau_GABA]
   
    #assigning variables and parameters- - -
    state_py=state[:6*n]
    state_fs=state[6*n:(6*n+4*m)]
    state_con=state[(6*n+4*m):]
    
    v_py= state_py[:n]
    m_py= state_py[n:n*2]
    n_py= state_py[n*2:n*3]
    h_py= state_py[n*3:n*4]
    Na= state_py[n*4:n*5]
    ATP= state_py[n*5:]
    
    v_fs= state_fs[:m]
    m_fs= state_fs[m:m*2]
    n_fs= state_fs[m*2:m*3]
    h_fs= state_fs[m*3:]
    
    x_AMPA = state_con[0:n]
    x_GABA = state_con[n:]
    
    E_Na= p[0]
    E_K= p[1]
    E_AMPA = p[2]
    E_GABA= p[3]
    g_Na= p[4]
    g_K = p[5]
    g_K_ATP = p[6]
    g_AMPA_py = p[7]
    g_AMPA_fs = p[8]
    g_GABA_py = p[9]
    g_GABA_fs= p[10]
    I_app_py= p[11]
    I_app_fs = p[12]
    J_ATP = p[13]
    ATP_max = p[14]
    K_m = p[15]
    F = p[16]
    tau_GABA = p[17] 
    
    
    #hudgkin-huxley equations - - -
    #
    alpha_m_py = ((0.32)*(v_py+54))/(1-np.exp(-(v_py+54)/4))
    beta_m_py = ((0.28)*(v_py+27))/(np.exp((v_py+27)/5)-1)
    alpha_h_py = (0.128)*np.exp(-(v_py+50)/18)
    beta_h_py = 4/(1+np.exp(-(v_py+27)/5))
    alpha_n_py = ((0.032)*(v_py+52))/(1-np.exp(-(v_py+52)/5))
    beta_n_py = 0.5*np.exp(-(v_py+57)/40)
    z = 1/(1+6*ATP)
    
    alpha_m_fs = ((0.32)*(v_fs+54))/(1-np.exp(-(v_fs+54)/4))
    beta_m_fs = ((0.28)*(v_fs+27))/(np.exp((v_fs+27)/5)-1)
    alpha_h_fs = (0.128)*np.exp(-(v_fs+50)/18)
    beta_h_fs = 4/(1+np.exp(-(v_fs+27)/5))
    alpha_n_fs = ((0.032)*(v_fs+52))/(1-np.exp(-(v_fs+52)/5))
    beta_n_fs = 0.5*np.exp(-(v_fs+57)/40)
    
    #gating fuctions 
    mdot_py = alpha_m_py*(1-m_py) - beta_m_py*m_py
    hdot_py= alpha_h_py*(1-h_py) - beta_h_py*h_py
    ndot_py = alpha_n_py*(1-n_py)-beta_n_py*n_py
    
    mdot_fs = alpha_m_fs*(1-m_fs) - beta_m_fs*m_fs
    hdot_fs= alpha_h_fs*(1-h_fs) - beta_h_fs*h_fs
    ndot_fs = alpha_n_fs*(1-n_fs)-beta_n_fs*n_fs
    
    #Connectivity Currents - - - - - - - - - 
    I_AMPA_py=np.zeros((n,n))    #preallocating space
    I_AMPA_fs=np.zeros((n,m))
    I_GABA_py=np.zeros((m,n))
    I_GABA_fs=np.zeros((m,m))
    
    #calculating input currents between all connections
    for j in np.arange(n):#post synaptic neuron recieving input (to)
        for i in np.arange(n): #presynaptic neuron giving input (from)
            #PY to PY
            I_AMPA_py[i,j]= g_AMPA_py*x_AMPA[i]*(v_py[j]-E_AMPA) 
        for i in np.arange(m): #from
            # FS to PY
            I_GABA_py[i,j]=g_GABA_py*x_GABA[i]*(v_py[j]-E_GABA)           
    for j in np.arange(m): #to
        for i in np.arange(n): #from i-->j
            #PY to FS
            I_AMPA_fs[i,j]= g_AMPA_fs*x_AMPA[i]*(v_fs[j]-E_AMPA) 
        for i in np.arange(m): #from
            # FS to FS
            I_GABA_fs[i,j]=g_GABA_fs*x_GABA[i]*(v_fs[j]-E_GABA) 
    
    #take out connections to self
    I_AMPA_py[np.eye(n)>0]=0
    I_GABA_fs[np.eye(m)>0]=0
    
    #summing total input from all other neurons
    I_AMPA_py=I_AMPA_py.sum(axis=0)
    I_AMPA_fs=I_AMPA_fs.sum(axis=0)
    I_GABA_py=I_GABA_py.sum(axis=0)
    I_GABA_fs=I_GABA_fs.sum(axis=0)
    
    #connection gating differentials
    xdot_AMPA = 5*(1+np.tanh(v_py/4))*(1-x_AMPA)-x_AMPA/2 #py cells
    xdot_GABA = 2*(1+np.tanh(v_fs/4))*(1-x_GABA)-x_GABA/tau_GABA  #fs cells
    
    #Other Currents - - - - -  - - - -
    #pyramidal cell currents - - -
    I_Na_py = g_Na*(m_py**3)*h_py*(v_py-E_Na) #sodium current 
    I_K_py = g_K*(n_py**4)*(v_py-E_K) #potassium current 
    I_leak_py = 0.1*(v_py+61) #leak current 
    I_K_ATP = g_K_ATP*z*(v_py-E_K) #potassium ATP current (metabolism)
    #pyramidal cell ATP dif eqs
    Nadot = F*I_Na_py-3*K_m*(Na**3)*ATP
    ATPdot = J_ATP*(ATP_max - ATP) - K_m*(Na**3)*ATP 
    
    #FS cell currents - - -
    I_Na_fs = g_Na*(m_fs**3)*h_fs*(v_fs-E_Na)
    I_K_fs = g_K*(n_fs**4)*(v_fs-E_K)
    I_leak_fs = 0.1*(v_fs+61)
 

    
    #voltage equations for both cells - - -
    vdot_py = I_app_py-I_Na_py-I_K_py-I_K_ATP-I_leak_py-I_GABA_py-I_AMPA_py
    vdot_fs = I_app_fs-I_Na_fs-I_K_fs-I_leak_fs -I_AMPA_fs-I_GABA_fs
    
    #state variable arrays - - - 
    dx_py= np.concatenate([vdot_py, mdot_py, ndot_py, hdot_py, Nadot, ATPdot])
    dx_fs= np.concatenate([vdot_fs, mdot_fs, ndot_fs, hdot_fs])
    dx_con= np.concatenate([xdot_AMPA, xdot_GABA])
    dx= np.concatenate([dx_py, dx_fs, dx_con])
    return dx

##Runge kutta order 4 integration of pyramidal cell-----------------------------------
for i in np.arange(T-1):
    k1= h*hodghuxATP(y[:,int(i)],n,m,p) #first step
    k2=h*hodghuxATP(y[:,int(i)]+k1,n,m,p) #second step
    k3=h*hodghuxATP(y[:,int(i)]+k2/2,n,m,p) #third step
    k4=h*hodghuxATP(y[:,int(i)],n,m,p) #fourth step
    y[:,int(i+1)]=y[:,int(i)]+(1/6)*(k1+2*k2+2*k3+k4) #RK order 4

y_py=y[0:6*n,:] #extracting py cell states
y_fs=y[6*n:(6*n+4*m),:] #extracting fs cell states

v_py=y_py[:n,:] #pyramidal cell membarane potentials
v_fs=y_fs[:m,:] #fast spiking interneuron cell membrane potentials


## plotting ---------------------------
plt.figure
for i in np.arange(n):
    plt.plot(t,v_py[i])
for i in np.arange(m):
    plt.plot(t,v_fs[i])
    plt.legend(['FS Cell '+ str(i+1)])
plt.show 

