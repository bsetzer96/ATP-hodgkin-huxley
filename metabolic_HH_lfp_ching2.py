# -*- coding: utf-8 -*-
"""
Spyder Editor
neurophysiological model for burst suppression
This is a temporary script file.
"""

import numpy as np
import matplotlib.pyplot as plt

#check constants
#look at time point right before it blows up
#which variable blows up first
#mess with the time step

## basic values for integration -------------------------------
t1=0            #start time
t2=800           #end time
delta=0.05      #step size
h=delta         #runge kutta step
T= (t2-t1)/delta+1 #length of time vector
t=np.linspace(int(t1),int(t2),int(T)) #time vector
n= 10               #number of pyramidal cells
m= 4                 #number of fs cells
ATP_scale=.3

## parameters   ------------------------------------------------------
y0_py=np.zeros(6*n)    #initial condition
y0_py[:n]=-70 #first n are voltage 
y0_py[n*4:n*5] = 10 #Na 0
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
J_ATP= 2*ATP_scale              #production rate of ATP (linked to metabolism)
ATP_max= 1.5         # ?
K_m= 6*(10**(-8))     # govern the NA-ATP pump dynamics
F= 2*3*.000168*1.8 #8.8*(10**(-5))     # help govern the NA-ATP pump dynamics
tau_GABA= 5

#Adding noise to parameters
E_Na=E_Na*np.random.uniform(.95,1.05,n+m)
E_k=E_k*np.random.uniform(.95,1.05,n+m)
E_AMPA_py=E_AMPA*np.random.uniform(.95,1.05,n)
E_GABA_py=E_GABA*np.random.uniform(.95,1.05,n)
E_AMPA_fs=E_AMPA*np.random.uniform(.95,1.05,m)
E_GABA_fs=E_GABA*np.random.uniform(.95,1.05,m)
g_Na=g_Na*np.random.uniform(.95,1.05,n+m)
g_K=g_K*np.random.uniform(.95,1.05,n+m)
g_K_ATP=g_K_ATP*np.random.uniform(.95,1.05,n)
g_AMPA_py=g_AMPA_py*np.random.uniform(.95,1.05,n)
g_AMPA_fs=g_AMPA_fs*np.random.uniform(.95,1.05,m)
g_GABA_py=g_GABA_py*np.random.uniform(.95,1.05,n)
g_GABA_fs=g_GABA_fs*np.random.uniform(.95,1.05,m)


I_app_py=I_app_py +np.random.normal(0, 0.1, n)
#I_app_py=np.array([0,0])
I_app_fs=I_app_fs+np.random.normal(0,.1, m)
J_ATP=J_ATP*np.ones(n) #should this also have added noise for each neuron?
ATP_max=ATP_max*np.ones(n)
K_m=K_m*np.ones(n)
F=F*np.ones(n)

#ATP_max=ATP_max*np.random.uniform(.95,1.05,n)
#K_m=K_m*np.random.uniform(.95,1.05,n)
#F=F*np.random.uniform(.95,1.05,n)
tau_GABA=tau_GABA*np.ones(1) #should this one also have added noise

p=np.concatenate([E_Na, E_k, E_GABA_py, E_GABA_fs, E_AMPA_py, E_AMPA_fs, g_Na, g_K, g_K_ATP, g_AMPA_py, g_AMPA_fs, \
            g_GABA_py, g_GABA_fs, I_app_py, I_app_fs, J_ATP, ATP_max, K_m, F, tau_GABA])
p=np.transpose(p)


##preallocate space
y=np.zeros((7*n+5*m,int(T)))
y[:6*n,0]=y0_py
y[6*n:(6*n+4*m),0]=y0_fs
y[(6*n+4*m):(7*n+5*m),0]=y0_x



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
    state_con=state[(6*n+4*m):(7*n+5*m)]
    
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
    
    E_Na_py= p[:n]
    E_Na_fs=p[n:(m+n)]
    E_K_py= p[(m+n):(m+2*n)]
    E_K_fs=p[(m+2*n):(2*m+2*n)]
    E_GABA_py = p[(2*m+2*n):(2*m+3*n)]
    E_GABA_fs= p[(2*m+3*n):(3*m+3*n)]
    E_AMPA_py=p[(3*m+3*n):(3*m+4*n)]
    E_AMPA_fs=p[(3*m+4*n):(4*m+4*n)]
    g_Na_py= p[(4*m+4*n):(4*m+5*n)]
    g_Na_fs=p[(4*m+5*n):(5*m+5*n)]
    g_K_py = p[(5*m+5*n):(5*m+6*n)]
    g_K_fs = p[(5*m+6*n):(6*m+6*n)]
    g_K_ATP = p[(6*m+6*n):(6*m+7*n)]
    g_AMPA_py = p[(6*m+7*n):(6*m+8*n)]
    g_AMPA_fs = p[(6*m+8*n):(7*m+8*n)]
    g_GABA_py = p[(7*m+8*n):(7*m+9*n)]
    g_GABA_fs= p[(7*m+9*n):(8*m+9*n)]
    I_app_py= p[(8*m+9*n):(8*m+10*n)]
    I_app_fs = p[(8*m+10*n):(9*m+10*n)]
    J_ATP = p[(9*m+10*n):(9*m+11*n)]
    ATP_max = p[(9*m+11*n):(9*m+12*n)]
    K_m = p[(9*m+12*n):(9*m+13*n)]
    F = p[(9*m+13*n):(9*m+14*n)]
    tau_GABA = p[(9*m+14*n):(9*m+15*n)]

    
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
            I_AMPA_py[i,j]= g_AMPA_py[i]*x_AMPA[i]*(v_py[j]-E_AMPA_py[i]) 
        for i in np.arange(m): #from
            # FS to PY
            I_GABA_py[i,j]=g_GABA_py[j]*x_GABA[i]*(v_py[j]-E_GABA_py[j])           
    for j in np.arange(m): #to
        for i in np.arange(n): #from i-->j
            #PY to FS
            I_AMPA_fs[i,j]= g_AMPA_fs[j]*x_AMPA[i]*(v_fs[j]-E_AMPA_fs[j]) 
        for i in np.arange(m): #from
            # FS to FS
            I_GABA_fs[i,j]=g_GABA_fs[j]*x_GABA[i]*(v_fs[j]-E_GABA_fs[j]) 
    
    #take out connections to self
    I_AMPA_py[np.eye(n)>0]=0
    I_GABA_fs[np.eye(m)>0]=0
    
    I_AMPA_test=I_AMPA_py[0,1]
    I_AMPA_test=np.ones((1))*I_AMPA_test
    
    #summing total input from all other neurons
    I_AMPA_py_sum=I_AMPA_py.sum(axis=0)
    I_AMPA_fs_sum=I_AMPA_fs.sum(axis=0)
    I_GABA_py_sum=I_GABA_py.sum(axis=0)
    I_GABA_fs_sum=I_GABA_fs.sum(axis=0)
    
    #connection gating differentials
    xdot_AMPA = 5*(1+np.tanh(v_py/4))*(1-x_AMPA)-x_AMPA/2 #py cells
    xdot_GABA = 2*(1+np.tanh(v_fs/4))*(1-x_GABA)-x_GABA/tau_GABA  #fs cells
    
    #Other Currents - - - - -  - - - -
    #pyramidal cell currents - - -
    I_Na_py = g_Na_py*(m_py**3)*h_py*(v_py-E_Na_py) #sodium current 
    I_K_py = g_K_py*(n_py**4)*(v_py-E_K_py) #potassium current 
    I_leak_py = 0.1*(v_py+61) #leak current 
    I_K_ATP = g_K_ATP*z*(v_py-E_K_py) #potassium ATP current (metabolism)
    #pyramidal cell ATP dif eqs
    Nadot = 2*3*((0.000168)*1.8*np.abs(I_Na_py)-(3*(0.00000006)*ATP*Na*Na*Na))#F*np.abs(I_Na_py)-3*K_m*(Na**3)*ATP
    ATPdot = 2*5*((J_ATP*0.0004)*(2.00-ATP)-(2*(0.00000006)*ATP*Na*Na*Na))#J_ATP*(ATP_max - ATP) - K_m*(Na**3)*ATP 
    
    
    #FS cell currents - - -
    I_Na_fs = g_Na_fs*(m_fs**3)*h_fs*(v_fs-E_Na_fs)
 
    I_K_fs = g_K_fs*(n_fs**4)*(v_fs-E_K_fs)
    I_leak_fs = 0.1*(v_fs+61)
    
    #voltage equations for both cells - - -
    vdot_py = I_app_py-I_Na_py-I_K_py-I_K_ATP-I_leak_py-I_GABA_py_sum-I_AMPA_py_sum
    vdot_fs = I_app_fs-I_Na_fs-I_K_fs-I_leak_fs -I_AMPA_fs_sum-I_GABA_fs_sum
    
    #state variable arrays - - - 
    dx_py= np.concatenate([vdot_py, mdot_py, ndot_py, hdot_py, Nadot, ATPdot])
    dx_fs= np.concatenate([vdot_fs, mdot_fs, ndot_fs, hdot_fs])
    dx_con= np.concatenate([xdot_AMPA, xdot_GABA])
    dx= np.concatenate([dx_py, dx_fs, dx_con])
    return dx

##Runge kutta order 4 integration of pyramidal cell-----------------------------------
for i in np.arange(T-1):
    k1= h*hodghuxATP(y[:,int(i)],n,m,p) #first step
    k2=h*hodghuxATP(y[:,int(i)]+k1/2,n,m,p) #second step
    k3=h*hodghuxATP(y[:,int(i)]+k2/2,n,m,p) #third step
    k4=h*hodghuxATP(y[:,int(i)]+k3,n,m,p) #fourth step
    y[:,int(i+1)]=y[:,int(i)]+(1/6)*(k1+2*k2+2*k3+k4) #RK order 4
    
y_py=y[0:6*n,:] #extracting py cell states
y_fs=y[6*n:(6*n+4*m),:] #extracting fs cell states

na = y_py[n*4:n*5,:]
atp = y_py[n*5:,:]

v_py=y_py[:n,:] #pyramidal cell membarane potentials
v_fs=y_fs[:m,:] #fast spiking interneuron cell membrane potentials

x_AMPA=y[(6*n+4*m):(7*n+4*m),:]

plt.plot(t,x_AMPA[0,:])
plt.title('x_ampa')
plt.show()

I_AMPA_py = np.zeros((int(n),int(T)))
for i in np.arange(n):
    for j in np.arange(n):
        I_AMPA_py[j,:]=g_AMPA_py[i]*x_AMPA[j,:]*(v_py[i,:]-E_AMPA_py[i])


lfs = I_AMPA_py.sum(axis=0)

## plotting ---------------------------
fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = 12
fig_size[1] = 2

#plt.plot(t, I_AMPA_py[0,:])
#plt.title('actual I ampa ?')
#plt.show()

plt.plot(t,lfs, linewidth=.5)
plt.xlabel("time")
plt.ylabel("lfp")
plt.title("lfp metabolism rate = " + str(ATP_scale) + " of origional value")
plt.show()

plt.figure
for i in np.arange(n):
    plt.plot(t,v_py[i])
plt.title('membrane potential')
plt.show() 

plt.rcParams["figure.figsize"] = fig_size

plt.plot(t, atp[1,:])
plt.plot(t,na[1,:])

plt.plot(t[:4000],v_py[1,:4000])
plt.title('Hodgkin-Huxley Action Potentials')
plt.xlabel('Time (mS)')
plt.ylabel('Membrane Potential (mV)')

#beeping when code is finished
import winsound
duration = 1000  # millisecond
freq = 440  # Hz
winsound.Beep(freq, duration)