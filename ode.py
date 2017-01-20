import numpy as np

def ode_euler(x0,t0,v0,tf,dt,m,k):
    # initialize counter 
    i = 0
    
    # make arrays for computation 
    time = np.arange(t0,tf,dt)
    pos = np.zeros(len(time))
    vel = np.zeros(len(time))
    
    # intial conditions
    pos[0] = x0
    vel[0] = v0
    while i < len(time)-1: 
        vel[i+1] = vel[i] - (k/m)*pos[i]*dt
        pos[i+1] = pos[i] + vel[i]*dt
        
        # update the counter
        i = i + 1
    return time,pos,vel


def ode_pc(x0,t0,v0,tf,dt,m,k):
    # initialize counter 
    i = 0
    
    # make arrays for computation 
    time = np.arange(t0,tf,dt)
    pos = np.zeros(len(time))
    vel = np.zeros(len(time))
        
    # intial conditions
    pos[0] = x0
    vel[0] = v0
    
    while i < len(time)-1: 
        # predict next values
        pos[i+1] = pos[i] + vel[i]*dt
        vel[i+1] = vel[i] - pos[i]*dt
        
        # correct them 
        pos[i+1] = pos[i] + (vel[i]+vel[i+1])*(dt/2.)
        vel[i+1] = vel[i] - (pos[i]+pos[i+1])*(dt/2.)
        
        # update the counter
        i = i + 1
    
    return time,pos,vel