# -*- coding: utf-8 -*-
"""
Created on Sun Mar 28 17:42:21 2021

@author: georg
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 19:14:57 2021

@author: georg
"""

# -*- coding: utf-8 -*-
"""
Created on Sun Mar 14 10:59:45 2021

@author: georg
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Mar 11 10:04:21 2021

@author: georg
"""

import math
# this imports the random number generation functions - vital for this project!
import random
import matplotlib.pyplot as plt
import numpy as np



J = 1.0
N = 10
dims = 1
tps = 30 #number of temperature points
eqsteps = 1000 #number of MC steps excluded to allow for equilibration
NTIMES = 1000 #number of MS steps for calculation
KT = 1
B = 1/KT


#starting state
def startingstate(N):
    xlattice = np.random.random(size=(N,N,N))
    ylattice = np.random.random(size=(N,N,N))
    zlattice = np.random.random(size=(N,N,N))
    for i in range(len(xlattice)):
        for j in range(len(ylattice)):
            for k in range(len(zlattice)):
                x = random.uniform(-1,1)
                y = random.uniform(-1,1)
                z = random.uniform(-1,1)
                mag = math.sqrt(x**2+y**2+z**2)
                x = x/mag
                y = y/mag
                z = z/mag
                xlattice[i,j,k]=x
                ylattice[i,j,k]=y
                zlattice[i,j,k]=z

    return xlattice,ylattice,zlattice
    
       

def calculate_energy(xlattice,ylattice,zlattice):
    energy = 0
    for i in range(len(xlattice)):
        for j in range(len(ylattice)):
            for k in range(len(zlattice)):
                sigmax = xlattice[i,j,k]
                sigmay = ylattice[i,j,k]
                sigmaz = zlattice[i,j,k]
                xenergy_from_neighbours = xlattice[(i+1)%N,j,k]+xlattice[i,(j+1)%N,k]+xlattice[i,j,(k+1)%N]+xlattice[(i-1)%N,j,k]+xlattice[i,(j-1)%N,k]+xlattice[i,j,(k-1)%N]
                yenergy_from_neighbours = ylattice[(i+1)%N,j,k]+ylattice[i,(j+1)%N,k]+ylattice[i,j,(k+1)%N]+ylattice[(i-1)%N,j,k]+ylattice[i,(j-1)%N,k]+ylattice[i,j,(k-1)%N]
                zenergy_from_neighbours = zlattice[(i+1)%N,j,k]+zlattice[i,(j+1)%N,k]+zlattice[i,j,(k+1)%N]+zlattice[(i-1)%N,j,k]+zlattice[i,(j-1)%N,k]+zlattice[i,j,(k-1)%N]
                energy += (sigmax*xenergy_from_neighbours+sigmay*yenergy_from_neighbours+sigmaz*zenergy_from_neighbours)*(-J)
                
    return energy

def calculate_magnetization(xlattice,ylattice,zlattice):
    magnetization = math.sqrt((np.sum(xlattice)**2)+(np.sum(ylattice)**2)+(np.sum(zlattice)**2))/(N)
    return magnetization

def MCstep(xlattice,ylattice,zlattice,B):
    for i in range(len(xlattice)):
        for j in range(len(ylattice)):
            for k in range(len(zlattice)):
                a,b,c=np.random.randint(0, N),np.random.randint(0, N),np.random.randint(0, N)
                sigmax=xlattice[a,b,c]
                sigmay=ylattice[a,b,c]
                sigmaz=zlattice[a,b,c]
                xenergy_from_neighbours = xlattice[(i+1)%N,j,k]+xlattice[i,(j+1)%N,k]+xlattice[i,j,(k+1)%N]+xlattice[(i-1)%N,j,k]+xlattice[i,(j-1)%N,k]+xlattice[i,j,(k-1)%N]
                yenergy_from_neighbours = ylattice[(i+1)%N,j,k]+ylattice[i,(j+1)%N,k]+ylattice[i,j,(k+1)%N]+ylattice[(i-1)%N,j,k]+ylattice[i,(j-1)%N,k]+ylattice[i,j,(k-1)%N]
                zenergy_from_neighbours = zlattice[(i+1)%N,j,k]+zlattice[i,(j+1)%N,k]+zlattice[i,j,(k+1)%N]+zlattice[(i-1)%N,j,k]+zlattice[i,(j-1)%N,k]+zlattice[i,j,(k-1)%N]
                deltaE = 2*(sigmax*xenergy_from_neighbours+sigmay*yenergy_from_neighbours+sigmaz*zenergy_from_neighbours)
                if deltaE <0:
                    sigmax*=-1
                    sigmay*=-1
                    sigmaz*=-1
                else:
                    r = random.uniform(0,1)
                    w = math.exp(-B*deltaE)
                    if r<w:
                        sigmax*=-1
                        sigmay*=-1
                        sigmaz*=-1
                xlattice[a,b,c]=sigmax
                ylattice[a,b,c]=sigmay
                zlattice[a,b,c]=sigmaz
    return xlattice,ylattice,zlattice



T = np.linspace(0,5,tps)  
E,M = [], []

for temp in range(tps):
    E1 = M1 = 0
    xlattice,ylattice,zlattice = startingstate(N)
    B = 1.0/T[temp]
    
    for i in range(eqsteps):
        MCstep(xlattice,ylattice,zlattice, B)
        
    for i in range(NTIMES):
        MCstep(xlattice,ylattice,zlattice, B)
        energy = calculate_energy(xlattice,ylattice,zlattice)
        magnetization = calculate_magnetization(xlattice,ylattice,zlattice)
        E1+=energy
        M1+=magnetization
        
        
    E.append( E1*(1.0/(NTIMES*N*N)))
    M.append(M1*(1.0/(NTIMES*N*N)))  

plt.figure(1)                    
plt.scatter(T, E, s=50, marker='o', color='Red')
plt.xlabel("Temperature (T)", fontsize=20)
plt.ylabel("Energy ", fontsize=20)        

plt.figure(2)
plt.scatter(T, M, s=50, marker='o', color='Blue')
plt.xlabel("Temperature (T)", fontsize=20); 
plt.ylabel("Magnetization ", fontsize=20);   plt.axis('tight')