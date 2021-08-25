import numpy as np
import matplotlib.pyplot as plt
import sys
import astropy
import random

def init_state(N,phi0,sigma):
    #Lattice = np.zeros((N,N))
    dist = np.random.normal(phi0,sigma,size=(N,N))
    #print(dist)
    #Lattice = Lattice + dist
    return dist

def initial_state(N):
    #PRIMARY STATE FOR A RANDOM CONFIG
    #Lattice = np.zeros((N,N))
    Lattice = (2*np.random.random_sample((N,N)) - 1)
    return Lattice

def neighbours(Lattice,N,i,j):
    neighbours_set = Lattice[(i+1)%N,j] + Lattice[i,(j+1)%N] +\
    Lattice[(i-1),j] + Lattice[i,(j-1)]

    return neighbours_set

def mu_func(Lattice,N,i,j):
        
    mu = -a*Lattice[i,j] + a*Lattice[i,j]**3 - k*(dt/dx)\
        *(neighbours(Lattice,N,i,j)-4*Lattice[i,j])  
    return mu

def phi_function(Lattice,N,i,j):
    return True

def step_func(Lattice,N,a,M):
    New_Lattice = Lattice.copy()
    f= np.zeros((N,N))
    for i in range (N):
        for j in range (N):
            #mu = mu_func(Lattice,N,i,j)

            New_Lattice[i,j] = Lattice[i,j] + M*(dt/dx)*\
                (mu_func(Lattice,N,(i+1)%N,j) + mu_func(Lattice,N,i-1,j) + \
                    mu_func(Lattice,N,i,(j+1)%N) + mu_func(Lattice,N,i,j-1)  - 4*mu_func(Lattice,N,i,j))

            interm = (1/4)*(Lattice[(i+1)%N,j] - Lattice[i-1,j] + Lattice[i,(j+1)%N] - Lattice[i,j-1])

            f[i,j] = -(a/2)*(Lattice[i,j]**2) + (a/4)*(Lattice[i,j]**4) + (k/2)* (interm**2)

    return New_Lattice,f

def main():
    global a
    global k
    global dt 
    dt = 1
    global dx
    dx = 1

    N = int(sys.argv[1])
    phi0 = float(sys.argv[2])
    timesteps = int(sys.argv[3])

    """
    INPUTS:
    1 = N
    2 = phi
    3 = timesteps
    """

    Lattice2 = initial_state(N)

    Lattice = init_state(N,phi0,0.25)
    
    
    a,M,k = 0.1,0.1,0.1

    energy_list = []

    for i in range(timesteps):    
        #plt.cla()
        #fig,im = plt.subplots(1,1)
        #im=plt.imshow(Lattice,cmap="pink")
        #fig.colorbar(im)
        #plt.draw()
        #plt.pause(0.00002)
        #plt.close(fig)
        
        if (i%100 == 0):
            print(i)
        
        Lattice,f = step_func(Lattice,N,a,M)
        
        Energy = np.sum(f)
        
        energy_list.append(Energy)
    #print(energy_list)

    plt.show()
    x_list = np.linspace(0,timesteps-1,timesteps)
    plt.plot(x_list,energy_list)
    plt.xlabel("Time")
    plt.ylabel("Free energy")
    plt.savefig("free_energy_vs_time.png")
    plt.show()

    stack = np.column_stack((x_list,energy_list))
    np.savetxt(f"Cahn_Hilliard_phi_{phi0}.csv",stack, fmt='%1.3f', delimiter=",", newline="\n")
    print("done")
        

main()
