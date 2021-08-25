import numpy as np
import scipy
from numpy.random import rand
from matplotlib import pyplot as plt
from matplotlib import animation
import sys

def primarystate_oscillator(N):
    #PRIMARY STATE FOR A OSCILLATOR
    Lattice = np.zeros([N,N])
    a = int(N/2)
    Lattice[a-2:a+1,a] = 1
    print(Lattice)
    return Lattice    

def primarystate_pentomino(N):
    #PRIMARY STATE FOR A PENTOMINO
    Lattice = np.zeros([N,N])
    Lattice[25,26:28] = 1
    Lattice[26,25:27] = 1
    Lattice[27,26] = 1
    return Lattice

def primarystate_glider(N):
    #PRIMARY STATE FOR A GLIDER
    
    Lattice = np.zeros([N,N])
    
    Lattice[1,0] = 1
    Lattice[1,2] = 1
    Lattice[0,2] = 1
    Lattice[2,1:3] = 1

    return Lattice

def primarystate_spaceship(N):
    Lattice = np.zeros([N,N])

    Lattice[1:4,0] = 1
    Lattice[0,1] = 1
    Lattice[3,1] = 1
    Lattice[3,2:4] = 1
    Lattice[0,4] = 1
    Lattice[2,4] = 1
    return Lattice

def primarystate(N):
    #PRIMARY STATE FOR A RANDOM CONFIG
    Lattice = np.random.randint(2, size=(N,N))
    return Lattice

def com_glider(Lattice,N):
    com = np.where(Lattice==1)
    print(com)
    a = np.mean(com[0])
    b = np.mean(com[1])

    dif1 = a - com[0][0]
    dif2 = b - com[1][0]


    if ((np.abs(dif1) > 3) or (np.abs(dif2) > 3)):
        return None

    else:
        return a,b

def neighbours(Lattice,i,j,N):

    alive_neighbours = Lattice[(i+1)%N,j] + Lattice[i,(j+1)%N]\
     + Lattice[(i-1),j] + Lattice[i,(j-1)] + Lattice[(i+1)%N,(j+1)%N] + Lattice[i-1,(j+1)%N]\
     + Lattice[(i-1)%N,(j-1)%N] + Lattice[(i+1)%N,(j-1)]

    return alive_neighbours

def total_alive(Lattice):
    return np.sum(Lattice)

def step_func(Lattice,N):
    New_Lattice = Lattice.copy()
    for i in range(N):
        for j in range(N):
            #print("latticeconfig",i,j)
            
            if (Lattice[i,j] == 0):
                a = neighbours(Lattice,i,j,N)
                #print(a)
                if (a ==3):
                    New_Lattice[i,j] = 1

            else:
                a = neighbours(Lattice,i,j,N)
                #print("a",a)
                if (a < 2) or (a > 3):
                    New_Lattice[i,j] = 0
    
    return New_Lattice

def write_to_file(name,x,y):
    with open(name, 'w') as f:
        f.write('{0:20},{1:20}\n'.format("Generation","Position",))
        for a,b in zip(x,y):
            f.write('{0:20},{1:20}\n'.format(a,b))

def plotting_func(plot_data,length):
    hist,binedges,p = plt.hist(plot_data, bins=20,range = (1,2000))
    plt.xlabel("Generations before equilibration")
    plt.ylabel("Counts")
    plt.savefig("histogram_gol.png")
    plt.show()

    plot_data = np.asarray(plot_data)

    with open('game_of_life_data.csv', 'w') as f:
        f.write('{0:20},{1:20}\n'.format("Run Number","Average length of run",))
        for a,b in zip(length,plot_data):
            f.write('{0:20},{1:20}\n'.format(a,b))
