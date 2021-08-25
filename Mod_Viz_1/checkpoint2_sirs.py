import numpy as np
import scipy
from numpy.random import rand
from matplotlib import pyplot as plt
from matplotlib import animation
import sys

def primary_state(N):
    Lattice = np.random.randint(3,size=(N,N)) -1
    return Lattice

def neighbours(Lattice,N,i,j):
    neighbours_set = np.array([Lattice[(i+1)%N,j], Lattice[i,(j+1)%N]\
    ,Lattice[(i-1),j],Lattice[i,(j-1)]])

    if 0 in neighbours_set:
        return True

    else:
        return False

def average_infected(Lattice,N):
    a = N*N - np.count_nonzero(np.abs(Lattice))
    return a

def step_func(Lattice,N,p1,p2,p3):

    rand1 = np.random.randint(N)
    rand2 = np.random.randint(N)

    # print(Lattice)

    if (Lattice[rand1,rand2] == -1) and neighbours(Lattice,N,rand1,rand2):
        if np.random.rand() <=  p1:
            Lattice[rand1,rand2] = 0

    elif Lattice[rand1,rand2] == 0:
        if np.random.rand() <=  p2:
            Lattice[rand1,rand2] = 1

    elif Lattice[rand1,rand2] == 1:
        if np.random.rand() <=  p3:
            Lattice[rand1,rand2] = -1

    return Lattice

def exp_1(N):
    p2 = 0.5
    p1 = np.linspace(0,1,21)
    p3 = np.linspace(0,1,21)

    aveInfected = np.zeros([21,21])
    aveVariance = np.zeros([21,21])
    
    sweeps = 10
    for i in range(21):       #TWO FOR LOOPS TO GO OVER FULL MATRIX OF P1,P3 VALUES
        #print(i)
        for j in range(21):
            Lattice = primary_state(N)
            #print(j,"--------------------------------")
            b=1
            Infected = []
            Variance = []

            for k in range(sweeps): #LOOP TO DO THE TOTAL NUMBER OF SWEEPS
                for l in range(N*N):  #LOOP THAT DOES A FULL LOOP
                    Lattice = step_func(Lattice,N,p1[j],p2,p3[i])

                if (k>1):
                    a = average_infected(Lattice,N)
                    #print(a)
                    Infected.append(a)

                    if  (a==0):
                        b=0
                        break

                
            if b == 0:
                aveVariance[i,j] = 0
                aveInfected[i,j] = 0
            
            Infected = np.asarray(Infected)
            
            meanInfect = np.mean(Infected)
            varInfect = np.mean(Infected**2) - meanInfect**2
            

            aveVariance[i,j] = varInfect/(N*N)
            aveInfected[i,j] = meanInfect/(N*N)

    X,Y = np.meshgrid(p1,p3)
    print(aveInfected)
    fig,ax=plt.subplots(1,1)
    cp = ax.contourf(X, Y, aveInfected)
    fig.colorbar(cp)
    plt.savefig("SIRS_PHASE_DIAGRAM.png")
    plt.show()

def main():
    if (str(sys.argv[1]) == "exp_1"):
        N = int(sys.argv[2])
        exp_1(N)
        exit(0)
    try:
        N = int(sys.argv[1])
        p1 = float(sys.argv[2])
        p2 = float(sys.argv[3])
        p3 = float(sys.argv[4])

    except IndexError:
        raise IndexError("Please hand 4 arguments at the command line: Size of Lattice, p1,p2,p3 ")

    except ValueError:
        raise ValueError(" Please give an integer then 3 floats")
    
    Lattice = primary_state(N)

    # p1,p2,p3 = 0.8,0.1,0.01     THIS IS FOR WAVES
    # p1,p2,p3 = 0.5,0.5,0.5
    # p1,p2,p3 = 0.1,0.7,0.1

    Sweeps = 1000

    for i in range(Sweeps):
        #print(Lattice)
        a = average_infected(Lattice,N)
        print(a)

        for i in range(N*N):
            Lattice = step_func(Lattice,N,p1,p2,p3)
        
        plt.cla()
        im=plt.imshow(Lattice,vmin=-1,vmax = +1)
        plt.draw()
        plt.pause(0.00002)

main()