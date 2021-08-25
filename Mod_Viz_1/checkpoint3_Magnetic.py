import math 
import numpy as np
import matplotlib.pyplot as plt
import sys
import astropy
import random
from mpl_toolkits.mplot3d import axes3d
#np.seterr(divide='ignore', invalid='ignore')

def init_state(N):
    a = np.zeros((N,N,N))
    return a

def init_wire(N):
    a = np.zeros((N,N,N))
    a[int(N/2),int(N/2),:] = 1
    return a

def normaliser(N,Ex,Ey):
    Ax = np.zeros((N,N,N))
    Ay = np.zeros((N,N,N))
    for i in range (1,N-1):
        for j in range(1,N-1):
            Ax[i,j,int(N/2)] = Ex[i,j,int(N/2)]/(np.sqrt((Ex[i,j,int(N/2)]**2+Ey[i,j,int(N/2)]**2)))
            Ay[i,j,int(N/2)] = Ey[i,j,int(N/2)]/(np.sqrt((Ex[i,j,int(N/2)]**2+Ey[i,j,int(N/2)]**2)))
    
    return Ax,Ay

def Field(Lattice,N,i,j,k):
    
     
    Ey = -(1/2)*(Lattice[(i+1)%N,j,k]-Lattice[(i-1),j,k]) 
    Ex = (1/2)*(Lattice[i,(j+1)%N,k]-Lattice[(i),j-1,k]) 
    
    
    Ez = (1/2)*(Lattice[(i),j,(k+1)%N]-Lattice[(i),j,k-1])

    return Ex,Ey,Ez

def neighbours(Lattice,N,i,j,k):
    neighbours_set = Lattice[(i+1)%N,j,k] + Lattice[i,(j+1)%N,k] + Lattice[i,j,(k+1)%N] +\
    Lattice[(i-1),j,k] + Lattice[i,(j-1),k]  + Lattice[i,j,k-1]

    return neighbours_set

def Jacobi(Lattice,rho_latt,N):
    Bx = np.zeros((N,N,N))
    By = np.zeros((N,N,N))
    Bz = np.zeros((N,N,N))
    
    New_Lattice = Lattice.copy()
    for i in range(1,N-1):
        for j in range(1,N-1):
            for k in range(1,N-1):
                neighbour_sum = neighbours(Lattice,N,i,j,k)
                New_Lattice[i,j,k] = (1/6)*(neighbour_sum + rho_latt[i,j,k])

                B = Field(Lattice,N,i,j,k)
                Bx[i,j,k] =B[0]
                By[i,j,k] = B[1]
                Bz[i,j,k] = B[2]

    return New_Lattice,Bx,By,Bz

def Gauss_seidel(Lattice,rho_latt,N):
    Ex = np.zeros((N,N,N))
    Ey = np.zeros((N,N,N))
    Ez = np.zeros((N,N,N))

    for i in range(1,N-1):
        for j in range(1,N-1):
            for k in range(1,N-1):
                neighbour_sum = neighbours(Lattice,N,i,j,k)
                Lattice[i,j,k] = (1/6)*(neighbour_sum + rho_latt[i,j,k])

                E = Field(Lattice,N,i,j,k)
                Ex[i,j,k] = E[0]
                Ey[i,j,k] = E[1]
                Ez[i,j,k] = E[2]

    return Lattice,Ex,Ey,Ez

def Run_Bfield_exp(N,accuracy,Method):
    Method_dict = {"Jacobi":Jacobi,"Gauss_seidel":Gauss_seidel}
    
    """
    Initialising Lattices
    """
    Lattice = init_state(N)  #potential Latticee
    Az_latt = init_wire(N)   #rho Lattice
    
    convergence = 10

    while convergence > accuracy:
        """
        updating
        """
        New_Lat = Lattice.copy()
        Lattice,Bx,By,Bz = Method_dict[Method](Lattice,Az_latt,N)
        
        difference = (New_Lat - Lattice)
        convergence = np.sum(np.abs(difference))
        #percen = (convergence/(np.mean(New_Lat))) * 100

        print(convergence)
   

    x,y = np.meshgrid(np.arange(0,N,1), np.arange(0,N,1)) #use to define x,y coords

    Bx,By = normaliser(N,Bx,By)
    plt.quiver(y,x,Bx[:,:,int(N/2)],By[:,:,int(N/2)])
    plt.title("B Field")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.savefig("B_Field_slice.png")
    plt.show()

    im=plt.imshow(Lattice[:,:,int(N/2)],cmap="pink")
    plt.title("Z component of Vector Potential")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.savefig("Vector_potential_slice.png")
    plt.show()


    np.savetxt("B_field_x_data.csv",Bx[:,:,int(N/2)], fmt='%1.3f', delimiter=",", newline="\n")
    np.savetxt("B_field_y_data.csv",By[:,:,int(N/2)], fmt='%1.3f', delimiter=",", newline="\n")
    np.savetxt("B_Potential_data.csv",Lattice[int(N/2),:,:], fmt='%1.3f', delimiter=",", newline="\n")


def main():
    global dx
    global e_0
    global dt
    x,e_0 = 1,1
    dt = 1/2
    """
    input 1 = N
    input 2 = accuracy
    input 3 = gauss or jacobi
    """
    if len(sys.argv)!=4:
        print("Number of arguments incorrect")
        print("input 1 = N \ninput 2 = accuracy  \ninput 3 = gauss or jacobi")
        quit()

    N = int(sys.argv[1]) 
    accuracy = float(sys.argv[2])
    Method = str(sys.argv[3])
    
    if (Method == ("J") or Method == ("j")):
        Method = "Jacobi"
    
    elif (Method == ("G") or Method == ("g")):
        Method = "Gauss_seidel"
    
    Run_Bfield_exp(N,accuracy,Method)


main()