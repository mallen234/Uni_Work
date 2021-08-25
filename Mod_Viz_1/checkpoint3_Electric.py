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

def init_charge_state(N):
    a = np.zeros((N,N,N))
    a[int(N/2),int(N/2),int(N/2)] = 1
    return a

def neighbours(Lattice,N,i,j,k):
    neighbours_set = Lattice[(i+1)%N,j,k] + Lattice[i,(j+1)%N,k] + Lattice[i,j,(k+1)%N] +\
    Lattice[(i-1),j,k] + Lattice[i,(j-1),k]  + Lattice[i,j,k-1]

    return neighbours_set

def Field(Lattice,N,i,j,k):

    Ex = (1/2)*(Lattice[(i+1)%N,j,k]-Lattice[(i-1),j,k])
    Ey = (1/2)*(Lattice[i,(j+1)%N,k]-Lattice[(i),j-1,k])
    Ez = (1/2)*(Lattice[(i),j,(k+1)%N]-Lattice[(i),j,k-1])

    return Ex,Ey,Ez

def normaliser(N,Ex,Ey):
    Ax = np.zeros((N,N,N))
    Ay = np.zeros((N,N,N))
    for i in range (1,N-1):
        for j in range(1,N-1):
            Ax[i,j,int(N/2)] = Ex[i,j,int(N/2)]/(np.sqrt((Ex[i,j,int(N/2)]**2+Ey[i,j,int(N/2)]**2)))
            Ay[i,j,int(N/2)] = Ey[i,j,int(N/2)]/(np.sqrt((Ex[i,j,int(N/2)]**2+Ey[i,j,int(N/2)]**2)))
    
    # norm = np.sqrt(Ex**2+Ey**2+Ez**2)
    # Ax =Ex/norm
    # Ay = Ey/norm
    # Az = Ez/norm
    
    
    return Ax,Ay

def Jacobi(Lattice,stat_latt,N):
    Fx = np.zeros((N,N,N))
    Fy = np.zeros((N,N,N))
    Fz = np.zeros((N,N,N))
    
    New_Lattice = Lattice.copy()
    for i in range(1,N-1):
        for j in range(1,N-1):
            for k in range(1,N-1):
                neighbour_sum = neighbours(Lattice,N,i,j,k)
                New_Lattice[i,j,k] = (1/6)*(neighbour_sum + stat_latt[i,j,k])

                F = Field(Lattice,N,i,j,k)
                Fx[i,j,k] =F[0]
                Fy[i,j,k] = F[1]
                Fz[i,j,k] = F[2]

    return New_Lattice,Fx,Fy,Fz

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

def SOR_Gauss_seidel(Lattice,rho_latt,N,w):
    # Ex = np.zeros((N,N,N))
    # Ey = np.zeros((N,N,N))
    # Ez = np.zeros((N,N,N))

    for i in range(1,N-1):
        for j in range(1,N-1):
            for k in range(1,N-1):
                previous = Lattice[i,j,k]
                neighbour_sum = neighbours(Lattice,N,i,j,k)

                
                Lattice[i,j,k] = (1/6)*(neighbour_sum + rho_latt[i,j,k])

                Lattice[i,j,k] = (1-w)*previous + w*Lattice[i,j,k]

                # E = E_field(Lattice,N,i,j,k)
                # Ex[i,j,k] = E[0]
                # Ey[i,j,k] = E[1]
                # Ez[i,j,k] = E[2]

    return Lattice

def Run_Efield_exp(N,accuracy,Method):
    Method_dict = {"Jacobi":Jacobi,"Gauss_seidel":Gauss_seidel}
    
    """
    Initialising Lattices
    """
    Lattice = init_state(N)  #potential Latticee
    rho_latt = init_charge_state(N)   #rho Lattice
    
    
    convergence = 10

    while convergence > accuracy:
        """
        updating
        """
        New_Lat = Lattice.copy()
        Lattice,Ex,Ey,Ez = Method_dict[Method](Lattice,rho_latt,N)

        
        difference = (New_Lat - Lattice)
        convergence = np.sum(np.abs(difference))
        #percen = (convergence/(np.mean(New_Lat))) * 100

        print(convergence)
   

    x,y = np.meshgrid(np.arange(0,N,1), np.arange(0,N,1)) #use to define x,y coords

    Ex,Ey = normaliser(N,Ex,Ey)

    #plot slice of B-field at the center of the lattice
    plt.quiver(y,x,Ex[:,:,int(N/2)],Ey[:,:,int(N/2)])
    plt.title("E Field")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.savefig("E_Field_slice.png")
    plt.show()

    print((Ex[:,:,int(N/2)]))

    #stack = np.column_stack((Ex[:,:,int(N/2)],Ey[:,:,int(N/2)]))
    
    np.savetxt("E_field_x_data.csv",Ex[:,:,int(N/2)], fmt='%1.3f', delimiter=",", newline="\n")
    np.savetxt("E_field_y_data.csv",Ey[:,:,int(N/2)], fmt='%1.3f', delimiter=",", newline="\n")
    np.savetxt("E_Potential_data.csv",Lattice[int(N/2),:,:], fmt='%1.3f', delimiter=",", newline="\n")

    #fig,im = plt.subplots(1,1)
    im=plt.imshow(Lattice[int(N/2),:,:],cmap="pink")
    plt.title("Electrostatic_Potential")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.savefig("Electrostatic_Potential.png")
    plt.show() 

def SOR_exp(N,accuracy):
    #rho Lattice

    x_list = np.linspace(1,1.9,10)
    w_list = []
    rho_latt = init_charge_state(N)

    for i in range(10,20):
            
        convergence = 10
        Lattice = init_state(N)  #potential Latticee
    
        
        print(i)
        w = (i/10)
        ticker = 0
        while convergence > accuracy:
            """
            updating
            """
            New_Lat = Lattice.copy()
            Lattice = SOR_Gauss_seidel(Lattice,rho_latt,N,w)

            
            difference = (New_Lat - Lattice)
            convergence = np.sum(np.abs(difference))
            #percen = (convergence/(np.mean(New_Lat))) * 100

            print(convergence)
            ticker += 1
        
        w_list.append(ticker)
    
    final_latt = np.column_stack((x_list,w_list))
    
    np.savetxt("SOR_data.csv",final_latt, fmt='%1.3f', delimiter=",", newline="\n")
    #np.savetxt("SOR_data.csv", x_list, fmt='%1.3f', delimiter=",", newline="\n")
    plt.plot(x_list,w_list)
    plt.xlabel("w")
    plt.ylabel("Steps to convergence")
    plt.savefig("SOR_plot.png")
    plt.show()


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
        print("In order to RELAX type SOR followed by N and accuracy")
        quit()
    
    elif sys.argv[1] == "SOR":
        N = int(sys.argv[2]) 
        accuracy = float(sys.argv[3])
        SOR_exp(N,accuracy)
        quit()

    elif sys.argv[1] != "SOR":
        N = int(sys.argv[1]) 
        accuracy = float(sys.argv[2])
        Method = str(sys.argv[3])
    
    
    if (Method == ("J") or Method == ("j")):
        Method = "Jacobi"
    
    elif (Method == ("G") or Method == ("g")):
        Method = "Gauss_seidel"
    
    Run_Efield_exp(N,accuracy,Method)


main()