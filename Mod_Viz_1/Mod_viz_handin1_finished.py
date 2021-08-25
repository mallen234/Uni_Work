import numpy as np
import scipy
from numpy.random import rand
from matplotlib import pyplot as plt
from matplotlib import animation
import sys
from astropy.stats import jackknife_resampling
from astropy.stats import jackknife_stats

def primarystate_glauber(N):
    Lattice = np.zeros([N,N])
    Lattice+=1

    return Lattice

def primarystate_kawasaki(N):
    Lattice = np.zeros([N,N])
    a = int(N/2)
    Lattice[:,:a]+=1
    Lattice[:,a:]-=1

    print(np.sum(Lattice))

    return Lattice    

def primarystate(N):
    Lattice = 2*np.random.randint(2, size=(N,N))-1

    return Lattice

def mag_func(Lattice,N):
    magnetisation = np.sum(Lattice)
    return magnetisation

def energy_func(Lattice,N):
    energy = 0
    for i in range (N):
        for j in range(N):
            S = Lattice[i,j]
            nb = Lattice[(i+1)%N, j] + Lattice[i,(j+1)%N] + Lattice[(i-1)%N, j] + Lattice[i,(j-1)%N]
            energy += -nb*S

    return energy

def plotting_func_Glauber(x,mag,sus,ene,shc,shc_error):
    with open('data_file_Glauber.csv', 'w') as f:
        f.write('{0:>20},{1:>20},{2:>20},{3:>20},{4:>20},{5:>20}\n'.format("T","MAG","SUS","ENERGY","SHC","SHC_error"))
        for a,b,c,d,e,g in zip(x,mag,sus,ene,shc,shc_error):
            f.write('{0:20},{1:20},{2:20},{3:20},{4:20},{5:20}\n'.format(a,b,c,d,e,g))
    
    plt.show()
    plt.scatter(x,sus)
    plt.xlabel("Temperature")
    plt.ylabel("Susceptibility")
    plt.title("Susceptibility against Temperature: Glauber")
    plt.savefig("fig_Glauber_1.png")
    plt.show()
    
    plt.scatter(x,mag)
    plt.xlabel("Temperature")
    plt.ylabel("Magnetisation")
    plt.title("Magnetisation against Temperature: Glauber")
    plt.savefig("fig_Glauber_2.png")
    plt.show()
    
    plt.scatter(x,ene)
    plt.xlabel("Temperature")
    plt.ylabel("energy")
    plt.title("energy against Temperature: Glauber")
    plt.savefig("fig_Glauber_3.png")
    plt.show()
    
    plt.scatter(x,shc)
    plt.xlabel("Temperature")
    plt.ylabel("shc")
    plt.title("shc against Temperature: Glauber")
    plt.errorbar(x,shc,yerr = shc_error)
    plt.savefig("fig_Glauber_4.png")
    plt.show()

def plotting_func_Kawasaki(x,ene,shc,shc_error):
    with open('data_file_Kawasaki.csv', 'w') as f:
        f.write('{0:15},    {1:15},    {2:15}    \n'.format("T","ENERGY","SHC","Error"))
        for a,b,c,d in zip(x,ene,shc,shc_error):
            f.write('{0:15},  {1:15},  {2:15},   {3:15}  \n'.format(a,b,c,d))        
    
    plt.scatter(x,ene)
    plt.xlabel("Temperature")
    plt.ylabel("energy")
    plt.title("energy against Temperature: Kawaski")
    plt.savefig("fig_Kawasaki_1.png")
    plt.show()
    
    plt.scatter(x,shc)
    plt.xlabel("Temperature")
    plt.ylabel("shc")
    plt.title("shc against Temperature: Kawaski")
    plt.errorbar(x,shc,yerr = shc_error)
    plt.savefig("fig_Kawasaki_2.png")
    plt.show()

def Glauber(Lattice,N,beta):
    rand1 = np.random.randint(N)
    rand2 = np.random.randint(N)
    
    state = Lattice[rand1,rand2]
    
    deltaE = Lattice[rand1-1,rand2] + Lattice[[0,rand1+1][rand1<N-1],rand2] \
        + Lattice[rand1,rand2-1] + Lattice[rand1,[0,rand2+1][rand2<N-1]]
    deltaE *= 2*state
    
    if (deltaE <= 0):
        state *= -1
    elif (np.random.rand() <= np.exp(-deltaE*beta)):
        state *= -1
          
    Lattice[rand1,rand2] = state
    
    return Lattice

def Kawasaki(Lattice,N,beta):
    rand1,rand2 = np.random.randint(N),np.random.randint(N)
    rand3,rand4 = np.random.randint(N),np.random.randint(N)
    
    state1 = Lattice[rand1,rand2]
    state2 = Lattice[rand3,rand4]
    
    if (state1==state2):
        pass
    else:
        
        deltaE1= Lattice[rand1-1,rand2] + Lattice[[0,rand1+1][rand1<N-1],rand2]\
        + Lattice[rand1,rand2-1] + Lattice[rand1,[0,rand2+1][rand2<N-1]]
        deltaE1 *= state1*2
        
        deltaE2 = Lattice[rand3-1,rand4] + Lattice[[0,rand3+1][rand3<N-1],rand4]\
        + Lattice[rand3,rand4-1] + Lattice[rand3,[0,rand4+1][rand4<N-1]]
        deltaE2 *= state2*2 
        
        delta E
        xdiff = abs(rand1-rand3) % N
        ydiff = abs(rand2-rand4) % N

        if ((rand1 == rand3) and (ydiff == 1)) or ((rand2 == rand4) and (xdiff ==1)):
            T_DE = (deltaE1 + deltaE2) +4
        

        if (abs(rand1 % (len(lattice)-1) - rand3 % (len(lattice)-1)) <= 1 and\
        abs(rand2 % (len(lattice)-1) - rand4 % (len(lattice)-1)) <= 1):
            T_DE = (deltaE1 + deltaE2) +4

        

        else:
            T_DE = (deltaE1 + deltaE2)
        
        
        if (T_DE <= 0):
            state1 *= -1
            state2 *= -1
        elif (np.random.rand() <= np.exp(-T_DE*beta)):
            state1 *= -1
            state2 *= -1

        Lattice[rand1,rand2] = state1
        Lattice[rand3,rand4] = state2
        
    return Lattice

def Run_exp(Method,N):
    """
    SETTING UP THE EXPERIMENT
    """
    kb = 1
    Temp_graph_xaxis = np.round(np.linspace(1,3,21),decimals=1)
    Temp_graph_var = []
    Temp_graph_mag = []
    Temp_graph_ene = []
    Temp_graph_shc = []
    ErrorShc =[]
    Method_dict = {"Glauber":Glauber,"Kawasaki":Kawasaki}
    

    """
    FOR LOOP TO INCREASE TEMPERATURE
    """
    if (Method == "Glauber"):
        Lattice = primarystate_glauber(N)

    elif(Method == "Kawasaki"):
        Lattice = primarystate_kawasaki(N)

    for k in range(10,31):
        
        T = np.round(k/10,decimals=1)
        beta = kb/T
        magnetisation = []
        energy = []

        print(T)

        """
        THE BIG LOOP TO UPDATE THE LATTICE
        """
        for i in range(10000):                   # NUMBER OF SWEEPS
            for j in range(N*N):
                Lattice = Method_dict[Method](Lattice,N,beta)
            
            if ((i >= 200) and (i%10 == 0)):
                if (Method == "Glauber"):
                    a = mag_func(Lattice,N)
                    magnetisation.append(a)
                
                b = energy_func(Lattice,N)
                energy.append(b)

        
        
        """
        ENERGY LISTS AND SUCH
        """
        energy = np.asarray(energy)
        ene_ave = (np.mean(energy))
        ene_sq = np.mean(energy**2)
        shc = (1/(kb*(T**2)))*(ene_sq - ene_ave**2)

        """
        JACKKNIFE -------------ERROR ANALYSIS ----------------------
        """
        resamples = jackknife_resampling(energy)
        resSquaredMean = np.mean(resamples**2,axis=1)
        meanResSquared = (np.mean(resamples,axis=1))**2
        
        Error_single = 0
        for l in range(len(meanResSquared)):
            single_shc = ((1/(kb*(T**2)))*(resSquaredMean[l] -meanResSquared[l]))

            Error_single+= (single_shc - shc)**2

        Error_single = np.sqrt(Error_single)
        """
        END OF ERROR ANALYSIS
        """

        ErrorShc.append(Error_single)
        Temp_graph_shc.append(shc)
        Temp_graph_ene.append(ene_ave)
        
        """
        MAGNETISATION LISTS AND SUCH
        """
        if (Method == "Glauber"):
            magnetisation = np.asarray(magnetisation)
            ave_sq = np.abs(np.mean(magnetisation))
            #print(ave_sq)
            mag_sq = np.mean(magnetisation**2)
            susceptibility = (1/(N*N*kb*T))*(mag_sq - ave_sq**2)
            
            Temp_graph_mag.append(ave_sq)
            Temp_graph_var.append(susceptibility)

    """
    SENDING DATA TO PLOT AND FILE
    """
    if (Method == "Glauber"):
        plotting_func_Glauber(Temp_graph_xaxis,Temp_graph_mag,Temp_graph_var,Temp_graph_ene,Temp_graph_shc,ErrorShc)

    elif(Method == "Kawasaki"):
        plotting_func_Kawasaki(Temp_graph_xaxis,Temp_graph_ene,Temp_graph_shc,ErrorShc)

def main():
    
    if (sys.argv[1] == "run_exp"):
        Method = str(sys.argv[2])
        N = int(sys.argv[3])
        Run_exp(Method,N)
        exit(0)

    else:
        Method = str(sys.argv[1])
        N = int(input("Enter N:"))
        T = int(input("Enter Temperature:"))
    
    kb = 1
    
    Method_dict = {"Glauber":Glauber,"Kawasaki":Kawasaki}
    
    
    Lattice = primarystate(N) 
    beta = kb/T

    for i in range(10000):
        for j in range(N**2):
            Lattice = Method_dict[Method](Lattice,N,beta)
            
            #if (j % 10000 == 0):
        plt.cla()
        im=plt.imshow(Lattice, animated=True,cmap = "Pastel1",vmin=-1,vmax=1)
        plt.draw()
        plt.pause(0.0001)
main()