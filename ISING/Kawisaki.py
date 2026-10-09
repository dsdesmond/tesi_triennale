import time
import sys
import numpy as np
from Lattice import Lattice


def main(T_low, N = 50, num_sweeps = 1000, plot_anim = False, T_high = None, stepsize = None):
    """
    Function to run Kawasaki dynamics on an ising lattice.

    Saves results to file named "KAW_ALL" where each column is an ising lattice evaluated at a different temperature, and
    each row is as follows: T  M  E  c  Chi  c_err


    Parameters
    ----------
    T_low : float
        Temperature to evaluate ising lattice. Lowest value of a range of values if T_high and stepsize are provided.
    N: int
        Size of lattice. Assumes square NxN lattice.
    num_sweeps: int
        total number of sweeps required in the simulation.
    plot_anim: bool
        Whether animation should be plotted every 10 sweeps or not.
    T_high: float
        Upper bound of temperature to evaluate if a range of temperatures are to be evaluated
    stepsize: int
        stepsize between temperatures if a range of temperatures are to be evaluated

    Returns
    -------


    """

    #check if a range or T are to be evaluated
    if T_high is not None and stepsize is not None:
        #nte that a range isnt wanted, to be used later
        rang = True
        #make required range
        T_arr = np.arange(T_low, T_high, stepsize)

    if T_high is None or stepsize is None:
        #range isnt required
        rang = False
        T_arr = np.asarray([int(T_low)])




    #define how many wait sweeps to do. Generally 100 was seen to be appropriate,
    #however if for some reason less total sweeps than this are required, it adjusts
    wait_sweeps = 100
    if wait_sweeps >= num_sweeps:
        wait_sweeps = num_sweeps//10
    
    #loop through temperatures
    lattices = []
    for count, T in enumerate(T_arr):

        #if temperature = 1, assume ground state
        if T == 1:
            L = Lattice(N,T, dynamics = "Kawasaki", lattice = "halved")
        else:
            #use starting lattice from end of last lattice if not first loop
            if count !=0:
                L = Lattice(N,T, dynamics = "Kawasaki", lattice = lattices[-1].lattice)
            #if first loop assume uniform distribution
            else:
                L = Lattice(N,T, dynamics = "Kawasaki", lattice = "uniform")

        #run simulation
        L.run(wait_sweeps = wait_sweeps , num_tot_sweeps = num_sweeps , plot_anim = plot_anim)
        
        lattices.append(L)
        

    
    
    susc_list = []
    hc_list = []
    M_list = []
    E_list = []
    M_all_list = []
    E_all_list = []
    hc_err_list = []

        
    for L in lattices:
        M_list.append(np.mean(L.magnetisation))
        E_list.append(np.mean(L.energies))
        hc_list.append(L.find_heat_capacity())
        hc_err_list.append(L.jacknife_c())
        M_all_list.append(L.magnetisation)
        E_all_list.append(L.energies)
        susc_list.append(L.susceptibility())
    
    
    np.savetxt("G_results/KAW_ALL",np.vstack((T_arr,M_list,E_list,hc_list,susc_list, hc_err_list)))#, delimiter = ",")
    #np.savetxt("G_results/KAW_M_ALL",M_all_list)#, delimiter = ",")
    #np.savetxt("G_results/KAW_E_ALL",E_all_list)#, delimiter = ",")
    
    

    return
  
 
 
 
 
 
 
 




if __name__ == "__main__":





    #if incorrect number of arguments are given, state what arguments are expected for which process
    if len(sys.argv) != 7 and len(sys.argv) != 5:
        print("Usage for a single T: \n T (float), N (int), num_sweeps (int), plot_anim (bool)")
        print("Usage for a range of T: \n T_low (float), N (int), num_sweeps (int), plot_anim (bool), T_high (float), stepsize (int)")
        print("For further usage instructions, see README or documentation")
        sys.exit(1)

    #if only one Temperature is to be evaluated
    elif len(sys.argv) == 5:
        main(float(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), eval(sys.argv[4]))
    #if a range of temperatures re to be evaluated
    elif len(sys.argv) == 7:
        main(float(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), eval(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6]))

