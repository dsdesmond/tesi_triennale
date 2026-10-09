import random
import numpy as np
import matplotlib.pyplot as plt






    
    
class Lattice(object):
    """
    Class for a Lattice object, which can be manipulated by either Glauber or Kawasaki dynamics after being initialised
    in a specified state.

    """
    def __init__(self, N, T, dynamics = "Glauber", lattice = None):
        """
        Initialisation function for Lattice class.

        Parameters
        ----------

        N : int
            determines size of square lattce (NxN)
        T : float
            Temperature of lattice. Tested between 1K and 3K
        dynamics : str
            Dynamics to use within the simulation, options are "Glauber" or "Kawasaki"

        lattice : numpy array or str
            An option to input starting lattice. Must be of size (NxN).
            Can also input string with options: "uniform" "ground" "halved" (see functions for specifics)

        Returns
        -------


        """

        #initialise variables
        self.N = N
        self.T = T
        self.beta = 1/self.T

        self.energies = []
        self.magnetisation = []
        self.sweep_list = []
        #N**2 is the required length of time between visualisations as in the lecture notes
        self.sweep_size = self.N**2

        #define dynamice from input
        if dynamics == "Glauber":
            self.dynamics = self.glauber
        elif dynamics == "Kawasaki":
            self.dynamics = self.kawasaki
        else:
            raise ValueError('The only options for dynamics are "Glauber" or "Kawasaki". Please input either of these')
        #print(type(lattice))
        
        if type(lattice) == np.ndarray:#'numpy.ndarray':
            if lattice.shape == (N,N):
                #print("yeye")
                self.lattice = lattice
            else:
                print("input lattice must be correct shape!!")
        
        elif lattice == "uniform":
            #assume uniform
            self.lattice = self.uniform_generate()
        elif lattice == "ground":
            self.lattice = self.ground_generate()
        elif lattice == "halved":
            self.lattice = self.halved_generate()
        else:
            print("your lattice input is wrong")
        #print(self.lattice)
            
        
        
       
    def halved_generate(self):
        """
        Generates "halved" lattice, which has half of the lattice 1 and the other half -1.
        This approximates a ground state for kawisaki dynamics

        Parameters
        ----------

        Returns
        -------
        a : numpy array
            halved (NxN) array
        """


        a = np.hstack((np.ones((self.N**2)//2),np.ones((self.N**2)//2)*-1))
        if (self.N**2)%2 !=0:
            a = np.hstack((a,[1]))
        a = np.reshape(a, (self.N,self.N))
        return a

        
        

    def uniform_generate(self):

        """
        Generates "uniform" lattice, where each point on the lattice has equal probability of being -1 or 1

        Parameters
        ----------

        Returns
        -------
        a : numpy array
            uiform (NxN) array
        """
        a = np.random.choice((1,-1), size=[self.N,self.N], replace=True, p=[0.5,0.5])
        return a
    
    def ground_generate(self):

        """
        Generates "ground" lattice, where every point on the lattice is set to 1,
        approximating the ground state for Grauber dynamics.

        Parameters
        ----------

        Returns
        -------
        a : numpy array
            uniform (NxN) array
        """


        #a = np.rint(np.random.uniform(0,1,(self.N,self.N))).astype(int)
        #a = np.where(a == 1,a,-1)
        a = np.ones((self.N,self.N))
        return a
    
    
    def glauber(self):

        """
        Function to update the lattice with glauber dynamics

        Parameters
        ----------

        Returns
        -------

        """

        #randomly choose a point to possibly flip
        flip_coords = self.rand_flip_coords()

        #find energy corresponding to this point wrt the nearest neighbours
        E_init = self.compute_E(flip_coords)
        
        #change in E  if flipped is defined by -2*init_E (proven in my notes)
        delta_E = -2*E_init

        # if energy efficient, change
        if delta_E<0:
            #flip point
            self.lattice[flip_coords[0],flip_coords[1]] = -1*self.lattice[flip_coords[0],flip_coords[1]]


        # if energy innefficient, change with probability related to temperature
        else:

            p = min(1,np.exp(-self.beta*delta_E))
            if np.random.uniform(0,1,1)<=p:
                self.lattice[flip_coords[0],flip_coords[1]] = -1*self.lattice[flip_coords[0],flip_coords[1]]

    
    
    
    def kawasaki(self):

        """
        Function to update the klattice with kawasaki dynamics

        Parameters
        ----------

        Returns
        -------

        """


        #get two different coordinates on grid
        #check not nearest neighbours (probabily a quicker way to do this
        #but will happen very rarely so doesn't matter much
        coords_same = True
        while coords_same:
            flip_coords_1 = self.rand_flip_coords()
            flip_coords_2 = self.rand_flip_coords()
            if not np.array_equal(flip_coords_1,flip_coords_2):
                coords_same = False
                
        #if lattice values are identical then theres no point in evaluating them
        if self.lattice[flip_coords_1[0],flip_coords_1[1]] == self.lattice[flip_coords_2[0],flip_coords_2[1]]:
            return

        #if nearest neighbours
        elif abs(np.sum(flip_coords_1 - flip_coords_2)) == 1 and abs(
                np.sum(flip_coords_1[0] - flip_coords_2[0])) <= 1 and abs(
                np.sum(flip_coords_1[1] - flip_coords_2[1])) <= 1:


            #compute initial energies of each point wrt the nearest neighbours
            E_init_1 = self.compute_E(flip_coords_1)
            E_init_2 = self.compute_E(flip_coords_2)

            # find the change in energy if points are swapped.
            # Same as glauber but with added factor of 4 (proven in my notes)
            delta_E = (-2 * E_init_1 - 2 * E_init_2) + 4


        
        else:
            #compute change in energy the same way as glauber dynamics (adding the two individual points)
            E_init_1 = self.compute_E(flip_coords_1)
            E_init_2 = self.compute_E(flip_coords_2)
            delta_E = -2*E_init_1 -2*E_init_2
            
            
        #if energy efficient, change
        if delta_E<0:
            self.lattice[flip_coords_1[0],flip_coords_1[1]] = -1*self.lattice[flip_coords_1[0],flip_coords_1[1]]
            self.lattice[flip_coords_2[0],flip_coords_2[1]] = -1*self.lattice[flip_coords_2[0],flip_coords_2[1]]
            return
        

        #if energy innefficient, change with probability related to temperature
        else:
            p = min(1,np.exp(-self.beta*delta_E))
            if np.random.uniform(0,1,1)<=p:
                self.lattice[flip_coords_1[0],flip_coords_1[1]] = -1*self.lattice[flip_coords_1[0],flip_coords_1[1]]
                self.lattice[flip_coords_2[0],flip_coords_2[1]] = -1*self.lattice[flip_coords_2[0],flip_coords_2[1]]
            return

        
        
        
    def find_magnetisation(self, return_val = False):
        """
        Function to update the magnetisation of the lattice

        Parameters
        ----------

        Returns
        -------

        """
        M = np.sum(self.lattice)
        self.magnetisation.append(M)
        if return_val:
            return M
        

        
    def rand_flip_coords(self):
        """
        Function to choose a random point on the lattice

        Parameters
        ----------

        Returns
        -------
        rand_coords : numpy array
            [x,y] random point on graph
        """
        #randint is very slow
        r = int(self.N * random.random())
        r_2 = int(self.N * random.random())
        #return np.random.randint(0,self.N,2)
        return np.asarray([r, r_2])


    
    def compute_tot_E(self, return_val = False):

        """
        Function to calculate the total energy of the lattice

        Parameters
        ----------

        return_val : bool
            Option to return energy of the current lattice if True.
            If False, it updates the self.energies list

        Returns
        -------
        energy : float
            optionally returns energy of the current lattice
        """

        #loop through lattice points and add energy of each.
        #should probably vectorize this but it will be ok
        E_list = []
        for i in range(self.N-1):
            for j in range(self.N-1):
                E_list.append(self.compute_E([i,j]))
        # sum energies but divide by 2 to negate overcounting
        self.energies.append(np.sum(E_list)/2)       
        if return_val:
            return self.energies[-1]  
    
    
    def find_heat_capacity(self, energies = None):
        """
        Function to calculate the heat capacity given energies.
        Can calculate from provided energies or from the self.energies.
        This dua; functionality is helpful when evaluating errors

        Parameters
        ----------
        energies : list
            optionally provide a list of energies to evaluate the heat capacity of

        Returns
        -------

        hc : float
            heat capacity

        """



        #ie if it is the full array of energies
        if energies is None:
            Es = self.energies
            
        #ie if it is a partial array of energies for evaluating errors
        else:
            Es = energies
            
        #calculate heat capacity
        hc = (1/(self.T**2*self.N**2)) *(np.mean(np.square(Es))-np.square(np.mean(Es)))
        return hc
        
        
    def susceptibility(self):
        """
        Function to calculate the susceptability

        Parameters
        ----------


        Returns
        -------

        susc : float
             susceptability
        """
        susc = 1/(self.N**2*self.T) *(np.mean(np.square(self.magnetisation))-np.square(np.mean(self.magnetisation)))       
        
        return susc
    
    
    
    def jacknife_c(self):
        """
        Function to calculate the energy in heat capacity by jacknife method
        Parameters
        ----------

        Returns
        -------

        err : float
            error in heat capacity given by jacknife method
        """
        
        c = self.find_heat_capacity()

        
        sums = 0
        for i in range(len(self.energies)-1):
            Es = np.hstack((self.energies[:i],self.energies[i:-1]))
            #print(Es)
            c_i = self.find_heat_capacity(energies = Es)
            sums += (c_i - c)**2
        err = np.sqrt(sums)
        return err
        
    
    
    
    def compute_E(self,flip_coords):

        """
        Function to compute energy with respect to the nearest neighbours of a point on the lattice

        Parameters
        ----------
        flip_coords : (x,y)
            lattice coordinates to evaluate energy of

        Returns
        -------

        E : float
            energy of point wrt. nearest neightbours

        """

        #E = -S11 * (S01+S10,S21,S12)
        #use mod (%) for periodic boundary conditions
        return  -self.lattice[flip_coords[0],flip_coords[1]]*(self.lattice[flip_coords[0]-1,flip_coords[1]]+self.lattice[flip_coords[0],flip_coords[1]-1]+self.lattice[(flip_coords[0]+1)%self.N,flip_coords[1]]+self.lattice[flip_coords[0],(flip_coords[1]+1)%self.N])
    
    
    def plot_lattice(self): 


        """
        Function to plot the lattice

        This function uses pyplot not gnuplot or other.
        It does not write to file to simplify the procedure.
        It plots directly from self.lattice.

        Parameters
        ----------

        Returns
        -------

        """

        plt.cla()
        plt.imshow(self.lattice, animated=True)
        plt.draw()
        plt.pause(0.01)
        
        
        
    def run(self,wait_sweeps = 100, num_tot_sweeps = 1000, plot_anim = True):

        """
         Function to run full dynamics given simulation inputs.
         It uses the specified dynamics of the system and simply loops through when required.

         Parameters
         ----------
         wait_sweeps : int
             number of sweeps to wait before starting measurements

        num_total_sweeps : int
            number of total sweeps to run
        plot_anim : bool
            whether or not to animate as it runs. Saves time on long simulations to not animate.

         Returns
         -------

         """


        #loop for required number of sweeps
        for i in range(num_tot_sweeps*self.sweep_size):
            #self.N**2 is the required length of time between visualisations as in the lecture notes
            if i%self.sweep_size == 0:
                #if it is time to take a measurement
                if i%(self.sweep_size*10) == 0 and i > wait_sweeps*self.sweep_size:


                    self.sweep_list.append(i/self.sweep_size )
                    self.find_magnetisation()
                    self.compute_tot_E()
                    #plot animation if required
                    if plot_anim:
                        self.plot_lattice()
                #print every 100 sweeps to check the sim progress
                if i%(self.sweep_size*100) == 0:
                    print(i/self.sweep_size)
            #run dynamics
            self.dynamics()

    
        


        

            
