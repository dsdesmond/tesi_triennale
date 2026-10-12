import random
import numpy as np
import matplotlib.pyplot as plt






    
    
class Lattice(object):

    """
    Class for a Lattice object, which can be manipulated by either Glauber or Kawasaki dynamics after being initialised
    in a specified state.
    """

    #def __init__(self, N, T, dynamics = "Glauber", lattice = None):
    def __init__(self, N, T, dynamics = "Glauber", lattice = None, J = -1.0):
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
        self.J = J

        self.energies = []
        self.magnetisation = []
        self.sweep_list = []
        #N**2 is the required length of time between visualisations as in the lecture notes
        self.sweep_size = self.N**2

        if type(lattice) == str:
            self.nome_config = lattice
        else:
            self.nome_config = "custom"

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


        #tipi di reticoli disponibili
        
        elif lattice == "uniform":
            #assume uniform
            self.lattice = self.uniform_generate()
        # elif lattice == "ground":  è una configurazione che non esiste nel caso della lega binaria
        #     self.lattice = self.ground_generate()
        elif lattice == "halved":
            self.lattice = self.halved_generate()
        elif lattice == "circle":
            self.lattice = self.circle_generate()
        elif lattice == "edges":
            self.lattice = self.edges_generate()

        else:
            print("your lattice input is wrong")
        #print(self.lattice)
            
        

    def calc_num_B(self):
        """
        Calcola il numero esatto (intero) di atomi di tipo B in base alla concentrazione richiesta.
        Garantisce che, in caso di concentrazioni estreme, sia sempre 
        presente almeno 1 atomo della minoranza per mantenere il senso fisico della lega.
        """
        N_tot = self.N ** 2
        N_B = int(round(N_tot * self.conc_B))
        
        # Se l'utente ha chiesto una concentrazione > 0, ma l'arrotondamento 
        # l'ha portata a 0 (es: 100 atomi al 0.1%), forziamo almeno 1 atomo B.
        if self.conc_B > 0.0 and N_B == 0:
            N_B = 1
            
        # Allo stesso modo, se la concentrazione è < 1 (es: 99.9%), ma l'arrotondamento
        # satura la griglia, togliamo un atomo B per garantire almeno 1 atomo A.
        elif self.conc_B < 1.0 and N_B == N_tot:
            N_B = N_tot - 1
            
        return N_B

    

    # def circle_generate(self):

    #     """
    #     Genera un reticolo con tutti spin SU (+1) e un cerchio centrale di spin GIU (-1).
    #     """

    #     # partiamo da un reticolo tutto a +1
    #     a = np.ones((self.N, self.N))
        
    #     # raggio è un quarto della grandezza totale
    #     centro = self.N // 2
    #     raggio = self.N // 4  
        
        
    #     for i in range(self.N):
    #         for j in range(self.N):
    #             # equazione del cerchio: (x-x0)^2 + (y-y0)^2 <= r^2
    #             if (i - centro)**2 + (j - centro)**2 <= raggio**2:
    #                 a[i, j] = -1
                    
    #     return a


    def circle_generate(self):

        """
        Genera un cerchio centrale di spin -1 (atomi B).
        Gli atomi vengono assegnati partendo dal centro esatto e riempiendo a raggiera in modo radiale, garantendo
        matematicamente che il numero di atomi B sia esattamente quello richiesto dalla concentrazione.
        """
        # reticolo di atomi A (+1)
        a = np.ones((self.N, self.N))
        
        # numero atomi B che ci servono esattamente
        N_B = self.calc_num_B()
        
        # Se la concentrazione è zero, non serve fare nulla
        if N_B == 0:
            return a
            
        # Calcoliamo il centro della griglia
        # Usiamo (N-1)/2 per centrare perfettamente sia su griglie pari che dispari
        centro = (self.N - 1) / 2.0
        
        # costruiamo una lista che conterrà (distanza_dal_centro, coordinata_x, coordinata_y)
        distanze = []
        for i in range(self.N):
            for j in range(self.N):
                dist_quadrata = (i - centro)**2 + (j - centro)**2
                distanze.append((dist_quadrata, i, j))
                
        # ordiniamo la lista dalla distanza minore (il centro) alla maggiore (i bordi)
        distanze.sort(key=lambda x: x[0])
        
        # Coloriamo di -1 esattamente i primi N_B punti della lista
        for k in range(N_B):
            _, i, j = distanze[k]  # Il _ serve a ignorare la distanza, ci servono solo i e j (conosco distanze[k])
            a[i, j] = -1
            
        return a



    # def edges_generate(self):

    #     """
    #     Genera un reticolo con spin SU (+1) all'interno e spin GIU (-1) solo sui bordi da spessore variabile.
    #     """

    #     spessore = 4  #spessore del bordo
        

    #     # partiamo da un reticolo tutto a +1
    #     a = np.ones((self.N, self.N))
        
    #     # q -1 prime e ultime righe in base allo spessore
    #     a[0:spessore, :] = -1
    #     a[-spessore:, :] = -1
        
    #     # a -1 prime e ultime colonne in base allo spessore
    #     a[:, 0:spessore] = -1
    #     a[:, -spessore:] = -1
        
    #     return a


    def edges_generate(self):
        """
        Genera atomi B (-1) sui bordi del reticolo, 
        procedendo per anelli quadrati concentrici dall'esterno verso l'interno. 
        Gli anelli vengono riempiti in senso orario. 
        Se gli atomi finiscono a metà, l'anello interno si interrompe, creando uno spessore asimmetrico continuo.
        """
        a = np.ones((self.N, self.N))
        N_B = self.calc_num_B()
        
        if N_B == 0:
            return a
            
        layer = 0
        
        # Continuiamo a sfogliare anelli finché abbiamo atomi B e non superiamo il centro
        while N_B > 0 and layer < self.N / 2.0:
            min_idx = int(layer)
            max_idx = int(self.N - 1 - layer)
            
            if min_idx == max_idx:
                # Caso limite: se il reticolo è dispari ed è rimasto solo il pixel centrale
                if N_B > 0:
                    a[min_idx, min_idx] = -1
                    N_B -= 1
            else:
                # Lato superiore (da sinistra a destra)
                for j in range(min_idx, max_idx + 1):
                    if N_B > 0:
                        a[min_idx, j] = -1
                        N_B -= 1
                        
                # Lato destro (dall'alto verso il basso, escludendo l'angolo già fatto)
                for i in range(min_idx + 1, max_idx + 1):
                    if N_B > 0:
                        a[i, max_idx] = -1
                        N_B -= 1
                        
                # Lato inferiore (da destra a sinistra, escludendo l'angolo già fatto)
                for j in range(max_idx - 1, min_idx - 1, -1):
                    if N_B > 0:
                        a[max_idx, j] = -1
                        N_B -= 1
                        
                # Lato sinistro (dal basso verso l'alto, escludendo gli angoli)
                for i in range(max_idx - 1, min_idx, -1):
                    if N_B > 0:
                        a[i, min_idx] = -1
                        N_B -= 1
                        
            layer += 1
            
        return a


       
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



    
    # def ground_generate(self):        configurazione che non esiste nel caso della lega binaria

    #     """
    #     Generates "ground" lattice, where every point on the lattice is set to 1,
    #     approximating the ground state for Grauber dynamics.

    #     Parameters
    #     ----------

    #     Returns
    #     -------
    #     a : numpy array
    #         uniform (NxN) array
    #     """

    #     #a = np.rint(np.random.uniform(0,1,(self.N,self.N))).astype(int)
    #     #a = np.where(a == 1,a,-1)
    #     a = np.ones((self.N,self.N))
    #     return a
    

    
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
        
        #change in E  if flipped is defined by -2*init_E (proven in my notes)  QUANTITA ENERGIA PER FLIPPARE UNO SPIN
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
        Function to update the lattice with kawasaki dynamics

        Parameters
        ----------

        Returns
        -------
        """


        # #get two different coordinates on grid
        # #check not nearest neighbours (probabily a quicker way to do this
        # #but will happen very rarely so doesn't matter much
        # coords_same = True
        # while coords_same:
        #     flip_coords_1 = self.rand_flip_coords()
        #     flip_coords_2 = self.rand_flip_coords()
        #     if not np.array_equal(flip_coords_1,flip_coords_2):
        #         coords_same = False
                
        # #if lattice values are identical then theres no point in evaluating them
        # if self.lattice[flip_coords_1[0],flip_coords_1[1]] == self.lattice[flip_coords_2[0],flip_coords_2[1]]:
        #     return

        # #if nearest neighbours
        # elif abs(np.sum(flip_coords_1 - flip_coords_2)) == 1 and abs(
        #         np.sum(flip_coords_1[0] - flip_coords_2[0])) <= 1 and abs(
        #         np.sum(flip_coords_1[1] - flip_coords_2[1])) <= 1:


        #     #compute initial energies of each point wrt the nearest neighbours
        #     E_init_1 = self.compute_E(flip_coords_1)
        #     E_init_2 = self.compute_E(flip_coords_2)

        #     # find the change in energy if points are swapped.
        #     # Same as glauber but with added factor of 4 (proven in my notes)
        #     delta_E = (-2 * E_init_1 - 2 * E_init_2) + 4


        
        # else:
        #     #compute change in energy the same way as glauber dynamics (adding the two individual points)
        #     E_init_1 = self.compute_E(flip_coords_1)
        #     E_init_2 = self.compute_E(flip_coords_2)
        #     delta_E = -2*E_init_1 -2*E_init_2
            
            
        # #if energy efficient, change
        # if delta_E<0:
        #     self.lattice[flip_coords_1[0],flip_coords_1[1]] = -1*self.lattice[flip_coords_1[0],flip_coords_1[1]]
        #     self.lattice[flip_coords_2[0],flip_coords_2[1]] = -1*self.lattice[flip_coords_2[0],flip_coords_2[1]]
        #     return
        

        # #if energy innefficient, change with probability related to temperature
        # else:
        #     p = min(1,np.exp(-self.beta*delta_E))
        #     if np.random.uniform(0,1,1)<=p:
        #         self.lattice[flip_coords_1[0],flip_coords_1[1]] = -1*self.lattice[flip_coords_1[0],flip_coords_1[1]]
        #         self.lattice[flip_coords_2[0],flip_coords_2[1]] = -1*self.lattice[flip_coords_2[0],flip_coords_2[1]]
        #     return
        #               questa funzione scambia 2 spin anche se non vicini, ergo c'è un teletrasporto che non ha senso fisico


        # seleziono il primo spin a caso
        coords_1 = self.rand_flip_coords()
        r1, c1 = coords_1[0], coords_1[1]
        
        # seleziono un vicino a caso (0: Su, 1: Giù, 2: Sinistra, 3: Destra), il 4 è ESCLUSO
        # Usiamo il modulo % self.N per i bordi periodici (effetto Pac-Man)
        direzione = np.random.randint(0, 4)
        if direzione == 0:
            r2, c2 = (r1 - 1) % self.N, c1
        elif direzione == 1:
            r2, c2 = (r1 + 1) % self.N, c1
        elif direzione == 2:
            r2, c2 = r1, (c1 - 1) % self.N
        else:
            r2, c2 = r1, (c1 + 1) % self.N

        # se spin sono identici, lo scambio non fa nulla
        if self.lattice[r1, c1] == self.lattice[r2, c2]:
            return

        # calcolo delle energie iniziali per i due punti scelti
        E_init_1 = self.compute_E([r1, c1])
        E_init_2 = self.compute_E([r2, c2])

        # Variazione di energia per lo scambio di due PRIMI VICINI
        #(include il fattore di correzione +4 che fa sì che il legame condiviso non venga contato siccome i 2 spin sono necessariamente diversi)
        # delta_E = (-2 * E_init_1 - 2 * E_init_2) + 4
        delta_E = (-2 * E_init_1 - 2 * E_init_2) + 4*self.J


        #if energy efficient, change
        if delta_E < 0:
            self.lattice[r1, c1] *= -1
            self.lattice[r2, c2] *= -1

        #if energy innefficient, change with probability related to temperature   
        else:
            p = min(1, np.exp(-self.beta * delta_E))
            if np.random.uniform(0, 1) <= p:
                self.lattice[r1, c1] *= -1
                self.lattice[r2, c2] *= -1

        
        
        
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

        #return  -self.lattice[flip_coords[0],flip_coords[1]]*(self.lattice[flip_coords[0]-1,flip_coords[1]]+self.lattice[flip_coords[0],flip_coords[1]-1]+self.lattice[(flip_coords[0]+1)%self.N,flip_coords[1]]+self.lattice[flip_coords[0],(flip_coords[1]+1)%self.N])
        return  -self.J*self.lattice[flip_coords[0],flip_coords[1]]*(self.lattice[flip_coords[0]-1,flip_coords[1]]+self.lattice[flip_coords[0],flip_coords[1]-1]+self.lattice[(flip_coords[0]+1)%self.N,flip_coords[1]]+self.lattice[flip_coords[0],(flip_coords[1]+1)%self.N])
        



    def save_photo(self, step_name):

        """
        Salva la foto in sottocartelle ordinate per: 
        Dinamica / Configurazione / Temperatura
        """
        import os
        
        # percorso file
        cartella = os.path.join("simulation_photo", self.dynamics.__name__, self.nome_config, f"T_{self.T}")
        
        # creo cartelle se non esistono
        if not os.path.exists(cartella):
            os.makedirs(cartella)
            
        plt.clf()
        plt.imshow(self.lattice, cmap='coolwarm', vmin=-1, vmax=1)
        
        plt.title(f"{self.dynamics.__name__} | {self.nome_config} | T={self.T} | Sweep: {step_name}")
        
        # Salvo i vari step tutti dentro la cartella della temperatura
        percorso_file = f"{cartella}/step_{step_name}.png"
        plt.savefig(percorso_file)


    def run(self, wait_sweeps=100, num_tot_sweeps=1000, plot_anim=True):

        """
        Gestisce l'evoluzione temporale.
        Scatta SEMPRE all'inizio (0) e SEMPRE alla fine (num_tot_sweeps).
        """

        # lista sweep intermedi, UNICA COSA DA MODIFICARE (non serve mettere 0 e l'ultimo sweep)
        sweeps_intermedi = [200, 500, 2000, 5000]

        # FOTO INIZIALE
        if plot_anim:
            self.save_photo("00000_iniziale")

        # loop for required number of sweeps
        for i in range(num_tot_sweeps * self.sweep_size):
            if i % self.sweep_size == 0:
                sweep_corrente = i // self.sweep_size
                
                # FOTO INTERMEDIE
                if plot_anim and sweep_corrente > 0 and sweep_corrente in sweeps_intermedi:
                    self.save_photo(f"{sweep_corrente:05d}")

                if i % (self.sweep_size * 10) == 0 and i > wait_sweeps * self.sweep_size:
                    self.sweep_list.append(i / self.sweep_size)
                    self.find_magnetisation()
                    self.compute_tot_E()
                
                if i % (self.sweep_size * 100) == 0:
                    print(f"Sweep completati: {sweep_corrente}/{num_tot_sweeps}")
                    
            # run dynamics
            self.dynamics()
            
        # 3. FOTO FINALE    
        if plot_anim:
            self.save_photo(f"{num_tot_sweeps:05d}_finale")



    # def plot_lattice(self): 


    #     """
    #     Function to plot the lattice

    #     This function uses pyplot not gnuplot or other.
    #     It does not write to file to simplify the procedure.
    #     It plots directly from self.lattice.

    #     Parameters
    #     ----------

    #     Returns
    #     -------

    #     """

    #     plt.cla()
    #     plt.imshow(self.lattice, animated=True)
    #     plt.draw()
    #     plt.pause(0.01)
        
        
        
    # def run(self,wait_sweeps = 100, num_tot_sweeps = 1000, plot_anim = True):

    #     """
    #      Function to run full dynamics given simulation inputs.
    #      It uses the specified dynamics of the system and simply loops through when required.

    #      Parameters
    #      ----------
    #      wait_sweeps : int
    #          number of sweeps to wait before starting measurements

    #     num_total_sweeps : int
    #         number of total sweeps to run
    #     plot_anim : bool
    #         whether or not to animate as it runs. Saves time on long simulations to not animate.

    #      Returns
    #      -------

    #      """


    #     #loop for required number of sweeps
    #     for i in range(num_tot_sweeps*self.sweep_size):
    #         #self.N**2 is the required length of time between visualisations as in the lecture notes
    #         if i%self.sweep_size == 0:
    #             #if it is time to take a measurement
    #             if i%(self.sweep_size*10) == 0 and i > wait_sweeps*self.sweep_size:


    #                 self.sweep_list.append(i/self.sweep_size )
    #                 self.find_magnetisation()
    #                 self.compute_tot_E()
    #                 #plot animation if required
    #                 if plot_anim:
    #                     self.plot_lattice()
    #             #print every 100 sweeps to check the sim progress
    #             if i%(self.sweep_size*100) == 0:
    #                 print(i/self.sweep_size)
    #         #run dynamics
    #         self.dynamics()