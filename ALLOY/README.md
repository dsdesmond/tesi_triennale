# Ising Model Simulation

This repository hold an Ising model simulation of Glauber and Kawisaki dynamics aswell as example results.

Example results are for a simulation of a (50x50) lattice, with Magnetisation, Energy, Heat capacity and Susceptability measurements every 10 sweeps for 1000 measurements.

# Running Ising simulations

Glauber.py and Kawisaki.py can be called from the commandline, with both requiring arguments that determine the parameters of their respective Ising simulations.

Visualisation is also optional, chosen through a commandline argument. It is normally switched off during long computation to save time, and can be turned on to qualitatively study behaviour.

Results are saved to respective files "GLA_ALL" for glauber and "KAW_ALL" for Kawasaki, with the form of the document being a seperate column for each temperature evaluated, thus rows are ordered Temperature, Magnetisation, Energy, Heat capacity, Susceptability, error in heat capactity.

Example usage from terminal:
To run Kawisaki at T = 1, N = 50, 10000 total sweeps, no animation

>python3 Kawisaki.py, 1, 50, 10000, False


To run Kawisaki at T_low = 1, N = 50, 10000 total sweeps, animation, T_high = 3 and stepsize = 0.1 (ie a range of Temperatures between 1 and 3 at 0.1 intervals will be evaluated)

> python3 Kawisaki.py, 1, 50, 10000, True, 3, 0.1


The same could be done with Glauber.py


To ru  analysis of Glauber and Kawisaki simulation where the results are stored in "GLA_ALL" and "KAW_ALL" in a seperate directory names "G_results" as they are currently:
> python3 analysis.py, "GLA_ALL", "KAW_ALL"    


# Analysing results

analysis.py can also be called straight from the command line, with the arguments being the filenames of the results (output from Glauber.py and Kawisaki.py)
It will output graphs detailing how magnetisation, energy, Heat capacity and susceptability vary with temperature in the range 0K to 3K (from the baseline experiment detailed above).
Simple fits have been employed to estimate the critical temperatures for Magnetisation, energy and heat capacity graphs, with results ranging between 2.25-2.29K as one would expect from literature.

Note that for low quality results (such as few different temperatures assessed or too few sweeps) the fitting algorithm will likely throw an error.



# Code structure and usage details

When initialised, the Lattice class (in Lattice.py) creates a Lattice object capable of evolution through Glauber or Kawisaki dynamics. Once initialised the run() function can be called in order to loop through the appropriate dynamics function as many times as required by the input. throughout the function, magnetisation and Energy are periodically (every 10 sweeps, where a sweep is equal to the number of lattice points) calculated and stores, aswell as an optional visualisation. At the end of the simulation, the heat capacity and susceptability is calculated. The error in heat capacity is calculated using the jacknife method

The Glauber.py and Kawisaki.py files employ the Lattice class multiple times (as many as required for the job at hand) before saving the results to file as above

code was written and tested on python 3.8 using pycharm IDE. 
