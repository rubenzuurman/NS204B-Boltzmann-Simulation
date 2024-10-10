import random as rnd

import matplotlib.pyplot as plt
import numpy as np

def attempt_exchange(sites):
    """
    sites: List of sites containing the number of quanta for each site.
    """
    # Get two random numbers representing the indices of the source and destination site.
    source_site_index      = rnd.randint(0, len(sites) - 1)
    destination_site_index = rnd.randint(0, len(sites) - 1)
    
    # Return if the source site has zero quanta.
    if sites[source_site_index] == 0:
        return
    
    # Remove quantum from the source and add it to the destination.
    sites[source_site_index] -= 1
    sites[destination_site_index] += 1

def execute_sweep(sites, num_exchanges):
    """
    sites:         List of sites containing the number of quanta for each site.
    num_exchanges: The number of exchanges to attempt.
    """
    # Execute num_exchanges attempts at an exchange.
    for _ in range(num_exchanges):
        attempt_exchange(sites)

def simulate(num_sites, num_quanta, num_sweeps):
    """
    num_sites (M):     Number of sites to use in the simulation.
    num_quanta (N):    Number of energy quanta to distribute over the sites.
    num_sweeps (N_sw): Total number of sweeps to perform. Every sweep makes num_sites attempts at an exchange.
    """
    # Create initial list of sites and distribute the available quanta over the sites.
    sites = [0] * num_sites
    for i in range(num_quanta):
        sites[i % num_sites] += 1
    
    # Initialize dictionary containing the bincount results for every sweep.
    sweep_bincounts = {}
    
    # Execute num_sweeps sweeps of num_sites swaps each.
    for sweep_number in range(num_sweeps):
        # Execute the sweep. No return value is needed since sites is passed by reference.
        execute_sweep(sites=sites, num_exchanges=num_sites)
        
        # Add bincount of this sweep of bincount dictionary.
        sweep_bincounts[sweep_number] = np.bincount(sites)
    
    # Return dictionary of bincounts per sweep.
    return sweep_bincounts

def analyse_results(sweep_bincounts, num_sites):
    """
    sweep_bincounts: Dictionary produced by the simulate() function.
    num_sites (M):   Number of sites to use in the simulation.
    """
    # Extract number of sweeps from bincounts dictionary.
    num_sweeps = len(sweep_bincounts)
    
    # Calculate maximum number of quanta at any site at the end of any sweep.
    max_quanta = max([len(bincount) for bincount in sweep_bincounts.values()])
    
    # Calculate average number of sites with n quanta over all sweeps (this is <N_n>).
    average_num_sites_with_n_quanta = {n: 0 for n in range(max_quanta)}
    for bincount in sweep_bincounts.values():
        for n, number_of_quanta in enumerate(bincount):
            average_num_sites_with_n_quanta[n] += number_of_quanta / num_sweeps
    
    # Calculate average number of sites with n quanta squared over all sweeps (this is <(N_n)^2>).
    average_num_sites_with_n_quanta_squared = {n: 0 for n in range(max_quanta)}
    for bincount in sweep_bincounts.values():
        for n, number_of_quanta in enumerate(bincount):
            average_num_sites_with_n_quanta_squared[n] += (number_of_quanta * number_of_quanta) / num_sweeps
    
    # Calculate probability for a site to have n quanta at the end of any sweep (this is P_n).
    prob_to_have_n_quanta = {n: avg_quanta / num_sites for n, avg_quanta in average_num_sites_with_n_quanta.items()}
    
    # Generate elementary plot of the probabilities.
    fig, ax = plt.subplots(nrows=1, ncols=1)
    
    ax.bar(prob_to_have_n_quanta.keys(), prob_to_have_n_quanta.values())
    
    fig.savefig("test.png")

def main():
    # Set random seed.
    RANDOM_SEED = 1234
    rnd.seed(RANDOM_SEED)
    
    # Simulation parameters (will be passed around in the functions).
    num_sites  = 5
    num_quanta = 10
    num_sweeps = 500
    
    # Execute simulation.
    sweep_bincounts = simulate(num_sites=num_sites, num_quanta=num_quanta, num_sweeps=num_sweeps)
    
    # Analyse simulation results.
    analyse_results(sweep_bincounts=sweep_bincounts, num_sites=num_sites)

if __name__ == "__main__":
    main()
