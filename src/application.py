import random as rnd

import matplotlib.pyplot as plt
import numpy as np

# Set font size for plots.
plt.rcParams.update({'font.size': 20})

BOLTZMANN_CONSTANT = 1.38e-23

PRIMARY_COLOR = "royalblue"
SECONDARY_COLOR = "indianred"

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

def analyse_results(simulation_data):
    """
    simulation_data: Dictionary containing for each simulation the following keys:
        num_sites (M):     Number of sites.
        num_quanta (N):    Number of quanta distributed over all the sites in the simulation.
        num_sweeps (N_sw): Total number of sweeps to perform.
        sweep_bincounts:   Dictionary containing the result of bincount every sweep with sweep number as the keys.
    Returns a dictionary containing the simulation_data keys and values, together with the following added keys and values (spacing for readability):
        avg_Nn (<N_n>):    Dictionary of length n containing the average number of sites with n quanta over all sweeps (aka the entire simulation). For example, the average number of sites with 1 quanta is avg_Nn[1].
        
        avg_Nn_squared (<(N_n)^2>): Dictionary of length n containing the average number of sites squared with n quanta over all sweeps (aka the entire simulation).
        
        num_quanta_probabilities_for_any_site (P_n): Dictionary of length n containing the probability, on average, that any site has n quanta. This is essentially just <N_n> normalized with the number of sites.
        
        sweep_entropies: Entropy calculated for every sweep in the simulation.
    """
    for simulation_index, params in simulation_data.items():
        # Define sweep_bincounts for readability.
        sweep_bincounts = params["sweep_bincounts"]
        
        # Calculate maximum number of quanta at any site at the end of any sweep.
        max_quanta = max([len(bincount) for bincount in sweep_bincounts.values()])
        
        # Calculate average number of sites with n quanta over all sweeps (this is <N_n>).
        average_num_sites_with_n_quanta = {n: 0 for n in range(max_quanta)}
        for bincount in sweep_bincounts.values():
            for n, number_of_quanta in enumerate(bincount):
                average_num_sites_with_n_quanta[n] += number_of_quanta / params["num_sweeps"]
        params["avg_Nn"] = average_num_sites_with_n_quanta
        
        # Calculate average number of sites with n quanta squared over all sweeps (this is <(N_n)^2>).
        average_num_sites_with_n_quanta_squared = {n: 0 for n in range(max_quanta)}
        for bincount in sweep_bincounts.values():
            for n, number_of_quanta in enumerate(bincount):
                average_num_sites_with_n_quanta_squared[n] += (number_of_quanta * number_of_quanta) / params["num_sweeps"]
        params["avg_Nn_squared"] = average_num_sites_with_n_quanta_squared
        
        # Calculate probability for a site to have n quanta at the end of any sweep (this is P_n).
        prob_to_have_n_quanta = {n: avg_quanta / params["num_sites"] for n, avg_quanta in average_num_sites_with_n_quanta.items()}
        params["num_quanta_probabilities_for_any_site"] = prob_to_have_n_quanta
        
        # Calculate entropy over 'time'. Initialize list to store entropy for every sweep.
        entropy_over_sweeps = []
        
        # Calculate probability to find n quanta at any site for every bincount (aka every sweep).
        for sweep_number, bincount in sweep_bincounts.items():
            # Normalize bincount to get probability for any site to have some number of quanta.
            # Form is [prob to have 0 quanta, prob to have 1 quanta, ...].
            probabilities = [entry / sum(bincount) for entry in bincount]
            
            # Calculate entropy for this sweep using Gibbs' expression for entropy.
            # Skip probabilities of zero, since lim(x->0+)xlnx=0, and ln(0+)->-inf.
            entropy = -BOLTZMANN_CONSTANT * sum([p * np.log(p) for p in probabilities if p != 0])
            
            # Add entropy to entropy list for this simulation.
            entropy_over_sweeps.append(entropy)
        
        # Add entropy data to simulation data.
        params["sweep_entropies"] = {index: entropy for index, entropy in enumerate(entropy_over_sweeps)}
    
    # Return simulation data dictionary.
    return simulation_data

def generate_entropy_over_time_plots():
    # Simulation parameters (will be passed around in the functions).
    num_sweeps = 200
    
    # Initialize simulation parameters.
    simulation_parameters = [
        {"num_sites": 10, "num_quanta": 10, "num_sweeps": num_sweeps}, 
        {"num_sites": 500, "num_quanta": 10, "num_sweeps": num_sweeps}, 
        {"num_sites": 10, "num_quanta": 500, "num_sweeps": num_sweeps}, 
        {"num_sites": 500, "num_quanta": 500, "num_sweeps": num_sweeps}, 
    ]
    
    # Copy simulation parameters to another dictionary. This dictionary will also contain the sweep_bincounts and entropy for every simulation.
    simulation_data = {index: {k: v for k, v in d.items()} for index, d in enumerate(simulation_parameters)}
    
    # Execute simulations.
    for index, params in simulation_data.items():
        sweep_bincounts = simulate(num_sites=params["num_sites"], num_quanta=params["num_quanta"], num_sweeps=params["num_sweeps"])
        simulation_data[index]["sweep_bincounts"] = sweep_bincounts
    
    # Analyse simulation results. This function adds analysis results to the keys already present in the simulation_data dictionary.
    simulation_data = analyse_results(simulation_data)
    
    # Plot entropy evolution of several simulation parameters (N and M).
    fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(12, 12))
    
    for index, params in simulation_data.items():
        # Plot data.
        axs[index % 2][index // 2].plot(params["sweep_entropies"].keys(), [x / BOLTZMANN_CONSTANT for x in params["sweep_entropies"].values()], label=f"M={params["num_sites"]} N={params["num_quanta"]}", color=PRIMARY_COLOR)
        
        # Set axis labels.
        axs[index % 2][index // 2].set_xlabel("Sweep Number")
        axs[index % 2][index // 2].set_ylabel("Entropy $\\frac{S}{k_B}$")
        
        # Set axis limits.
        #min_ylim = float(f"{min(params["sweep_entropies"].values()):.1g}")
        #max_ylim = float(f"{max(params["sweep_entropies"].values()):.1g}")
        axs[index % 2][index // 2].set_xlim([0, 200])
        
        # Add minimal legend indicating simulation parameters.
        legend = axs[index % 2][index // 2].legend(handlelength=0, handletextpad=0, fancybox=True, loc="lower right")
        for item in legend.legend_handles:
            item.set_visible(False)
    
    fig.tight_layout()
    fig.savefig("entropy_over_time.png")

def generate_boltzmann_bar_chart():
    # Simulation parameters (will be passed around in the functions).
    num_sweeps = 200
    
    # Initialize simulation parameters.
    simulation_parameters = [
        {"num_sites": 20, "num_quanta": 20, "num_sweeps": num_sweeps}
    ]
    
    # Copy simulation parameters to another dictionary. This dictionary will also contain the sweep_bincounts and entropy for every simulation.
    simulation_data = {index: {k: v for k, v in d.items()} for index, d in enumerate(simulation_parameters)}
    
    # Execute simulations.
    for index, params in simulation_data.items():
        sweep_bincounts = simulate(num_sites=params["num_sites"], num_quanta=params["num_quanta"], num_sweeps=params["num_sweeps"])
        simulation_data[index]["sweep_bincounts"] = sweep_bincounts
    
    # Analyse simulation results. This function adds analysis results to the keys already present in the simulation_data dictionary.
    simulation_data = analyse_results(simulation_data)
    
    num_squanta = simulation_data[0]["num_quanta"]
    
    # Plot bar chart of last sweep against Boltzmann bar chart on the left, and entropy over sweeps on the right to verify thermalization.
    fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(12, 6))
    
    # Calculate probability distribution from bincounts over the last 100 sweeps. Aka calculate the total number of occurrences of 0 quanta at any site at any sweep, the total number of occurrences of 1 quanta at any site at any sweep, etc. Then normalize this distribution to make it a valid probability distribution.
    max_quanta_at_any_site = max([len(sweep_bincounts) for index, sweep_bincounts in simulation_data[0]["sweep_bincounts"].items() if index >= 100])
    bincounts_total = [0 for _ in range(max_quanta_at_any_site)]
    for index, sweep_bincount in simulation_data[0]["sweep_bincounts"].items():
        if index < 100:
            continue
        
        for num_quanta, occurrences in enumerate(sweep_bincount):
            bincounts_total[num_quanta] += occurrences
    
    # Normalize the distribution.
    simulation_probability_dist = list(bincounts_total / sum(bincounts_total))
    
    # Appends zeros to end of prob dist to make sure all possible numbers of quanta are present.
    simulation_probability_dist.extend([0] * (num_squanta + 1 - len(simulation_probability_dist)))
    
    # Calculate the theoretical Boltzmann distribution from the Boltzmann factor for every number of quanta and the partition function. (I don't know what beta is so I'm setting it to 1 in case it matters. Maybe we fit beta?)
    beta = 0.66
    partition_function = sum([np.exp(-beta * n) for n in range(num_squanta + 1)])
    theoretical_probability_dist = [np.exp(-beta * n) / partition_function for n in range(num_squanta + 1)]
    
    #print(simulation_data[0]["sweep_bincounts"][num_sweeps - 1])
    bar_chart_x = list(range(0, num_squanta + 1))
    #bar_chart_height = simulation_data[0]["sweep_bincounts"][num_sweeps - 1] / sum(simulation_data[0]["sweep_bincounts"][num_sweeps - 1])
    axs[0].bar([x - 0.225 for x in bar_chart_x], theoretical_probability_dist, width=0.4, label="Theory", color=PRIMARY_COLOR)
    axs[0].bar([x + 0.225 for x in bar_chart_x], simulation_probability_dist, width=0.4, label="Simulation", color=SECONDARY_COLOR)
    axs[0].set_xlabel("Number of Quanta")
    axs[0].set_ylabel("Probability")
    axs[0].legend()
    
    axs[1].plot(simulation_data[0]["sweep_entropies"].keys(), [x / BOLTZMANN_CONSTANT for x in simulation_data[0]["sweep_entropies"].values()], color=PRIMARY_COLOR)
    axs[1].set_xlabel("Sweep Number")
    axs[1].set_ylabel("Entropy $\\frac{S}{k_B}$")
    axs[1].set_xlim([0, 200])
    
    fig.tight_layout()
    fig.savefig("boltzmann_bar_chart.png")

def main():
    # Set random seed.
    RANDOM_SEED = 1234
    rnd.seed(RANDOM_SEED)
    
    # Generate entropy over time plots.
    generate_entropy_over_time_plots()
    
    # Generate bar chart of Boltzmann distribution vs simulation result.
    generate_boltzmann_bar_chart()

if __name__ == "__main__":
    main()
