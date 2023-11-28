import numpy as np
import numpy.random as rnd
import itertools
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.animation import FuncAnimation
from tqdm import tqdm

from Class import *
        
            
#create swarm of birds and do an animation of the evolution

v = np.load('name_file_v.npy')
eta = np.load('name_file_eta.npy')

import numpy as np
import matplotlib.pyplot as plt

def plot_v_vs_eta_subplot(v, eta, eta_critical_values):
    num_plots = len(eta_critical_values)
    
    # Create subplots
    fig, axes = plt.subplots(num_plots, 1, figsize=(8, 4 * num_plots), sharex=True)

    for i, eta_critical in enumerate(eta_critical_values):
        # Calculate the x-values: (eta_critical - eta) / eta_critical
        x_values = (eta_critical - eta) / eta_critical

        # Plot v in function of (eta_critical - eta) / eta_critical in each subplot
        axes[i].plot(x_values, v, label=f'eta_critical = {eta_critical}')
        
        # Add labels and title to the first subplot
        if i == 0:
            axes[i].set_xlabel(f'({eta_critical} - eta) / {eta_critical}')
            axes[i].set_ylabel('v')
            axes[i].set_title('v vs. (eta_critical - eta) / eta_critical')

        # Add a legend to each subplot
        axes[i].legend()

    # Adjust layout
    plt.tight_layout()

    # Show the plot
    plt.show()

# Example usage:
# Generate some example data
eta_critical_values = [2, 3, 4, 5]

# Call the function to plot the data with different eta_critical values
plot_v_vs_eta_subplot(v, eta, eta_critical_values)

