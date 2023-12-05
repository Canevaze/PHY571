import numpy as np
import numpy.random as rnd
import itertools
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.animation import FuncAnimation
from tqdm import tqdm
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error

from Class import *
        
            
#create swarm of birds and do an animation of the evolution

#v = np.load('https://github.com/Canevaze/PHY571/blob/21a57986093a7471c2b700f9360e4145c796fa35/v_it_200_N_100_V_0.03_L_5_R_1_avg_20_nb_pts_50.npy')
#eta = np.load('https://github.com/Canevaze/PHY571/blob/21a57986093a7471c2b700f9360e4145c796fa35/eta_it_200_N_100_V_0.03_L_5_R_1_avg_20_nb_pts_50.npy')

v= np.load('v_it_200_N_100_V_0.03_L_5_R_1_avg_20_nb_pts_50.npy')
eta = np.load('eta_it_200_N_100_V_0.03_L_5_R_1_avg_20_nb_pts_50.npy')

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
        axes[i].scatter(x_values, v, label=f'eta_critical = {eta_critical}')
        
        axes[i].set_ylim(0.6, 1)
        axes[i].set_yscale('log')
        axes[i].set_xlim(0.01, 1)
        axes[i].set_xscale('log')
        # Add labels and title to the first subplot
        if i == 0:
            axes[i].set_xlabel(f'({eta_critical} - eta) / {eta_critical}')
            axes[i].set_ylabel('v')
            axes[i].set_title('v vs. (eta_critical - eta) / eta_critical')

        # Add a legend to each subplot
        axes[i].legend()
        axes[i].grid(True)


    # Adjust layout
    plt.tight_layout()

    # Show the plot
    plt.show()

# Example usage:
# Generate some example data
eta_critical_values = [3,4]

# Call the function to plot the data with different eta_critical values
plot_v_vs_eta_subplot(v, eta, eta_critical_values)

# i want to build a function that finds the parameter eta_critical for which the plot is the straightest

def find_eta_critical():
    eta_critical_values = np.linspace(0.1, 5, 100)
    v = np.load('v_it_200_N_100_V_0.03_L_5_R_1_avg_20_nb_pts_50.npy')
    eta = np.load('eta_it_200_N_100_V_0.03_L_5_R_1_avg_20_nb_pts_50.npy')
    
    v_eta_critical = []
    for eta_critical in eta_critical_values:
        x_values = (eta_critical - eta) / eta_critical
        v_eta_critical.append(np.polyfit(x_values, v, 1)[0])
    return eta_critical_values[np.argmax(v_eta_critical)]

def find_eta_critical_v2(eta_range, generate_plot_data):
    min_error = float('inf')
    best_eta = None

    for eta in eta_range:
        x, y = generate_plot_data(eta)
        model = LinearRegression().fit(x.reshape(-1, 1), y)
        y_pred = model.predict(x.reshape(-1, 1))
        error = mean_squared_error(y, y_pred)

        if error < min_error:
            min_error = error
            best_eta = eta

    return best_eta

