import numpy as np
import numpy.random as rnd
import itertools
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.animation import FuncAnimation
from tqdm import tqdm
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error

import sys
sys.path.append('.')

from our_library import Complex_model as cm
        
            
#create swarm of birds and do an animation of the evolution

#v = np.load('https://github.com/Canevaze/PHY571/blob/21a57986093a7471c2b700f9360e4145c796fa35/v_it_200_N_100_V_0.03_L_5_R_1_avg_20_nb_pts_50.npy')
#eta = np.load('https://github.com/Canevaze/PHY571/blob/21a57986093a7471c2b700f9360e4145c796fa35/eta_it_200_N_100_V_0.03_L_5_R_1_avg_20_nb_pts_50.npy')

v = np.load('v_it_100_N_40_V_0.03_L_6.2_R_1_avg_50_nb_pts_50.npy')
eta = np.load('eta_it_100_N_40_V_0.03_L_6.2_R_1_avg_50_nb_pts_50.npy')

import numpy as np
import matplotlib.pyplot as plt

def plot_v_vs_eta_subplot(v, eta, eta_critical_values):
    num_plots = 1
    
    # Create subplots
    fig, axes = plt.subplots(num_plots, 1, figsize=(8, 4 * num_plots), sharex=True)

    # Ensure eta_critical_values is iterable
    if not isinstance(eta_critical_values, (list, np.ndarray)):
        eta_critical_values = [eta_critical_values]

    for i, eta_critical in enumerate(eta_critical_values):
        # Calculate the x-values: (eta_critical - eta) / eta_critical
        x_values = (eta_critical - eta) / eta_critical

        # Plot v in function of (eta_critical - eta) / eta_critical in each subplot
        axes.scatter(x_values, v, label=f'eta_critical = {eta_critical}')
        
        axes.set_ylim(0.6, 1)
        axes.set_yscale('log')
        axes.set_xlim(0.01, 1)
        axes.set_xscale('log')
        
        # Add labels and title to the first subplot
        if i == 0:
            axes.set_xlabel(f'({eta_critical} - eta) / {eta_critical}')
            axes.set_ylabel('v')
            axes.set_title('v vs. (eta_critical - eta) / eta_critical')

        # Add a legend to each subplot
        axes.legend()
        axes.grid(True)

    # Adjust layout
    plt.tight_layout()
    plt.title(f'v vs. (eta_critical - eta) / eta_critical with eta_critical = {eta_critical}')
    # Show the plot
    plt.show()




# i want to build a function that finds the parameter eta_critical for which the plot is the straightest

from sklearn.metrics import mean_squared_error

def find_eta_critical():
    eta_critical_values = np.linspace(0.1, 5, 100)
    v = np.load('v_it_100_N_40_V_0.03_L_6.2_R_1_avg_50_nb_pts_50.npy')
    eta = np.load('eta_it_100_N_40_V_0.03_L_6.2_R_1_avg_50_nb_pts_50.npy')
    
    mse_values = []

    for eta_critical in eta_critical_values:
        x_values = (eta_critical - eta) / eta_critical

        # Use log-transformed v for straightness on a log scale
        log_v = np.log(v)

        # Filter out NaN and infinity values
        valid_indices = np.isfinite(log_v) & np.isfinite(x_values)
        log_v_valid = log_v[valid_indices]
        x_values_valid = x_values[valid_indices]

        # Calculate the mean squared error between log_v and x_values
        mse = mean_squared_error(log_v_valid, x_values_valid)
        mse_values.append(mse)

    return eta_critical_values[np.argmin(mse_values)]

def find_eta_critical_straightest():
    eta_critical_values = np.linspace(0.1, 5, 100)
    v = np.load('v_it_100_N_40_V_0.03_L_6.2_R_1_avg_50_nb_pts_50.npy')
    eta = np.load('eta_it_100_N_40_V_0.03_L_6.2_R_1_avg_50_nb_pts_50.npy')
    
    min_error = float('inf')
    best_eta_critical = None

    for eta_critical in eta_critical_values:
        x_values = (eta_critical - eta) / eta_critical

        # Filter out NaN and infinity values
        valid_indices = np.isfinite(v) & np.isfinite(x_values)
        v_valid = v[valid_indices]
        x_values_valid = x_values[valid_indices]

        # Reshape for Linear Regression
        x_values_valid = x_values_valid.reshape(-1, 1)

        # Fit a line using Linear Regression
        model = LinearRegression().fit(x_values_valid, np.log(v_valid))

        # Get the mean squared error between the data and the fitted line
        y_pred = model.predict(x_values_valid)
        error = mean_squared_error(np.log(v_valid), y_pred)

        # Update the best eta_critical if the error is lower
        if error < min_error:
            min_error = error
            best_eta_critical = eta_critical

    return best_eta_critical

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

# Generate some example data
eta_critical_values = find_eta_critical()
print(f'eta_critical_value = {eta_critical_values}')

# Call the function to plot the data with different eta_critical values
plot_v_vs_eta_subplot(v, eta, eta_critical_values)

#plot v vs eta
plt.scatter(eta, v)
plt.xlabel("eta")
plt.ylabel("v")
plt.title("v vs eta")
plt.show()