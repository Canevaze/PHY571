import numpy as np
import numpy.random as rnd
import itertools
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.animation import FuncAnimation
from tqdm import tqdm

import sys
sys.path.append('D:\Github\PHY571')

from our_library import Simple_model as sm
        
            
#create swarm of birds and do an animation of the evolution
avg = 10
nb_pts = 40

def retrieve_param_eta(eta,it,N,V, L, R):
    swarm = sm.Swarm(L, N, V, eta, R)
    swarm.initialize()

    for i in range(it):
        swarm.evolve()
    
    return swarm.get_swarm_mean_velocity()

def calculate_avg_and_uncertainty(data):
    avg = np.mean(data)
    uncertainty = np.std(data, ddof=1) / np.sqrt(len(data))  # ddof=1 for sample standard deviation
    return avg, uncertainty

def calculate_r_squared(y, y_pred):
    ss_residual = np.sum((y - y_pred)**2)
    ss_total = np.sum((y - np.mean(y))**2)
    r_squared = 1 - (ss_residual / ss_total)
    return r_squared

def plot(it=200,V=0.03, R=1):
    eta = 0.75
    N = np.linspace(100,400,nb_pts).astype(int)
    L = np.sqrt(N/40)*6.2
    # v = []
    # for e in tqdm(eta):
    #     average = []
    #     for i in range(avg):
    #         average.append(retrieve_param_eta(e,it,N,V, L, R))
    #     v.append(np.mean(average))

    v_avg = []
    v_uncertainty = []

    for n, l in tqdm(zip(N,L)):
        velocities = []
        for i in range(avg):
            velocities.append(retrieve_param_eta(eta, it, n, V, l, R))

        average, uncertainty = calculate_avg_and_uncertainty(velocities)
        v_avg.append(average)
        v_uncertainty.append(uncertainty)

    # Linear regression
    coeffs = np.polyfit(N, v_uncertainty, 1)
    linear_fit = np.polyval(coeffs, N)

    #save the values
    np.save(f'N_it_{it}_V_{V}_eta_{eta}_R_{R}_avg_{avg}_nb_pts_{nb_pts}_STD_Class1.npy', N)
    np.save(f'v_uncertainty_it_{it}_V_{V}_eta_{eta}_R_{R}_avg_{avg}_nb_pts_{nb_pts}_STD_Class1.npy', v_uncertainty)

    # Calculate R-squared
    r_squared = calculate_r_squared(v_uncertainty, linear_fit)

    #plot it all
    
    plt.scatter(N,v_uncertainty)
    plt.plot(N, linear_fit, color='red', label="Linear Fit")
    plt.xlabel("N")
    plt.ylabel("v_uncertainty")

    # Annotate linear regression parameters next to the red line
    slope, intercept = coeffs
    annotation_text = f'Slope: {slope:.4f}\nIntercept: {intercept:.4f}\n$R^2$: {r_squared:.4f}'
    plt.annotate(annotation_text, xy=(N[-1], linear_fit[-1]), color='red', fontsize=8)

    plt.legend()
    #display the parameter values
    plt.title("V = " + str(V) + ", eta = " + str(eta) + ", R = " + str(R) + ", avg = " + str(avg))
    plt.show()



plot()