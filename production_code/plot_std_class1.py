import numpy as np
import numpy.random as rnd
import itertools
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.animation import FuncAnimation
from tqdm import tqdm

import sys
sys.path.append('.')

from our_library import Simple_model as sm
        
            
#create swarm of birds and do an animation of the evolution
avg = 5
nb_pts = 50

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

def plot(it=200,V=0.03, R=1):
    eta = 0.5
    N = np.linspace(100,500,nb_pts)
    L = np.sqrt(N/40)*6.2
    
    # v = []
    # for e in tqdm(eta):
    #     average = []
    #     for i in range(avg):
    #         average.append(retrieve_param_eta(e,it,N,V, L, R))
    #     v.append(np.mean(average))

    v_avg = []
    v_uncertainty = []

    for n in tqdm(N):
        velocities = []
        for n, l in zip(N,L):
            velocities.append(retrieve_param_eta(eta, it, n, V, l, R))

        average, uncertainty = calculate_avg_and_uncertainty(velocities)
        v_avg.append(average)
        v_uncertainty.append(uncertainty)

    #save the values
    np.save(f'N_it_{it}_V_{V}_eta_{eta}_R_{R}_avg_{avg}_nb_pts_{nb_pts}_STD_Class1.npy', N)
    np.save(f'v_uncertainty_it_{it}_V_{V}_eta_{eta}_R_{R}_avg_{avg}_nb_pts_{nb_pts}_STD_Class1.npy', v_uncertainty)


    #plot it all
    
    plt.scatter(N,v_uncertainty)
    plt.xlabel("N")
    plt.ylabel("v_uncertainty")
    #display the parameter values
    plt.title("V = " + str(V) + ", eta = " + str(eta) + ", R = " + str(R))
    plt.show()



plot()