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

def plot_2_a(it=200,N=40,V=0.03, L=6.2, R=1):
    eta = np.linspace(0.0,5,nb_pts)
    # v = []
    # for e in tqdm(eta):
    #     average = []
    #     for i in range(avg):
    #         average.append(retrieve_param_eta(e,it,N,V, L, R))
    #     v.append(np.mean(average))

    v_avg = []
    v_uncertainty = []

    for e in tqdm(eta):
        velocities = []
        for i in range(avg):
            velocities.append(retrieve_param_eta(e, it, N, V, L, R))

        average, uncertainty = calculate_avg_and_uncertainty(velocities)
        v_avg.append(average)
        v_uncertainty.append(uncertainty)

    #save the values
    np.save(f'v_avg_it_{it}_N_{N}_V_{V}_L_{L}_R_{R}_avg_{avg}_nb_pts_{nb_pts}_Class1.npy', v_avg)
    np.save(f'eta_it_{it}_N_{N}_V_{V}_L_{L}_R_{R}_avg_{avg}_nb_pts_{nb_pts}_Class1.npy', eta)
    np.save(f'v_uncertainty_it_{it}_N_{N}_V_{V}_L_{L}_R_{R}_avg_{avg}_nb_pts_{nb_pts}_Class2.npy', v_uncertainty)


    #plot it all
    
    plt.errorbar(eta, v_avg, yerr=v_uncertainty, fmt='o', label='Average with Uncertainty')
    plt.scatter(eta,v)
    plt.xlabel("eta")
    plt.ylabel("v")
    #display the parameter values
    plt.title("N = " + str(N) + ", V = " + str(V) + ", L = " + str(L) + ", R = " + str(R))
    plt.show()

def plot_2_b(it=200,N=40,V=0.03, eta=2 , R=1):
    L = np.linspace(2,25,nb_pts)
    #make density equal to N/L^2
    density = N/(L**2)


    v = []
    for l in tqdm(L):
        average = []
        for i in range(avg):
            average.append(retrieve_param_eta(eta,it,N,V, l, R))
        v.append(np.mean(average))

    #save the values for Class1
    np.save(f'v_it_{it}_N_{N}_V_{V}_eta_{eta}_R_{R}_avg_{avg}_nb_pts_{nb_pts}_Class1.npy', v)
    np.save(f'density_it_{it}_N_{N}_V_{V}_eta_{eta}_R_{R}_avg_{avg}_nb_pts_{nb_pts}_Class1.npy', density)

    #plot it all

    plt.scatter(density,v)
    plt.xlabel("density")
    plt.ylabel("v")
    #display the parameter values
    plt.title("N = " + str(N) + ", V = " + str(V) + ", eta = " + str(eta) + ", R = " + str(R))
    # plt.title("N = " + str(N) + ", V = " + str(V) + ", eta = " + str(eta) + ", R1 = " + str(interaction_radius_1)
    #            + ", R2 = " + str(interaction_radius_2) + ", R3 = " + str(interaction_radius_3))
    plt.show()

plot_2_a()