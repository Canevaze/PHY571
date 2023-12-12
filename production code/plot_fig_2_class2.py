import numpy as np
import numpy.random as rnd
import itertools
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.animation import FuncAnimation
from tqdm import tqdm

import sys
sys.path.append('.')

from our_library import Complex_model as cm
        
            
#create swarm of birds and do an animation of the evolution
avg = 5
nb_pts = 50

def retrieve_param_eta(eta,it,N,V, L, interaction_radius_1, interaction_radius_2, interaction_radius_3):
    swarm = cm.Swarm(L, N, V, eta, interaction_radius_1, interaction_radius_2, interaction_radius_3)
    swarm.initialize()

    for i in range(it):
        swarm.evolve()
    
    return swarm.get_swarm_mean_velocity()

def calculate_avg_and_uncertainty(data):
    avg = np.mean(data)
    uncertainty = np.std(data, ddof=1) / np.sqrt(len(data))  # ddof=1 for sample standard deviation
    return avg, uncertainty

# def plot_2_aa(it=100,N=40,V=0.03, L=6.2, interaction_radius_1=1, interaction_radius_2=2, interaction_radius_3=3):
#     eta = np.linspace(0.0,5,nb_pts)
#     v = []
#     for e in tqdm(eta):
#         average = []
#         for i in range(avg):
#             average.append(retrieve_param_eta(e,it,N,V, L, interaction_radius_1, interaction_radius_2, interaction_radius_3))
#         v.append(np.mean(average))

#     #save the values
#     np.save(f'v_it_{it}_N_{N}_V_{V}_L_{L}_R1_{interaction_radius_1}_R2_{interaction_radius_2}_R3_{interaction_radius_3}_avg_{avg}_nb_pts_{nb_pts}_Class2.npy', v)
#     np.save(f'eta_it_{it}_N_{N}_V_{V}_L_{L}_R1_{interaction_radius_1}_R2_{interaction_radius_2}_R3_{interaction_radius_3}_avg_{avg}_nb_pts_{nb_pts}_Class2.npy', eta)


#     #plot it all

#     plt.scatter(eta,v)
#     plt.xlabel("eta")
#     plt.ylabel("v")
#     #display the parameter values
#     plt.title("N = " + str(N) + ", V = " + str(V) + ", L = " + str(L) + ", R1 = " + str(interaction_radius_1)
#                + ", R2 = " + str(interaction_radius_2) + ", R3 = " + str(interaction_radius_3))
#     plt.show()

def plot_2_a(it=200, N=40, V=0.03, L=6.2, interaction_radius_1=1, interaction_radius_2=2, interaction_radius_3=3):
    eta = np.linspace(0.0, 5, nb_pts)
    v_avg = []
    v_uncertainty = []

    for e in tqdm(eta):
        velocities = []
        for i in range(avg):
            velocities.append(retrieve_param_eta(e, it, N, V, L, interaction_radius_1, interaction_radius_2, interaction_radius_3))

        avgrage, uncertainty = calculate_avg_and_uncertainty(velocities)
        v_avg.append(avgrage)
        v_uncertainty.append(uncertainty)

    # Save the values
    np.save(f'v_avg_it_{it}_N_{N}_V_{V}_L_{L}_R1_{interaction_radius_1}_R2_{interaction_radius_2}_R3_{interaction_radius_3}_avg_{avg}_nb_pts_{nb_pts}_Class2.npy', v_avg)
    np.save(f'v_uncertainty_it_{it}_N_{N}_V_{V}_L_{L}_R1_{interaction_radius_1}_R2_{interaction_radius_2}_R3_{interaction_radius_3}_avg_{avg}_nb_pts_{nb_pts}_Class2.npy', v_uncertainty)
    np.save(f'eta_it_{it}_N_{N}_V_{V}_L_{L}_R1_{interaction_radius_1}_R2_{interaction_radius_2}_R3_{interaction_radius_3}_avg_{avg}_nb_pts_{nb_pts}_Class2.npy', eta)

    # Plot the data
    plt.errorbar(eta, v_avg, yerr=v_uncertainty, fmt='o', label='Average with Uncertainty')
    plt.xlabel("eta")
    plt.ylabel("v")
    plt.title(f"N = {N}, V = {V}, L = {L}, R1 = {interaction_radius_1}, R2 = {interaction_radius_2}, R3 = {interaction_radius_3}")
    plt.legend()
    plt.show()

def plot_2_b(it=200,N=40,V=0.03, eta=2 , interaction_radius_1=0.5, interaction_radius_2=1, interaction_radius_3=2):
    L = np.linspace(2,25,nb_pts)
    #make density equal to N/L^2
    density = N/(L**2)


    v = []
    for l in tqdm(L):
        average = []
        for i in range(avg):
            average.append(retrieve_param_eta(eta,it,N,V, l, interaction_radius_1, interaction_radius_2, interaction_radius_3))
        v.append(np.mean(average))

    #save the values for Class1
    np.save(f'v_it_{it}_N_{N}_V_{V}_eta_{eta}_R1_{interaction_radius_1}_R2_{interaction_radius_2}_R3_{interaction_radius_3}_avg_{avg}_nb_pts_{nb_pts}_Class2.npy', v)
    np.save(f'density_it_{it}_N_{N}_V_{V}_eta_{eta}_R1_{interaction_radius_1}_R2_{interaction_radius_2}_R3_{interaction_radius_3}_avg_{avg}_nb_pts_{nb_pts}_Class2.npy', density)

    #plot it all

    plt.scatter(density,v)
    plt.xlabel("density")
    plt.ylabel("v")
    #display the parameter values
    plt.title("N = " + str(N) + ", V = " + str(V) + ", eta = " + str(eta) + ", R1 = " + str(interaction_radius_1)
               + ", R2 = " + str(interaction_radius_2) + ", R3 = " + str(interaction_radius_3))
    plt.show()

plot_2_a()