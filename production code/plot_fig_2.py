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

def retrieve_param_eta(eta,it,N,V, L, R, interaction_radius_1, interaction_radius_2, interaction_radius_3):
    #swarm = Swarm(L, N, V, eta, R)
    swarm = cm.Swarm(L, N, V, eta, interaction_radius_1, interaction_radius_2, interaction_radius_3)
    swarm.initialize()

    for i in range(it):
        swarm.evolve()
    
    return swarm.get_swarm_mean_velocity()

def plot_2_a(it=200,N=40,V=0.03, L=5, R=1, interaction_radius_1=1, interaction_radius_2=2, interaction_radius_3=3):
    eta = np.linspace(0.0,5,nb_pts)
    v = []
    for e in tqdm(eta):
        average = []
        for i in range(avg):
            average.append(retrieve_param_eta(e,it,N,V, L, R, interaction_radius_1, interaction_radius_2, interaction_radius_3))
        v.append(np.mean(average))

    #save the values
    np.save(f'v_it_{it}_N_{N}_V_{V}_L_{L}_R_{R}_avg_{avg}_nb_pts_{nb_pts}_Class2.npy', v)
    np.save(f'eta_it_{it}_N_{N}_V_{V}_L_{L}_R_{R}_avg_{avg}_nb_pts_{nb_pts}_Class2.npy', eta)


    #plot it all

    plt.scatter(eta,v)
    plt.xlabel("eta")
    plt.ylabel("v")
    #display the parameter values
    plt.title("N = " + str(N) + ", V = " + str(V) + ", L = " + str(L) + ", R1 = " + str(interaction_radius_1)
               + ", R2 = " + str(interaction_radius_2) + ", R3 = " + str(interaction_radius_3))
    plt.show()

def plot_2_b(it=200,N=40,V=0.03, eta=2 , R=1, interaction_radius_1=0.5, interaction_radius_2=1, interaction_radius_3=2):
    L = np.linspace(2,25,nb_pts)
    #make density equal to N/L^2
    density = N/(L**2)


    v = []
    for l in tqdm(L):
        average = []
        for i in range(avg):
            average.append(retrieve_param_eta(eta,it,N,V, l, R, interaction_radius_1, interaction_radius_2, interaction_radius_3))
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