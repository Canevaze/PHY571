import numpy as np
import numpy.random as rnd
import itertools
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.animation import FuncAnimation
from tqdm import tqdm

from Class2 import *
        
            
#create swarm of birds and do an animation of the evolution


def retrieve_param_eta(eta,it,N,V, L, R, interaction_radius_1, interaction_radius_2, interaction_radius_3):
    #swarm = Swarm(L, N, V, eta, R)
    swarm = Swarm(L, N, V, eta, interaction_radius_1, interaction_radius_2, interaction_radius_3)
    swarm.initialize()

    for i in range(it):
        swarm.evolve()
    
    return swarm.get_swarm_mean_velocity()

def plot_2_a(it=200,N=40,V=0.03, L=5, R=1, interaction_radius_1=0.5, interaction_radius_2=1, interaction_radius_3=2):
    eta = np.linspace(0.0,5,50)
    v = []
    for e in tqdm(eta):
        average = []
        for i in range(20):
            average.append(retrieve_param_eta(e,it,N,V, L, R, interaction_radius_1, interaction_radius_2, interaction_radius_3))
        v.append(np.mean(average))



    #plot it all

    plt.scatter(eta,v)
    plt.xlabel("eta")
    plt.ylabel("v")
    #display the parameter values
    plt.title("N = " + str(N) + ", V = " + str(V) + ", L = " + str(L) + ", R1 = " + str(interaction_radius_1)
               + ", R2 = " + str(interaction_radius_2) + ", R3 = " + str(interaction_radius_3))
    plt.show()

def plot_2_b(it=200,N=40,V=0.03, eta=2 , R=1, interaction_radius_1=0.5, interaction_radius_2=1, interaction_radius_3=2):
    L = np.linspace(2,25,50)
    #make density equal to N/L^2
    density = N/(L**2)


    v = []
    for l in tqdm(L):
        average = []
        for i in range(5):
            average.append(retrieve_param_eta(eta,it,N,V, l, R, interaction_radius_1, interaction_radius_2, interaction_radius_3))
        v.append(np.mean(average))



    #plot it all

    plt.scatter(density,v)
    plt.xlabel("density")
    plt.ylabel("v")
    #display the parameter values
    plt.title("N = " + str(N) + ", V = " + str(V) + ", eta = " + str(eta) + ", R1 = " + str(interaction_radius_1)
               + ", R2 = " + str(interaction_radius_2) + ", R3 = " + str(interaction_radius_3))
    plt.show()

plot_2_b()