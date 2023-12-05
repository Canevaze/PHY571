import numpy as np
import numpy.random as rnd
import itertools
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.animation import FuncAnimation
from tqdm import tqdm

from Class2 import *
        
            
#create swarm of birds and do an animation of the evolution

interaction_radius_1 = 1
interaction_radius_2 = 2
interaction_radius_3 = 3

def retrieve_param_eta(eta,it,N,V, L, R, interaction_radius_1, interaction_radius_2, interaction_radius_3):
    #swarm = Swarm(L, N, V, eta, R)
    swarm = Swarm(L, N, V, eta, interaction_radius_1, interaction_radius_2, interaction_radius_3)
    swarm.initialize()

    for i in range(it):
        swarm.evolve()
    
    return swarm.get_swarm_mean_velocity()

def plot(it=200,N=40,V=0.03, L=5, R=1, interaction_radius_1=1, interaction_radius_2=2, interaction_radius_3=3):
    eta = np.linspace(0.0,5,50)
    v = []
    for e in tqdm(eta):
        average = []
        for i in range(5):
            average.append(retrieve_param_eta(e,it,N,V, L, R, interaction_radius_1, interaction_radius_2, interaction_radius_3))
        v.append(np.mean(average))



    #plot it all

    plt.scatter(eta,v)
    plt.xlabel("eta")
    plt.ylabel("v")
    #display the parameter values
    plt.title("N = " + str(N) + ", V = " + str(V) + ", L = " + str(L) + ", R = " + str(R))
    plt.show()


plot()