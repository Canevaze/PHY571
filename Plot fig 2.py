import numpy as np
import numpy.random as rnd
import itertools
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.animation import FuncAnimation
from tqdm import tqdm

from Class import *
        
            
#create swarm of birds and do an animation of the evolution


def retrieve_param_eta(eta,it,N,V, L, R):
    swarm = Swarm(L, N, V, eta, R)
    swarm.initialize()

    for i in range(it):
        swarm.evolve()
    
    return swarm.get_swarm_mean_velocity()

def plot(it=200,N=40,V=0.03, L=3.1, R=1):
    eta = np.linspace(0.0,5,50)
    v = []
    for e in tqdm(eta):
        average = []
        for k in range(10):
            average.append(retrieve_param_eta(e,it,N,V, L, R))
        v.append(np.mean(average))
    plt.scatter(eta,v)
    plt.xlabel("eta")
    plt.ylabel("v")
    plt.show()

    
plot()