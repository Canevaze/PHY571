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

def plot(it=200,N=100,V=0.03, L=5, R=1, avg=20, nb_pts=50):
    eta = np.linspace(0,5,nb_pts)
    v = []
    for e in tqdm(eta):
        average = []
        for i in range(avg):
            average.append(retrieve_param_eta(e,it,N,V, L, R))
        v.append(np.mean(average))

    #save the values
    np.save(f'v_it_{it}_N_{N}_V_{V}_L_{L}_R_{R}_avg_{avg}_nb_pts_{nb_pts}.npy', v)
    np.save(f'eta_it_{it}_N_{N}_V_{V}_L_{L}_R_{R}_avg_{avg}_nb_pts_{nb_pts}.npy', eta)
    
    
    plt.ylim(0, 1)
    plt.xlim(0, 5)
    
    #plot it all

    plt.scatter(eta,v)
    plt.xlabel("eta")
    plt.ylabel("v")
    #display the parameter values
    plt.title("N = " + str(N) + ", V = " + str(V) + ", L = " + str(L) + ", R = " + str(R))
    plt.show()

    
plot()
