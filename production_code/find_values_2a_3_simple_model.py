import numpy as np
import sys
sys.path.append('.')
from our_library import Simple_model as sm
        
            
#create swarm of birds and do an animation of the evolution
avg = 20
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
    v_avg = []
    v_uncertainty = []

    for e in eta:
        velocities = []
        for i in range(avg):
            velocities.append(retrieve_param_eta(e, it, N, V, L, R))

        average, uncertainty = calculate_avg_and_uncertainty(velocities)
        v_avg.append(average)
        v_uncertainty.append(uncertainty)

    #save the values
    np.save(f'v_avg_it_{it}_N_{N}_V_{V}_L_{L}_R_{R}_avg_{avg}_nb_pts_{nb_pts}_Class1.npy', v_avg)
    np.save(f'eta_it_{it}_N_{N}_V_{V}_L_{L}_R_{R}_avg_{avg}_nb_pts_{nb_pts}_Class1.npy', eta)
    np.save(f'v_uncertainty_it_{it}_N_{N}_V_{V}_L_{L}_R_{R}_avg_{avg}_nb_pts_{nb_pts}_Class1.npy', v_uncertainty)

plot_2_a()