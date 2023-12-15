import numpy as np
import sys
sys.path.append('.')
from our_library import Complex_model as cm
        
            
#create swarm of birds and do an animation of the evolution
avg = 20
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

def plot_2_a(it=200, N=40, V=0.03, L=6.2, interaction_radius_1=1, interaction_radius_2=2, interaction_radius_3=3):
    eta = np.linspace(0.0, 5, nb_pts)
    v_avg = []
    v_uncertainty = []

    for e in eta:
        velocities = []
        for i in range(avg):
            velocities.append(retrieve_param_eta(e, it, N, V, L, interaction_radius_1, interaction_radius_2, interaction_radius_3))

        average, uncertainty = calculate_avg_and_uncertainty(velocities)
        v_avg.append(average)
        v_uncertainty.append(uncertainty)

    # Save the values
    np.save(f'v_avg_it_{it}_N_{N}_V_{V}_L_{L}_R1_{interaction_radius_1}_R2_{interaction_radius_2}_R3_{interaction_radius_3}_avg_{avg}_nb_pts_{nb_pts}_Class2.npy', v_avg)
    np.save(f'v_uncertainty_it_{it}_N_{N}_V_{V}_L_{L}_R1_{interaction_radius_1}_R2_{interaction_radius_2}_R3_{interaction_radius_3}_avg_{avg}_nb_pts_{nb_pts}_Class2.npy', v_uncertainty)
    np.save(f'eta_it_{it}_N_{N}_V_{V}_L_{L}_R1_{interaction_radius_1}_R2_{interaction_radius_2}_R3_{interaction_radius_3}_avg_{avg}_nb_pts_{nb_pts}_Class2.npy', eta)

plot_2_a()