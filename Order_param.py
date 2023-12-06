import numpy as np
import matplotlib.pyplot as plt
from Predator import *

import numpy as np
import numpy.random as rnd
import itertools
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.animation import FuncAnimation

from Predator import *
        
            
#create swarm of birds and do an animation of the evolution
L = 15
N = 50
V = 0.03
eta = 0.1
interaction_radius_1 = 1
interaction_radius_2 = 2
interaction_radius_3 = 3
birds_awarness = 0.8   #close to 1 = very aware, close to 0 = not aware
birds_acceleration = 1.5  #how many times the base speed

# Create a predator
predator = Predator(X=7, Y=7, velocity=0.07, detection_radius=3)

swarm = Swarm(L, N, V, eta, interaction_radius_1, interaction_radius_2, interaction_radius_3, birds_awarness, birds_acceleration)
swarm.add_predator(predator)
#swarm = Swarm(L, N, V, eta, interaction_radius_1)
swarm.initialize()
# Number of iterations
num_iterations = 1000

# Lists to store mean velocities
mean_velocities = []

# Run iterations
for _ in range(num_iterations):
    swarm.evolve()
    mean_velocities.append(swarm.get_mean_norm_velocity())

# Plot the mean velocity over iterations
plt.plot(range(num_iterations), mean_velocities, label='Mean Velocity')
plt.xlabel('Iteration')
plt.ylabel('Mean Velocity')
plt.legend()
plt.show()
