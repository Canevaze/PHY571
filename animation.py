import numpy as np
import numpy.random as rnd
import itertools
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.animation import FuncAnimation

from Predator import *
        
            
#create swarm of birds and do an animation of the evolution
L = 15
N = 80
V = 0.03
eta = 0.1
interaction_radius_1 = 1
interaction_radius_2 = 2
interaction_radius_3 = 3
birds_awarness = 0.8   #close to 1 = very aware, close to 0 = not aware
birds_acceleration = 2  #how many times the base speed

# Create a predator
predator = Predator(X=7, Y=7, velocity=0.07, detection_radius=3)

swarm = Swarm(L, N, V, eta, interaction_radius_1, interaction_radius_2, interaction_radius_3, birds_awarness, birds_acceleration)
swarm.add_predator(predator)
#swarm = Swarm(L, N, V, eta, interaction_radius_1)
swarm.initialize()


# Create the figure and axis for the plot
fig, ax = plt.subplots()




# Set the axis limits
ax.set_xlim(0, L)
ax.set_ylim(0, L)



# Initialize the plot objects
points, = ax.plot([bird.X for bird in swarm.birds], [bird.Y for bird in swarm.birds], 'bo')
predator_point, = ax.plot([], [], 'ro')

# Update function
def update(q):
    swarm.evolve()
    points.set_data([bird.X for bird in swarm.birds], [bird.Y for bird in swarm.birds])
    predator_point.set_data(swarm.predator.X, swarm.predator.Y)

    return points, predator_point


# Create the animation
animation = FuncAnimation(fig, update, frames=range(100), interval=1)

# Show the animation
plt.show()





