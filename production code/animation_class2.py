import numpy as np
import numpy.random as rnd
import itertools
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.animation import FuncAnimation

from our_library import Complex_model as cm
        
            
#create swarm of birds and do an animation of the evolution
L = 15
N = 50
V = 0.03
eta = 0.1
interaction_radius_1 = 1
interaction_radius_2 = 2
interaction_radius_3 = 3


swarm = cm.Swarm(L, N, V, eta, interaction_radius_1, interaction_radius_2, interaction_radius_3)
swarm.initialize()


# Create the figure and axis for the plot
fig, ax = plt.subplots()




# Set the axis limits
ax.set_xlim(0, L)
ax.set_ylim(0, L)



# Initialize the plot objects
points, = ax.plot([bird.X for bird in swarm.birds], [bird.Y for bird in swarm.birds], 'bo')

# Update function
def update(q):
    swarm.evolve()
    points.set_data([bird.X for bird in swarm.birds], [bird.Y for bird in swarm.birds])

    return points


# Create the animation
animation = FuncAnimation(fig, update, frames=range(100), interval=1)

# Show the animation
plt.show()





