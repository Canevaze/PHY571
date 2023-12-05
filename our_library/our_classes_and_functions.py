import numpy as np
import numpy.random as rnd
import itertools
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.animation import FuncAnimation

def haversine(lon1: float, lat1: float, lon2: float, lat2: float) -> float:
    """
    Calculate the great circle distance between two points on the 
    earth (specified in decimal degrees), returns the distance in
    kilometers.
    All arguments must be of equal length.
    :param lon1: longitude of first place
    :param lat1: latitude of first place
    :param lon2: longitude of second place
    :param lat2: latitude of second place
    :return: distance in kilometers between the two sets of coordinates
    """
    # Convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])

    # Haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a))
    r = 6371 # Radius of earth in kilometers
    return c * r
class Bird:
    """
    
    A configuration of bird parameters

    
    """

    def __init__(self, X, Y, theta, V):
        self.X = X
        self.Y = Y
        self.theta = theta
        self.velocity = V
        self.all_thetas = [theta]
    




    def get_neighbors(self, swarm, R, L) -> list:
        """
        
        Get the neighbors of a bird. This includes birds on the other side of the grid because there are periodic boundary conditions.
        
        :param self: The bird itself.
        :type self: Class

        :param swarm: The swarm the bird belongs to.
        :type swarm: Class
    
        :param R: The radius of interaction, i.e. the distance in which birds are considered as neighbors.
        :type R: int

        :param L: The size of the box, with peridic conditions.
        :type L: int
        
        :return: The neighbors of the chosen bird.
        :rtype: list
        
        """

        neighbors = [self]
        for bird in swarm:
            if bird != self:
                # Calculate the distance between birds with periodic boundary conditions
                dx = abs(self.X - bird.X)
                dy = abs(self.Y - bird.Y)

                # Apply periodic boundary conditions
                dx = min(dx, L - dx)
                dy = min(dy, L - dy)

                distance = np.sqrt(dx**2 + dy**2)

                if distance <= R:
                    neighbors.append(bird)

        return neighbors

    
    def get_mean_theta(self, neighbors):
        "get the average theta of the neighbors, according to the formula in the paper"
        thetas = [bird.theta for bird in neighbors]

        mean = np.mean(np.sin(thetas))/np.mean(np.cos(thetas))
        if np.mean(np.cos(thetas))<0 :
            return np.arctan(mean) + np.pi
        else :
            return np.arctan(mean)

    # def evolve(self, dt):
    #     "get the new position of the bird"

    #     self.X += self.velocity * np.cos(self.theta) * dt
    #     self.Y += self.velocity * np.sin(self.theta) * dt
    #     self.theta += self.get_mean_theta(self.get_neighbors(swarm, R)) * dt
    
class Swarm :
    "Creates the swarm state with many birds"

    def __init__(self, L, N, V, eta, radius):
        self.length = L
        self.number = N
        self.velocity_norm = V
        self.eta = eta
        self.interaction_radius = radius

        self.dt = 1
        self.rho = N/(L**2)
        self.birds = []

    def initialize(self):
        "initialize the swarm with random positions and velocities"

        for i in range(self.number):
            X = rnd.uniform(0, self.length)
            Y = rnd.uniform(0, self.length)
            theta = rnd.uniform(0, 2*np.pi)

            self.birds.append(Bird(X, Y, theta, self.velocity_norm))

    def get_swarm_mean_velocity(self):
        "get the mean vectorial velocity of the swarm, which is the order parameter"

        mean_vx = np.mean([bird.velocity * np.cos(np.mean(bird.all_thetas[-1])) for bird in self.birds])
        mean_vy = np.mean([bird.velocity * np.sin(np.mean(bird.all_thetas[-1])) for bird in self.birds])

        return np.sqrt(mean_vx**2 + mean_vy**2)/self.velocity_norm

    def evolve(self):
        "evolve the swarm"
        
        updated_birds = []

        for bird in self.birds:
            new_X = bird.X + bird.velocity * np.cos(bird.theta) * self.dt    # birds are indistinguishable, so it doesn't matter which way the list goes
            new_Y = bird.Y + bird.velocity * np.sin(bird.theta) * self.dt

            #if the bird is out of the box, it comes back from the other side
            if new_X > self.length:
                new_X -= self.length
            elif new_X < 0:
                new_X += self.length

            if new_Y > self.length:
                new_Y -= self.length
            elif new_Y < 0:
                new_Y += self.length

            random_theta = np.random.uniform(-self.eta/2, self.eta/2)
            new_theta = random_theta + bird.get_mean_theta(bird.get_neighbors(self.birds, self.interaction_radius,self.length))
            updated_birds.append([new_X, new_Y, new_theta])

        for i, bird in enumerate(self.birds):
            bird.X = updated_birds[i][0]
            bird.Y = updated_birds[i][1]
            bird.theta = updated_birds[i][2]
            bird.all_thetas.append(bird.theta)
        