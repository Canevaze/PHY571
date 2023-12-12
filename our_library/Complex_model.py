import numpy as np
import numpy.random as rnd
import itertools
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.animation import FuncAnimation
import scipy.stats as sc
class Bird:
    """
    
    A configuration of bird parameters, within the class Complex_model.
    The bird only interact with its neighbors, and the neighbors are defined by a radius of interaction.
    Its neighbors are divided into three groups.
    The ones within R1 are considered short range neighbors, and the bird tries to avoid them in order to avoid collisions.
    The ones within R2 are considered medium range neighbors, and the bird tries to align its velocity with the mean velocity of its neighbors.
    The ones within R3 are considered long range neighbors, and the bird tries to go straight towards them.

    """
    def __init__(self, X, Y, theta, V):
        self.X = X
        self.Y = Y
        self.theta = theta
        self.velocity = V
        self.all_thetas = [theta]
    
    def get_neighbors(self, swarm, R1, R2, R3, L):
        """
        
        Get the neighbors of a bird. This includes birds on the other side of the grid because there are periodic boundary conditions.
        
        :param self: The bird itself.
        :type self: Class

        :param swarm: The swarm the bird belongs to.
        :type swarm: Class
    
        :param R1: The short range radius of interaction, i.e. the distance in which birds are considered as short range neighbors. The bird tries to avoid them.
        :type R1: int

        :param R2: The mdeium range radius of interaction, i.e. the distance in which birds are considered as mdeium range neighbors. The bird tries to align its velocity with their mean velocity.
        :type R2: int

        :param R3: The long range radius of interaction, i.e. the distance in which birds are considered as long range neighbors. The bird tries to go straight towards them.
        :type R3: int

        :param L: The size of the box, with periodic conditions.
        :type L: int
        
        :return: The three groups of neighbors of the bird.
        :rtype: tuple of lists
        
        """

        neighbors_short = []
        neighbors_long = []
        neighbors_medium = []
        for bird in swarm:
            if bird != self:
                # Calculate the distance between birds with periodic boundary conditions
                dx = abs(self.X - bird.X)
                dy = abs(self.Y - bird.Y)
                # Apply periodic boundary conditions
                dx = min(dx, L - dx)
                dy = min(dy, L - dy)
                # #other way of calculating the distance
                # dx = (bird.X - self.X + L / 2) % L - L / 2
                # dy = (bird.Y - self.Y + L / 2) % L - L / 2
                distance = np.sqrt(dx**2 + dy**2)
                if distance <= R1:
                    neighbors_short.append(bird)
                elif distance <= R2:
                    neighbors_medium.append(bird)
                elif distance <= R3:
                    neighbors_long.append(bird)
        return neighbors_short, neighbors_medium, neighbors_long, distance
    
    def get_theta_medium(self, neighbors_medium):
        """"
        
        Get the average velocity direction of the medium range neighbors, according to the formula in the paper.
        
        :param self: The bird itself.
        :type self: Class

        :param neighbors: The medium range neighbors of the chosen bird.
        :type neighbors: list

        :return: The average velocity direction of the medium range neighbors.
        :rtype: float
        
        """
        if not neighbors_medium:   # so it is verified when calculating the new theta
            return np.nan
        
        thetas = [bird.theta for bird in neighbors_medium]
        mean = np.mean(np.sin(thetas))/np.mean(np.cos(thetas))
        if np.mean(np.cos(thetas))<0 :
            return np.arctan(mean) + np.pi
        else :
            return np.arctan(mean)
        
        
    def get_theta_long(self,neighbors_long,L):       
        """"
        
        Get the mean angle of all the angles that point to long range neighbors, so that the bird goes straight towards them.
        
        :param self: The bird itself.
        :type self: Class

        :param neighbors: The long range neighbors of the chosen bird.
        :type neighbors: list

        :param L: The size of the box, with peroidic conditions.
        :type L: int

        :return: The mean angle of all the angles that point to long range neighbors.
        :rtype: float
        
        """
        if not neighbors_long:   # so it is verified when calculating the new theta
            return np.nan
        # Calculate the differences in coordinates for each neighbor
        delta_x = np.array([(neighbor.X - self.X + L/2) % L - L/2 for neighbor in neighbors_long])
        delta_y = np.array([(neighbor.Y - self.Y + L/2) % L - L/2 for neighbor in neighbors_long])
        #delta_x = np.array([neighbor.X - self.X for neighbor in neighbors])
        #delta_y = np.array([neighbor.Y - self.Y for neighbor in neighbors])
        # Calculate the angles using arctan2
        angles = np.arctan2(delta_y, delta_x)

        return sc.circmean(angles)

    def get_theta_short(self,neighbors_short, L):
        """"
        
        Get the mean angle of all the angles that point to the opposite of the short range neighbors, so that the bird avoids them.
        
        :param neighbors: The short range neighbors of the chosen bird.
        :type neighbors: list

        :param L: The size of the box, with periodic conditions.
        :type L: int

        :return: The mean angle of all the angles that point to the opposite of the short range neighbors.
        :rtype: float
        
        """
        if not neighbors_short:   # so it is verified when calculating the new theta
            return np.nan
        
        # Calculate the differences in coordinates for each neighbor
        delta_x = np.array([(neighbor.X - self.X + L/2) % L - L/2 for neighbor in neighbors_short])
        delta_y = np.array([(neighbor.Y - self.Y + L/2) % L - L/2 for neighbor in neighbors_short])
        #delta_x = np.array([neighbor.X - self.X for neighbor in neighbors])
        #delta_y = np.array([neighbor.Y - self.Y for neighbor in neighbors])

        # Calculate the angles using arctan2
        angles = np.arctan2(delta_y, delta_x)

        return sc.circmean(angles) + np.pi

        # Convert angles to radians
        angles = np.radians(angles)

        weights = 1 / np.sqrt(delta_x**2 + delta_y**2) ## we tried to ponderate the angles by the distance but it made the behavior less realistic

        # Calculate weighted sums
        sum_sin = np.sum(weights * np.sin(angles))
        sum_cos = np.sum(weights * np.cos(angles))

        # Calculate the circular mean
        mean_angle = np.arctan2(sum_sin, sum_cos)

        # Convert mean angle to degrees
        mean_angle_deg = np.degrees(mean_angle)

        # Ensure the result is between 0 and 360 degrees
        mean_angle_deg = (mean_angle_deg + 360) % 360

        #return sc.circmean(angles) + np.pi


        #do the circmean but by weighting by the distance
        return mean_angle_deg + np.pi




    # def evolve(self, dt):
    #     "get the new position of the bird"
    #     self.X += self.velocity * np.cos(self.theta) * dt
    #     self.Y += self.velocity * np.sin(self.theta) * dt
    #     self.theta += self.get_mean_theta(self.get_neighbors(swarm, R)) * dt
    
class Swarm :
    """"
    
    Creates the swarm state with many birds.

    """
    def __init__(self, L, N, V, eta, radius1, radius2, radius3):
        self.length = L
        self.number = N
        self.velocity_norm = V
        self.eta = eta
        self.interaction_radius_1 = radius1
        self.interaction_radius_2 = radius2
        self.interaction_radius_3 = radius3
        self.dt = 1
        self.rho = N/(L**2)
        self.birds = []
    
    def initialize(self):
        """
        
        Initialize the swarm with random positions and velocities.

        :param self: The swarm itself.
        :type self: Class
        
        """
        for i in range(self.number):
            X = rnd.uniform(0, self.length)
            Y = rnd.uniform(0, self.length)
            theta = rnd.uniform(0, 2*np.pi)
            
            self.birds.append(Bird(X, Y, theta, self.velocity_norm))
    
    def get_swarm_mean_velocity(self):
        """
        
        Get the mean vectorial velocity of the swarm, which is the order parameter.

        :param self: The swarm itself.
        :type self: Class
        
        """
        mean_vx = np.mean([bird.velocity * np.cos(sc.circmean(bird.all_thetas[-1])) for bird in self.birds])
        mean_vy = np.mean([bird.velocity * np.sin(sc.circmean(bird.all_thetas[-1])) for bird in self.birds])
        return np.sqrt(mean_vx**2 + mean_vy**2)/self.velocity_norm
    
    def evolve(self):
        """

        Evolve the swarm to the next time step.

        :param self: The swarm itself.
        :type self: Class

        """
        
        updated_birds = []
        for bird in self.birds:
            #assure the angle is between 0 and 2pi
            bird.theta = (bird.theta + 2*np.pi) % (2*np.pi)
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
            neigh = bird.get_neighbors(self.birds, self.interaction_radius_1,
                                                                        self.interaction_radius_2, self.interaction_radius_3
                                                                        ,self.length)
            new_theta_long = bird.get_theta_long(neigh[2],self.length) if neigh[2] else bird.theta
            new_theta_short = bird.get_theta_short(neigh[0],self.length) if neigh[0] else bird.theta
            new_theta_medium = bird.get_theta_medium(neigh[1]) if neigh[1] else bird.theta
            
            new_theta = sc.circmean([new_theta_medium, new_theta_long, new_theta_short]) + random_theta
            new_theta = (new_theta + 2*np.pi) % (2*np.pi)
            updated_birds.append([new_X, new_Y, new_theta])
        for i, bird in enumerate(self.birds):
            bird.X = updated_birds[i][0]
            bird.Y = updated_birds[i][1]
            bird.theta = updated_birds[i][2]
            bird.all_thetas.append(bird.theta)