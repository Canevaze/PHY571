import numpy as np
import numpy.random as rnd
import itertools
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.animation import FuncAnimation
import scipy.stats as sc

class Bird:
    """
    
    A configuration of bird parameters, within the class Predator_model.
    The bird only interact with its neighbors, and the neighbors are defined by a radius of interaction.
    Its neighbors are divided into three groups.
    The ones within R1 are considered short range neighbors, and the bird tries to avoid them in order to avoid collisions.
    The ones within R2 are considered medium range neighbors, and the bird tries to align its velocity with the mean velocity of its neighbors.
    The ones within R3 are considered long range neighbors, and the bird tries to go straight towards them.

    In addition, a predator is added.
    The bird tries to avoid the predator, and the predator tries to catch the bird.
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

            number_of_endangered_birds = 0

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

                if bird.velocity > self.velocity:
                    number_of_endangered_birds += 1
                    

        return neighbors_short, neighbors_medium, neighbors_long, distance, number_of_endangered_birds

    
    def get_theta_medium(self, neighbors_medium):
        """"
        
        Get the average velocity direction of the mdeium range neighbors, according to the formula in the paper.
        
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
        
        :param self: The bird itself.
        :type self: Class
        
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

        

        #do the circmean but by weighting by the distance
        return mean_angle_deg + np.pi
        
    def theta_predator(self, predator, L):
        """"
        
        Get the angle that points to the predator
        
        :param self: The bird itself.
        :type self: Class
        
        :param predator: The predator.
        :type predator: Class

        :param L: The size of the box, with periodic conditions.
        :type L: int

        :return: The angle that points to the predator.
        :rtype: float

        """
        delta_x = (predator.X - self.X + L/2) % L - L/2
        delta_y = (predator.Y - self.Y + L/2) % L - L/2

        # Calculate the angles using arctan2
        angles = np.arctan2(delta_y, delta_x)

        return angles + np.pi


class Predator:
    """
    
    A class for the predator bird
    
    """

    def __init__(self, X, Y, velocity, detection_radius):
        self.X = X
        self.Y = Y
        self.velocity = velocity
        self.detection_radius = detection_radius
        self.all_positions = [(X, Y)]

    def update_position(self, prey_positions, dt, L):
        """

        Update the position of the predator, according to the position of the preys.
        
        :param self: The predator itself.
        :type self: Class

        :param prey_positions: The positions of the preys inside the detection radius. If no bird is inside the detection raidus, the predator aims for the mean position of the swarm.
        :type prey_positions: list

        :param dt: The time step.
        :type dt: float
        
        :param L: The size of the box, with periodic conditions.
        :type L: int

        :return: the list new_positions, updated with the new position of the predator at the end of the list. 
        :rtype: list

        """      
        # For example, you can make the predator move towards the average position of the prey
        if prey_positions:
            direction_x = np.array([(pos[0] - self.X + L/2) % L - L/2 for pos in prey_positions])
            direction_y = np.array([(pos[1] - self.Y + L/2) % L - L/2 for pos in prey_positions])

            #calculate average direction
            direction_x = np.mean(direction_x)
            direction_y = np.mean(direction_y)

            distance = np.sqrt(direction_x**2 + direction_y**2)


            # mean_prey_x = np.mean([pos[0] for pos in prey_positions])
            # mean_prey_y = np.mean([pos[1] for pos in prey_positions])
            # direction_x = (mean_prey_x - self.X + length/2) % length - length/2
            # direction_y = (mean_prey_y - self.Y + length/2) % length - length/2
            # distance = np.sqrt(direction_x**2 + direction_y**2)
            
            if distance > 0:
                direction_x /= distance
                direction_y /= distance
                
                self.X += direction_x * self.velocity * dt
                self.Y += direction_y * self.velocity * dt

                # Apply periodic boundary conditions  (other way of doing it, compared to what we did with the birds)
                self.X = (self.X + L) % L
                self.Y = (self.Y + L) % L

        # Save the new position
        self.all_positions.append((self.X, self.Y))
    
    def get_velocity_vector(self, prey_positions, L):
        """

        Get the components of the velocity vector of the predator, according to the position of the preys.
        
        :param self: The predator itself.
        :type self: Class
        
        :param prey_positions: The positions of the preys inside the detection radius. If no bird is inside the detection raidus, the predator aims for the mean position of the swarm.
        :type prey_positions: list

        :param L: The size of the box, with periodic conditions.
        :type L: int

        :return: The components of the velocity vector of the predator.
        :rtype: tuple of floats

        """ 
        # For example, you can make the predator move towards the average position of the prey
        if prey_positions:
            direction_x = np.array([(pos[0] - self.X + L/2) % L - L/2 for pos in prey_positions])
            direction_y = np.array([(pos[1] - self.Y + L/2) % L - L/2 for pos in prey_positions])

            #calculate average direction
            direction_x = np.mean(direction_x)
            direction_y = np.mean(direction_y)

        return direction_x, direction_y


    def get_predator_neighbors(self, swarm, R, L):
        """
        
        Get the neighbors of the predator. This includes birds on the other side of the grid because there are periodic boundary conditions.
        
        :param self: The predator itself.
        :type self: Class

        :param swarm: The swarm the predator belongs to.
        :type swarm: Class
    
        :param R: The radius of interaction, i.e. the distance in which birds are considered as neighbors.
        :type R: int

        :param L: The size of the box, with periodic conditions.
        :type L: int
        
        :return: The neighbors of the predator.
        :rtype: list
        
        """

        neighbors = [self]
        for bird in swarm:
            if bird != self:
                # Calculate the distance between birds with periodic boundary conditions
                dx = self.X - bird.X
                dy = self.Y - bird.Y

                # Apply periodic boundary conditions
                dx = (dx + L / 2) % L - L / 2
                dy = (dy + L / 2) % L - L / 2

                distance = np.sqrt(dx**2 + dy**2)

                if 0 < distance <= R:
                    neighbors.append(bird)

        return neighbors
    




    
class Swarm :
    """
    
    Creates the swarm state with many birds.
    
    """
    def __init__(self, L, N, V, eta, radius1, radius2, radius3, awarness, birds_acceleration):
        self.length = L
        self.number = N
        self.velocity_norm = V
        self.eta = eta
        self.interaction_radius_1 = radius1
        self.interaction_radius_2 = radius2
        self.interaction_radius_3 = radius3
        self.birds_awarness = awarness
        self.birds_acceleration = birds_acceleration

        self.dt = 1
        self.rho = N/(L**2)
        self.birds = []
        self.predator = None

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

        :return: The mean vectorial velocity of the swarm.  
        :rtype: float
        
        """
        mean_vx = np.mean([bird.velocity * np.cos(sc.circmean(bird.all_thetas[-1])) for bird in self.birds])
        mean_vy = np.mean([bird.velocity * np.sin(sc.circmean(bird.all_thetas[-1])) for bird in self.birds])

        return np.sqrt(mean_vx**2 + mean_vy**2)/self.velocity_norm

    def get_mean_norm_velocity(self):
        """
        
        Get the mean norm velocity of the swarm.

        :param self: The swarm itself.
        :type self: Class

        :return: The mean norm velocity of the swarm.
        :rtype: float
        
        """

        mean_v = np.mean([bird.velocity for bird in self.birds])

        return mean_v

    # ------------------ Predator methods ------------------ #

    def add_predator(self, predator):
        """

        Add the predator to the swarm.

        :param self: The swarm itself.
        :type self: Class

        :param predator: The predator.
        :type predator: Class

        """
        self.predator = predator


    def get_prey_positions(self):
        """
        
        Get the positions of the preys inside the detection radius of the predator. If no bird is inside the detection radius, gives the positions of all the birds.

        :param self: The swarm itself.
        :type self: Class

        :return: The positions of the preys inside the detection radius of the predator. If no bird is inside the detection radius, gives the positions of all the birds.
        :rtype: list
        
        """
        if self.predator:
            neighbors = self.predator.get_predator_neighbors(self.birds, self.predator.detection_radius, self.length)
            if len(neighbors) > 1:  # Check if there are more than just the predator in the list
                return [(bird.X, bird.Y) for bird in neighbors]
        # If no prey or neighbors, return positions of all birds
        return [(bird.X, bird.Y) for bird in self.birds]



    # ------------------ Evolution methods ------------------ #

    def evolve(self):
        """"
        
        Evolve the swarm to the next time step.

        :param self: The swarm itself.
        :type self: Class
        
        """


        #predator part to begin with
        prey_positions = self.get_prey_positions()

        # Update predator position
        self.predator.update_position(prey_positions, self.dt, self.length)
        
        updated_birds = []

        for bird in self.birds:

            #assure the angle is between 0 and 2pi
            bird.theta = (bird.theta + 2*np.pi) % (2*np.pi)

            random_theta = np.random.uniform(-self.eta/2, self.eta/2)

            disX = abs(bird.X - self.predator.X)
            disY = abs(bird.Y - self.predator.Y)

            disX = min(disX, self.length - disX)
            disY = min(disY, self.length - disY)

            dis = np.sqrt(disX**2 + disY**2)

            #if the bird is to close from the predator, it tries to escape
            if dis < self.interaction_radius_3:

                #get the velocity vector of the predator
                predator_vel = self.predator.get_velocity_vector(prey_positions, self.length)
                predator_velocity_x = predator_vel[0]
                predator_velocity_y = predator_vel[1]

                #get the vector from the predator to the bird
                vector_to_bird_x = bird.X - self.predator.X
                vector_to_bird_y = bird.Y - self.predator.Y

                #get the vetorial product of the two vectors
                vetorial_product = predator_velocity_x*vector_to_bird_y - predator_velocity_y*vector_to_bird_x

                #if the vetorial product is positive, the bird is on the right of the predator, so it goes even more to the right
                if vetorial_product > 0:

                    new_theta = bird.theta_predator(self.predator, self.length) + random_theta + np.pi/4
                    new_theta = (new_theta + 2*np.pi) % (2*np.pi)
                
                #if the vetorial product is negative, the bird is on the left of the predator, so it goes even more to the left
                else:
                    new_theta = bird.theta_predator(self.predator, self.length) + random_theta - np.pi/4
                    new_theta = (new_theta + 2*np.pi) % (2*np.pi)

                new_X = bird.X + (bird.velocity*self.birds_acceleration) * np.cos(bird.theta) * self.dt    # birds are indistinguishable, so it doesn't matter which way the list goes
                new_Y = bird.Y + (bird.velocity*self.birds_acceleration) * np.sin(bird.theta) * self.dt

            else:

                neigh = bird.get_neighbors(self.birds, self.interaction_radius_1,
                                                                            self.interaction_radius_2, self.interaction_radius_3
                                                                            ,self.length)



                new_theta_long = bird.get_theta_long(neigh[2],self.length) if neigh[2] else bird.theta
                new_theta_short = bird.get_theta_short(neigh[0],self.length) if neigh[0] else bird.theta
                new_theta_medium = bird.get_theta_medium(neigh[1]) if neigh[1] else bird.theta




                
                new_theta = sc.circmean([new_theta_medium, new_theta_long, new_theta_short]) + random_theta


                new_theta = (new_theta + 2*np.pi) % (2*np.pi)
                
                #if to many endangered birds, bird velocity increases
                if neigh[4] > (1 - self.birds_awarness)*(len(neigh[2])+len(neigh[1])+len(neigh[0])) + 1:
                    new_X = bird.X + (bird.velocity*self.birds_acceleration) * np.cos(bird.theta) * self.dt
                    new_Y = bird.Y + (bird.velocity*self.birds_acceleration) * np.sin(bird.theta) * self.dt

                else:

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


            updated_birds.append([new_X, new_Y, new_theta])

        for i, bird in enumerate(self.birds):
            bird.X = updated_birds[i][0]
            bird.Y = updated_birds[i][1]
            bird.theta = updated_birds[i][2]
            bird.all_thetas.append(bird.theta)
        