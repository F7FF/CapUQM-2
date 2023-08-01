#File for all of the simulation tools used. Analysis tools are stored seperately.

import networkx as nx 
import numpy as np
import matplotlib.pyplot as plt
import random

class QMgraph: #a little bit like Qmatrix, but a couple more tricks. mutable!
    def __init__(self, dim):
        self.graph = nx.Graph()
        self.dim = dim
        #create nodes in a grid of size dim, and assign each one a position ("pos") using **attr. of course, pos isn't actually passed to the QPU, but it's tidy to keep it all together.
        for x in range(dim):
            for y in range(dim):
                self.graph.add_node(dim*y + x, pos=(x,y))

    def distance_between(self, a, b): #returns the distance between node number A and node number B
        positions = nx.get_node_attributes(self.graph,'pos')
        ax, ay = positions[a]
        bx, by = positions[b]
        return ((ax-bx)**2 + (ay-by)**2)**0.5
    
    def visualize(self, filename = "QMgraph"): #saves a diagram of the graph
        positions = nx.get_node_attributes(self.graph,'pos')
        nx.draw(self.graph, 
                positions,
                with_labels=True, #label each node
                )
        plt.savefig(filename)
    
    def assign_chemical(self, chemical_potential, temperature_beta): #mutates self.graph, adding the chemical potential nodes.
        for x in self.graph.nodes:
            self.graph.add_edge(x,x, potential = temperature_beta * chemical_potential)
    
    def assign_potentials(self, potential_function, cutoff=999.0): #mutates self.graph, adding in edges with potentials from potential_function. Any edges longer than cutoff will be ignored.
        for a in range(self.dim**2): 
            for b in range(a+1, self.dim**2):
                distance = self.distance_between(a,b)
                if(distance <= cutoff):
                    self.graph.add_edge(a, b, potential=potential_function.potential(distance))

    def stagger(self, amount): #moves each point a little bit, randomly, constrained by "amount". for example, an amount of 0.1 will take each coordinate of each node and will add some value between 0.1 and -0.1.
        positions = nx.get_node_attributes(self.graph,'pos')
        for a in range(self.dim**2):
            newposition = (positions[a][0] + random.uniform(amount, -amount), positions[a][1] + random.uniform(amount, -amount))
            nx.set_node_attributes(self.graph, {a:newposition}, name="pos")


class LennardJones():
    def __init__(self, sigma, epsilon):
        self.sigma = sigma
        self.epsilon = epsilon
    def potential(self, distance):
        dist_term = (distance / self.sigma)**-6
        return 4 * self.epsilon * (dist_term * (dist_term - 1))