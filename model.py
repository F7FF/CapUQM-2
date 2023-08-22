"""
>           PROGRAMMING MODEL:

1: Make a fluid_problem, either by using one of the helper methods, or by manually assigning potentials using the class functions.

(optional step: save/load your fluid_problem)

2: Use a sampler to sample fluid_problem and generate a fluid_results object.

(again, fluid_results can be saved/loaded freely. the sampler should automatically save the results, too)

3: Use whatever analytical tool you please to analyze fluid_results.

>           WHY DID YOU DO IT LIKE THIS:

-The current version is heavily limited by its qmatrix creating system (networkx is not needed, and we have a whole load of classes that are incredibly hard to read.
-This program is much less flexible but significantly shorter (200 lines vs 1500 lines). Although this will need modification to run 3D simulations, this can be much more easily peer-reviewed (I still don't 100% understand how q_matrix works after three months).
-I've left room in the code to potentially stagger the cells in the future - Always use the positions dict to get the position of a cell!

g(r) works by:
1: pick a bin size - smaller bins will mean more data points but also more noise.
2: consider all possible "rings" of that width
3: number of cells inside that ring divided by the area is g(r)
"""

import json #used to save/load fluid_problems
import datetime #used to timestamp samples/problems
import random


#do NOT use "from x import *" - importing normally will make module management so much easier!
import dimod, dwave.samplers, dwave.system

pi = 3.1415926 #let's just define this right away

#import dwave.system
#import dwave.samplers

class fluid_problem: #EXCLUSIVELY FOR 2D PROBLEMS - 3D PROBLEMS WILL REQUIRE A REWRITE
    def __init__(self, width, name=""):
        self.linearconnections = {} #Contains all self-loops (think: chemical potential) in form cellID:strength. USED TO CREATE BQMS!
        self.quadraticconnections = {} #Contains all cross-cell connections in form (cellID_1, cellID_2):strength. USED TO CREATE BQMS!
        self.positions = {} #maps coordinates to cell IDs, using form cellID:(x,y). Usually they will be arranged in a simple grid, but they might be staggered in future.
        self.width = width #used for iterating later - PLEASE DON'T CHANGE ME!
        self.notes = [name + " - created with width=" + str(width)] #Whenever mutating this class, leave a little note explaining what you did! This can be read after from the file. Handy for keeping track.
        for y in range(width):
            for x in range(width):
                self.positions[y*width+x] = (x,y)

    def assign_linear_potentials(self, strength): #Assigns a linear connection of a certain strength to all nodes. I'm making this a class function because it's pretty obvious and we don't need to worry about varying between versions
        self.notes.append("Added constant linear potential of strength "+str(strength))
        for cellID in self.positions:
            self.linearconnections[cellID] = float(strength)

    def assign_quadratic_potentials(self, function, note=None): #Assigns potentials to all pairs of nodes - includes PBC. Ask Max if you'd like an explanation!
        if note: self.notes.append(note)
        #notekeeping done - now let's iterate over pairs of cells
        for cell_A, cell_B, distance in all_cell_pairs_tiled(self.width, self.positions):
            #now we find the potential for that distance and add it on to the quadratic connection for these two cells
            self.quadraticconnections[(cell_A, cell_B)] = self.quadraticconnections.get((cell_A, cell_B), 0) + function(distance)
            #self.quadraticconnections.get is simply getting the current potential, or 0 if it doesn't exist (just like defaultdict)

def simple_fluidproblem(width, sigma, epsilon, chem_potential, name=""): #This constructs a fluid_problem from the simple parameters - this should be identical to the older version's model.
    def LJ(s, e, r): #sigma, epsilon, radius (distance)
        dist_term = (r/s)**-6
        return 4 * e * (dist_term * (dist_term - 1)) #unreduced!
    output = fluid_problem(width, name=name)
    output.assign_linear_potentials(chem_potential)
    output.assign_quadratic_potentials(lambda x:LJ(sigma, epsilon,x), note="Assigning quadratic potentials using LJ - sigma="+str(sigma)+", epsilon="+str(epsilon))
    return output

def all_cell_pairs_tiled(width, positions): #this yields (cell_A, cell_B, distance) for all pairs of cells INCLUDING TILING! ask Max for a full explanation
    for cell_A in positions:
        cell_A_position = positions[cell_A]
        for cell_B in positions:
            if cell_A >= cell_B:
                continue #do NOT self connect, and do not accidentally connect nodes twice! By only connecting if cell B's ID is higher than cell A's ID, nodes will never be connected twice.
            cell_B_true_position = positions[cell_B]
            #now: we iterate over the nine possible microstate tiles (eight neighbors plus the "centre" tile). We could increase this for more accuracy (ie. iterating over 25 tiles), but at the cost of speed, and LJ drops off fast anyway.
            for dy in range(-1, 2):
                for dx in range(-1, 2):
                    #we now calculate where cell B should be in this tile
                    cell_B_position = (cell_B_true_position[0] + width*dx, cell_B_true_position[1] + width*dy)
                    #now we calculate distance between cell_A_position and cell_B_position (the one that is in a different tile)
                    distance = ((cell_B_position[0]-cell_A_position[0])**2 + (cell_B_position[1]-cell_A_position[1])**2)**0.5
                    yield (cell_A, cell_B, distance)

"""
VALID SAMPLER NAMES:

"hybrid"                -hybrid solver (part QPU, part CPU)
"gradient_descent"      -CPU solver that simply performs a gradient descent algorithm. Very fast, gives reasonable results.
"random"                -CPU solver that simply returns a random microstate (50% chance any one cell is filled). Useless results.
"QPU"                   -the pure QPU solver from the older versions (only works up to 13x13).

POSSIBLE FUTURE IDEAS:
"tabu"                  -use the dwave-system tabu solver (seems to be a good quantum simulator)
"tree"                  -uses dwave-samplers tree decomposition solver... might be fast? maybe?

part of the reason this codebase exists is because different solvers give different datatypes and require different parameters - hence why I'm using custom classes for everything, so I can actually know what's going on
"""


def get_fluid_results(current_fluid_problem, number_of_samples, sampler_name): #Accepts a fluid_problem and returns a fluid_results. This is a black box in order to control the nightmare that is dimod/dwave-system/dwave-samplers.
    assert isinstance(current_fluid_problem, fluid_problem) #didya put a BQM in here?

    output = fluid_results() #prepare the object to return information to
    output.width = current_fluid_problem.width
    output.notes = current_fluid_problem.notes #pass the notes down the chain
    output.notes.append("Running " +str(number_of_samples)+ " samples using " + sampler_name)
    output.positions = current_fluid_problem.positions

    BQM = dimod.BinaryQuadraticModel(current_fluid_problem.linearconnections,
                                     current_fluid_problem.quadraticconnections,
                                     vartype='BINARY') #create a BQM that reflects the fluid_problem

    #now - check what sampler we want to use
    if(sampler_name == "gradient_descent"):
        sampler = dwave.samplers.SteepestDescentSolver()
        sampleresults = sampler.sample(BQM, num_reads=number_of_samples)
        output.actual_microstates = [tuple(cell for cell in result.values()) for result in sampleresults]
    elif(sampler_name == "random"):
        sampler = dwave.samplers.RandomSampler()
        sampleresults = sampler.sample(BQM, num_reads=number_of_samples)
        output.actual_microstates = [tuple(cell for cell in result.values()) for result in sampleresults]
    else: #clearly the sampler name was invalid
        raise Exception(sampler_name + " is not a valid sampler name!")

    #after all the sampling and preparing is done, return the final fluid_results object.
    return output

def print_microstate(microstate, dim): #Makes a pretty little thing to visualize a single 2d microstate. This ignores stagger.
    assert isinstance(microstate, (list, tuple)) #if this fails, you're passing print_microstate the wrong datatype!
    for y in range(dim):
        for x in range(dim):
            print(("░", "█")[microstate[y*dim+x]], end='')
        print("")


class fluid_results: #Note that this is immutable and should not be altered after being created by get_samples
    def __init__(self): #You shouldn't be initializing this on its own - this class should only be created by get_samples() or loading from file!
        self.actual_microstates = [] #THIS COULD BE A LIST OR A TUPLE! Contains tuples, each tuple is a microstate that was returned by the sampler
        self.notes = [] #same functionality as fluid_problem.notes. Just a list of strings. Whenever something is changed, leave a little note to be kind to others.
        self.width = 0 #width of the simulation
        self.positions = {} #maps coordinates to cell IDs, using form cellID:(x,y).

    def save_as(self, filename): #Saves to a filename so it can be loaded later. UNTESTED!
        output = ""
        for microstate in self.actual_microstates:
            for cell in microstate:
                output += str(cell)
            output += "\n"
        parameters = {"positions":self.positions, "notes":self.notes, "width":self.width}
        #now, add the parameters
        output += repr(parameters)
        with open(filename, "w") as savefile:
            savefile.write(output)
    def load_from(self, filename): #this should be called as fluid_results.load_from("thing.results")
        with open(filename, "r") as f:
            for line in f:
                if(line[0] == '{'): #if this is the parameter line
                    parameters = eval(line) #this could probably be done with JSON, but I'm a little lazy and this works great
                else: #if this is a microstate line
                    self.actual_microstates.append([int(letter) for letter in line[:-1]]) #trim off the last letter (the newline)
        self.positions, self.notes, self.width = parameters["positions"], parameters["notes"], parameters["width"]

    def prettyprint(self): #prints the contained microstates in a pretty manner, for viewing.
        print("\n".join(self.notes)) #print the notes
        for microstate in self.actual_microstates:
            print("- "*20)
            print_microstate(microstate, self.width)
    def get_g_of_r(self): #returns a dict of edge_length:g(r), just like the older version
        pairs_total = {} #dict of every possible pair length and the number of times it exists in the microstates, in form dsquared:count
        pairs_filled = {} #just like pairs_total, except only includes the actual filled pairs, in form dsquared:count

        for microstate in self.actual_microstates:
            for cell_A, cell_B, distance in all_cell_pairs_tiled(self.width, self.positions):
                if(microstate[cell_A] and microstate[cell_B]):
                    pairs_filled[distance] = pairs_filled.get(distance, 0) + 1
                pairs_total[distance] = pairs_total.get(distance, 0) + 1

        output = {}
        for distance in pairs_total:
            output[distance] = pairs_filled.get(distance,0) / pairs_total[distance]

        output = {x:output[x] for x in sorted(output)} #sort output by distances (for nicer displaying - really, this step can be skipped)

        return output

    def average_population(self): #returns the average population of the microstates
        return sum([sum(microstate) for microstate in self.actual_microstates]) / len(self.actual_microstates)
    
    def to_csv(filename): #saves a CSV of the results (including g(r)!) for easy exporting into excel
        output = ""
        #MAKE THIS WORK!

def load_fluid_results(filename):
    output = fluid_results()
    output.load_from(filename)
    return output

def fluid_results_example(): #returns an example fluid_results object that's roughly like our normal simulations
    test = simple_fluidproblem(20, 4, 1, 0.01, name="fluid_results Example Problem")
    output = get_fluid_results(test, 20, "gradient_descent")
    return output

def dict_to_csv(x, filename, label=None): #saves a dict to a CSV file that can be opened in excel
    rows = []
    if(label != None): rows.append(label)
    for y in x:
        rows.append(str(y) + ',' + str(x[y]))
    output = '\n'.join(rows)
    with open(filename + ".csv", "w") as savefile:
        savefile.write(output)