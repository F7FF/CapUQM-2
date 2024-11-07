#DO NOT WRITE YOUR EXPERIMENT IN THIS FILE! Make a new file for your experiment, and use "import model" to use all of these classes.

import json #used to save/load things
import datetime #used to timestamp samples/problems
import random


#do NOT use "from x import *" - importing normally will make module management so much easier!
import dimod, dwave.samplers, dwave.system

pi = 3.1415926 #let's just define this right away, just in case g(r) calculations need it

#import dwave.system
#import dwave.samplers


"""------------- DEFINING fluid_problem CLASS HERE, AND ITS LOADING FUNCTION-------------"""

class fluid_problem: #EXCLUSIVELY FOR 2D PROBLEMS - 3D PROBLEMS WILL REQUIRE A REWRITE
    def __init__(self, width):
        self.linearconnections = {}         #Contains all self-loops (think: chemical potential) in form cellID:strength. USED TO CREATE BQMS!
        self.quadraticconnections = {}      #Contains all cross-cell connections in form (cellID_1, cellID_2):strength. USED TO CREATE BQMS!
        self.positions = {}                 #maps coordinates to cell IDs, using form cellID:(x,y).
        self.metadata = {}                  #keeps track of all parameters used to generate me (width, sigma, epsilon, temp, etc) - whenever a change is made, add any variables here!
        #now, set up the positions and width metadata
        for y in range(width):
            for x in range(width):
                self.positions[y*width+x] = (x,y)
        self.metadata["width"] = width

    def assign_linear_potentials(self, strength): #Assigns a linear connection of a certain strength to all nodes. I'm making this a class function because it's pretty obvious and we don't need to worry about varying between versions
        for cellID in self.positions:
            self.linearconnections[cellID] = float(strength)

    def assign_quadratic_potentials(self, function): #Assigns potentials to all pairs of nodes - includes PBC. Ask Max if you'd like an explanation!
        #let's iterate over pairs of cells
        for cell_A, cell_B, distance in all_cell_pairs_tiled(self.metadata["width"], self.positions):
            #now we find the potential for that distance and add it on to the quadratic connection for these two cells
            self.quadraticconnections[(cell_A, cell_B)] = self.quadraticconnections.get((cell_A, cell_B), 0) + function(distance)
            #self.quadraticconnections.get is simply returning the current potential, or 0 if it doesn't exist (just like defaultdict)

    def save_as(self, filename): #this is essentially just a JSON file, but it works nicely. It's a little bloated, though.
        data = {
            "linearconnections":self.linearconnections,
            "quadraticconnections":self.quadraticconnections,
            "positions":self.positions,
            "metadata":self.metadata
            }
        with open(filename + ".problem", "w") as savefile:
            json.dump(data, savefile)

def load_fluid_problem(filename):
    output = fluid_problem(0) #create the object we'll return later
    with open(filename + ".problem", "r") as f:
        param = json.load(f)
    #load all the other data in
    output.linearconnections = param["linearconnections"]
    output.quadraticconnections = param["quadraticconnections"]
    output.positions = param["positions"]
    output.metadata = param["metadata"]
    return output

"""
------------- DEFINING get_fluid_results FUNCTION HERE-------------

VALID SAMPLER NAMES:

"hybrid"                -hybrid solver (part QPU, part CPU)
"gradient_descent"      -CPU solver that simply performs a gradient descent algorithm. Very fast, gives reasonable results.
"sim"                   -CPU solver that tries to simulate QPU-ish performance. Takes a little while to load, but gives good results.

"random"                -CPU solver that simply returns a random microstate (50% chance any one cell is filled). Useless results.
"exact"                 -exact CPU solver. impractical for widths greater than 5


POSSIBLE FUTURE IDEAS:
"tabu"                  -use the dwave-system tabu solver (seems to be a good quantum simulator)
"tree"                  -uses dwave-samplers tree decomposition solver... might be fast? maybe?
"QPU"                   -the pure QPU solver from the older versions (only works up to 13x13). currently broken and I don't know why lmao
"fastannealer"          -???? no idea, let's see what it does

part of the reason this codebase exists is because different solvers give different datatypes and require different parameters - hence why I'm using custom classes for everything, so I can actually know what's going on
"""

def get_fluid_results(current_fluid_problem, number_of_samples, sampler_name): #Accepts a fluid_problem and returns a fluid_results. This is a black box in order to control the nightmare that is dimod/dwave-system/dwave-samplers.
    assert isinstance(current_fluid_problem, fluid_problem) #didya put a BQM in here? you should be giving this a fluid_problem!

    output = fluid_results() #prepare the object to return information to

    output.metadata = current_fluid_problem.metadata #pass variables down
    output.metadata["number of samples"] = number_of_samples
    output.metadata["sampler"] = sampler_name

    output.positions = current_fluid_problem.positions

    BQM = dimod.BinaryQuadraticModel(current_fluid_problem.linearconnections,
                                     current_fluid_problem.quadraticconnections,
                                     vartype='BINARY')


    #now - check what sampler we want to use
    if(sampler_name == "gradient_descent"):
        sampler = dwave.samplers.SteepestDescentSolver()
        sampleresults = sampler.sample(BQM, num_reads=number_of_samples)
    elif(sampler_name == "random"):
        sampler = dwave.samplers.RandomSampler()
        sampleresults = sampler.sample(BQM, num_reads=number_of_samples)
    elif(sampler_name == "sim"): #simulated annealer
        sampler = dwave.samplers.SimulatedAnnealingSampler()
        sampleresults = sampler.sample(BQM, num_reads=number_of_samples)
        output.actual_microstates = [tuple(cell for cell in result.values()) for result in sampleresults]
    elif(sampler_name == "QPU"):
        sampler = dwave.system.LazyFixedEmbeddingComposite(dwave.system.DWaveSampler())
        sampleresults = sampler.sample(BQM, num_reads=number_of_samples) #might need chain strength?
        output.actual_microstates = [tuple(cell for cell in result.values()) for result in sampleresults]
    else: #if the sampler name was invalid
        raise Exception(sampler_name + " is not a valid sampler name!")

    #IF THE SAMPLER RETURNED A DIMOD-ISH OBJECT, PROCESS IT INTO fluid_results HERE! (keeping this all together for neatness)
    if(sampler_name in ("gradient_descent", "random", "sim")):
        output.actual_microstates = [tuple(cell for cell in result.values()) for result in sampleresults]
        output.energies = [float(x.energy) for x in sampleresults.data(fields=['energy'])]

    #after all the sampling and preparing is done, return the final fluid_results object.
    return output

def print_microstate(microstate, dim): #Makes a pretty little thing to visualize a single 2d microstate. This ignores stagger.
    assert isinstance(microstate, (list, tuple)) #if this fails, you're passing print_microstate the wrong datatype!
    for y in range(dim):
        for x in range(dim):
            print(("░", "█")[microstate[y*dim+x]], end='')
        print("")


"""------------- DEFINING fluid_results CLASS HERE, AND ITS LOADING FUNCTION-------------"""

class fluid_results:                        #Note that this is immutable and should not be altered after being created by get_samples
    def __init__(self):                     #You shouldn't be initializing this on its own - this class should only be created by get_samples() or loading from file! It's a bit of a nightmare to initialize
        self.actual_microstates = []        #THIS COULD BE A LIST OR A TUPLE! Contains tuples, each tuple is a microstate that was returned by the sampler
        self.energies = []                  # Contains the energies of each microstates, zipped to acutal_microstates (ie. energies[4] and actual_microstates[4] will refer to the same state)
        self.metadata = {}                  #same functionality as fluid_problem.metadata
        self.positions = {}                 #maps coordinates to cell IDs, using form cellID:(x,y).
    
    def save_as(self, filename): #Saves to a filename so it can be loaded later. UNTESTED!
        output = ""
        for microstate in self.actual_microstates:
            for cell in microstate:
                output += str(cell)
            output += "\n"
        #now, add the parameters line
        parameters = {"positions":self.positions, "metadata":self.metadata, "energies":self.energies}
        output += json.dumps(parameters)
        with open(filename + ".results", "w") as savefile:
            savefile.write(output)

    def prettyprint(self): #prints the contained microstates in a pretty manner, for viewing.
        print("\n".join([x + "; " + str(self.metadata[x]) for x in self.metadata])) #print the metadata too
        for x in range(len(self.actual_microstates)):
            print("- "*20)
            print("Population:", sum(self.actual_microstates[x]), ", Energy: ", self.energies[x])
            print_microstate(self.actual_microstates[x], self.metadata["width"])
    
    def get_pair_lengths(self): #returns a dict of edge_length:(number of filled pairs/number of possible pairs) for each length in the sample set
        #In future, this might be rewritten in either C or numpy arrays, to save time. This is the step that takes the longest!
        try: #this caching step will save a LOT of time. but note - the cache is NOT saved with the file! DO NOT PULL DIRECTLY FROM CACHE
            return self._secret_pairlength_cache
        except:
            pairs_total = {} #dict of every possible pair length and the number of times it exists in the microstates, in form dsquared:count
            pairs_filled = {} #just like pairs_total, except only includes the actual filled pairs, in form dsquared:count

            for microstate in self.actual_microstates:
                for cell_A, cell_B, distance in all_cell_pairs_tiled(self.metadata["width"], self.positions):
                    if(microstate[cell_A] and microstate[cell_B]):
                        pairs_filled[distance] = pairs_filled.get(distance, 0) + 1
                    pairs_total[distance] = pairs_total.get(distance, 0) + 1
            output = {}
            for distance in pairs_total:
                output[distance] = pairs_filled.get(distance,0) / pairs_total[distance]

            output = {x:output[x] for x in sorted(output)} #sort output by distances (for nicer reading - really, this step could be skipped)

            #finally - return values and exit
            self._secret_pairlength_cache = output
            return output
    
    def get_g_of_r(self):
        pairs = self.get_pair_lengths()

        above_5_sigma = [pairs[r] for r in pairs if r>5* self.metadata["sigma"]]
        average_above_5_sigma = sum(above_5_sigma) / len(above_5_sigma)

        output = {r:pairs[r]/average_above_5_sigma for r in pairs} #normalize all values using the average calculated
        self.metadata["g(r) normalizing method"] = "normalized by dividing all g(r) by the average >5*sigma, " + str(average_above_5_sigma)

        return output

    def average_population(self): #returns the average population of the microstates
        return sum([sum(microstate) for microstate in self.actual_microstates]) / len(self.actual_microstates)

    def to_csv(self, filename): #saves a CSV of the useful results (including g(r)!) for easy exporting into excel. this may take a while, as it calculates g(r)!
        output = ""
        #make the g(r) table
        output += "r,g(r)\n"
        gr = self.get_g_of_r()
        for r in gr:
            output += str(r) + ',' + str(gr[r]) + '\n'
        
        #now add the metadata header

        header = ''
        header += '\n'.join([str(data) + "," + str(self.metadata[data]) for data in self.metadata]) #take all metadata points - make them into neat lines
        header += '\nAVERAGE POPULATION:,' + str(self.average_population()) + '\n'
        
        with open(filename + ".csv", "w") as csvfile:
            csvfile.write(header + output)

def load_fluid_results(filename):
    output = fluid_results()
    with open(filename + ".results", "r") as f:
        for line in f:
            if(line[0] == '{'): #if this is the parameter line
                param = json.loads(line)
            else: #if this is a microstate line
                output.actual_microstates.append([int(letter) for letter in line[:-1]]) #trim off the last letter (the newline)
    #load all the other data too
    output.positions = param["positions"]
    output.metadata = param["metadata"]
    output.energies = param["energies"]
    return output


"""------------- DEFINING GLOBAL FUNCTIONS HERE -------------"""

def standard_fluidproblem(width, sigma, epsilon, chem_potential): #This constructs a fluid_problem from the basic parameters - this should be identical to the older version's model. TO CHANGE MODELS - MAKE A NEW FUNCTION!
    #create the fluid problem
    output = fluid_problem(width)
    #leave some metadata so we know what we did
    output.metadata["function"] = "Lennard Jones: s=" + str(sigma) + ", e=" + str(epsilon)
    output.metadata["sigma"] = sigma
    output.metadata["epsilon"] = epsilon
    output.metadata["chemical potential"] = chem_potential
    #assign potentials
    def LJ(s, e, r): #sigma, epsilon, radius (distance)
        dist_term = (r/s)**-6
        return 4 * e * (dist_term * (dist_term - 1)) #unreduced!
    output.assign_linear_potentials(chem_potential)
    output.assign_quadratic_potentials(lambda x:LJ(sigma, epsilon,x))
    return output

def all_cell_pairs_tiled(width, positions): #this yields (cell_A, cell_B, distance) for all pairs of cells INCLUDING TILING! ask Max for a full explanation
    positions_list = tuple([x for x in positions.values()])
    for cell_A in positions:
        cell_A_position = positions_list[cell_A]
        for cell_B in positions:
            if cell_A >= cell_B:
                continue #do NOT self connect, and do not accidentally connect nodes twice! By only connecting if cell B's ID is higher than cell A's ID, nodes will never be connected twice.
            cell_B_true_position = positions_list[cell_B]
            #now: we iterate over the nine possible microstate tiles (eight neighbors plus the "centre" tile). We could increase this for more accuracy (ie. iterating over 25 tiles), but at the cost of speed, and LJ drops off fast anyway.
            for dy in range(-1, 2):
                for dx in range(-1, 2):
                    #we now calculate where cell B should be in this tile
                    cell_B_position = (cell_B_true_position[0] + width*dx, cell_B_true_position[1] + width*dy)
                    #now we calculate distance between cell_A_position and cell_B_position (the one that is in a different tile)
                    distance = ((cell_B_position[0]-cell_A_position[0])**2 + (cell_B_position[1]-cell_A_position[1])**2)**0.5
                    yield (cell_A, cell_B, distance)
