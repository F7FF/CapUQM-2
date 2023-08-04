#Let's see if we can fit all of this into one file!

import dimod
import dwave.system
import dwave.samplers
import time

#PART 1 - SOME UTILITY FUNCTIONS

def LJ(s, e, r): #sigma, epsilon, radius (distance)
    dist_term = (r/s)**-6
    return 4 * e * (dist_term * (dist_term - 1)) #unreduced!

def point_distance(dim, a, b): #returns the distance (radius) between grid cell A and grid cell B
    return ((a//dim - b//dim)**2 + (a%dim - b%dim)**2)**0.5

def print_sample(sampledata): #Makes a pretty little thing to visualize a single SampleView object.
    dim = int((max(sampledata.keys()) + 1)**0.5)
    for y in range(dim):
        for x in range(dim):
            print(("░", "█")[sampledata[y*dim+x]], end='')
        print("")

def print_sampleset(sampleset): #same as above, but for an entire set.
    for sample in sampleset.samples():
        print_sample(sample)
        print("\n")

def runlabel(dim, **kwargs): #creates a nice little string useful for labelling QPU runs
    output = [str(dim) + "^2"]
    for x in kwargs:
        output.append(x + "=" + str(kwargs[x]))
    return ", ".join(output)

"""
class Stopwatch: #Used to time different sections of code. Works just like a real stopwatch. Note that this is global, though.
    starttime = 0.0
    def start():
        Stopwatch.starttime = time.time()
    def read(): 
        return time.time() - Stopwatch.starttime
    def lap(): #returns result and resets timer
        output = Stopwatch.read()
        Stopwatch.start()
        return output
"""

#PART 2 - DEFINING HELPER CLASSES (for easy creation of BQMs)

class problem_2D: #A wrapper that holds a potential graph and various misc. data, all in one trenchcoat. Useful for making BQMs.
    def __init__(self, dim):
        self.connections = {} #All the links between nodes, in form (nodeID1, nodeID2):strength
        self.linearconnections = {} #All of the self-linking nodes, in form nodeID:strength
        self.positions = {} #all of the "physical" positions of the nodes, in form nodeID:(x,y). Usually this will simply be aligned with the grid, but might be staggered in future.
        self.dim = dim
        self.cellcount = dim**2
        for x in range(dim):
            for y in range(dim):
                self.positions[dim*y + x] = (x,y)

    def assignPotentials(self, function, cutoff):
        maximum_offset = int(self.dim * (cutoff+1)) #this is the maximum number of nodes above A that you could feasibly add and stay inside the cutoff
        for a in range(self.cellcount):
            for b in range(a+1, min(a+maximum_offset, self.cellcount)):
                distance = point_distance(self.dim,a,b)
                if(distance <= cutoff):
                    self.connections[(a,b)] = function(distance)

    def assignPotentials_chemical(self, chemicalpotential):
        for x in range(self.cellcount):
            self.linearconnections[x] = chemicalpotential

    def getBQM(self): #attempts to return a dimod.BQM object from this problem_2D.
        return dimod.BinaryQuadraticModel(self.linearconnections, self.connections, vartype='BINARY')

#PART 3: BQM SOLVERS

def solve_BQM(BQM, method, num_reads=1, label=None):
    if method == "hybrid":
        esttime = BQM.num_variables / 206 #rough estimate of how many seconds it will take
        print("BQM ready for QPU. Est.", esttime, " seconds per run. Continue? Y/N")
        if("n" in input("").lower()): return
        sampler = dwave.system.LeapHybridSampler()
        return sampler.sample(BQM, label=label)
    elif method == "sim":
        sampler = dwave.samplers.SteepestDescentSolver()
        return sampler.sample(BQM, num_reads=num_reads)
    elif method == "exact":
        sampler = dimod.ExactSolver().sample
        return sampler(BQM)
    elif method == "random":
        sampler = dimod.RandomSampler().sample
        return sampler(BQM, num_reads=num_reads)
    else:
        raise Exception("Unknown solver, '", method, "'!")

#PART 4: STATISTICS AND ANALYSIS

def getlengths(sample): #Takes a single sample, returns a dict with each possible bond length and the number of times it occured, in form length:count
    dim = int((max(sample.keys()) + 1)**0.5) #get the dimensions of the sample
    output = {}
    s = list(sample.values()) #extract dict values and make a list representing the microstate
    for a in range(len(s)):
        if(s[a] == 0): continue
        for b in range(a+1, len(s)):
          if(s[b] == 0): continue
          distance = point_distance(dim, a, b)
          output[distance] = output.get(distance, 0) + 1
    return output

def getlengths_sampleset(sampleset): #same as getlengths, except it works over a large sampleset instead.
    output = {}
    for sample in sampleset.samples():
        result = getlengths(sample)
        for x in result:
            output[x] = output.get(x, 0) + result[x] #add up all values in the two dicts
    return output

def print_dict_excel(x): #print in a way that can be copied into excel
    for item in x:
        print(item, x[item])

#todo: make histogram work