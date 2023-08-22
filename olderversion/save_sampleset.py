#A program to make a fluidmodel and save its results to a large JSON

print("Importing...")
import fluidmodel as fm 
import json

#EDIT THE SETTINGS HERE
dim = 4
sigma = 5
epsilon = 1.0
chem_pot = 0.01
cutoff = 20
sampler = "sim" # "sim", "hybrid", "QPU", "exact", "random"

print("Creating problem...")
testproblem = fm.problem_2D(dim)

print("Assigning chemical potentials...")
testproblem.assignPotentials_chemical(chem_pot)

print("Assigning LJ potentials...")
testpotentials = fm.LennardJonesCached(sigma, epsilon)

testproblem.assignPotentials(testpotentials, cutoff) #DON'T WRAP

print("Generating BQM...")
testBQM = testproblem.getBQM()

print("BQM has", testBQM.num_interactions, "interactions and", testBQM.num_variables, "variables.")
print("Running sampler...")

testresults = fm.solve_BQM(
    testBQM, 
    sampler,
    num_reads=2,
    label=fm.runlabel(dim, s=sigma, c=cutoff)
    )

print(type(testresults)) #<class 'dimod.sampleset.SampleSet'>

fm.print_sampleset(testresults)

def save_sampleset(filename, sampleset, dim, sigma, epsilon, **kwargs): #kwargs is for any other arguments you'd like to save
    output = {}
    output["dim"] = dim
    output["sigma"] = sigma
    output["epsilon"] = epsilon

    rawsamples = []
    for sample in sampleset:
        rawsamples.append(tuple(sample.values()))
        rawsamples[-1] = tuple([int(x) for x in rawsamples[-1]]) #convert all the numpy.int8s to builtin ints, so JSON can serialize them
    
    output["data"] = rawsamples

    with open(filename, "w") as write_file:
        json.dump(output, write_file) #, indent=4)

print(save_sampleset("test", testresults, dim, sigma, epsilon))