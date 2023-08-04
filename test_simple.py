#A quick file to test fluidmodel and confirm it works.

import fluidmodel as fm #please don't use "from fluidmodel import *", let's try and keep namespaces seperate

#EDIT THE SETTINGS HERE
dim = 14
sigma = 1.4314
epsilon = 1.0
chem_pot = 0.01
cutoff = 4*sigma

print("Creating problem...")
testproblem = fm.problem_2D(dim)

print("Assigning chemical potentials...")
testproblem.assignPotentials_chemical(chem_pot)

print("Assigning LJ potentials...")
testproblem.assignPotentials_LJ(sigma, epsilon, cutoff)

print("Generating BQM...")
testBQM = testproblem.getBQM()

print("BQM has", testBQM.num_interactions, "interactions.")
print("Running sampler...")

testresults = fm.solve_BQM(
    testBQM, 
    "sim", 
    label=fm.runlabel(dim, s=sigma, c=cutoff)
    )

fm.print_sampleset(testresults)
