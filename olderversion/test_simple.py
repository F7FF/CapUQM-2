#A quick file to test fluidmodel and confirm it works.

import fluidmodel as fm #please don't use "from fluidmodel import *", let's try and keep namespaces seperate

#EDIT THE SETTINGS HERE
dim = 50
sigma = 5
epsilon = 1.0
chem_pot = 0.01
cutoff = 20

print("Creating problem...")
testproblem = fm.problem_2D(dim)

print("Assigning chemical potentials...")
testproblem.assignPotentials_chemical(chem_pot)

print("Assigning LJ potentials...")
testpotentials = fm.LennardJonesCached(sigma, epsilon)

testproblem.assignPotentials_wrap(testpotentials, cutoff)

"""
   8 6 7 8 6
   2 0 1 2 0
   5 3 4 5 3
   8 6 7 8 6
   2 0 1 2 0
"""

print("Generating BQM...")
testBQM = testproblem.getBQM()

print("BQM has", testBQM.num_interactions, "interactions and", testBQM.num_variables, "variables.")
print("Running sampler...")

testresults = fm.solve_BQM(
    testBQM, 
    "hybrid",
    label=fm.runlabel(dim, s=sigma, c=cutoff)
    )

fm.print_sampleset(testresults)
