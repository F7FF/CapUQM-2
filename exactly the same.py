#This file does exactly the same thing as the first CapUQM of the project.

import model

width = 30
sigma = 4
epsilon = 1
chemical_potential = 0.01
num_runs = 2
sampler_name = "gradient_descent" #see model.py for a list of valid solvers (near line 95)

#Simulation begins here!

print("Generating problem...")
problem = model.standard_fluidproblem(width, sigma, epsilon, chemical_potential) #This creates a fluid_problem object, with parameters identical to our current version.

print("Converting to BQM and sampling...")
results = model.get_fluid_results(problem, num_runs, sampler_name)
#results is a fluid_results object

del problem #frees up some RAM when working with wide simuations (as much as 1.5 gb - sheesh)

print("Done sampling!")
#at this point, do whatever you like with the result object. have fun!

print("Finding g(r)... (this may take a while)")

if(width == 70):
    print("This should take around",num_runs*2,"minutes...")

gofr = results.get_g_of_r()

print("r - g(r)")
for r in gofr:
    if(r < width): #if we don't do this, it'll completely flood the terminal
        print(r, ',', gofr[r])

print("Done!")