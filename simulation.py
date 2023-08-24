#this file uses model.py to generate a 2D BQM and sample it, pretty much identically to what the current version of CapUQM does.
import model

width = 20
sigma = 4
epsilon = 1
chemical_potential = 0.01

num_runs = 3
sampler_name = "gradient_descent" #see model.py for a list of valid solvers (near line 95)

#Simulation begins here!

def lennardjones(s, e, r): #sigma, epsilon, radius (distance). note: this is NOT a class!
    dist_term = (r/s)**-6
    return 4 * e * (dist_term * (dist_term - 1)) #unreduced!

print("Generating problem...")
problem = model.standard_fluidproblem(width, sigma, epsilon, chemical_potential) #This creates a fluid_problem object, with parameters identical to our current version.

print("Converting to BQM and sampling...")
results = model.get_fluid_results(problem, num_runs, sampler_name)

print("Done sampling!")
#at this point, do whatever you like with the result object. have fun!

#results.prettyprint() #for a quick peek to see if it came out right

if(False): #set me to True to save the microstates to a .results file for loading later
    print("Saving results as file...")
    results.save_as("output")
    #exit() #closes this entire program

if(True): #export g(r) to a csv file (along with the notes)
    print("Saving g(r) to CSV...")
    results.to_csv("output")
    print("Done saving!")
    exit()
