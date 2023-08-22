#this file uses model.py to generate a 2D BQM and sample it, pretty much identically to what the current version of CapUQM does.
import model

width = 70
sigma = 4
epsilon = 1
chemical_potential = 0.01

num_runs = 50
sampler_name = "gradient_descent" #see model.py for a list of valid solvers (near line 95)

#Simulation begins here!

def lennardjones(s, e, r): #sigma, epsilon, radius (distance). note: this is NOT a class!
    dist_term = (r/s)**-6
    return 4 * e * (dist_term * (dist_term - 1)) #unreduced!

print("Generating problem...")
problem = model.fluid_problem(width, name="Simple Problem") #name is optional

print("Assigning chemical potentials...")
problem.assign_linear_potentials(chemical_potential)

print("Assigning lennard jones potentials (this may take a while)...")
problem.assign_quadratic_potentials(lambda x:lennardjones(sigma, epsilon,x), note="Assigning lennard-jones potentials - s="+str(sigma))

print("Converting to BQM and sampling...")
results = model.get_fluid_results(problem, num_runs, sampler_name)

#print("Saving results as file...")
#results.save_as("simple results.results")
#exit()

print("Finding g(r) (this may take a while)...")
gr = results.get_g_of_r() #this is a dict

print("Done sampling!")
#at this point, have fun! 

#find and print g(r) in a nice way
print("r                    g(r)")
for r in gr:
    if(r > sigma*4): break #this will cut off the print statement, but the full results will still be in the CSV file
    print(str(r).ljust(20) + str(gr[r]))

#export g(r) to a csv file (along with the notes)
model.dict_to_csv(gr, "simple g(r)", label = '; '.join(results.notes) + " AVERAGE POPULATION:" + str(results.average_population()))
# todo: turn this nightmare of a line into a fluid_results function