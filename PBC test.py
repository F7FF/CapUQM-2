#this should test whether or not PBC is working

import model

sigma = 1.2
epsilon = 1

width = 4

def lennardjones(s, e, r): 
    dist_term = (r/s)**-6
    return 4 * e * (dist_term * (dist_term - 1)) 

problem = model.fluid_problem(width, name="PBC Test") #name is optional
problem.assign_quadratic_potentials(lambda x:lennardjones(sigma, epsilon,x))

#this sets up a problem that should generate a perfect chessboard pattern

results = model.get_fluid_results(problem, 1, "exact")

results.prettyprint()


"""
if this looks like a perfect chessboard, PBC is not working! 
if the chessboard is missing squares and messed up, PBC is working

the ideal pattern with these settings is a perfect chessboard, but because width is an odd number, a perfect chessboard pattern would cause two cells to border along the edges, incurring a huge energy penalty.
if PBC is working, it will generate a messed up chessboard.
"""