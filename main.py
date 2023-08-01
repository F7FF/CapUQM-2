#Quick rewrite of CapUQM - let's see if this works!

from simulation import *

test = QMgraph(5)

test.stagger(0.2)

potential_function = LennardJones(1.3, 1.0)
test.assign_potentials(potential_function, cutoff=2.5)
test.assign_chemical(1,1)

test.visualize()