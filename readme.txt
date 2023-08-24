(...I figured I'd make a good documentation file, for everyone looking to read this)


>           WORKFLOW MODEL:

    1: Import model, the Python file that contains all classes and functions needed to perform experiments.
        -You shouldn't be editing model.py! You should only be making your own python files and importing model, so we can keep the classes seperate from the experiments.

    2: Make a fluid_problem, either by using , or by manually assigning potentials using the class functions.
        (optional step: save/load your fluid_problem)

    3: Use a sampler to sample fluid_problem and generate a fluid_results object.
        (again, fluid_results can be saved/loaded freely. the sampler should automatically save the results, too)

    4: Use whatever analytical tool you please to analyze fluid_results. 
        -You can find g(r) by calling XXXX.get_g_of_r(), where XXXX is the results object.
        -You can

>           THINGS TO AVOID:
    -Wildcard imports (from XXXX import *);             They are painful to keep track of
    -Modifying model.py in any way;                     We should try and keep all of our copies identical!
    -Running model.py;                                  Since it only defines classes and things, it shouldn't do anything!


>           WHY DID YOU DO IT LIKE THIS:

-The current version is heavily limited by its qmatrix creating system (networkx is not needed, and we have a whole load of classes that are incredibly hard to read with a tangled web of inheritences.)
-This program is much less flexible but significantly shorter (200 lines vs 1500 lines). Although this will need modification to run 3D simulations, this can be much more easily peer-reviewed.
-I've left room in the code to potentially stagger the cells in the future - Always use the positions dict to get the position of a cell! Don't just assume position from cell_ID!

By the way - let's try and keep this all on one branch! Just make your pull requests and I'll add them in, so we all stay in sync. 
The reason I've moved all of the definitions to "model.py" is so we can all work off of the same simulation model while doing our own experiments.

>           WHAT SHOULD I DO 
