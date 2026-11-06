#================================
#|     Luke Bain                |
#|     Undergraduate Thesis     |
#|     November 2026            |
#================================

This file contains a quick start guide and a short description of the files present
in this tandem-wing design optimisation program.

===========
Quick Start
===========
AVL should be installed in the same directory as the .py files.

train.py should be run first

opt.py should be run after train.py has been run

opt.py will create a text file in the same directory, containing the output

test.py can be run to solve one geometry at a time, after manually transferring 
the geometry definition from opt.py's output text file into test.py

changes to the bounds of the geometry should be edited into train.py 
(run train.py again)

changes to the mission specifications should be edited into opt.py and test.py

to specify a starting point for optimisation, uncomment line 75 in opt.py
and comment the 'while' loop at line 202



========
train.py
========
train.py creates the lists of minimums and maximums (bounds) for the state vector.


========
opt.py
========
opt.py calls the optimisation methods, and writes the optimised output to a text file.
to choose a different optimisation method, and change the solver options, see lines 141-176
to change the number of optimisations to complete, change the definition of the while loop at line 191

at every function call, opt.py will print a set of 5 numbers, then a set of 2 numbers

the line of 5 numbers represent:
the solve number, 
number of times AVL has crashed,
the total number of function calls,
the total number of mass evaluations in the solve,
the total number of stability evaluations in the solve

the line of 2 numbers represent:
the mass at this function call,
the stability at this function call

========
base.py
========
base.py is the solver that takes an aircraft geometry, payload, range, and endurance
it is not intended to be run stand-alone (see test.py)
to see the aircraft geometry every time base.py is run, uncomment line 363
to print the final mass breakdown and performance, uncomment lines 911-930


========
test.py
========
test.py will call base.py using the aircraft geometry and mission requirements specified