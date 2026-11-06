#================================
#|     Luke Bain                |
#|     Undergraduate Thesis     |
#|     November 2026            |
#================================
import sys

# import functions
from scipy.optimize import minimize
from scipy.optimize import differential_evolution
# from smt.surrogate_models import RBF
# from smt.sampling_methods import Random
# from smt.surrogate_models import KRG, QP
from base import get_mass
from pathlib import Path
import numpy as np
import avlwrapper as avl

import matplotlib.pyplot as plt

mass_payload = 1000 #kg
aircraft_range = 3000*1e3 #m
aircraft_endurance = 20*3600 #s

maximums = np.loadtxt("input_maximums")
minimums = np.loadtxt("input_minimums")
x = 0

if True: # proteus geometry
    fwng_xpos0	=-6
    fwng_zpos0	=0
    fwng_chrd0	=1.32
    fwng_inci0	=0
    fwng_span0	=3.78
    fwng_dihe0	=5
    fwng_swep0	=0
    fwng_chrd1	=1.132
    fwng_inci1	=0
    fwng_span1	=2.63
    fwng_dihe1	=5
    fwng_swep1	=0
    fwng_chrd2	=0.943
    fwng_inci2	=0
    fwng_span2	=2.08
    fwng_dihe2	=5
    fwng_swep2	=0
    fwng_chrd3	=0.755
    fwng_inci3	=0
    rwng_xpos0	= 1
    rwng_zpos0	=0
    rwng_chrd0	=1.51
    rwng_inci0	=0
    rwng_span0	=3.9
    rwng_dihe0	=0
    rwng_swep0	=0
    rwng_chrd1	=1.51
    rwng_inci1	=0
    rwng_span1	=4.63
    rwng_dihe1	=8
    rwng_swep1	=0
    rwng_chrd2	=0.943
    rwng_inci2	=0
    rwng_span2	=3.45
    rwng_dihe2	=-8
    rwng_swep2	=0
    rwng_chrd3	=0.755
    rwng_inci3	= 0
    tail_xpos	= 8
    tail_zpos	= -2
    tail_span	= 4.15
    tail_chrd0	= 1.89
    tail_chrd1  = 1.8
    tail_swep   = 0
#uncomment next line to use for first guess (and uncomment random generation below)
# x0 = [fwng_xpos0,	fwng_zpos0,	fwng_chrd0,	fwng_inci0,	fwng_span0,	fwng_dihe0,	fwng_swep0,	fwng_chrd1,	fwng_inci1,	fwng_span1,	fwng_dihe1,	fwng_swep1,	fwng_chrd2,	fwng_inci2,	fwng_span2,	fwng_dihe2,	fwng_swep2,	fwng_chrd3,	fwng_inci3,	rwng_xpos0,	rwng_zpos0,	rwng_chrd0,	rwng_inci0,	rwng_span0,	rwng_dihe0,	rwng_swep0,	rwng_chrd1,	rwng_inci1,	rwng_span1,	rwng_dihe1,	rwng_swep1,	rwng_chrd2,	rwng_inci2,	rwng_span2,	rwng_dihe2,	rwng_swep2,	rwng_chrd3,	rwng_inci3,	tail_xpos,	tail_zpos,	tail_span,	tail_chrd0, tail_chrd1, tail_swep]


n = 1
a = 1

def pred_func2(x0): #main function call
    global solves
    global fails
    global calls
    global n
    global a
    print(solves, fails, calls, n, a)
    calls = calls + 1
    n = n + 1
    x0 = np.array(x0)
    # print(list(x0))
    mass,stability = get_mass(x0, mass_payload, aircraft_range, aircraft_endurance)
    print(round(mass,1), round(stability,2))
    print("===")
    return mass

def output_constraint(x0): #stability constraint
    global solves
    global fails
    global calls
    global a
    global n
    print(solves, fails, calls, n, a)
    calls = calls + 1
    a = a + 1
    x0 = np.array(x0)
    mass,stability = get_mass(x0, mass_payload, aircraft_range, aircraft_endurance)
    print(round(mass,1), round(stability,2))
    print("===")
    return -stability  # SciPy keeps this >= 0

def output_constraint2(x0): #wing relative angle constraint
    fwng_xpos0	=	x0	[	0	]#default no touchie
    fwng_zpos0	=	x0	[	1	]
    rwng_xpos0	=	x0	[	19	]
    rwng_zpos0	=	x0	[	20	]

    zdiff = abs(fwng_zpos0-rwng_zpos0)
    xdiff = abs(fwng_xpos0-rwng_xpos0)
    theta = 30 #line between wings should be within this angle (deg) of horizontal
    return np.tan(np.deg2rad(theta))-(zdiff/xdiff)

def output_constraint3(x0): #wing relative x-pos constraint
    fwng_xpos0	=	x0	[	0	] #default no touchie
    fwng_chrd0	=	x0	[	2	]
    rwng_xpos0	=	x0	[	19	]
    rwng_chrd0	=	x0	[	21	]

    diff1 = abs(fwng_xpos0-rwng_xpos0)/fwng_chrd0
    diff2 = abs(fwng_xpos0-rwng_xpos0)/rwng_chrd0
    return min(diff1,diff2)-2.5 #keep wings more than 2.5 root chords apart


constraints = [
    {'type':'ineq','fun':output_constraint},          # 'ineq' means >= 0
    {'type':'ineq','fun':output_constraint2},
    {'type':'ineq','fun':output_constraint3}
]


#select method
# method = 'SLSQP'
# method = 'COBYQA'
method = 'COBYLA'
# method = 'L-BFGS-B'
# method = 'trust-constr'

#select method options

# solver_options = { # SLSQP
#     'maxiter': 30,
#     'ftol': 0.03,
#     'disp': True
# }

# solver_options = { # COBYQA
#     'maxiter': 30,
#     'disp': True,
#     'final_tr_radius': 0.001
# }

solver_options = { # COBYLA
    'maxiter': 50*42,
    'rhobeg': 2,
    # 'rhoend':0.01,
    'tol':0.001,
    }

# solver_options = { # L-BFGS-B
#     'ftol': 0.005,
#     'maxiter':30
#     }

# solver_options = { # trust-constr
#     'xtol': 0.005,
#     }

length = 44 #number of design vars

bnds = list(np.zeros(length))
# print(bnds)
y = 0
while y < length:
    bnds[y] = [minimums[y],maximums[y]]
    y = y + 1

random_id = (str(np.random.random_sample())[2:])#create a random number for output file name
solves = 1
fails = 0
calls = 0
while solves <= 5: #number of solves to complete
    # solves = solves + 1
    print("solve = ",solves)
    x = 0
    y = 0
    n = 1
    a = 1

    mass = 100001
    stability = 1

    while y < 0 and (mass > 20000 or stability > 0):
        x = 0
        x0 = list(np.zeros(44))
        x0 = np.array(x0)
        while x < 44:
            x0[x] = np.random.random_sample()*(maximums[x]-minimums[x]) + minimums[x]
            x = x + 1
        # print("y = ",y)
        try:
            mass, stability = get_mass(x0, mass_payload, aircraft_range, aircraft_endurance)
        except:#FileNotFoundError:
            print("AVL crashed!")
        else:
            print("mass = ",round(mass,1))
            print("stability = ",round(stability,2))
        print("===")
        calls = calls + 1

    print("Initial x0 = ", x0)

    try:
        res = minimize(pred_func2, x0,method=method, bounds=bnds, constraints=constraints, options=solver_options)#local, direct optimisation
    except:
        print("AVL crashed!")
        fails = fails + 1
    else:
        # print(res)
        solves = solves + 1
        xref2 = res.x
        fref2 = res.fun
        neval2 = res.nfev
        print(
            "The reference solution with "+method+" is ",
            list(xref2),
            "resulting in ",
            fref2,
            " obtained with ",
            neval2,
            "function evaluations",
        )
        print(res.message)

        # find current working directory
        output_dir = Path.cwd()
        file_name = "output_"+random_id+".txt"
        #print(output_dir)

        # make filepath
        file_path = output_dir / file_name

        # create directory
        output_dir.mkdir(parents=True, exist_ok=True)

        # save to file
        with open(file_path, "a", encoding="utf-8") as file:
            file.write("# payload = "+str(mass_payload)+" kg")
            file.write("\n")
            file.write("# range = "+str(aircraft_range)+" m")
            file.write("\n")
            file.write("# endurance = "+str(aircraft_endurance)+" s")
            file.write("\n")
            file.write("# mass = "+str(fref2)+" kg")
            file.write("\n")
            file.write("# "+str(res.message))
            file.write("\n")
            # file.write("# By random start") #used for extra note
            # file.write("\n")
            file.write("# With method "+method)
            file.write("\n")
            file.write("# With n="+str(n)+" a="+str(a))
            file.write("\n")
            file.write("# geom = "+str(list(xref2)))
            file.write("\n")
            file.write("#=============================")
            file.write("\n")
            file.close()

#================================
#|     Luke Bain                |
#|     Undergraduate Thesis     |
#|     November 2026            |
#================================