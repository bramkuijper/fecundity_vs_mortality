#!/usr/bin/env python3
# generate parameter combinations for iterations


import os, re, sys, math, itertools, copy

import numpy as np

import pandas as pd

# auxiliary expand grid function
def expand_grid(data_dict):

    for key, value in data_dict.items():
        assert(type(value) is type([]))


    # calculate total product of data frame
    rows = itertools.product(*data_dict.values())
    
    return(pd.DataFrame([row for row in rows],
        columns=data_dict.keys(),
        dtype=object))

# make the dataframe with x,y values ourselves
# here x is the frequency of environment 2
# here y is sbar :w

col_dict = {}

freq_nstep = 50
freq_min = 0.00001
freq_max = 0.99999
# generate all values of x2 (frequency of envt 2)
col_dict["x2"] = list(np.linspace(freq_min, freq_max, freq_nstep))

sbar_nstep = 200
sbar_min = -1.5 
sbar_max = 0.5 

# generate all values of y2 (sbar)
col_dict["y2"] = list(np.linspace(sbar_min, sbar_max, sbar_nstep))

col_dict["y2"] = [-1.4]

dat = expand_grid(col_dict)

d = [0.1 ]

exe = "./xnumsolve"


i = 0.0
t = 0.0

# constant reflecting genetic variation in 
# evolving traits
eul = 0.01

# constant that reflects the constant with which
# iterations of reproductive values get updated
eulrv = 0.0001

# fecundity vs mortality combinations
# first element is fecundity selection, second mortality
fecmort_combis = [ [2.0, 1.0], [1.0,2.0] ]

# trembling hand error in evolving traits
epsilon = 0.05

# initial values for p1 and p2
p1_init = [ 0.5 ]
p2_init = [ 0.5 ]

# maximum number of iterations
max_iter = 1e08


# initial values for patch freqs and relatedness coeffs
fval_init = " ".join([ str(1.0/6) for xi in range(0,6) ])

# initial values for the reproductive values
vval_init = " ".join([ str(1.0) for xi in range(0,8) ])

# initial values for the relatedness coefficients 
relval_init = " ".join([str(xi) for xi in [0.5]*6])

ctr = 0

for combi_i in fecmort_combis:

    # mortality selection
    cm1 = combi_i[1]
    cm2 = combi_i[1]
    ca1 = 1
    ca2 = 1

    # fecundity selection
    Fa1 = combi_i[0]
    Fa2 = combi_i[0]
    Fm1 = 1
    Fm2 = 1

    for p1_init_i in p1_init:
        for p2_init_i in p2_init:
            for index, row in dat.iterrows():

                f2 = float(row["x2"])
                sbar_i = float(row["y2"])

                s2 = math.sqrt(((1.0-f2)/f2) * 10**(2*sbar_i))
                s1 = 10**(2*sbar_i) / s2
                for d_i in d:
                    print("echo " + str(ctr))
                    print(exe + " " 
                            + str(max_iter) + " " 
                            + str(ca1) + " " 
                            + str(ca2) + " " 
                            + str(cm1) + " " 
                            + str(cm2) + " " 
                            + str(Fa1) + " " 
                            + str(Fa2) + " "  
                            + str(Fm1) + " "  
                            + str(Fm2) + " "  
                            + str(d_i) + " " 
                            + str(s1) + " " 
                            + str(s2) + " " 
                            + str(eulrv) + " " 
                            + str(epsilon) + " " 
                            + fval_init + " " 
                            + vval_init + " " 
                            + relval_init + " "
                            + str(p1_init_i) + " " 
                            + str(p2_init_i) + " " 
                            )
                    ctr+=1
                
