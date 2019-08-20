#!/usr/bin/env python3

import os, re, sys,math

from numpy import *

# frequency of envt 2
step = 0.02
#freq_patch_2 = list(arange(0+step,1.0,step))

freq_patch_2 = [ 0.4 ]

step = 0.02
#freq_patch_2 = list(arange(0.3,0.7 + step,step))

# avarage switch rate
#sbar = list(arange(-1.5, 0.5, 2.0/50))
sbar = [ -1.0, -1.42 ]

#k = [ 0.01, 0.25, 0.35, 0.5, 0.65, 0.95 ]
k = [ 0 ]

# first is fecundity, second is mortality
fecmortcombis = [[ 2.0, 1.0 ], [1.0, 2.0 ] ]

exe = "./xp1p2n"

ctr = 0

step = 0.1

# initial values along the edges of the square
pi_init = list(arange(step,1.0 - step + 0.05,step))

initvals_p1p2 = [[0,0],[0,1.0],[1.0,0],[1.0,1.0]]

for pi_init_i in pi_init:
    initvals_p1p2 += [[ 0, pi_init_i]] 
    initvals_p1p2 += [[ 1.0, pi_init_i]]
    initvals_p1p2 += [[ pi_init_i, 0]] 
    initvals_p1p2 += [[ pi_init_i, 1.0]]


d = [ 0.25, 1.0 ]

error_prob = 0.1

nrep = 5

for rep_i in range(0,nrep):
    for d_i in d:
        for combi_i in fecmortcombis:

            ca1 = 1.0
            cm1 = combi_i[1]
            ca2 = 1.0
            cm2 = combi_i[1]

            Fa1 = combi_i[0]
            Fm1 = 1.0
            Fa2 = combi_i[0]
            Fm2 = 1.0

            for pi_init_i in initvals_p1p2:

                p1_init_i = pi_init_i[0]
                p2_init_i = pi_init_i[1]


                for f2 in freq_patch_2:
                    for sbar_i in sbar:
                        for k_i in k:
                            s2 = sqrt(((1.0-f2)/f2) * 10**(2*sbar_i))
                            s1 = 10**(2*sbar_i) / s2

                            ctr += 1
                            print("echo " + str(ctr))

                            print(exe + " " 
                                    + str(s1) + " " 
                                    + str(s2) + " " 
                                    + str(cm1) + " " 
                                    + str(cm2) + " " 
                                    + str(ca1) + " " 
                                    + str(ca2) + " " 
                                    + str(Fm1) + " "  
                                    + str(Fm2) + " "  
                                    + str(Fa1) + " " 
                                    + str(Fa2) + " "  
                                    + str(k_i) + " " 
                                    + " 0.02 0.01 " 
                                    + str(d_i) + " " 
                                    + str(p1_init_i) + " " 
                                    + str(error_prob) + " " 
                                    + str(p1_init_i) + " " 
                                    + str(p2_init_i) + " " 
                                    )
