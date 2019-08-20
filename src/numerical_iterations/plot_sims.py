#!/usr/bin/env python3

import sys, re, os.path, itertools
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import pandas as pd

from pylab import *

params={'axes.linewidth' : .5}
rcParams.update(params)

# get the filename from the command line
filename = sys.argv[1]

# first figure out last line of dataset
f = open(filename)
fl = f.readlines()
f.close()


# find the parameters
parline = -1

for idx, line in enumerate(fl):
    if re.match("^Ca1.*",line) != None:
        parline = idx - 1;
        break;

if parline == -1:
    nrow = len(fl)
else:
    nrow = parline -2  


data = pd.read_csv(filename, sep=";", nrows=nrow)


# initialize and specify size 
fig = plt.figure(figsize=(10,10))

num_rows = 5

# names of columns in dataframe
colnames = list(data.columns.values)

# add first subplot depicting % type 1 offspring
plt.subplot(num_rows,1,1)

plt.plot(
        data["time"],data["p1"],'b',
        data["time"],data["p2"],'r',
        data["time"],data["mt"],'c',
        )

plt.ylabel(r'Prob. offspring is hawk')
plt.legend((r'$p_{1}$',
                r'$p_{2}$',
                r'$m$'
                ),
                bbox_to_anchor=(1.3,1.0))
plt.ylim(-0.05,1.05)

# add 2nd subplot depicting patch frequencies
plt.subplot(num_rows,1,2)

plt.plot(
        data["time"],data["fa0"],'y',
        data["time"],data["fa1"],'g',
        data["time"],data["fa2"],'magenta',
        data["time"],data["fb0"],'b',
        data["time"],data["fb1"],'r',
        data["time"],data["fb2"],'c',
        linewidth=1)

plt.legend((
                r'$f_{a0}$',
                r'$f_{a1}$',
                r'$f_{a2}$',
                r'$f_{b0}$',
                r'$f_{b1}$',
                r'$f_{b2}$')
                ,bbox_to_anchor=(1.3,1.0)
                )

plt.ylabel(r'Patch freqs')
plt.ylim(-0.05,1.05)


# add 3rd subplot depicting relatedness
plt.subplot(num_rows,1,3)

plt.plot(
        data["time"],data["raa12"],'r',
        data["time"],data["raa22"],'g',
        data["time"],data["ram11"],'b',
        data["time"],data["ram21"],'y',
        data["time"],data["rmm10"],'magenta',
        data["time"],data["rmm20"],'black',
        linewidth=1)

plt.legend((
                r'$r_{aa,e_{1}}$',
                r'$r_{aa,e_{2}}$',
                r'$r_{am,e_{1}}$',
                r'$r_{am,e_{2}}$',
                r'$r_{mm,e_{1}}$',
                r'$r_{mm,e_{2}}$')
                ,bbox_to_anchor=(1.3,1.0)
                )

plt.ylabel(r'Relatedness')
plt.ylim(-0.05,1.05)


# add 4th subplot depicting reprovals
plt.subplot(num_rows,1,4)

plt.plot(
        data["time"],data["v100"],'y',
        data["time"],data["v101"],'g',
        data["time"],data["v110"],'k',
        data["time"],data["v111"],'m',
        linewidth=1)

plt.legend((
                r'$v_{mm,e_{1}}$',
                r'$v_{ma,e_{1}}$',
                r'$v_{am,e_{1}}$',
                r'$v_{aa,e_{1}}$')
                ,bbox_to_anchor=(1.3,1.0)
        )

plt.ylabel(r'Reproductive value')

# add 5th subplot depicting reprovals
plt.subplot(num_rows,1,5)

plt.plot(
        data["time"],data["v200"],'y',
        data["time"],data["v201"],'g',
        data["time"],data["v210"],'k',
        data["time"],data["v211"],'m',
        linewidth=1)

plt.legend((
                r'$v_{mm,e_{2}}$',
                r'$v_{ma,e_{2}}$',
                r'$v_{am,e_{2}}$',
                r'$v_{aa,e_{2}}$')
                ,bbox_to_anchor=(1.3,1.0)
        )

plt.ylabel(r'Reproductive value')

graphname = os.path.dirname(filename)
if graphname != '':
    graphname += "/"
graphname += "graph_" + os.path.basename(filename) + ".pdf"

plt.tight_layout()

plt.savefig(graphname,format="pdf",bbox_inches="tight")
