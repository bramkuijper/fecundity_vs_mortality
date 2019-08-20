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
    if re.match("^d0.*",line) != None:
        parline = idx - 1;
        break;

nrow = len(fl)
skiprows = 0

if parline < 4:
    skiprows = 22
elif parline > 0:
    nrow = parline -2  

data = pd.read_csv(filename, sep=";", nrows=nrow,skiprows=skiprows)


# initialize and specify size 
fig = plt.figure(figsize=(10,10))

num_rows = 4

# names of columns in dataframe
colnames = list(data.columns.values)

# add first subplot depicting % type 1 offspring
plt.subplot(num_rows,1,1)

plt.plot(
        data["generation"],data["meanp1"],'b',
        data["generation"],data["meanp2"],'r',
        data["generation"],data["meand"],'darkgreen',
        )

plt.ylabel(r'$\mathrm{Pr}\left(z_{1}\right)$')
plt.legend((r'$p_{1}$',
                r'$p_{2}$',
                r'$d$'
                ),bbox_to_anchor=(1.1,1.0))
plt.ylim(-0.05,1.05)

# add subplot depicting patch frequencies environment e1
plt.subplot(num_rows,1,2)

plt.plot(data["generation"],data["f_1_0"],'y',
        data["generation"],data["f_1_1"],'g',
        data["generation"],data["f_1_2"],'magenta',
        linewidth=1)
plt.legend((r'$f_{e_{1},0}$',r'$f_{e_{1},1}$',r'$f_{e_{1},2}$'))

plt.ylabel(r'Patch freqs')
plt.ylim(-0.05,1.05)

# add subplot depicting patch frequencies environment e2
plt.subplot(num_rows,1,3)

plt.plot(data["generation"],data["f_2_0"],'y',
        data["generation"],data["f_2_1"],'g',
        data["generation"],data["f_2_2"],'magenta',
        linewidth=1)
plt.legend((r'$f_{e_{2},0}$',r'$f_{e_{2},1}$',r'$f_{e_{2},2}$'))

plt.ylabel(r'Patch freqs')
plt.ylim(-0.05,1.05)


# add subplot depicting variances
plt.subplot(num_rows,1,4)

plt.plot(data["generation"],data["varp1"],'y',
        data["generation"],data["varp2"],'g',
        data["generation"],data["vard"],'k',
        linewidth=1)
plt.legend((r'$\sigma_{p_{1}}^{2}$',r'$\sigma_{p_{2}}^{2}$',r'$\sigma_{d}^{2}$'))

plt.ylabel(r'Variances')

graphname = os.path.dirname(filename)
if graphname != '':
    graphname += "/"
graphname += "graph_" + os.path.basename(filename) + ".pdf"

plt.savefig(graphname,format="pdf")
