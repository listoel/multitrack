import csv
import io
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from os import listdir

pylossfile = './out/mt/losses.dat'
madxlossfile = './out/madx/losses.tfs'

fulltrack = True

pytrackfile = './out/mt/tracks.dat'
madxtrackfolder = './out/madx/tracks/'

# Read python losses
pyloss = pd.read_csv(pylossfile, delim_whitespace = True,
                     skipinitialspace = True,
                     names = ['NUMBER', 'TURN', 'X', 'PX', 'PT'], index_col = 0)

# Read madx losses
nskip=1
with open(madxlossfile, 'r') as datafile:
    for line in datafile:
        nskip += 1
        if line.startswith('*'):
            colnames = line.strip().split()[1:]
            break

madxloss = pd.read_csv(madxlossfile, delim_whitespace = True,
                       skipinitialspace = True, skiprows = nskip,
                       names = colnames, index_col = 0)

if fulltrack:
    # Read python tracks
    pytracks = pd.read_csv(pytrackfile, delim_whitespace = True,
                           skipinitialspace = True,
                           names = ['NUMBER', 'TURN', 'X', 'PX', 'PT'])

    # Read madx tracks
    colnames = ['NUMBER', 'TURN', 'X', 'PX', 'Y', 'PY', 'T', 'PT', 'S', 'E']
    nskip = 8
    with io.BytesIO() as output:
        csv_writer = csv.writer(output)
        for trackfile in listdir(madxtrackfolder):
            with open(madxtrackfolder+trackfile) as f:
                i = 1
                for line in f:
                    if i>nskip:
                        csv_writer.writerow(line.split())
                    i = i+1
        output.seek(0)
        madxtracks = pd.read_csv(output, names = colnames)

# Make some graphs (also plot diffs?)
# Turn number extracted
f, (ax1, ax2) = plt.subplots(1, 2)

madxloss['TURN'].hist(ax=ax1, bins=100)
pyloss['TURN'].hist(ax=ax2, bins=100)

plt.show()

# Extracted beam
xlim = (0.023,0.041)
ylim = (-120E-6,-80E-6)
f, (ax1, ax2) = plt.subplots(1, 2, sharex=True, sharey=True)

if not madxloss.empty:
    madxloss.plot.scatter('X','PX', ax=ax1, xlim=xlim, ylim=ylim)
if not pyloss.empty:
    pyloss.plot.scatter('X','PX', ax=ax2, xlim=xlim, ylim=ylim)

plt.show()

# Circulating beam
if fulltrack:
    f, (ax1, ax2) = plt.subplots(1, 2, sharex=True, sharey=True)

    madxtracks.plot.scatter('X','PX', ax=ax1)
    pytracks.plot.scatter('X','PX', ax=ax2)

    plt.show()

# Initial conditions
if fulltrack:
    f, (ax1, ax2) = plt.subplots(1, 2, sharex=True, sharey=True)

    madxtracks[madxtracks['TURN']==0].plot.scatter('X','PX', ax=ax1)
    pytracks[pytracks['TURN']==0].plot.scatter('X','PX', ax=ax2)

    plt.show()

# Brag about stats
if not (madxloss.empty or pyloss.empty):
    x = madxloss['X'] - pyloss['X']
    y = madxloss['PX'] - pyloss['PX']
    pd.concat([x, y], axis=1).plot.scatter('X','PX')
    plt.show()

if fulltrack:
    x = (madxtracks.set_index(['NUMBER','TURN']))['X'] - (pytracks.set_index(['NUMBER','TURN']))['X']
    y = (madxtracks.set_index(['NUMBER','TURN']))['PX'] - (pytracks.set_index(['NUMBER','TURN']))['PX']
    pd.concat([x, y], axis=1).plot.scatter('X','PX')
    plt.show()
    



