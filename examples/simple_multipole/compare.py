import csv
import io
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np
from os import listdir
import pandas as pd

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

# Make some graphs
cmap = cm.get_cmap('viridis')
vmin = -0.0016
vmax = 0.0016

# Turn number extracted
f, (ax1, ax2) = plt.subplots(1, 2)

madxloss['TURN'].hist(ax=ax1, bins=100)
pyloss['TURN'].hist(ax=ax2, bins=100)

plt.show()

# Extracted beam
xlim = (0.023,0.041)
ylim = (-175E-6,-35E-6)
f, (ax1, ax2) = plt.subplots(1, 2, sharex=True, sharey=True)

if not madxloss.empty:
    madxloss.plot.scatter('X', 'PX', c='PT', cmap=cmap, ax=ax1, xlim=xlim, ylim=ylim, vmin=vmin, vmax=vmax)
if not pyloss.empty:
    pyloss.plot.scatter('X', 'PX', c='PT', cmap=cmap, ax=ax2, xlim=xlim, ylim=ylim, vmin=vmin, vmax=vmax)

plt.show()

# Circulating beam
if fulltrack:
    xlim = (-0.03,0.035)
    ylim = (-3.8E-4,1.7E-4)
    f, (ax1, ax2) = plt.subplots(1, 2, sharex=True, sharey=True)

    madxtracks.plot.scatter('X', 'PX', c='PT', cmap=cmap, ax=ax1, xlim=xlim, ylim=ylim, vmin=vmin, vmax=vmax)
    pytracks.plot.scatter('X', 'PX', c='PT', cmap=cmap, ax=ax2, xlim=xlim, ylim=ylim, vmin=vmin, vmax=vmax)

    plt.show()

# Initial conditions
if fulltrack:
    xlim = (-0.01,0.01)
    ylim = (-1.5E-4,1.5E-4)
    f, (ax1, ax2) = plt.subplots(1, 2, sharex=True, sharey=True)

    madxtracks[madxtracks['TURN']==0].plot.scatter('X', 'PX', c='PT', cmap=cmap, ax=ax1, xlim=xlim, ylim=ylim, vmin=vmin, vmax=vmax)
    pytracks[pytracks['TURN']==0].plot.scatter('X', 'PX', c='PT', cmap=cmap, ax=ax2, xlim=xlim, ylim=ylim, vmin=vmin, vmax=vmax)

    plt.show()

# Brag about stats
if not (madxloss.empty or pyloss.empty):
    xlim = (-0.017,0.017)
    ylim = (-3E-5,3E-5)
    x = madxloss['X'] - pyloss['X']
    y = madxloss['PX'] - pyloss['PX']
    pd.concat([x, y, pyloss['PT']], axis=1).plot.scatter('X', 'PX', c='PT', cmap=cmap, xlim=xlim, ylim=ylim, vmin=vmin, vmax=vmax)
    plt.show()

if fulltrack:
    xlim = (-0.03,0.03)
    ylim = (-4E-4,4E-4)
    madxtracks.set_index(['NUMBER','TURN'], inplace=True)
    pytracks.set_index(['NUMBER','TURN'], inplace=True)
    x = madxtracks['X'] - pytracks['X']
    y = madxtracks['PX'] - pytracks['PX']
    pd.concat([x, y, pytracks['PT']], axis=1).plot.scatter('X', 'PX', c='PT', cmap=cmap, xlim=xlim, ylim=ylim, vmin=vmin, vmax=vmax)
    plt.show()
    



