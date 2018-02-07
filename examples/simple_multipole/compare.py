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
                           names = ['NUMBER', 'TURN', 'X', 'PX', 'PT'], index_col = [0,1])

    # Read madx tracks
    madxtracks = pd.DataFrame()
    colnames = ['NUMBER', 'TURN', 'X', 'PX', 'Y', 'PY', 'T', 'PT', 'S', 'E']
    nskip = 8
    for i, trackfile in enumerate(listdir(madxtrackfolder)):
        temptrack = pd.read_csv(madxtrackfolder+trackfile, delim_whitespace = True,
                                skipinitialspace = True, skiprows = nskip,
                                names = colnames, index_col = [0,1])
        madxtracks = madxtracks.append(temptrack)
        print 'waitforit '+str(i)

# Make some graphs (also plot diffs?)
# Turn number extracted
f, (ax1, ax2) = plt.subplots(1, 2)

madxloss['TURN'].hist(ax=ax1, bins=100)
pyloss['TURN'].hist(ax=ax2, bins=100)

plt.show()

# Extracted beam
xlim = (0.023,0.041)
ylim = (-120E-6,-80E-6)
f, (ax1, ax2) = plt.subplots(1, 2, sharey=True)

madxloss.plot.scatter('X','PX', ax=ax1, xlim=xlim, ylim=ylim)
pyloss.plot.scatter('X','PX', ax=ax2, xlim=xlim, ylim=ylim)

plt.show()

# Circulating beam
if fulltrack:
    f, (ax1, ax2) = plt.subplots(1, 2)

    madxtracks.plot.scatter('X','PX', ax=ax1)
    pytracks.plot.scatter('X','PX', ax=ax2)

    plt.show()


