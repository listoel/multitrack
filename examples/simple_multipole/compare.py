import pandas as pd
import matplotlib.pyplot as plt

madxfile = './out/madx/losses.tfs'
pyfile = './out/mt/losses.dat'

# Read madx losses
nskip=1
with open(madxfile, 'r') as datafile:
    for line in datafile:
        nskip += 1
        if line.startswith('*'):
            colnames = line.strip().split()[1:]
            break

madxloss = pd.read_csv(madxfile, delim_whitespace = True,
                       skipinitialspace = True, skiprows = nskip,
                       names = colnames, index_col = 0)

# Read python losses
pyloss = pd.read_csv(pyfile, delim_whitespace = True,
                     skipinitialspace = True, skiprows = nskip,
                     names = ['NUMBER', 'TURN', 'X', 'PX', 'PT'], index_col = 0)

# Make some graphs
xlim = (0.023,0.041)
ylim = (-120E-6,-80E-6)

madxloss.plot.scatter('X','PX')
plt.xlim(xlim)
plt.ylim(ylim)
plt.show()

pyloss.plot.scatter('X','PX')
plt.xlim(xlim)
plt.ylim(ylim)
plt.show()


