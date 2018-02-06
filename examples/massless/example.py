# -*- coding: utf-8 -*-
import numpy as np
import os
import sys

codedir = os.path.abspath('../../..')
sys.path.insert(1, codedir)
import multitrack as mt
import multitrack.dataprocessing as datproc

# Kristof's parameters
xwire = 0.02427
alphazs=2.358
betazs=96.78
muzs = 2*np.pi*0.141#0.106#

mag6 = 0.09#0.03#
beta6 = betazs

strms = -46#-17#
xcircms = 0.0177#0.0245#
xextrms = xcircms+0.0085#0.0055#
betams = betazs

qres = 1.0/3.0
qstart = 0.33
qend = 0.334

beamvar = 0.4
nturns = 5000
npart = 500
# Emittance??


# Slow extraction simulation

extraction = mt.Extraction(alpha=alphazs, beta=betazs, mu=muzs,
                           dx=0, dp=0, xbump=0, pbump=0,
                           xwire=xwire, wirethickness=0.0002)

k2 = 1000*mag6*np.sqrt(beta6)
kick = strms/1E3/betams
septum = mt.simplems(kick, xcircms, xextrms, beta=betams)

ring = mt.Ring([[0,{2:k2, -1:septum}]], tune=qres, chroma=0)

sigma = np.sqrt(beamvar)/1000/np.sqrt(beta6)
init = mt.get_init(ring, btype='gaussian', scale=sigma, dpp=0.0000,
                   npart=npart, seed=0)

tracks, extractt = mt.track(ring, init, extraction=extraction,
                            epsilonstart=6*np.pi*(qstart-qres), epsilonend=6*np.pi*(qend-qres),
                            nturns=nturns, fulltrack=False)


# Some plots

import matplotlib.pyplot as plt
npart = len(extractt)
circ = [i for i in range(npart) if extractt[i]==-1]
extr = [i for i in range(npart) if extractt[i]!=-1]
xcirc = [x for i,x in enumerate(tracks[:,0]) if i in circ]
pcirc = [p for i,p in enumerate(tracks[:,1]) if i in circ]
xextr = [x for i,x in enumerate(tracks[:,0]) if i in extr]
pextr = [p for i,p in enumerate(tracks[:,1]) if i in extr]
print len(extr)
plt.hist(extractt)
plt.show()
plt.scatter(xcirc,pcirc)
plt.show()
plt.scatter(xextr,pextr)
plt.show()

