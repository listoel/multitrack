import math
import os
import sys

codedir = os.path.abspath('../../..')
sys.path.insert(1, codedir)
import multitrack as mt

extraction = mt.Extraction(alpha=0, beta=100, mu=11*math.pi/60,
                           dx=0, dp=0, xbump=0, pbump=0,
                           xwire=0.02, wirethickness=0.0002)

ring = mt.Ring([[0,{2:400}]], tune=2.0/3, chroma=0)

init = mt.get_init(ring, btype="gaussian", scale=100E-6, dpp=0.0015,
                   npart=10000, seed=0)

tracks, extractt = mt.track(ring, init, extraction=extraction,
                            epsilonstart=0.0, epsilonend=0.0,
                            nturns=1000, fulltrack=False)

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

report = mt.datproc.get_report(tracks, extractt, extraction)
mt.datproc.print_report(report)
