# -*- coding: utf-8 -*-
import numpy as np
import os
import sys

codedir = os.path.abspath('../../..')
sys.path.insert(1, codedir)
import multitrack as mt
import multitrack.dataprocessing as datproc

extraction = mt.Extraction(alpha=0, beta=100, mu=np.pi/3,
                           dx=0, dp=0, xbump=0, pbump=0,
                           xwire=0.024, wirethickness=0.0002)

ring = mt.Ring([[0,{2:210}]], tune=2.0/3, chroma=0)

init = mt.get_init(ring, btype='gaussian', scale=2E-4, dpp=0.0000,
                   npart=9999, seed=0)
                   #npart=100, seed=0)

datproc.init_to_madx(ring, init, './out/mt/init.madx', alpha=0, beta=100)

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

report = datproc.get_report(tracks, extractt, extraction)
datproc.print_report(report)

with open('./out/mt/losses.dat', 'w') as f:
    for i in extr:
	t = extractt[i]
	x = tracks[i,0]
	p = tracks[i,1]
	dpp = init['dpp'].iloc[i]
	f.write(str(i+1)+'\t'+str(t)+'\t'+str(x)+'\t'+str(p)+'\t'+str(dpp/(1+dpp))+'\n')
