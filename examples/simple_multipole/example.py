# -*- coding: utf-8 -*-
import numpy as np
import os
import sys

codedir = os.path.abspath('../../..')
sys.path.insert(1, codedir)
import multitrack as mt
import multitrack.dataprocessing as datproc

#Settings

npart = 9999
dpp = 0.0015
chroma = -0.7351072162

use_fulltrack = True
nturns = 1000

pyk2l = 210

alpha_init = 0
beta_init = 100


#Simulation

extraction = mt.Extraction(alpha=0, beta=100, mu=np.pi/3,
                           dx=0, dp=0, xbump=0, pbump=0,
                           xwire=0.024, wirethickness=0.0002)

ring = mt.Ring([[0,{2:pyk2l}]], tune=2.0/3, chroma=chroma)

init = mt.get_init(ring, btype='gaussian', scale=2E-4, dpp=dpp,
                   npart=npart, seed=0)

datproc.init_to_madx(ring, init, './out/mt/init.madx', alpha=0, beta=100)

tracks, extractt = mt.track(ring, init, extraction=extraction,
                            dqstart=0.0, dqend=0.0,
                            nturns=nturns, fulltrack=use_fulltrack)

# Hacky plotting

nturn = tracks.shape[1]
fulltrack = (nturn!=1)

if fulltrack:
    turnind = extractt
else:
    turnind = [0 for i in extractt]

import matplotlib.pyplot as plt
npart = len(extractt)
circ = [i for i in range(npart) if extractt[i]==-1]
extr = [i for i in range(npart) if extractt[i]!=-1]
xcirc = [tracks[i, turnind[i], 0] for i in circ]
pcirc = [tracks[i, turnind[i], 1] for i in circ]
xextr = [tracks[i, turnind[i], 0] for i in extr]
pextr = [tracks[i, turnind[i], 1] for i in extr]
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
	x = tracks[i,turnind[i],0]
	p = tracks[i,turnind[i],1]
	dpp = init['dpp'].iloc[i]
	f.write(str(i+1)+'\t'+str(t)+'\t'+str(x)+'\t'+str(p)+'\t'+str(dpp/(1+dpp))+'\n')

if fulltrack:
    normalization = ring.get_normalization()
    with open('./out/mt/tracks.dat', 'w') as f:
        for i in range(npart):
            maxt = nturn if extractt[i]==-1 else extractt[i]
            for t in range(maxt):
                x = tracks[i,t,0]
                p = tracks[i,t,1]
                p = (p - alpha_init*x)/normalization/np.sqrt(beta_init)
                x = x/normalization*np.sqrt(beta_init)
                dpp = init['dpp'].iloc[i]
                f.write(str(i+1)+'\t'+str(t)+'\t'+str(x)+'\t'+str(p)+'\t'+str(dpp/(1+dpp))+'\n')
