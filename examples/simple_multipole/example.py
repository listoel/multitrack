# -*- coding: utf-8 -*-
import numpy as np
import os
import sys

codedir = os.path.abspath('../../..')
sys.path.insert(1, codedir)
import multitrack as mt
import multitrack.dataprocessing as datproc

#Settings

npart = 5000
chromatic = True
thin = False
dpp_offset = 0.0

use_fulltrack = True
dispersive = True
nturns = 1000

pyk2l = 350
pyk3l = -0.2*pyk2l**2
pyk4l = 0

alpha_init = 0
beta_init = 100
alpha_ex = 0
beta_ex = 100

if thin:
    dpp = 0.0 
else:
    dpp = 0.0015

if chromatic:
    if dispersive:
        chroma = -0.7351188573
    else:
        chroma = -0.7351072162
else:
    chroma = 0.0

if dispersive:
    dx_m = 1.417037058/np.sqrt(beta_init)
    dp_m = 0
    dx_ex = 1.417037058
    dp_ex = 0
else:
    dx_m = 0.0
    dp_m = 0.0
    dx_ex = 0.0
    dp_ex = 0.0


#Simulation

extraction = mt.Extraction(alpha=0, beta=100, mu=np.pi/3,
                           dx=dx_ex, dp=dp_ex, xbump=0, pbump=0,
                           xwire=0.024, wirethickness=0.0002)

ring = mt.Ring([[0,{2:pyk2l, 3:pyk3l, 4:pyk4l},[dx_m,dp_m]]], tune=2.0/3, chroma=chroma)

init = mt.get_init(ring, btype='gaussian', scale=2E-4, dpp=dpp,
                   dpp_offset=dpp_offset, npart=npart, seed=0)

datproc.init_to_madx(ring, init, './out/mt/init.madx', alpha=0, beta=100, dx=dx_m, dp=dp_m)

tracks, extractt = mt.track(ring, init, extraction=extraction,
                            dqstart=0.0, dqend=0.0,
                            nturns=nturns, fulltrack=use_fulltrack)

# Hacky plotting

nturn = tracks.shape[1]
fulltrack = (nturn!=1)
normalization = ring.get_normalization()

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
plt.xlim(0.023,0.041)
plt.ylim(-175E-6,-35E-6)
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
    with open('./out/mt/tracks.dat', 'w') as f:
        for i in range(npart):
            maxt = nturn if extractt[i]==-1 else extractt[i]
            for t in range(maxt):
                x = tracks[i,t,0]/normalization+dx_m*init['dpp'].iloc[i]
                p = tracks[i,t,1]/normalization+dp_m*init['dpp'].iloc[i]
                p = (p - alpha_init*x)/np.sqrt(beta_init)
                x = x*np.sqrt(beta_init)
                dpp = init['dpp'].iloc[i]
                f.write(str(i+1)+'\t'+str(t)+'\t'+str(x)+'\t'+str(p)+'\t'+str(dpp/(1+dpp))+'\n')
