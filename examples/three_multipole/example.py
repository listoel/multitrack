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
chromaticelements = True
thin = False
dpp_offset = 0.0
disphack = False

use_fulltrack = True
dispersive = True
nturns = 1000

pyk2l = 510
pyk3l = 0
pyk4l = -0.24*pyk2l**3

alpha_init = 0 #dispersion assumes alpha_init==0
beta_init = 100

if thin:
    dpp = 0.0 
else:
    dpp = 0.0015

if chromatic:
    if dispersive:
        chroma = -1.837512547
    else:
        chroma = -1.837768041
else:
    chroma = 0.0

if dispersive:
    alpha_ex = 1.352459616
    beta_ex = 70.73881986
    mu_ex = 1.018000344
    dx_m1 = 1.558571225/np.sqrt(beta_init)
    dp_m1 = 0.004089154318*np.sqrt(beta_init)
    dx_m2 = 1.70022488/np.sqrt(beta_init)
    dp_m2 = 0.003271265846*np.sqrt(beta_init)
    dx_m3 = 0.9920265427/np.sqrt(beta_init)
    dp_m3 = 0.0008177727138*np.sqrt(beta_init)
    dx_ex = 1.417037058
    dp_ex = -0.02290287252
else:
    alpha_ex = 1.352410808
    beta_ex = 70.72537481
    mu_ex = 1.017999925
    dx_m1 = 0.0
    dp_m1 = 0.0
    dx_m2 = 0.0
    dp_m2 = 0.0
    dx_m3 = 0.0
    dp_m3 = 0.0
    dx_ex = 0.0
    dp_ex = 0.0

if disphack:
    dx_m1 = 0.0
    dp_m1 = 0.0
    dx_m2 = 0.0
    dp_m2 = 0.0
    dx_m3 = 0.0
    dp_m3 = 0.0


#Simulation

extraction = mt.Extraction(alpha=alpha_ex, beta=beta_ex, mu=2*np.pi*mu_ex,
                           dx=dx_ex, dp=dp_ex, xbump=0, pbump=0,
                           xwire=0.024269015, wirethickness=0.0002)

ring = mt.Ring([[0,{2:pyk2l/3, 3:pyk3l/3, 4:pyk4l/3},[dx_m1,dp_m1]],
                [2*np.pi/3,{2:pyk2l/3, 3:pyk3l/3, 4:pyk4l/3},[dx_m2,dp_m2]],
                [4*np.pi/3,{2:pyk2l/3, 3:pyk3l/3, 4:pyk4l/3},[dx_m3,dp_m3]]],
               tune=5.0/3, chroma=chroma)

init = mt.get_init(ring, btype='gaussian', scale=2E-4, dpp=dpp,
                   dpp_offset=dpp_offset, npart=npart, seed=0)

datproc.init_to_madx(ring, init, './out/mt/init.madx', alpha=0, beta=100, dx=dx_m1, dp=dp_m1)

tracks, extractt = mt.track(ring, init, extraction=extraction,
                            dqstart=0.0, dqend=0.0,
                            nturns=nturns, fulltrack=use_fulltrack,
                            chromaticelements=chromaticelements)

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
plt.xlim(0.023,0.045)
#plt.ylim(-50E-6,100E-6)
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
                x = tracks[i,t,0]/normalization+dx_m1*init['dpp'].iloc[i]
                p = tracks[i,t,1]/normalization+dp_m1*init['dpp'].iloc[i]
                p = (p - alpha_init*x)/np.sqrt(beta_init)
                x = x*np.sqrt(beta_init)
                dpp = init['dpp'].iloc[i]
                f.write(str(i+1)+'\t'+str(t)+'\t'+str(x)+'\t'+str(p)+'\t'+str(dpp/(1+dpp))+'\n')
