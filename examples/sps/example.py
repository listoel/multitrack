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

use_fulltrack = True
nturns = 1000


#SPS optics
dpp = 0.0 if thin else 0.0015
chroma = -31.87652221 if chromatic else 0.0
tune = 80.0/3

madk2l = -0.0887408

alpha_init = -2.179766259
beta_init = 98.19309432
dx_init = 1.585632486/np.sqrt(beta_init)
dp_init = alpha_init*dx_init-0.03840948179*np.sqrt(beta_init)

alpha_ex = 2.364235034
beta_ex = 97.13807829
mu_ex = 6.429077311
dx_ex = -0.1640607634/np.sqrt(beta_ex)
dp_ex = alpha_ex*dx_ex-0.001175087036*np.sqrt(beta_ex)
xbump_ex = 0.04366340344
pbump_ex = -0.0009612258572

alpha_m1 = -2.424574211
beta_m1 = 100.9844863
mu_m1 = 0.7476548715
dx_m1 = 2.578448815/np.sqrt(beta_m1)
dp_m1 = alpha_m1*dx_m1+0.04662179247*np.sqrt(beta_m1)
k2l_m1 = -0.5*madk2l*beta_m1**1.5

alpha_m2 = -2.397587252
beta_m2 = 99.37554046
mu_m2 = 7.405406312
dx_m2 = 4.123013489/np.sqrt(beta_m2)
dp_m2 = alpha_m2*dx_m2+0.09465540758*np.sqrt(beta_m2)
k2l_m2 = -0.5*madk2l*beta_m2**1.5

alpha_m3 = -2.27071997
beta_m3 = 99.7330477
mu_m3 = 14.06643485
dx_m3 = 2.864494031/np.sqrt(beta_m3)
dp_m3 = alpha_m3*dx_m3+0.04923058376*np.sqrt(beta_m3)
k2l_m3 = -0.5*madk2l*beta_m3**1.5

alpha_m4 = -2.107756887
beta_m4 = 89.12601545
mu_m4 = 20.74061979
dx_m4 = 4.505432689/np.sqrt(beta_m4)
dp_m4 = alpha_m4*dx_m4+0.1054524659*np.sqrt(beta_m4)
k2l_m4 = -0.5*madk2l*beta_m4**1.5


#Simulation

extraction = mt.Extraction(alpha=alpha_ex, beta=beta_ex, mu=2*np.pi*mu_ex,
                           dx=dx_ex, dp=dp_ex, xbump=xbump_ex, pbump=pbump_ex,
                           xwire=0.06795, wirethickness=0.0002)

ring = mt.Ring([[2*np.pi*mu_m1,{2:k2l_m1},[dx_m1,dp_m1]],
                [2*np.pi*mu_m2,{2:k2l_m2},[dx_m2,dp_m2]],
                [2*np.pi*mu_m3,{2:k2l_m3},[dx_m3,dp_m3]],
                [2*np.pi*mu_m4,{2:k2l_m4},[dx_m4,dp_m4]]],
               tune=tune, chroma=chroma)

init = mt.get_init(ring, btype='gaussian', scale=np.sqrt(12E-6/426.3156), dpp=dpp,
                   dpp_offset=dpp_offset, npart=npart, seed=0)

datproc.init_to_madx(ring, init, './out/mt/init.madx', alpha=alpha_init, beta=beta_init, dx=dx_m1, dp=dp_m1)

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
cm = plt.cm.get_cmap('viridis')
npart = len(extractt)
circ = [i for i in range(npart) if extractt[i]==-1]
extr = [i for i in range(npart) if extractt[i]!=-1]
xcirc = [tracks[i, turnind[i], 0] for i in circ]
pcirc = [tracks[i, turnind[i], 1] for i in circ]
xextr = [tracks[i, turnind[i], 0] for i in extr]
pextr = [tracks[i, turnind[i], 1] for i in extr]
dppextr = [init['dpp'].iloc[i] for i in extr]
print len(extr)
plt.hist(extractt)
plt.show()
plt.scatter(xcirc,pcirc)
plt.show()
s = plt.scatter(xextr,pextr, c=dppextr, cmap=cm, vmin=-0.0025, vmax=0.0020)
plt.xlim(0.0675, 0.084)
plt.ylim(-0.0018, -0.0013)
plt.colorbar(s)
plt.show()

report = datproc.get_report(tracks, extractt, extraction)
datproc.print_report(report)
