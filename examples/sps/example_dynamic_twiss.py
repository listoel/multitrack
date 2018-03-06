# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
import os
import sys

import matplotlib.pyplot as plt
cm = plt.cm.get_cmap('viridis')

codedir = os.path.abspath('../../..')
sys.path.insert(1, codedir)
import multitrack as mt
import multitrack.dataprocessing as datproc

sys.path.insert(1, '/afs/cern.ch/project/sloex/code/madxBatch')
from python.dataprocessing import getsettings, readtfs

# Settings for multitrack
npart = 500#0
chromatic = True
thin = False
dpp_offsets = list(np.arange(-0.0015,0.00150001,0.00005))

use_fulltrack = True
nturns = 300

disphack = False

dpp = 0.0 if thin else 5E-5

fig, ax = plt.subplots()

# For finding the twiss files
twissstudy = '/afs/cern.ch/project/sloex/test/multitwiss_nodb'
twissdpps = eval(getsettings(twissstudy)['slices'])
def gettwiss(dpp):
    twissnum = min(range(len(twissdpps)), key=lambda i: abs(twissdpps[i]-dpp))
    header, twiss = readtfs(twissstudy+'/output/twiss_'+str(twissnum)+'.tfs')
    return header, twiss

# Run simulations for slices
for i, dpp_offset in enumerate(dpp_offsets):
    # Load optics
    header, twiss = gettwiss(dpp_offset)
    chroma = header['DQ1'] if chromatic else 0.0
    tune = header['Q1']
    twissdpp = header['DELTAP']

    madk2l = -0.0887408

    alpha_ex = twiss.get_value('AP.UP.ZS21633','ALFX')
    beta_ex = twiss.get_value('AP.UP.ZS21633','BETX')
    mu_ex = twiss.get_value('AP.UP.ZS21633','MUX')
    dx_ex = twiss.get_value('AP.UP.ZS21633','DX')/np.sqrt(beta_ex)
    dp_ex = alpha_ex*dx_ex+twiss.get_value('AP.UP.ZS21633','DPX')*np.sqrt(beta_ex)
    xbump_ex = twiss.get_value('AP.UP.ZS21633','X')
    pbump_ex = twiss.get_value('AP.UP.ZS21633','PX')

    alpha_m1 = twiss.get_value('LSE.10602','ALFX')
    beta_m1 = twiss.get_value('LSE.10602','BETX')
    mu_m1 = twiss.get_value('LSE.10602','MUX')
    dx_m1 = twiss.get_value('LSE.10602','DX')/np.sqrt(beta_ex)
    dp_m1 = alpha_m1*dx_m1+twiss.get_value('LSE.10602','DPX')*np.sqrt(beta_m1)
    k2l_m1 = -0.5*madk2l*beta_m1**1.5

    alpha_m2 = twiss.get_value('LSE.22402','ALFX')
    beta_m2 = twiss.get_value('LSE.22402','BETX')
    mu_m2 = twiss.get_value('LSE.22402','MUX')
    dx_m2 = twiss.get_value('LSE.22402','DX')/np.sqrt(beta_ex)
    dp_m2 = alpha_m1*dx_m2+twiss.get_value('LSE.10602','DPX')*np.sqrt(beta_m2)
    k2l_m2 = -0.5*madk2l*beta_m2**1.5

    alpha_m3 = twiss.get_value('LSE.40602','ALFX')
    beta_m3 = twiss.get_value('LSE.40602','BETX')
    mu_m3 = twiss.get_value('LSE.40602','MUX')
    dx_m3 = twiss.get_value('LSE.40602','DX')/np.sqrt(beta_ex)
    dp_m3 = alpha_m3*dx_m3+twiss.get_value('LSE.10602','DPX')*np.sqrt(beta_m3)
    k2l_m3 = -0.5*madk2l*beta_m3**1.5

    alpha_m4 = twiss.get_value('LSE.52402','ALFX')
    beta_m4 = twiss.get_value('LSE.52402','BETX')
    mu_m4 = twiss.get_value('LSE.52402','MUX')
    dx_m4 = twiss.get_value('LSE.52402','DX')/np.sqrt(beta_ex)
    dp_m4 = alpha_m4*dx_m1+twiss.get_value('LSE.52402','DPX')*np.sqrt(beta_m4)
    k2l_m4 = -0.5*madk2l*beta_m4**1.5

    if disphack:
        dx_ex = 0
        dp_ex = 0
        dx_m1 = 0
        dp_m1 = 0
        dx_m2 = 0
        dp_m2 = 0
        dx_m3 = 0
        dp_m3 = 0
        dx_m4 = 0
        dp_m4 = 0


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
                       dpp_offset=dpp_offset-twissdpp, npart=npart, seed=None)

    tracks, extractt = mt.track(ring, init, extraction=extraction,
                                dqstart=0, dqend=0,
                                nturns=nturns, fulltrack=use_fulltrack)

    # Hacky plotting
    nturn = tracks.shape[1]
    fulltrack = (nturn!=1)
    normalization = ring.get_normalization()

    if fulltrack:
        turnind = extractt
    else:
        turnind = [0 for j in extractt]

    npart = len(extractt)
    extr = [j for j in range(npart) if extractt[j]!=-1]
    print i, len(extr)*1.0/npart
    xextr = [tracks[i, turnind[i], 0] for i in extr]
    pextr = [tracks[i, turnind[i], 1] for i in extr]
    dppextr = [init['dpp'].iloc[i]+twissdpp for i in extr]
    s = ax.scatter(xextr,pextr, c=dppextr, cmap=cm, vmin=-0.0025, vmax=0.0020, edgecolor='')
    ax.set_xlim(0.0675, 0.084)
    ax.set_ylim(-0.0018, -0.0013)

plt.colorbar(s)
plt.show()
