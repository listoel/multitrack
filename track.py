#import input
import math
import numpy as np

from plotting import histscatter, plotbeam

class ExtractIt(Exception): pass

# Code probably fails with negative x_wire, check!

# TODO!!!! Add functionality from sextdectrack small!! (fultrack is true, change track var)

# For extracted particles non-normalized coords at zs are given.
# For non-extracted particles final normalized coorda at s=0 are given.

def track(ring, beam, extraction=None, epsilonstart=0.0,
          epsilonend=0.0, nturns=10000, npart=10000, seed=0, plots=None,
          fulltrack=False, gaussinit=True,
          xbin=None, ybin=None, xlim=None, ylim=None, clim=[None,None]):
    np.random.seed(seed)    

    # Normalize based on K2 (virtual sextupole)
    resvirt = 0.0
    imsvirt = 0.0
    for multipole in ring.multipoles:
        if 2 in multipole[1]:
            resvirt += multipole[1][2]*math.cos(3.0*multipole[0])
            imsvirt += multipole[1][2]*math.sin(3.0*multipole[0])
    normalization = math.sqrt(resvirt**2+imsvirt**2)

    if normalization==0.0:
        print "ERROR: No sextupoles found"
        exit()

    multipoles = [[m[0],{i:m[1][i]/normalization**(i-1) for i in m[1]}]
                  for m in ring.multipoles]

    # Cos and sin of phase advance between multipoles
    rotcos = [math.cos(multipoles[0][0])]
    rotcos = rotcos+[math.cos(multipoles[i][0]-multipoles[i-1][0])
                     for i in range(1,len(multipoles))]
    rotsin = [math.sin(multipoles[0][0])]
    rotsin = rotsin+[math.sin(multipoles[i][0]-multipoles[i-1][0])
                     for i in range(1,len(multipoles))]

    # Final rotation to complete the ring at resonance
    phimul = multipoles[-1][0] % (2*math.pi)
    phifin = (2*math.pi*(1+ring.tune) - phimul) % (2*math.pi)

    # Initialize particle momenta and positions
    if beam.dpp==0.0:
        dpps = [0.0 for n in range(npart)]
        color = None
    else:
        dpps = np.random.uniform(-beam.dpp, beam.dpp, npart)
        color = dpps

    if gaussinit:
        init = np.random.normal(0, beam.sigma*abs(normalization), (npart, 2))
    else:
        init1 = np.array([[0, 0.1*i/zoom] for i in range(-15,16)])
        init2 = np.array([np.dot([[c,s],[-s,c]],v) for v in init1 if v[1]!=0])
        init = np.concatenate((init1,init2))
        npart = len(init)

    if fulltrack:
        tracks = np.zeros((npart, nturns+1, 2))
        tracks[:, 0, :] = np.copy(init)
    else:
        tracks = np.copy(init)

    # Initialize extraction related stuff
    iextr = -2
    if extraction is not None:
        iextr = -1
        for multipole in multipoles:
            if extraction.mu > multipole[0]:
                iextr+=1

        if iextr==-1:
            cmuextr = math.cos(extraction.mu)
            smuextr = math.sin(extraction.mu)
        else:
            cmuextr = math.cos(extraction.mu-multipoles[iextr][0])
            smuextr = math.sin(extraction.mu-multipoles[iextr][0])

        alpha = extraction.alpha
        beta = extraction.beta
        denorm = np.array(((math.sqrt(beta), 0),
                           (-alpha/math.sqrt(beta), 1/math.sqrt(beta))))
        rotate = np.array(((cmuextr, smuextr), (-smuextr, cmuextr)))
        zstrans = np.dot(denorm, rotate/normalization)
        dispersion = np.array((extraction.dx,extraction.dp))

        xbump = extraction.xbump
        pbump = extraction.pbump
        xwire = extraction.xwire

        wiretests = [extraction.wiretest(normalization,dppi) for dppi in dpps]
        extractt = [-1 for i in range(npart)]

    # Track particles
    depsilon = (epsilonend-epsilonstart)/nturns
    for part in range(npart):
        epsilon = epsilonstart
        if ring.chroma is not None:
            epsilon += 6*math.pi*ring.chroma*dpps[part]
        x, p = tracks[part, 0], tracks[part, 1]
        try:
            for turn in range(nturns):
                # If extraction point is before multipoles, check for extraction first
                if iextr==-1:
                    if np.sign(normalization)*(x*cmuextr+p*smuextr) > wiretests[part]:
                        extractt[part] = turn+1
                        tracks[part, 0] = x
                        tracks[part, 1] = p
                        raise ExtractIt
                # Then do the multipoles
                for i, multipole in enumerate(multipoles):
                    # Rotate to the multipole
                    xn = rotcos[i]*x + rotsin[i]*p
                    pn = rotcos[i]*p - rotsin[i]*x
                    x, p = xn, pn

                    # Give multipole kick
                    kick = 0
                    for j, strength in multipole[1].iteritems():
                        kick += strength*x**j
                    p = p+kick

                    # If extraction point is between here and the next multipole, check for extraction
                    if iextr==i:
                        if np.sign(normalization)*(x*cmuextr+p*smuextr) > wiretests[part]:
                            extractt[part] = turn+1
                            tracks[part, 0] = x
                            tracks[part, 1] = p
                            raise ExtractIt
                # Then rotate until the end of the ring
                cosfin = math.cos(phifin+epsilon/3)
                sinfin = math.sin(phifin+epsilon/3)
                xn = cosfin*x + sinfin*p
                pn = cosfin*p - sinfin*x
                x, p = xn, pn
                epsilon = epsilon+depsilon
            if fulltrack:
                tracks[part, turn+1, 0] = x
                tracks[part, turn+1, 1] = p
            else:
                tracks[part, 0] = x
                tracks[part, 1] = p
        except ExtractIt:
            if fulltrack:
                tracks[part, turn+1] = (np.dot(zstrans, tracks[part])
                                        +dispersion*dpps[part]+[xbump, pbump])
            else:
                tracks[part] = (np.dot(zstrans, tracks[part])
                                +dispersion*dpps[part]+[xbump, pbump])

    # Start analysis of results
    # ejected=no longer circulating, extracted=also did not hit septum
    remaining = [i for i in range(npart) if extractt[i] == -1]
    ejected = [i for i in range(npart) if extractt[i] != -1]
    extracted = [i for i in ejected
                 if (tracks[i, 0] > extraction.xwire+extraction.wirethickness)]

    rex = [tracks[i, 0] for i in remaining]
    rep = [tracks[i, 1] for i in remaining]
    ejx = [tracks[i, 0] for i in ejected]
    ejp = [tracks[i, 1] for i in ejected]
    exx = [tracks[i, 0] for i in extracted]
    exp = [tracks[i, 1] for i in extracted]

    # Calculate figures of merit
    ejectperc = 100.0*len(ejected)/npart
    try:
        hitperc = 100.0*(len(ejected)-len(extracted))/len(ejected)
        minexx = min(exx)
        avgexx = np.mean(exx)
        maxexx = max(exx)
        minexp = min(exp)
        avgexp = np.mean(exp)
        maxexp = max(exp)
        stdexx = np.std(exx)
        stdexp = np.std(exp)
        statemit = math.sqrt(np.linalg.det(np.cov(exx,exp)))
    except (ZeroDivisionError,ValueError) as e:
        print e
        hitperc = np.nan
        minexx = np.nan
        avgexx = np.nan
        maxexx = np.nan
        minexp = np.nan
        avgexp = np.nan
        maxexp = np.nan
        stdexx = np.nan
        stdexp = np.nan
        statemit = np.nan

    report = {'ejectperc': ejectperc, 'hitperc': hitperc, 'statemit': statemit,
              'minexx': minexx, 'avgexx': avgexx, 'maxexx':maxexx,
              'minexp': minexp, 'avgexp': avgexp, 'maxexp': maxexp,
              'stdexx': stdexx, 'stdexp':stdexp}

    # Plotting
    if plots is not None:
        if color is None:
            ejcolor = None
            recolor = None
        else:
            ejcolor = [color[i] for i in ejected]
            recolor = [color[i] for i in remaining]
        if len(ejected)>0 and "histscatter" in plots:
            histscatter([x*1000 for x in ejx], [p* 1000 for p in ejp],
                        xlabel=r"$x\ [\mathrm{mm}]$",
                        ylabel=r"$x'\ [\mathrm{mrad}]$",
                        xstart=extraction.xwire*1000,
                        xbin=extraction.wirethickness*1000, ybin=0.01,
                        save=plots["histscatter"], color=ejcolor, xlim=xlim, ylim=ylim,clim=clim)
        if "remaining" in plots:
            plotbeam(rex, rep, save=plots["remaining"], color=recolor, clim=clim)
        elif "fulltracks" in plots:
            #TODO replace with propper plot command from simple track?
            print "todo: implement missing plot"

    return {'figsofmerit': report, 'initzs': init, 'finalzs': tracks,
            'extractt': extractt, 'dpps': dpps}

def printreport(report):
    out = ("Percentage of particles extracted: "+str(report['ejectperc'])+"\n"+
           "Percentage of extracted particles hitting wire: "+str(report['hitperc'])+"\n"+
           "Statistical emittance of the new beam: "+str(report['statemit'])+"\n"+
           "Range of extracted x: ["+str(report['minexx'])+", "+str(report['maxexx'])+"]\n"+
           "Range of extracted p: ["+str(report['minexp'])+", "+str(report['maxexp'])+"]\n"+
           "Standard deviation extracted x: "+str(report['stdexx'])+"\n"+
           "Standard deviation extracted p: "+str(report['stdexp'])+"\n\n")
    print out

