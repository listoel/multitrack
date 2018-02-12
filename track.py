#import input
import math
import numpy as np

class ExtractIt(Exception):
    """Internal exception raised during tracking when the particle is extracted."""
    pass

# TODO Code probably fails with negative x_wire, check!

# TODO get access to all this stuff
#    return {'figsofmerit': report, 'init': init, 'finalzs': tracks,
#            'extractt': extractt}

# TODO split observation point normalized coord tracks and ZS coords

# For extracted particles non-normalized coords at zs are given.
# For non-extracted particles final normalized coords at s=0 are given.

def track(ring, init, extraction=None, dqstart=0.0, dqend=0.0,
          nturns=10000, fulltrack=False):
    """Generate particles from the beam and track them through the ring.

    Parameters
    ----------
    ring : multitrack.Ring
        Ring object to define the machine in which to track.
    init : dataframe
        Dataframe with initial dpp, X and P of particles to be tracked.
    extraction : multitrack.Extraction, optional
        Extraction object to define the characteristics of the
        extraction point, or None to avoid extracting.
    dqstart : float, optional
        :math:`delta Q` at the start of the tune sweep.
    dqend : float, optional
        :math:`delta Q` at the end of the tune sweep.
    nturns : int, optional
        Maximum number of turns for which to track particles. Indirectly
        determines the tune sweep speed.
    fulltrack : bool, optional
        True to store full tracks of all particles, or False to store
        only the most recent coordinate.

    Returns
    -------
    :obj:`array`
        npart x (nturns+1) x 2 array with full tracks of each particle
        if fulltrack=True, else npart x 2 array with latest coordinates.
        TODO coords at injection/extraction?
    :obj:`array`
        turnnumber on which extraction occured for each particle;
        defaults to -1 if it was not extracted
    """  
    npart = init.shape[0]

    phifin = (2*math.pi*ring.tune - ring.elements[-1][0])
    if phifin<0.0:
        print "ERROR: Final multipole beyond end of ring."
        exit()

    # Normalize based on K2 (virtual sextupole)
    normalization = ring.get_normalization()
    if normalization==0.0:
        print "ERROR: No sextupoles found."
        exit()

    def normalize(i,k):
        if i<0: # Custom function
            return k
        elif i==0: # Dipole kick
            return k*normalization
        else: # Multipole kick
            return k/normalization**(i-1)

    elements = [[m[0],{i:normalize(i, m[1][i]) for i in m[1]}]
                 for m in ring.elements]

    # Initialize particle momenta and positions
    if fulltrack:
        tracks = np.zeros((npart, nturns+1, 2))
    else:
        tracks = np.empty((npart, 1, 2))
    tracks[:, 0, 0] = np.copy(init['X'].values)
    tracks[:, 0, 1] = np.copy(init['P'].values)

    # Initialize extraction related stuff
    iextr = -2
    if extraction is not None:
        iextr = -1
        for element in elements:
            if extraction.mu > element[0]:
                iextr+=1

        alpha = extraction.alpha
        beta = extraction.beta
        denorm = np.array(((math.sqrt(beta), 0),
                           (-alpha/math.sqrt(beta), 1/math.sqrt(beta))))
        dispersion = np.array((extraction.dx,extraction.dp))

        xbump = extraction.xbump
        pbump = extraction.pbump
        xwire = extraction.xwire

        wiretests = [extraction.wiretest(normalization,dppi) for dppi in init['dpp']]
        extractt = [-1 for i in range(npart)]

    # Track particles
    tpdq_step = 2*math.pi*(dqend-dqstart)/nturns
    for part in range(npart):

        mydpp = init['dpp'].loc[part]
        myinvdpp = 1/(1+mydpp)

        tpdq = 2*math.pi*dqstart
        if ring.chroma is not None:
            tpdq += 2*math.pi*ring.chroma*mydpp

        turnind = 0
        x, p = tracks[part, turnind, 0], tracks[part, turnind, 1]

        phasemod = 1.0
        if ring.chroma is not None:
            phasemod = 1.0 + ring.chroma/ring.tune*mydpp

        # Cos and sin of phase advance between elements
        rotcos = [math.cos(elements[0][0]*phasemod)]
        rotcos = rotcos+[math.cos((elements[i][0]-elements[i-1][0])*phasemod)
                         for i in range(1,len(elements))]
        rotsin = [math.sin(elements[0][0]*phasemod)]
        rotsin = rotsin+[math.sin((elements[i][0]-elements[i-1][0])*phasemod)
                         for i in range(1,len(elements))]

        # And for the extraction
        if extraction is not None:
            if iextr==-1:
                cmuextr = math.cos(extraction.mu*phasemod)
                smuextr = math.sin(extraction.mu*phasemod)
            else:
                cmuextr = math.cos((extraction.mu-elements[iextr][0])*phasemod)
                smuextr = math.sin((extraction.mu-elements[iextr][0])*phasemod)

        # And for final rotation to complete the ring at resonance
        cosfin = math.cos(phifin*phasemod)
        sinfin = math.sin(phifin*phasemod)

        try:
            for turn in np.arange(nturns)+1:
                if fulltrack:
                    turnind = turn
                # If extraction point is before elements, check for extraction first
                if iextr==-1:
                    if (x*cmuextr+p*smuextr) > wiretests[part]:
                        extractt[part] = turn
                        raise ExtractIt
                # Then do the elements
                for i, element in enumerate(elements):
                    # Rotate to the element
                    xn = rotcos[i]*x + rotsin[i]*p
                    pn = rotcos[i]*p - rotsin[i]*x
                    x, p = xn, pn

                    # Give multipole/custom kick
                    kick = 0
                    for j, strength in element[1].iteritems():
                        if j<0: # custom kick
                            kick += strength(x, normalization)
                        else: #dipole/multipole kick
                            kick += strength*x**j
                    p = p+kick*myinvdpp

                    # If extraction point is between here and the next element, check for extraction
                    if iextr==i:
                        if (x*cmuextr+p*smuextr) > wiretests[part]:
                            extractt[part] = turn
                            raise ExtractIt
                # Then rotate until the end of the ring
                xn = cosfin*x + sinfin*p
                pn = cosfin*p - sinfin*x
                x, p = xn, pn
                tpdq += tpdq_step
                tracks[part, turnind, 0] = x
                tracks[part, turnind, 1] = p
        except ExtractIt:
            rotate = np.array(((cmuextr, smuextr), (-smuextr, cmuextr)))
            zstrans = np.dot(denorm, rotate/normalization)
            tracks[part, turnind] = (np.dot(zstrans, [x, p])
                                     +dispersion*mydpp
                                     +[xbump, pbump])

    return tracks, extractt

