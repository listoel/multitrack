#import input
import math
import numpy as np

class ExtractIt(Exception):
    """Internal exception raised during tracking when the particle is extracted."""
    pass

# TODO Code probably fails with negative x_wire, check!

# TODO!!!! Add functionality from sextdectrack small!! (fulltrack is true, change track var)
# (What did this mean?)

# TODO get access to all this stuff
#    return {'figsofmerit': report, 'init': init, 'finalzs': tracks,
#            'extractt': extractt}

# TODO split observation point normalized coord tracks and ZS coords

# For extracted particles non-normalized coords at zs are given.
# For non-extracted particles final normalized coords at s=0 are given.

def track(ring, init, extraction=None, epsilonstart=0.0, epsilonend=0.0,
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
    epsilonstart : float, optional
        :math:`6\\pi\\cdot\\delta Q` at the start of the tune sweep.
    epsilonend : float, optional
        :math:`6\\pi\\cdot\\delta Q` at the end of the tune sweep.
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

    # Cos and sin of phase advance between elements
    rotcos = [math.cos(elements[0][0])]
    rotcos = rotcos+[math.cos(elements[i][0]-elements[i-1][0])
                     for i in range(1,len(elements))]
    rotsin = [math.sin(elements[0][0])]
    rotsin = rotsin+[math.sin(elements[i][0]-elements[i-1][0])
                     for i in range(1,len(elements))]

    # Final rotation to complete the ring at resonance
    phimul = elements[-1][0] % (2*math.pi)
    phifin = (2*math.pi*(1+ring.tune) - phimul) % (2*math.pi)

    # Initialize particle momenta and positions
    if fulltrack:
        tracks = np.zeros((npart, nturns+1, 2))
        tracks[:, 0, 0] = np.copy(init['X'].values)
        tracks[:, 0, 1] = np.copy(init['P'].values)
    else:
        tracks = np.empty((npart, 2))
        tracks[:, 0] = np.copy(init['X'].values)
        tracks[:, 1] = np.copy(init['P'].values)

    # Initialize extraction related stuff
    iextr = -2
    if extraction is not None:
        iextr = -1
        for element in elements:
            if extraction.mu > element[0]:
                iextr+=1

        if iextr==-1:
            cmuextr = math.cos(extraction.mu)
            smuextr = math.sin(extraction.mu)
        else:
            cmuextr = math.cos(extraction.mu-elements[iextr][0])
            smuextr = math.sin(extraction.mu-elements[iextr][0])

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

        wiretests = [extraction.wiretest(normalization,dppi) for dppi in init['dpp']]
        extractt = [-1 for i in range(npart)]

    # Track particles
    depsilon = (epsilonend-epsilonstart)/nturns
    for part in range(npart):
        epsilon = epsilonstart
        if ring.chroma is not None:
            epsilon += 6*math.pi*ring.chroma*init['dpp'].loc[part]
        if fulltrack:
            x, p = tracks[part, 0, 0], tracks[part, 0, 1]
        else:
            x, p = tracks[part, 0], tracks[part, 1]
        try:
            for turn in np.arange(nturns)+1:
                # If extraction point is before elements, check for extraction first
                if iextr==-1:
                    if (x*cmuextr+p*smuextr) > wiretests[part]:
                        extractt[part] = turn
                        if fulltrack:
                            tracks[part, turn, 0] = x
                            tracks[part, turn, 1] = p
                        else:
                            tracks[part, 0] = x
                            tracks[part, 1] = p
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
                    p = p+kick

                    # If extraction point is between here and the next element, check for extraction
                    if iextr==i:
                        if (x*cmuextr+p*smuextr) > wiretests[part]:
                            extractt[part] = turn
                            #save coords before extraction
                            if fulltrack:
                                tracks[part, turn, 0] = x
                                tracks[part, turn, 1] = p
                            else:
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
                    tracks[part, turn, 0] = x
                    tracks[part, turn, 1] = p
                else:
                    tracks[part, 0] = x
                    tracks[part, 1] = p
        except ExtractIt:
            if fulltrack:
                tracks[part, turn] = (np.dot(zstrans, tracks[part, turn])
                                        +dispersion*init['dpp'].loc[part]+[xbump, pbump])
            else:
                tracks[part] = (np.dot(zstrans, tracks[part])
                                +dispersion*init['dpp'].loc[part]+[xbump, pbump])

    return tracks, extractt

