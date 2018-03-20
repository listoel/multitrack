from __future__ import division
import math
import numpy as np
import pandas as pd

class Ring:
    """Ring properties for tracking.

    Parameters
    ----------
    elements: :obj:`list` of :obj:`list` of [int, :obj:`dict`]
        A list of elements, each defined by a phase advance in radians
        since the start of the ring and a dictionary of multipole
        components, optionally folowed by the normalized dispersion
        vector and the normalized closed orbit position (if missing, these
        are presumed to be [0,0] and 0) e.g. [[0,{2:K_2a,4:K_4a}],
        [1.12,{2:K_2b,3:K_3b},[0.001,-7E-4]],...]
    tune: float
        Fractional part of the machine tune. (In radians over
        :math:`2\\pi`.)
    chroma: float
        Chromaticity of the machine. (TODO Units?)

    Attributes
    ----------
    elements: :obj:`list` of :obj:`list` like [int, :obj:`dict`]
        A list of elements, each defined by a phase advance in radians
        since the start of the ring and a dictionary of multipole
        components, e.g. [[0,{2:K_2a,4:K_4a}],[1.12,{2:K_2b,3:K_3b}],...]
    tune: float
        Fractional part of the machine tune. (In radians over
        :math:`2\\pi`.)
    chroma: float
        Chromaticity of the machine. (TODO Units?)
    """
    def __init__(self, elements, tune=2.0/3, chroma=None):
        self.elements = sorted(elements, key=lambda m: m[0])
        self.tune = float(tune)
        self.chroma = float(chroma)


    def add_element(self, element):
        """Adds an element to the ring.

        Parameters
        ----------
        element: :obj:`list` of [int, :obj:`dict`]
            An element, defined by a phase advance in radians since the
            start of the ring and a dictionary of multipole components,
            e.g. [1.57,{2:K_2c,3:K_3c}].
            Custom (non-multipole) elements can also be added, in this
            case the dictonary entry should specify a negative integer
            and a function that takes xhat and the normalizatin as input
            and the kick in p-hat as output.
        """
        self.elements = sorted(self.elements+[element],
                                 key=lambda m: m[0])

    def get_normalization(self):
        """Calculates normalized sextupole strength to normalise phase space.

        Returns
        -------
        float
            The normalized sextupole strength for normalisation of
            phase space as described in TODO.
        """
        resvirt = 0.0
        imsvirt = 0.0
        for element in self.elements:
            if 2 in element[1]:
                resvirt += element[1][2]*math.cos(3.0*element[0])
                imsvirt += element[1][2]*math.sin(3.0*element[0])
        return math.sqrt(resvirt**2+imsvirt**2)


class Extraction:
    """Properties of the extraction point. (Defaults for SPS LSS2 ZS.)
        
    Parameters
    ----------
    alpha : float, optional
        Horizontal Courant-Snyder alpha function at the extraction point. (TODO Units?)
    beta : float, optional
        Horizontal Courant-Snyder beta function at the extraction point. (TODO Units?)
    mu : float, optional
        Horizontal phase advance in radians from the start of the ring to the extraction point. (TODO Units?)
    dx : float, optional
        Horizontal linear spatial dispersion function at the extraction point. (TODO Units?)
    dp : float, optional
        Horizontal linear angular dispersion function at the extraction point. (TODO Units?)
    xbump : float, optional
        Horizontal position of the orbit with respect to the machine
        center at the extraction point. (In meters.)
    pbump : float, optional
        Horizontal angle of the orbit with respect to the machine center
        at the extraction point. (In radians.)
    xwire : float, optional
        Horizontal position of the closest point of the septum wire/blade
        with respect to the machine center at the extraction point. (In meters.)
    wirethickness : float, optional
        Thickness of the septum wire/blade at the extraction point. (In meters.)
    
    Attributes
    ----------
    alpha : float, optional
        Horizontal Courant-Snyder alpha function at the extraction point. (TODO Units?)
    beta : float, optional
        Horizontal Courant-Snyder beta function at the extraction point. (TODO Units?)
    mu : float, optional
        Horizontal phase advance in radians from the start of the ring to the extraction point. (TODO Units?)
    dx : float, optional
        Horizontal linear spatial dispersion function at the extraction point. (TODO Units?)
    dp : float, optional
        Horizontal linear angular dispersion function at the extraction point. (TODO Units?)
    xbump : float, optional
        Horizontal position of the orbit with respect to the machine
        center at the extraction point. (In meters.)
    pbump : float, optional
        Horizontal angle of the orbit with respect to the machine center
        at the extraction point. (In radians.)
    xwire : float, optional
        Horizontal position of the closest point of the septum wire/blade
        with respect to the machine center at the extraction point. (In meters.)
    wirethickness : float, optional
        Thickness of the septum wire/blade at the extraction point. (In meters.)
    dppbump: bool
        TODO to be removed??
    """
    def __init__(self, alpha=2.358057539, beta=96.78034795, mu=0,
                 dx=-0.1673632571, dp=-0.001240955694,
                 xbump=0.04368098538, pbump=-0.0009615184052,
                 xwire=0.06795, wirethickness=0.0002, dppbump=False):
        self.alpha = float(alpha)
        self.beta = float(beta)
        self.mu = float(mu)
        self.dx = float(dx)
        self.dp = float(dp)
        self.xbump = float(xbump)
        self.pbump = float(pbump)
        self.xwire = float(xwire)
        self.wirethickness = float(wirethickness)

        self._xbump = float(xbump)
        self._pbump = float(pbump)

    def adjustd(self, d):
        """
        Adjusts `xbump` and `pbump` proportionally from nominal, such
        that the distance between `xbump` and `xwire` becomes `d`.

        Parameters
        ----------
        d : float
            Desired distance between 'xbump` and `xwire`. (In meters.)
        """
        # Go back to nominal
        self.xbump = self._xbump
        self.pbump = self._pbump
        # Then scale the bump
        if d is not None:
            xbump = self.xwire-d
            self.pbump = self.pbump*(xbump/self.xbump)
            self.xbump = xbump

    def wiretest(self, normalization, pt):
        """Find the normalized distance beyond which a particle is extracted.

        Parameters
        ----------
        normalization : float
            Normalization constant as explained in TODO. (TODO Units?)
        pt : float
            Momentum deviation as in MAD-X for the
            particle to be tested for extraction.

        Returns
        -------
        float
            The distance beyond which a particle is lost on the septum
            wires/blade or blade. (In meters.)
        """
        return (normalization
                *(self.xwire-self.xbump-self.dx*pt)
                /math.sqrt(self.beta))


def get_init(ring, btype="gaussian", scale=math.sqrt(12E-6/426.3156),
             dpp=0.0015, dpp_offset=0, npart=10000, seed=0):
    """Beam parameters and particle initialization for the simulation.
        
    Parameters
    ----------
    ring : multitrack.Ring
        Ring object to define the machine in which to track.
    btype : string, optional
        String representing the type of beam distribution. Allowed
        options are \"gaussian\" or \"explore\".
    scale : float, optional
        If btype=gaussian: Standard deviation of the beam to be
        simulated, in standard normalized phase-space. (TODO Units?)
        If btype=explore: (TODO remember scaling rule for this?)
    dpp : float, optional
        Momentum spread of the beam to be simulated. Particles will be
        initialized with :math:`\\frac{\Delta p}{p_0}` uniformly in
        :math:`[\\texttt{dpp\\_offset}-\\texttt{dpp}, \\texttt{dpp\\_offset}+\\texttt{dpp}]`.
    dpp_offset : float, optional
        Momentum offset of the beam to be simulated. Centre of the uniform
        distribution from which :math:`\\frac{\Delta p}{p_0}` is initialized.
    npart : int, optional
        Number of particles to track. Overridden in case beam.btype=\"explore\".
    seed : int, optional
        Seed with which to initialize the numpy RNGs. (Specify None to
        use a random seed.)
    
    Returns
    ----------
    dataframe
        Initial particle distribution
    """
    np.random.seed(seed)

    normalization = ring.get_normalization()

    if dpp==0.0:
        dpps = [dpp_offset for n in range(npart)]
    else:
        dpps = np.random.uniform(dpp_offset-dpp, dpp_offset+dpp, npart)

    if btype=="gaussian":
        init = np.random.normal(0, scale*normalization, (npart, 2))
    elif btype=="explore":
        init1 = np.array([[0, 0.1*i/scale] for i in range(-15,16)])
        init2 = np.array([np.dot([[c,s],[-s,c]],v) for v in init1 if v[1]!=0])
        init = np.concatenate((init1,init2))
        npart = len(init)
    else:
        print 'ERROR: Desired beam type "'+btype+'" not recognized.'
        exit()

    return pd.DataFrame(data={'X': init[:,0], 'P': init[:,1], 'dpp': dpps})


def simplems(maxkick, xcirc, xextr, beta=100):
    """Generate a thin simple linear massless septum.
        
    Parameters
    ----------
    maxkick : float
        Septum kick at maximum field strength.
    xcirc : float
        Real space x position at the zero-field edge of the field rise.
    xextr : float
        Real space x position at the full-field edge of the field rise.
    beta : float
        Courant-Snyder beta function at the massless septum.
    
    Returns
    ----------
    function
        Takes xhat (float) and normalization (float) as input and returns
        the thin massless septum kick in p-hat(float).
    """
    if xextr>=xcirc:
        def myms(xhat, normalization):
            x = xhat/normalization*np.sqrt(beta)
            if x<=xcirc:
                return 0.0
            elif x>=xextr:
                return normalization*np.sqrt(beta)*maxkick
            else:
                return normalization*np.sqrt(beta)*maxkick*((xhat/normalization*np.sqrt(beta))-xcirc)/(xextr-xcirc)
    else:
        def myms(xhat, normalization):
            x = xhat/normalization*np.sqrt(beta)
            if x>=xcirc:
                return 0.0
            elif x<=xextr:
                return normalization*np.sqrt(beta)*maxkick
            else:
                return normalization*np.sqrt(beta)*maxkick*((xhat/normalization*np.sqrt(beta))-xcirc)/(xextr-xcirc)
    return myms

