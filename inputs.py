from __future__ import division
import math

class Ring:
    """Ring properties for tracking.

    Parameters
    ----------
    elements: :obj:`list` of :obj:`list` of [int, :obj:`dict`]
        A list of elements, each defined by a phase advance in radians
        since the start of the ring and a dictionary of multipole
        components, e.g. [[0,{2:K_2a,4:K_4a}],[1.12,{2:K_2b,3:K_3b}],...]
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
        """
        self.elements = sorted(self.elements+[element],
                                 key=lambda m: m[0])


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

    def wiretest(self, k2, dpp):
        """Find the normalized distance beyond which a particle is extracted.

        Parameters
        ----------
        k2 : float
            Normalization constant as explained in TODO. (TODO Units?)
        dpp : float
            Momentum deviation :math:`\\frac{\Delta p}{p_0}` for the
            particle to be tested for extraction.

        Returns
        -------
        float
            The distance beyond which a particle is lost on the septum
            wires/blade or blade. (In meters.)
        """
        return (abs(k2)
                *(self.xwire-self.xbump-self.dx*dpp)
                /math.sqrt(self.beta))


class Beam:
    """Beam parameters for the simulation.
        
    Parameters
    ----------
    sigma : float, optional
        Standard deviation of the beam to be simulated. (TODO Units?)
    dpp : float, optional
        Momentum spread of the beam to be simulated. Particles will be
        initialized with :math:`\\frac{\Delta p}{p_0}` uniformly in
        :math:`[-\\texttt{dpp}, \\texttt{dpp}]`.
    
    Attributes
    ----------
    sigma : float
        Standard deviation of the beam to be simulated. (TODO Units?)
    dpp : float
        Momentum spread of the beam to be simulated. Particles will be
        initialized with :math:`\\frac{\Delta p}{p_0}` uniformly in
        :math:`[-\\texttt{dpp}, \\texttt{dpp}]`.
    type : string
        String representing the type of beam distribution. Hardcoded to
        'gaussian' is allowed at the moment.
    """
    def __init__(self, sigma=math.sqrt(12E-6/426.3156), dpp=0.0015):
        self.sigma = float(sigma)
        self.dpp = float(dpp)
        self.type = "gaussian"
