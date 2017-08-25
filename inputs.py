import math

class Ring:
    """Ring settings for multipole tracking.

    Parameters
    ----------
    multipoles: list
        A list of multipoles, each defined by a phase advance and a
        dictionary of multipole components, i.e. [[0,{2:K_2a,4:K_4a}],[1.12,{2:K_2b,3:K_3b}],...]
    tune: float
        Fractional part of the machine tune.
    """
    def __init__(self, multipoles, tune=2.0/3, chroma=None):
        self.multipoles = sorted(multipoles, key=lambda m: m[0])
        self.tune = tune
        self.chroma = chroma


    def addmultipole(self, multipole):
        self.multipoles = sorted(self.multipoles+[multipole],
                                 key=lambda m: m[0])


class Extraction:
    def __init__(self, alpha=2.358057539, beta=96.78034795, mu=0,
                 dx=-0.1673632571, dp=-0.001240955694,
                 xbump=0.04368098538, pbump=-0.0009615184052,
                 xwire=0.06795, wirethickness=0.0002, dppbump=False):
        self.alpha = alpha
        self.beta = beta
        self.mu = mu
        self.dx = dx
        self.dp = dp
        self.xbump = xbump
        self.pbump = pbump
        self.xwire = xwire
        self.wirethickness = wirethickness

        self._xbump = xbump
        self._pbump = pbump

    def adjustd(self, d):
        # Go back to nominal
        self.xbump = self._xbump
        self.pbump = self._pbump
        # Then scale the bump
        if d is not None:
            xbump = self.xwire-d
            self.pbump = self.pbump*(xbump/self.xbump)
            self.xbump = xbump

    def wiretest(self, k2, dpp):
        return (abs(k2)
                *(self.xwire-self.xbump-self.dx*dpp)
                /math.sqrt(self.beta))


class Beam:
    def __init__(self, sigma=math.sqrt(12E-6/426.3156), dpp=0.0015):
        self.sigma = sigma
        self.dpp = dpp
        self.type = "gaussian"
