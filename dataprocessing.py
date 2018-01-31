import math
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
import numpy as np
import os

def get_report(tracks, extractt, extraction):
    """Get figures of merit.

    Parameters
    ----------
    Returns tracks, extractt
    -------
    tracks : array
        [npart x (nturns+1) x 2] array or [npart x 2] array of tracks
        TODO coords at injection/extraction?
    extractt : array
        turnnumber on which extraction occured for each particle;
        defaults to -1 if it was not extracted.
    extraction : multitrack.Extraction
        Extraction object that defined the characteristics of the
        extraction point simulated.

    Returns
    -------
    :obj:`dict`
        TODO It's complicated
    """
    npart = len(extractt)
    fulltrack = (len(tracks.shape)==3)

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

    return report


def print_report(report):
    """Fancy printing of figures of merit.
    """
    out = ("Percentage of particles extracted: "+str(report['ejectperc'])+"\n"+
           "Percentage of extracted particles hitting wire: "+str(report['hitperc'])+"\n"+
           "Statistical emittance of the new beam: "+str(report['statemit'])+"\n"+
           "Range of extracted x: ["+str(report['minexx'])+", "+str(report['maxexx'])+"]\n"+
           "Range of extracted p: ["+str(report['minexp'])+", "+str(report['maxexp'])+"]\n"+
           "Standard deviation extracted x: "+str(report['stdexx'])+"\n"+
           "Standard deviation extracted p: "+str(report['stdexp'])+"\n\n")
    print out


def make_plots(plots=None,
               xbin=None, ybin=None, xlim=None, ylim=None, clim=[None,None]):
    # TODO get color array
#        if beam.dpp==0.0:
#            dpps = [0.0 for n in range(npart)]
#            color = None
#        else:
#            dpps = np.random.uniform(-beam.dpp, beam.dpp, npart)
#            color = dpps

    # Plotting
    if plots is not None:
        if color is None:
            ejcolor = None
            recolor = None
        else:
            ejcolor = [color[i] for i in ejected]
            recolor = [color[i] for i in remaining]
        if len(ejected)>0 and "histscatter" in plots:
            histscatter([x*1000 for x in ejx], [p*1000 for p in ejp],
                        xlabel=r"$x\ [\mathrm{mm}]$",
                        ylabel=r"$x'\ [\mathrm{mrad}]$",
                        xstart=extraction.xwire*1000,
                        xbin=extraction.wirethickness*1000, ybin=0.01,
                        save=plots["histscatter"], color=ejcolor, xlim=xlim, ylim=ylim, clim=clim)
        if "remaining" in plots:
            plotbeam(rex, rep, save=plots["remaining"], color=recolor, clim=clim)
        elif "fulltracks" in plots:
            #TODO replace with propper plot command from simple track?
            print "todo: implement missing plot"


def histscatter(xdata, ydata, save="", xbin=None, ybin=None,
                xstart=0.068, xlabel=r'$x\ [\mathrm{m}]$',
                ylabel=r'$p\ [\mathrm{rad}]$', extra=None, color=None,
                xlim=None, ylim=None, clim=[None,None]):
    """Make a scatterplot with projected histograms

    Code based on
    http://matplotlib.org/examples/pylab_examples/scatter_hist.html

    TODO documentation
    """

    nullfmt = NullFormatter()         # no labels

    left, width = 0.1, 0.65
    left_h = left + width + 0.02
    bottom, height = 0.1, 0.65
    bottom_h = bottom + height + 0.02

    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2, height]

    plt.figure(1, figsize=(8, 8))

    axScatter = plt.axes(rect_scatter)

    axHistx = plt.axes(rect_histx)
    axHisty = plt.axes(rect_histy)

    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)
    axScatter.scatter(xdata, ydata, c=color, cmap='viridis', edgecolor='', vmin=clim[0], vmax=clim[1])

    if extra is not None:
        for line in extra:
            linex = [v[0] for v in line]
            liney = [v[1] for v in line]
            axScatter.plot(linex, liney, 'g--')

    weights = 100*np.ones_like(xdata)/len(xdata)

    if xbin is None:
        xbin = (max(xdata)-min(xdata))/100
    if ybin is None:
        ybin = (max(xdata)-min(xdata))/100

    if xlim is None:
        binsx = np.arange(min(xdata)-xbin, max(xdata)+2*xbin, xbin)
    else:
        binsx = np.arange(xlim[0]-xbin, xlim[1]+xbin, xbin)

    if ylim is None:
        binsy = np.arange(min(ydata)-ybin, max(ydata)+2*ybin, ybin)
    else:
        binsy = np.arange(ylim[0]-ybin, ylim[1]+ybin, ybin)




    axScatter.set_xlim((binsx[0], binsx[-1]))
    axScatter.set_ylim((binsy[0], binsy[-1]))

    axHistx.hist(xdata, bins=binsx, weights=weights)
    axHisty.hist(ydata, bins=binsy, weights=weights,
                 orientation='horizontal')

    axHistx.set_xlim(axScatter.get_xlim())
    axHisty.set_ylim(axScatter.get_ylim())

    axScatter.set_xlabel(xlabel)
    axScatter.set_ylabel(ylabel)
    axHistx.set_ylabel(r'$\%$')
    axHisty.set_xlabel(r'$\%$')

    plt.axis('tight')

    if os.path.exists(os.path.dirname(save)):
        plt.savefig(save)
        plt.clf()
    else:
        plt.show()


def plotbeam(xdata, ydata, xstart=None, save="", color=None,
             xlabel=r'$\hat{X}\ [1]$', ylabel=r'$\hat{P}\ [\mathrm{rad}]$', clim=[None,None]):
    """TODO documentation
    """
    fig, ax = plt.subplots()

    if xstart is not None:
        ax.set_xlim(left=xstart)

    ax.scatter(xdata, ydata, c=color, cmap='viridis', edgecolor='', vmin=clim[0], vmax=clim[1])

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    fig.tight_layout()

    if os.path.exists(os.path.dirname(save)):
        plt.savefig(save)
        plt.clf()
    else:
        plt.show()


def evolplot(tracks, lim, save="", xlabel=r'$\hat{X}\ [1]$',
                ylabel=r'$\hat{P}\ [\mathrm{rad}]$'):
    """TODO documentation
    """
    npart = tracks.shape[0]

    fig, ax = plt.subplots()
    ax.set_aspect(1)
    ax.set_xlim([-lim, lim])
    ax.set_ylim([-lim, lim])

    colors = iter(cm.nipy_spectral(np.linspace(0, 1, npart)))
    for part in range(npart):
        ax.scatter(tracks[part, :, 0], tracks[part, :, 1],
                   color=next(colors))

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    fig.tight_layout()

    if os.path.exists(os.path.dirname(save)):
        plt.savefig(save)
        plt.clf()
    else:
        plt.show()

