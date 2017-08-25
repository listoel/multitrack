import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
import numpy as np
import os

def histscatter(xdata, ydata, save="", xbin=None, ybin=None,
                xstart=0.068, xlabel=r'$x\ [\mathrm{m}]$',
                ylabel=r'$p\ [\mathrm{rad}]$', extra=None, color=None,
                xlim=None, ylim=None, clim=[None,None]):
    """Make a scatterplot with projected histograms

    Code based on
    http://matplotlib.org/examples/pylab_examples/scatter_hist.html
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

