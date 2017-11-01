#-----------------------------------------------------------------------------#
# Some functions for dealing with colourbars
# FJC 01/11/17
#-----------------------------------------------------------------------------#

import numpy as _np
import matplotlib.colors as _mcolors

def cmap_discretize(N, cmap):
    """Return a discrete colormap from the continuous colormap cmap.

    Arguments:
        cmap: colormap instance, eg. cm.jet.
        N: number of colors.

    Example:
        x = resize(arange(100), (5,100))
        djet = cmap_discretize(cm.jet, 5)
        imshow(x, cmap=djet)
    """

    if type(cmap) == str:
        cmap = _plt.get_cmap(cmap)
    colors_i = _np.concatenate((_np.linspace(0, 1., N), (0.,0.,0.,0.)))
    colors_rgba = cmap(colors_i)
    indices = _np.linspace(0, 1., N+1)
    cdict = {}
    for ki,key in enumerate(('red','green','blue')):
        cdict[key] = [ (indices[i], colors_rgba[i-1,ki], colors_rgba[i,ki])
                       for i in range(N+1) ]
    # Return colormap object.
    return _mcolors.LinearSegmentedColormap(cmap.name + "_%d"%N, cdict, 1024)
