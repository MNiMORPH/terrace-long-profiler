#-----------------------------------------------------------------------------#
# Some functions for dealing with colourbars
# FJC 01/11/17
#-----------------------------------------------------------------------------#

import numpy as _np
import matplotlib.colors as _mcolors
from matplotlib import ticker

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

def fix_colourbar_ticks(cbar,n_colours, cbar_type=float, min_value = 0, max_value = 0, colourbar_orientation='vertical'):
    """
    This function takes a discrete colourbar and fixes the ticks so they are
    in the middle of each colour

    Update: SMM, made this more flexible so you can set the minimum and maximum
     values. Required becase basin plotting and point plotting doesn't use
     raster values
    Update 2: BG, somehow, writing "update" in upper cases make my text editors (Atom and brackets) completely panic.
    so I changed it, sorry.

    Args:
        n_colours (int): number of colours in discrete colourbar.
        cbar_type (type): Sets the type of the colourbar (if you want int labels, set to int).
        min_value: the minimum value on the colourbar
        max_value: the maximum value on the colourbar
        colourbar_orientation: orientation, either horizontal or vertical. Default vertical.

    Returns:
        None but fixes ticks

    Author: FJC
    """

    vmin=min_value
    vmax=max_value

    # get the additional end spacing for colourbar
    tick_spacing = (vmax-vmin)/(n_colours)
    print tick_spacing
    new_vmin = vmin-(tick_spacing/2)
    new_vmax = vmax+(tick_spacing/2)+tick_spacing

    #get list of tick locations
    tick_locs = _np.arange(new_vmin, new_vmax, step=tick_spacing)
    print tick_locs

    # update ticks
    tick_locator = ticker.FixedLocator(tick_locs)
    cbar.locator = tick_locator
    cbar.update_ticks()

    # get tick labels
    tick_labels = _np.linspace(vmin, vmax, n_colours)
    if cbar_type == int:
        tick_labels = [str(int(x)) for x in tick_labels]
    else:
        tick_labels = [str(x) for x in tick_labels]
    print tick_labels

    if colourbar_orientation == "horizontal":
        cbar.ax.set_xticklabels(tick_labels)
    else:
        cbar.ax.set_yticklabels(tick_labels)
