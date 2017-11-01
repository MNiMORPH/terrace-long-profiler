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

def rand_cmap(nlabels, type='bright', first_color_black=True, last_color_black=False, verbose=True):
    """
    Creates a random colormap to be used together with matplotlib. Useful for segmentation tasks

    Args:
        nlabels: Number of labels (size of colormap)
        type: 'bright' for strong colors, 'soft' for pastel colors
        first_color_black: Option to use first color as black, True or False
        last_color_black: Option to use last color as black, True or False
        verbose: Prints the number of labels and shows the colormap. True or False

    Returns:
        colormap for matplotlib

    Author: FJC
    
    from https://github.com/delestro/rand_cmap
    """
    from matplotlib.colors import LinearSegmentedColormap
    import colorsys
    import numpy as np

    if type not in ('bright', 'soft'):
        print ('Please choose "bright" or "soft" for type')
        return

    if verbose:
        print('Number of labels: ' + str(nlabels))

    # Generate color map for bright colors, based on hsv
    if type == 'bright':
        randHSVcolors = [(np.random.uniform(low=0.0, high=1),
                          np.random.uniform(low=0.2, high=1),
                          np.random.uniform(low=0.9, high=1)) for i in xrange(nlabels)]

        # Convert HSV list to RGB
        randRGBcolors = []
        for HSVcolor in randHSVcolors:
            randRGBcolors.append(colorsys.hsv_to_rgb(HSVcolor[0], HSVcolor[1], HSVcolor[2]))

        if first_color_black:
            randRGBcolors[0] = [0, 0, 0]

        if last_color_black:
            randRGBcolors[-1] = [0, 0, 0]

        random_colormap = LinearSegmentedColormap.from_list('new_map', randRGBcolors, N=nlabels)

    # Generate soft pastel colors, by limiting the RGB spectrum
    if type == 'soft':
        low = 0.6
        high = 0.95
        randRGBcolors = [(np.random.uniform(low=low, high=high),
                          np.random.uniform(low=low, high=high),
                          np.random.uniform(low=low, high=high)) for i in xrange(nlabels)]

        if first_color_black:
            randRGBcolors[0] = [0, 0, 0]

        if last_color_black:
            randRGBcolors[-1] = [0, 0, 0]
        random_colormap = LinearSegmentedColormap.from_list('new_map', randRGBcolors, N=nlabels)

    # Display colorbar
    if verbose:
        from matplotlib import colors, colorbar
        from matplotlib import pyplot as plt
        fig, ax = plt.subplots(1, 1, figsize=(15, 0.5))

        bounds = np.linspace(0, nlabels, nlabels + 1)
        norm = colors.BoundaryNorm(bounds, nlabels)

        cb = colorbar.ColorbarBase(ax, cmap=random_colormap, norm=norm, spacing='proportional', ticks=None,
                                   boundaries=bounds, format='%1i', orientation=u'horizontal')

    return random_colormap
