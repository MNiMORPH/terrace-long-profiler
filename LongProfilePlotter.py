#------------------------------------------------------------------------------#
# terrace-long-profiler
# Scripts to plot the long profiles of river terraces compared to the main
# channel.
# Authors: A. Wickert
#          F. Clubb
#------------------------------------------------------------------------------#

# import modules
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from LSDPlottingTools import colours
import matplotlib.cm as cm
from matplotlib import rcParams

def read_terrace_csv(DataDirectory,fname_prefix):
    """
    This function reads in the csv file with the extension "_terrace_info.csv"
    and returns it as a pandas dataframe

    Args:
        DataDirectory (str): the data directory
        fname_prefix (str): the name of the DEM

    Returns:
        pandas dataframe with the terrace info

    Author: FJC
    """
    csv_suffix = '_terrace_info.csv'
    fname = DataDirectory+fname_prefix+csv_suffix

    df = pd.read_csv(fname)

    return df

def read_channel_csv(DataDirectory,fname_prefix):
    """
    This function reads in the csv file with the extension "_baseline_channel_info.csv"
    and returns it as a pandas dataframe

    Args:
        DataDirectory (str): the data directory
        fname_prefix (str): the name of the DEM

    Returns:
        pandas dataframe with the channel info

    Author: FJC
    """
    csv_suffix = '_baseline_channel_info.csv'
    fname = DataDirectory+fname_prefix+csv_suffix

    df = pd.read_csv(fname)

    return df

def long_profiler(DataDirectory,fname_prefix, FigFormat='png', size_format='ESURF'):
    """
    This function creates a plot of the terraces with distance
    downstream along the main channel.

    Args:
        DataDirectory (str): the data directory
        fname_prefix (str): the name of the DEM
        FigFormat: the format of the figure, default = png
        size_format (str): Can be "big" (16 inches wide), "geomorphology" (6.25 inches wide), or "ESURF" (4.92 inches wide) (defualt esurf).

    Returns:
        terrace long profile plot

    Author: AW, FJC
    """
    # make a figure
    if size_format == "geomorphology":
        fig = plt.figure(1, facecolor='white',figsize=(6.25,3.5))
        #l_pad = -40
    elif size_format == "big":
        fig = plt.figure(1, facecolor='white',figsize=(16,9))
        #l_pad = -50
    else:
        fig = plt.figure(1, facecolor='white',figsize=(4.92126,3.2))
        #l_pad = -35

    # make the plot
    gs = plt.GridSpec(100,100,bottom=0.15,left=0.05,right=0.85,top=1.0)
    ax = fig.add_subplot(gs[5:100,10:95])

    # read in the terrace csv
    terraces = read_terrace_csv(DataDirectory,fname_prefix)

    # read in the baseline channel csv
    lp = read_channel_csv(DataDirectory,fname_prefix)
    lp = lp[lp['Elevation'] != -9999]

    terraceIDs = sorted(list(set(list(terraces.TerraceID))))
    xTerraces = []
    zTerraces = []
    yTerraces = []
    newIDs = []

    # loop through the terrace IDs and get the x, y, and z values
    for terraceID in terraceIDs:
        _terrace_subset = (terraces.TerraceID.values == terraceID)
        _x = terraces['DistAlongBaseline'].values[_terrace_subset]
        _y = terraces['DistToBaseline'].values[_terrace_subset]
        _z = terraces['Elevation'].values[_terrace_subset]
        _x_unique = sorted(list(set(list(_x))))
        _z_unique = []
        # Filter
        if len(_x) > 50 and len(_x_unique) > 1 and len(_x_unique) < 1000:
            #print np.max(np.diff(_x_unique))
            #if len(_x_unique) > 10 and len(_x_unique) < 1000 \
            #and :
            for _x_unique_i in _x_unique:
                #_y_unique_i = np.min(np.array(_y)[_x == _x_unique_i])
                #_z_unique.append(np.min(_z[_y == _y_unique_i]))
                _z_unique.append(np.min(_z[_x == _x_unique_i]))
            if np.mean(np.diff(_z_unique)/np.diff(_x_unique)) < 10:
                xTerraces.append(_x_unique)
                zTerraces.append(_z_unique)
                newIDs.append(terraceID)

    # get discrete colours so that each terrace is a different colour
    this_cmap = cm.rainbow
    this_cmap = colours.cmap_discretize(len(newIDs),this_cmap)
    print "N COLOURS: ", len(newIDs)
    print newIDs
    colors = iter(this_cmap(np.linspace(0, 1, len(newIDs))))
    # plot the terraces
    for i in range(len(xTerraces)):
        plt.scatter(xTerraces[i], zTerraces[i], s=2, c=next(colors))

    # plot the main stem channel in black
    plt.plot(lp['DistAlongBaseline'],lp['Elevation'], c='k', lw=2)

    # add a colourbar
    cax = fig.add_axes([0.86,0.15,0.03,0.8])
    sm = plt.cm.ScalarMappable(cmap=this_cmap, norm=plt.Normalize(vmin=min(newIDs), vmax=max(newIDs)))
    sm._A = []
    cbar = plt.colorbar(sm,cmap=this_cmap,spacing='uniform',cax=cax, label='Terrace ID', orientation='vertical')
    colours.fix_colourbar_ticks(cbar,len(newIDs),cbar_type=int,min_value=min(newIDs),max_value=max(newIDs),labels=newIDs)

    # set axis params and save
    ax.set_xlabel('Distance upstream (m)')
    ax.set_ylabel('Elevation (m)')
    plt.savefig(DataDirectory+fname_prefix+'_terrace_plot.'+FigFormat,format=FigFormat,dpi=300)

    # x_terraces = sorted(list(set(list(lp.DistAlongBaseline))))
    # z_terraces = []
    # for x_i in x_terraces:
    #     z_terraces.append(np.mean(lp.Elevation.values[x_terraces == x_i]))

def raster_plotter(DataDirectory,fname_prefix, FigFormat='png', size_format='ESURF'):
    """
    This function makes a hillshade of the DEM with the terraces
    plotted onto it.

    Args:
        DataDirectory (str): the data directory
        fname_prefix (str): the name of the DEM without extension.
        FigFormat (str): the figure format, default='png'
        size_format (str): Can be "big" (16 inches wide), "geomorphology" (6.25 inches wide), or "ESURF" (4.92 inches wide) (defualt esurf).

    Returns:
        Raster plot of terrace IDs

    Author: FJC

    """
    from LSDMapFigure.PlottingRaster import BaseRaster
    from LSDMapFigure.PlottingRaster import MapFigure

    # Set up fonts for plots
    label_size = 10
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size

    # make a figure
    if size_format == "geomorphology":
        #fig = plt.figure(1, facecolor='white',figsize=(6.25,3.5))
        fig_width_inches=6.25
        #l_pad = -40
    elif size_format == "big":
        #fig = plt.figure(1, facecolor='white',figsize=(16,9))
        fig_width_inches=16
        #l_pad = -50
    else:
        fig_width_inches = 4.92126
        #fig = plt.figure(1, facecolor='white',figsize=(4.92126,3.2))
        #l_pad = -35

    # going to make the terrace plots - need to have bil extensions.
    print("I'm going to make the terrace raster plot. Your topographic data must be in ENVI bil format or I'll break!!")

    # get the rasters
    raster_ext = '.bil'
    BackgroundRasterName = fname_prefix+raster_ext
    HillshadeName = fname_prefix+'_hs'+raster_ext
    TerraceIDName = fname_prefix+'_terrace_IDs'+raster_ext

    # get the terrace csv
    terraces = read_terrace_csv(DataDirectory,fname_prefix)
    terraceIDs = sorted(list(set(list(terraces.TerraceID))))
    xTerraces = []
    zTerraces = []
    yTerraces = []
    newIDs = []

    # loop through the terrace IDs and get the x, y, and z values
    for terraceID in terraceIDs:
        _terrace_subset = (terraces.TerraceID.values == terraceID)
        _x = terraces['DistAlongBaseline'].values[_terrace_subset]
        _y = terraces['DistToBaseline'].values[_terrace_subset]
        _z = terraces['Elevation'].values[_terrace_subset]
        _x_unique = sorted(list(set(list(_x))))
        _z_unique = []
        # Filter
        if len(_x) > 50 and len(_x_unique) > 1 and len(_x_unique) < 1000:
            #print np.max(np.diff(_x_unique))
            #if len(_x_unique) > 10 and len(_x_unique) < 1000 \
            #and :
            for _x_unique_i in _x_unique:
                #_y_unique_i = np.min(np.array(_y)[_x == _x_unique_i])
                #_z_unique.append(np.min(_z[_y == _y_unique_i]))
                _z_unique.append(np.min(_z[_x == _x_unique_i]))
            if np.mean(np.diff(_z_unique)/np.diff(_x_unique)) < 10:
                xTerraces.append(_x_unique)
                zTerraces.append(_z_unique)
                newIDs.append(terraceID)

    n_colours=len(newIDs)

    # create the map figure
    MF = MapFigure(HillshadeName, DataDirectory, coord_type='UTM_km', colourbar_location='right')
    # add the terrace drape
    terrace_cmap = plt.cm.rainbow
    #terrace_cmap = colours.cmap_discretize(n_colours,terrace_cmap)
    MF.add_drape_image(TerraceIDName, DataDirectory, colourmap = terrace_cmap, discrete_cmap=True, cbar_type=int, n_colours=n_colours, colorbarlabel="Terrace ID", alpha=0.8)

    ImageName = DataDirectory+fname_prefix+'_terrace_IDs_raster_plot.'+FigFormat
    MF.save_fig(fig_width_inches = fig_width_inches, FigFileName = ImageName, FigFormat=FigFormat, Fig_dpi = 300) # Save the figure
