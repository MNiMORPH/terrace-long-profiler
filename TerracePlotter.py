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

#---------------------------------------------------------------------------------------------#
# CSV READERS
# Functions to read in the csv files
#---------------------------------------------------------------------------------------------#

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

#---------------------------------------------------------------------------------------------#
# ANALYSIS FUNCTIONS
# Functions to analyse the terrace info
#---------------------------------------------------------------------------------------------#

def filter_terraces(terrace_df,min_size=5000, max_size=1000000):
    """
    This function takes the initial terrace dataframe and sorts it to remove terraces
    that are too small.

    Args:
        terrace_df: pandas dataframe with the terrace info
        min_size (int): minimum n of pixels in each terrace
        max_size (int): max n of pixels in each terrace

    Returns:
        dataframe with filtered terrace info.

    Author: FJC
    """
    # first get the unique terrace IDs
    terraceIDs = terrace_df.TerraceID.unique()

    # loop through unique IDs and check how many rows correspond to this ID, then
    # remove any that are too small
    for i in terraceIDs:
        n_pixels = len(terrace_df[terrace_df['TerraceID'] == i])
        if n_pixels < min_size or n_pixels > max_size:
            terrace_df = terrace_df[terrace_df['TerraceID'] != i]

    return terrace_df

def get_terrace_dip_and_dipdir(terrace_df):
    """
    This function takes the initial terrace dataframe and calculates the dip and
    strike of the terrace surfaces. Fits a polynomial surface to the distribution
    of terrace elevations and then gets the dip and dip directions of this surface.

    Args:
        terrace_df: pandas dataframe with the terrace info

    Returns:
        dataframe with terrace dip and dip directions

    Author: AW and FJC
    """
    from scipy import linalg
    import math
    from mpl_toolkits.mplot3d import Axes3D

    # get the unique terrace IDs
    terraceIDs = terrace_df.TerraceID.unique()
    print terraceIDs

    dips = []
    dip_dirs = []
    XbarTerraces = []
    YbarTerraces = []

    # get the info for each terrace ID
    for terraceID in terraceIDs:
        _terrace_subset = (terrace_df.TerraceID.values == terraceID)
        _x = terrace_df['DistAlongBaseline'].values[_terrace_subset]
        _y = terrace_df['DistToBaseline'].values[_terrace_subset]
        _z = terrace_df['Elevation'].values[_terrace_subset]
        _X = terrace_df['X'].values[_terrace_subset]
        _Y = terrace_df['Y'].values[_terrace_subset]

        # fit a plane to these points
        # form: Z = C[0]*X + C[1]*Y + C[2]
        _XY = np.vstack((_X, _Y, np.ones(len(_Y)))).transpose()
        C,_,_,_ = linalg.lstsq(_XY, _z)

        # going to get the dip and dip direction using the unit normal vector
        # to the plane.
        a = -C[0]
        b = -C[1]
        c = 1
        n_vec = np.array([a,b,c])
        print n_vec
        # n vector projected onto the xy plane (multiply n_vec by (1,1,0))
        n_xy = np.array([a,b,0])

        #print n_vec

        # now get the dip = angle between n_vec and n_xy. angle between 2 vectors:
        # cos theta = (alpha . beta) / (|alpha| |beta|)
        # This gives the dip in radians
        dip = np.arccos((np.dot(n_vec,n_xy))/(np.linalg.norm(n_vec)*np.linalg.norm(n_xy)))
        # convert to degrees
        dip = 90 - math.degrees(dip)

        # get the dip direction = angle between n_proj and the due north vector y = (0,1,0)
        y_vec = np.array([0,1,0])
        theta = np.arccos((np.dot(n_xy,y_vec))/(np.linalg.norm(n_xy)*np.linalg.norm(y_vec)))
        theta = math.degrees(theta)
        print ("Theta", theta)

        # work out strike depending on orientation
        if a > 0 and b > 0: # x is positive so dip dir is just theta
            strike = 270 + theta
        elif a < 0:
            strike = 270 - theta
        elif a > 0 and b < 0: # x is negative so dip dir is 180+theta?
            strike = theta - 90

        # now get the dip dir using the right hand rule
        dip_dir = strike+90
        if dip_dir > 359:
            dip_dir = dip_dir - 360

        dips.append(dip)
        dip_dirs.append(dip_dir)

        # get the mean x and y for plotting
        XbarTerraces.append(np.mean(_X))
        YbarTerraces.append(np.mean(_Y))

    outarray = np.vstack((XbarTerraces, YbarTerraces, dips, dip_dirs)).transpose()
    _column_names = ('X', 'Y', 'dip', 'dip_azimuth')
    _index = np.arange(len(dips))+1
    output_pd = pd.DataFrame(data = outarray, index=_index, columns=_column_names)
    # output_pd.to_csv(DataDirectory+fname_prefix+'_Dip_DipDirection.csv')
    return output_pd

def get_terrace_areas(terrace_df, min_size=5000):
    """
    This function takes the initial terrace dataframe and calculates the
    area of each terrace.

    Args:
        terrace_df: pandas dataframe with the terrace info


    Returns:
        dict where key is the terrace ID and value is the terrace area in m^2

    Author: FJC
    """
    # filter the terraces to remove ones that are too small


#---------------------------------------------------------------------------------------------#
# XZ PLOTS
# Functions to make XZ plots of terraces
#---------------------------------------------------------------------------------------------#

def long_profiler(DataDirectory,fname_prefix, min_size=5000, FigFormat='png', size_format='ESURF'):
    """
    This function creates a plot of the terraces with distance
    downstream along the main channel.

    Args:
        DataDirectory (str): the data directory
        fname_prefix (str): the name of the DEM
        min_size (int): the minimum number of pixels for a terrace. Any smaller ones will be excluded.
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
    filter_terraces(terraces, min_size)

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
    cax = fig.add_axes([0.83,0.15,0.03,0.8])
    sm = plt.cm.ScalarMappable(cmap=this_cmap, norm=plt.Normalize(vmin=min(newIDs), vmax=max(newIDs)))
    sm._A = []
    cbar = plt.colorbar(sm,cmap=this_cmap,spacing='uniform',cax=cax, label='Terrace ID', orientation='vertical')
    colours.fix_colourbar_ticks(cbar,len(newIDs),cbar_type=int,min_value=min(newIDs),max_value=max(newIDs),labels=newIDs)

    # set axis params and save
    ax.set_xlabel('Distance upstream (m)')
    ax.set_ylabel('Elevation (m)')
    plt.savefig(DataDirectory+fname_prefix+'_terrace_plot.'+FigFormat,format=FigFormat,dpi=300)

#---------------------------------------------------------------------------------------------#
# RASTER PLOTS
# Functions to make raster plots of terrace attributes
#---------------------------------------------------------------------------------------------#

def MakeRasterPlotTerraceIDs(DataDirectory,fname_prefix, FigFormat='png', size_format='ESURF'):
    """
    This function makes a hillshade of the DEM with the terraces
    plotted onto it coloured by their ID

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

def MakeRasterPlotTerraceElev(DataDirectory,fname_prefix, FigFormat='png', size_format='ESURF'):
    """
    This function makes a hillshade of the DEM with the terraces
    plotted onto it coloured by their elevation

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
    print("I'm going to make a raster plot of terrace elevations. Your topographic data must be in ENVI bil format or I'll break!!")

    # get the rasters
    raster_ext = '.bil'
    BackgroundRasterName = fname_prefix+raster_ext
    HillshadeName = fname_prefix+'_hs'+raster_ext
    TerraceElevName = fname_prefix+'_terrace_relief_final'+raster_ext

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
    terrace_cmap = plt.cm.Reds
    #terrace_cmap = colours.cmap_discretize(n_colours,terrace_cmap)
    MF.add_drape_image(TerraceElevName, DataDirectory, colourmap = terrace_cmap, colorbarlabel="Elevation above channel (m)", alpha=0.8)

    ImageName = DataDirectory+fname_prefix+'_terrace_elev_raster_plot.'+FigFormat
    MF.save_fig(fig_width_inches = fig_width_inches, FigFileName = ImageName, FigFormat=FigFormat, Fig_dpi = 300) # Save the figure

def MakeRasterPlotTerraceDips(DataDirectory,fname_prefix,min_size=5000,FigFormat='png',size_format='ESURF'):
    """
    This function makes a raster plot of terrace locations with arrows showing the terrace
    dip and dip directions.
    Dip and dip direction are calculated by fitting a plane to each terrace using least-squares
    regression.

    Args:
        DataDirectory (str): the data directory
        fname_prefix (str): the name of the DEM without extension.
        min_size (int): minimum number of pixels for a terrace, smaller ones will be removed
        FigFormat (str): the figure format, default='png'
        size_format (str): Can be "big" (16 inches wide), "geomorphology" (6.25 inches wide), or "ESURF" (4.92 inches wide) (defualt esurf).

    Returns:
        plot of terrace locations and dip/dip directions

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
    print("I'm going to make a raster plot of terrace elevations. Your topographic data must be in ENVI bil format or I'll break!!")

    # get the rasters
    raster_ext = '.bil'
    BackgroundRasterName = fname_prefix+raster_ext
    HillshadeName = fname_prefix+'_hs'+raster_ext
    TerraceElevName = fname_prefix+'_terrace_relief_final'+raster_ext

    # get the terrace csv
    terraces = read_terrace_csv(DataDirectory,fname_prefix)
    filter_terraces(terraces)

    # get the terrace IDs
    terraceIDs = terraces.TerraceID.unique()
    n_colours=len(terraceIDs)

    # get the terrace dip and dip dirs
    terrace_dips = get_terrace_dip_and_dipdir(terraces)

    # create the map figure
    MF = MapFigure(HillshadeName, DataDirectory, coord_type='UTM_km', colourbar_location='right')
    # add the terrace drape
    terrace_cmap = plt.cm.Reds
    #terrace_cmap = colours.cmap_discretize(n_colours,terrace_cmap)
    MF.add_drape_image(TerraceElevName, DataDirectory, colourmap = terrace_cmap, colorbarlabel="Elevation above channel (m)", alpha=0.8)

    # add arrows oriented in the direction of dip. We might want to colour these by the dip angle?
    MF.add_arrows_from_points(terrace_dips,azimuth_header='dip_azimuth', arrow_length=100)


    ImageName = DataDirectory+fname_prefix+'_terrace_dips_raster_plot.'+FigFormat
    MF.save_fig(fig_width_inches = fig_width_inches, FigFileName = ImageName, FigFormat=FigFormat, Fig_dpi = 300) # Save the figure
