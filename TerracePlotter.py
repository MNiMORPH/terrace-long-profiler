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
import matplotlib.cm as cm
from matplotlib import rcParams
from matplotlib import colors as colors
from shapely.geometry import shape, Polygon, Point, LineString
import fiona
import os
from scipy.interpolate import UnivariateSpline

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
        cmap = plt.get_cmap(cmap)
    colors_i = np.concatenate((np.linspace(0, 1., N), (0.,0.,0.,0.)))
    colors_rgba = cmap(colors_i)
    indices = np.linspace(0, 1., N+1)
    cdict = {}
    for ki,key in enumerate(('red','green','blue')):
        cdict[key] = [ (indices[i], colors_rgba[i-1,ki], colors_rgba[i,ki])
                       for i in range(N+1) ]
    # Return colormap object.
    return colors.LinearSegmentedColormap(cmap.name + "_%d"%N, cdict, 1024)

#---------------------------------------------------------------------------------------------#
# Set up figure
#---------------------------------------------------------------------------------------------#
def CreateFigure(FigSizeFormat="default", AspectRatio=16./9.):
    """
    This function creates a default matplotlib figure object

    Args:
        FigSizeFormat: the figure size format according to journal for which the figure is intended
            values are geomorphology,ESURF, ESPL, EPSL, JGR, big
            ddefault is ESURF

        AspectRatio: The shape of the figure determined by the aspect ratio, default is 16./9.

    Returns:
        matplotlib figure object

    Author: FJC
    """
    # set figure sizes (in inches) based on format
    if FigSizeFormat == "geomorphology":
        FigWidth_Inches = 6.25
    elif FigSizeFormat == "big":
        FigWidth_Inches = 16
    elif FigSizeFormat == "ESURF":
        FigWidth_Inches = 4.92
    elif FigSizeFormat == "ESPL":
        FigWidth_Inches = 7.08
    elif FigSizeFormat == "EPSL":
        FigWidth_Inches = 7.48
    elif FigSizeFormat == "JGR":
        FigWidth_Inches = 6.6

    else:
        FigWidth_Inches = 4.92126

    # Set up fonts for plots
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = 10
    #rcParams['text.usetex'] = True

    Fig = plt.figure(figsize=(FigWidth_Inches,FigWidth_Inches/AspectRatio))

    return Fig
# #---------------------------------------------------------------------------------------------#
# ANALYSIS FUNCTIONS
# Functions to analyse the terrace info
#---------------------------------------------------------------------------------------------#

def write_dip_and_dipdir_to_csv(DataDirectory,fname_prefix):
    """
    Wrapper for dip and dipdir function

    Args:
        DataDirectory (str): the data directory
        fname_prefix (str): name of the DEM

    Author: FJC
    """
    # read in the terrace csv
    terraces = pd.read_csv(DataDirectory+fname_prefix+'_terrace_info_filtered.csv')

    # get the terrace dip and dip dirs
    terrace_dips = get_terrace_dip_and_dipdir(terraces)

    # write to csv
    terrace_dips.to_csv(DataDirectory+fname_prefix+'_Dip_DipDirection.csv')


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
    print(terraceIDs)

    dips = []
    dip_dirs = []
    strikes = []
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
        print(n_vec)
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
        print("Theta", theta)

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
        strikes.append(strike)

        # get the mean x and y for plotting
        XbarTerraces.append(np.mean(_X))
        YbarTerraces.append(np.mean(_Y))

    outarray = np.vstack((XbarTerraces, YbarTerraces, dips, dip_dirs, strikes)).transpose()
    _column_names = ('X', 'Y', 'dip', 'dip_azimuth', 'strike')
    _index = np.arange(len(dips))+1
    output_pd = pd.DataFrame(data = outarray, index=_index, columns=_column_names)
    return output_pd

def get_terrace_areas(terrace_df, fname_prefix):
    """
    This function takes the initial terrace dataframe and calculates the
    area of each terrace.

    Args:
        terrace_df: pandas dataframe with the terrace info
        fname_prefix: name of the DEM (to get data res)

    Returns:
        dict where key is the terrace ID and value is the terrace area in m^2

    Author: FJC
    """
    # get unique IDs
    terraceIDs = terrace_df.terraceID.unique()

    area_dict = {}

    for terraceID in terraceIDs:
        # get the n rows with this ID
        masked_df = terrace_df[terrace_df['terraceID'] == terraceID]
        n_pixels = len(masked_df.index)

        # get the data resolution of the DEM
        Cell_area = IO.GetPixelArea(fname_prefix)
        terrace_area = n_pixels * Cell_area

        area_dict[terraceID] = terrace_area

    return area_dict

#---------------------------------------------------------------------------------------------#
# XZ PLOTS
# Functions to make XZ plots of terraces
#---------------------------------------------------------------------------------------------#

def long_profiler(terraces, lp, FigFormat='png'):
    """
    This function creates a plot of the terraces with distance
    downstream along the main channel.

    Args:
        terraces: the dataframe with the terrace info
        lp: the dataframe with the baseline profile info
        FigFormat: the format of the figure, default = png

    Returns:
        terrace long profile plot

    Author: AW, FJC
    """

    # get a pandas groupby object to group terraces by IDs
    ax = terraces.plot.scatter(x='DistAlongBaseline', y='Elevation', c='new_ID', colormap='viridis', s=0.2)

    # plot the main stem channel in black
    plt.plot(lp['DistAlongBaseline'],lp['Elevation'], c='k', lw=2)

    # set axis params and save
    ax.set_xlabel('Distance downstream (m)')
    ax.set_ylabel('Elevation (m)')
    plt.savefig(DataDirectory+fname_prefix+'_terrace_plot.'+FigFormat,format=FigFormat,dpi=300)
    plt.clf()

def long_profiler_spline(terraces, lp, FigFormat='png'):
    """
    Fit a spline to each terrace surface and plot against the long profile of the
    main channel. NOTE - this is not working yet.

    Args:
        terraces: the dataframe with the terrace info
        lp: the dataframe with the baseline profile info
        FigFormat: the format of the figure, default = png

    Returns:
        terrace long profile plot

    Author: FJC
    """
    # make a figure
    fig = CreateFigure()
    ax = plt.subplot(111)

    # create a groupby object from the terrace dataframe
    terrace_ids = terraces.new_ID.unique()

    for id in terrace_ids:
        # get the x and z data for this terrace id
        this_df = terraces[terraces.new_ID == id]
        _xTerraces = this_df['DistAlongBaseline']
        _zTerraces = this_df['Elevation']
        print(_xTerraces)

        # fit a univariate spline through these points
        spl = UnivariateSpline(_xTerraces, _zTerraces)
        xs = np.linspace(_xTerraces.min(), _xTerraces.max(), 1000)
        plt.plot(xs, spl(xs), lw=3)

    # plot the main stem channel
    plt.plot(lp['DistAlongBaseline'],lp['Elevation'], c='k', lw=2)
    # set axis params and save
    ax.set_xlabel('Distance downstream (m)')
    ax.set_ylabel('Elevation (m)')
    plt.savefig(DataDirectory+fname_prefix+'_terrace_plot_spline.'+FigFormat,format=FigFormat,dpi=300)
    plt.clf()


def long_profiler_dist(DataDirectory,fname_prefix, min_size=5000, FigFormat='png', size_format='ESURF'):
    """
    Make long profile plot where terrace points are binned by
    distance along the channel
    """
    # make a figure
    fig = CreateFigure()
    ax = plt.subplot(111)

    # read in the terrace csv
    terraces = pd.read_csv(DataDirectory+fname_prefix+'_terrace_info_filtered.csv')

    # read in the baseline channel csv
    lp = pd.read_csv(DataDirectory+fname_prefix+'_baseline_channel_info.csv')
    lp = lp[lp['Elevation'] != -9999]

    # get the distance from outlet along the baseline for each terrace pixels
    new_terraces = terraces.merge(lp, left_on = "BaselineNode", right_on = "node")
    print(new_terraces)

    xTerraces = np.array(new_terraces['DistFromOutlet'])
    yTerraces = np.array(new_terraces['DistToBaseline'])
    zTerraces = np.array(new_terraces['Elevation_x'])

    MaximumDistance = xTerraces.max()

    # now bin by distance along the baseline
    bins = np.unique(xTerraces)
    nbins = len(np.unique(xTerraces))
    n, _ = np.histogram(xTerraces, bins=nbins)
    s_zTerraces, _ = np.histogram(xTerraces, bins=nbins, weights=zTerraces)
    s_zTerraces2, _ = np.histogram(xTerraces, bins=nbins, weights=zTerraces*zTerraces)
    mean = s_zTerraces / n
    std = np.sqrt(s_zTerraces2/n - mean*mean)

    # # invert to get distance from outlet
    # MS_DistAlongBaseline = np.array(lp['DistAlongBaseline'])[::-1]
    MS_Dist = np.array(lp['DistFromOutlet'])
    MS_Elevation = np.array(lp['Elevation'])
    Terrace_Elevation = mean

    print(MS_Dist)
    print(MS_Elevation)
    print(Terrace_Elevation)

    # plot the main stem channel in black
    plt.plot(MS_Dist/1000,MS_Elevation, c='k', lw=1)
    plt.scatter((_/1000)[:-1], Terrace_Elevation, s=2, zorder=2, c='r')

    # set axis params and save
    ax.set_xlabel('Distance from outlet (km)')
    ax.set_ylabel('Elevation (m)')
    ax.set_xlim(0,80)
    plt.tight_layout()
    plt.savefig(DataDirectory+fname_prefix+'_terrace_plot_binned.'+FigFormat,format=FigFormat,dpi=300)

    plt.clf()


def MakeTerraceHeatMap(DataDirectory,fname_prefix, prec=100, bw_method=0.03, FigFormat='png', ages=""):
    """
    Function to make a heat map of the terrace pixels using Gaussian KDE.
    see https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.stats.gaussian_kde.html
    for more details.

    Args:
        DataDirectory(str): the data directory
        fname_prefix(str): prefix of your DEM
        prec(int): the resolution for the KDE. Increase this to get a finer resolution, decrease for coarser.
        bw_method: the method for determining the bandwidth of the KDE.  This is apparently quite sensitive to this.
        Can either be "scott", "silverman" (where the bandwidth will be determined automatically), or a scalar. Default = 0.03
        FigFormat(str): figure format, default = png
        ages (str): Can pass in the name of a csv file with terrace ages which will be plotted on the profile. Must be in the same directory

    FJC 26/03/18
    """
    import scipy.stats as st

    # make a figure
    fig = CreateFigure()
    ax = plt.subplot(111)

    # read in the terrace DataFrame
    terrace_df = pd.read_csv(DataDirectory+fname_prefix+'_terrace_info_filtered.csv')
    terrace_df = terrace_df[terrace_df['BaselineNode'] != -9999]

    # read in the baseline channel csv
    lp = pd.read_csv(DataDirectory+fname_prefix+'_baseline_channel_info.csv')
    lp = lp[lp['Elevation'] != -9999]

    # get the distance from outlet along the baseline for each terrace pixels
    terrace_df = terrace_df.merge(lp, left_on = "BaselineNode", right_on = "node")
    print(terrace_df.columns)
    flow_dist = terrace_df['DistAlongBaseline_x']/1000
    print(terrace_df)

	## Getting the extent of our dataset
    xmin = 0
    xmax = flow_dist.max()
    ymin = 0
    ymax = terrace_df["Elevation_x"].max()

    ## formatting the data in a meshgrid
    X,Y = np.meshgrid(np.linspace(0,xmax,num = prec),np.linspace(0,ymax, num = prec))
    positions = np.vstack([X.ravel(), Y.ravel()[::-1]]) # inverted Y to get the axis in the bottom left
    values = np.vstack([flow_dist, terrace_df['Elevation_x']])
    if len(values) == 0:
        print("You don't have any terraces, I'm going to quit now.")
    else:
        # get the kernel density estimation
        KDE = st.gaussian_kde(values, bw_method = bw_method)
        Z = np.reshape(KDE(positions).T,X.shape)

        # plot the density on the profile
        cmap = cm.gist_heat_r
        cmap.set_bad(alpha=0)
        cb = ax.imshow(Z, interpolation = "None",  extent=[xmin, xmax, ymin, ymax], cmap=cmap, aspect = "auto")

        # plot the main stem channel
        ax.plot(lp['DistFromOutlet']/1000,lp['Elevation_y'],'k',lw=1)

        # if present, plot the ages on the profile
        if ages:
            # read in the ages csv
            ages_df = pd.read_csv(DataDirectory+ages)
            upstream_dist = list(ages_df['upstream_dist'])
            elevation = list(ages_df['elevation'])
            ax.scatter(upstream_dist, elevation, s=8, c="w", edgecolors="k", label="$^{14}$C age (cal years B.P.)")
            ax.legend(loc='upper left', fontsize=8, numpoints=1)

        # set some plot lims
        ax.set_xlim(xmin,xmax)
        ax.set_ylim(ymin,ymax)
        ax.set_xlabel('Flow distance (km)')
        ax.set_ylabel('Elevation (m)')

        # add a colourbar
        cbar = plt.colorbar(cb,cmap=cmap,orientation='vertical')
        cbar.set_label('Density')

        # save the figure
        plt.tight_layout()
        plt.savefig(T_directory+fname_prefix+'_terrace_plot_heat_map.png',format=FigFormat,dpi=300)
        plt.clf()
