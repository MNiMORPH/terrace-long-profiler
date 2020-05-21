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
from mpl_toolkits import mplot3d
from shapely.geometry import shape, Polygon, Point, LineString, mapping
import fiona
import os
import sys
from scipy.signal import savgol_filter
from scipy import linalg
from scipy import stats
import math

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
    #from mpl_toolkits.mplot3d import Axes3D

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

def get_terrace_areas(terrace_df, res=5):
    """
    This function takes the initial terrace dataframe and calculates the
    area of each terrace.

    Args:
        terrace_df: pandas dataframe with the terrace info
        res: DEM resolution (default =5m)

    Returns:
        dict where key is the terrace ID and value is the terrace area in m^2

    Author: FJC
    """
    # get unique IDs
    terraceIDs = terrace_df.new_ID.unique()

    area_dict = {}

    for terraceID in terraceIDs:
        # get the n rows with this ID
        masked_df = terrace_df[terrace_df['new_ID'] == terraceID]
        n_pixels = len(masked_df.index)

        # get the data resolution of the DEM
        #Cell_area = IO.GetPixelArea(fname_prefix)
        terrace_area = n_pixels * res * res

        area_dict[terraceID] = terrace_area

    return area_dict

def dist_along_line(X, Y, line):
    """
    find the distance of the point X,Y along the line shapefile
    FJC
    """
    dist = line.project(Point(X,Y))
    return dist

def get_distance_along_baseline(terraces, lp):
    """
    This function gets the distance along the baseline for each of the terrace
    points. This gives continuous distances along the baseline, compared to the
    DistAlongBaseline column in the CSV which is just the nearest point along the
    baseline.

    Args:
        terraces: the dataframe with the terrace info
        lp_shp: the shapefile of the points along the baseline

    Returns:
        terrace dataframe with additional column - 'DistAlongBaseline_new'.

    FJC
    """
    baseline_x = lp['X']
    baseline_y = lp['Y']
    points_list = zip(baseline_x, baseline_y)

    line = LineString(points_list)
    terraces['DistAlongBaseline_new'] = terraces.apply(lambda x: dist_along_line(x['X'], x['Y'], line), axis=1)

    return terraces

#---------------------------------------------------------------------------------------------#
# XZ PLOTS
# Functions to make XZ plots of terraces
#---------------------------------------------------------------------------------------------#

def long_profiler(DataDirectory, fname_prefix, terraces, lp, FigFormat='png'):
    """
    This function creates plots of the terraces with distance
    downstream along the main channel.
    It also attempts to fit a spline through each profile but this is not working
    very well yet.

    Args:
        terraces: the dataframe with the terrace info
        lp: the dataframe with the baseline profile info
        FigFormat: the format of the figure, default = png

    Returns:
        terrace long profile plot

    Author: FJC
    """

    fig = plt.figure()
    ax = plt.subplot(111)
    # now plot
    #ax = terraces.plot.scatter(x='DistAlongBaseline_new', y='Elevation', c='new_ID', colormap='viridis', s=0.2)

    terrace_ids = terraces.new_ID.unique()
    for id in terrace_ids:
        sys.stdout.write("This id is %d       \r" %(id))
        this_df = terraces[terraces.new_ID == id]
        this_df['DistAlongBaseline_new'] = this_df['DistAlongBaseline_new']/1000
        sorted_df = this_df.sort_values(by='DistAlongBaseline_new')
        _xTerraces = sorted_df['DistAlongBaseline_new'].values
        _zTerraces = sorted_df['Elevation'].values
        plt.scatter(_xTerraces, _zTerraces, s = 0.2)

        # bin the points into 100 m bins and find the median elevation within that bin.
        #bin_width = 0.05
        #n_bins = int((_xTerraces.max() - _xTerraces.min())/bin_width)
        n_bins = 50
        print("Number of bins:", n_bins)
        bin_medians, bin_edges, binnumber = stats.binned_statistic(_xTerraces, _zTerraces, statistic='median', bins=n_bins)
        bin_width = (bin_edges[1] - bin_edges[0])
        bin_centres = bin_edges[1:] - bin_width/2
        plt.plot(bin_centres, bin_medians, c='red')

        #print("Fitting spline...")
        #fit a spline through the terrace points
        #print(_xTerraces)
        #yhat = savgol_filter(_zTerraces, 101, 2)   # window size 51, polynomial order 3
        #plt.plot(_xTerraces, yhat, color='red')

        plt.xlabel('Distance downstream (km)')
        plt.ylabel('Elevation (m)')
        plt.savefig(DataDirectory+fname_prefix+'_terrace_plot_'+str(id)+'.'+FigFormat,format=FigFormat,dpi=300)
        plt.clf()

def long_profiler_all_terraces(DataDirectory, fname_prefix, terraces, lp, FigFormat='png'):
    """
    Plot each terrace surface against the long profile of the
    main channel.

    Args:
        terraces: the dataframe with the terrace info
        lp: the dataframe with the baseline profile info
        FigFormat: the format of the figure, default = png

    Returns:
        terrace long profile plot

    Author: FJC
    """

    fig = plt.figure()
    ax = plt.subplot(111)
    # now plot
    #ax = terraces.plot.scatter(x='DistAlongBaseline_new', y='Elevation', c='new_ID', colormap='viridis', s=0.2)

    # plot the main stem channel in black
    plt.plot(lp['DistAlongBaseline']/1000,lp['Elevation'], c='k', lw=2)

    # now plot each terrace individually
    terrace_ids = terraces.new_ID.unique()
    # normalize colours by relief above channel
    norm = colors.Normalize(vmin=terraces.ChannelRelief.min(),vmax=terraces.ChannelRelief.max())
    # get area of terraces: size of marker is scaled by area
    areas = get_terrace_areas(terraces, res=5)
    print(areas)

    for i, id in enumerate(terrace_ids):
        sys.stdout.write("This id is %d       \r" %(id))
        this_df = terraces[terraces.new_ID == id]
        this_df['DistAlongBaseline_new'] = this_df['DistAlongBaseline_new']/1000
        sorted_df = this_df.sort_values(by='DistAlongBaseline_new')
        _xTerraces = sorted_df['DistAlongBaseline_new'].values
        _zTerraces = sorted_df['Elevation'].values
        ChannelRelief = sorted_df['ChannelRelief'].values
        #plt.scatter(_xTerraces, _zTerraces, s = 0.2)

        # # bin the points into 100 m bins and find the median elevation within that bin.
        # bin_width = 0.1
        # n_bins = int((_xTerraces.max() - _xTerraces.min())/bin_width)
        # if not n_bins == 0:
        # #n_bins = 20
        #     print("Number of bins:", n_bins)
        #     bin_medians, bin_edges, binnumber = stats.binned_statistic(_xTerraces, _zTerraces, statistic='median', bins=n_bins)
        #     #bin_width = (bin_edges[1] - bin_edges[0])
        #     bin_centres = bin_edges[1:] - bin_width/2
        #     plt.plot(bin_centres, bin_medians, lw=2, c=colours[i])
        mean_elevation = np.mean(_zTerraces)
        std_elevation = np.std(_zTerraces)
        distance = np.take(_xTerraces, _xTerraces.size // 2)
        mean_relief = np.mean(ChannelRelief)

        #print(mean_relief)
        plt.scatter(distance, mean_elevation, c=mean_relief, s=areas[id]/100000, edgecolors='k', cmap=cm.Reds, norm=norm, zorder=1)
        plt.errorbar(distance, mean_elevation, yerr=std_elevation, zorder=0.1, c='0.5', lw=1, capsize=2, alpha=0.5)

        # to do -save the aggregated terrace data to csv and make a plot for the entire Mississippi.


    #ax.set_ylim(200,260)
    # set axis params and save
    ax.set_xlabel('Distance downstream (m)')
    ax.set_ylabel('Elevation (m)')
    plt.colorbar(cmap=cm.Reds,norm=norm, label="Elevation above modern channel (m)")
    #plt.legend(loc='upper right')
    plt.tight_layout()
    plt.savefig(DataDirectory+fname_prefix+'_terrace_plot.'+FigFormat,format=FigFormat,dpi=300)
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
    terrace_df = pd.read_csv(DataDirectory+fname_prefix+'_terrace_info_filtered_dist.csv')
    terrace_df = terrace_df[terrace_df['BaselineNode'] != -9999]

    # read in the baseline channel csv
    lp = pd.read_csv(DataDirectory+fname_prefix+'_baseline_channel_info.csv')
    lp = lp[lp['Elevation'] != -9999]

    # get the distance from outlet along the baseline for each terrace pixels
    #terrace_df = terrace_df.merge(lp, left_on = "BaselineNode", right_on = "node")
    #print(terrace_df.columns)
    flow_dist = terrace_df['DistAlongBaseline_new']/1000
    #print(terrace_df)

	## Getting the extent of our dataset
    xmin = 0
    xmax = flow_dist.max()
    ymin = 0
    ymax = terrace_df["Elevation"].max()

    ## formatting the data in a meshgrid
    X,Y = np.meshgrid(np.linspace(0,xmax,num = prec),np.linspace(0,ymax, num = prec))
    positions = np.vstack([X.ravel(), Y.ravel()[::-1]]) # inverted Y to get the axis in the bottom left
    values = np.vstack([flow_dist, terrace_df['Elevation']])
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
        ax.plot(lp['DistFromOutlet']/1000,lp['Elevation'],'k',lw=1)

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
        plt.clf

#--------------------------------------------------------------------------------------------------#
# 3D plots
#--------------------------------------------------------------------------------------------------#
def plane(x,y,C):
    return C[0]*x + C[1]*y + C[2]

def PlotTerraceSurfaces(DataDirectory, fname_prefix, terraces):
    """
    Make 3d plot of each terrace surface
    """
    # make a figure
    fig = plt.figure()

    # create a groupby object from the terrace dataframe
    terrace_ids = terraces.new_ID.unique()

    for id in terrace_ids:
        sys.stdout.write("This id is %d       \r" %(id))
        # get the x and z data for this terrace id
        this_df = terraces[terraces.new_ID == id]
        #_x = terrace_df['DistAlongBaseline'].values[_terrace_subset]
        #_y = terrace_df['DistToBaseline'].values[_terrace_subset]
        _z = this_df['Elevation'].values
        _X = this_df['X'].values
        _Y = this_df['Y'].values

        # fit a plane to these points
        # form: Z = C[0]*X + C[1]*Y + C[2]
        _XY = np.vstack((_X, _Y, np.ones(len(_Y)))).transpose()
        C,_,_,_ = linalg.lstsq(_XY, _z)
        plane_x = np.linspace(_X.min()-100, _X.max()+100, 1000)
        plane_y = np.linspace(_Y.min()-100, _Y.max()+100, 1000)
        plane_X, plane_Y = np.meshgrid(plane_x, plane_y)
        plane_Z = plane(plane_X,plane_Y,C)

        # plot the plane
        ax = plt.axes(projection='3d')
        ax.plot_wireframe(plane_X, plane_Y, plane_Z, color='black', lw=0.7, alpha=0.5, zorder=0)

        # plot the terrace points
        ax.scatter(_X, _Y, _z, c=_z, s=0.2, edgecolors=None, zorder=2)

        ax.set_xlabel('X (m)')
        ax.set_ylabel('Y (m)')
        ax.set_zlabel('Elevation (m)')
        plt.savefig(DataDirectory+fname_prefix+'_3d_plot_'+str(id)+'.png',format='png',dpi=300)
        plt.clf()
