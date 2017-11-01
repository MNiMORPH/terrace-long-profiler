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
import colours
import matplotlib.cm as cm

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

def long_profiler(DataDirectory,fname_prefix):
    """
    This function creates a plot of the terraces with distance
    downstream along the main channel.

    Args:
        DataDirectory (str): the data directory
        fname_prefix (str): the name of the DEM

    Returns:
        terrace long profile plot

    Author: AW, FJC
    """
    """
    lp = pd.read_csv('Rio_Toro_baseline_channel_info.csv')

    x_lp = sorted(list(set(list(lp.DistAlongBaseline))))
    z_lp = []
    for x_i in x_lp:
        z_lp.append(np.mean(lp.Elevation.values[x_lp == x_i]))
    # /usr/bin/ipython:2: VisibleDeprecationWarning: boolean index did not match indexed array along dimension 0; dimension is 3790331 but corresponding boolean dimension is 6490

    plt.plot(x_lp, z_lp, 'k-', linewidth=2)
    """
    # read in the terrace csv
    terraces = read_terrace_csv(DataDirectory,fname_prefix)

    # read in the baseline channel csv
    lp = read_channel_csv(DataDirectory,fname_prefix)

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

    all_x = []
    all_z = []
    all_ids = []
    #flatten into arrays for plotting
    for n in range(len(xTerraces)):
        for x in range(len(xTerraces[n])):
            all_x.append(xTerraces[n][x])
            all_z.append(zTerraces[n][x])
            all_ids.append(newIDs[n])

    # make the plot
    fig = plt.figure()
    gs = plt.GridSpec(100,100,bottom=0.15,left=0.05,right=0.85,top=1.0)
    ax = fig.add_subplot(gs[5:100,10:95])

    # make a random colormap
    this_cmap = cm.rainbow
    this_cmap = colours.cmap_discretize(len(newIDs),this_cmap)
    colors = iter(this_cmap(np.linspace(0, 1, len(newIDs))))

    for i in range(len(xTerraces)):
        plt.plot(xTerraces[i], zTerraces[i], 'o', c=next(colors))

    #norm=plt.Normalize(vmin=min(all_ids), vmax=max(all_ids))
    #ax.scatter(all_x,all_z,c=all_ids, cmap=new_cmap)

    # add a colourbar
    #this_cmap = colours.cmap_discretize(cm.rainbow,len(newIDs))
    cax = fig.add_axes([0.86,0.15,0.03,0.8])
    sm = plt.cm.ScalarMappable(cmap=this_cmap, norm=plt.Normalize(vmin=min(newIDs), vmax=max(newIDs)))
    sm._A = []
    cbar = plt.colorbar(sm,cmap=this_cmap,spacing='uniform',cax=cax, label='Terrace ID', orientation='vertical')
    colours.fix_colourbar_ticks(cbar,len(newIDs),cbar_type=int,min_value=min(all_ids),max_value=max(all_ids))

    #plt.ion()
    #plt.show()
    ax.set_xlabel('Distance upstream (m)')
    ax.set_ylabel('Elevation (m)')
    plt.savefig(DataDirectory+fname_prefix+'_terrace_plot_test.png',format='png',dpi=300)


    x_terraces = sorted(list(set(list(lp.DistAlongBaseline))))
    z_terraces = []
    for x_i in x_terraces:
        z_terraces.append(np.mean(lp.Elevation.values[x_terraces == x_i]))
