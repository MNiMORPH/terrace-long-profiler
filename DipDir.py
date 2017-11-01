import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from scipy import linalg

"""
lp = pd.read_csv('Rio_Toro_baseline_channel_info.csv')

x_lp = sorted(list(set(list(lp.DistAlongBaseline))))
z_lp = []
for x_i in x_lp:
    z_lp.append(np.mean(lp.Elevation.values[x_lp == x_i]))
# /usr/bin/ipython:2: VisibleDeprecationWarning: boolean index did not match indexed array along dimension 0; dimension is 3790331 but corresponding boolean dimension is 6490

plt.plot(x_lp, z_lp, 'k-', linewidth=2)
"""

#terraces = pd.read_csv('Rio_Toro_terrace_info.csv')
terraces = pd.read_csv('Eel_DEM_clip_terrace_info.csv')
#channel = pd.read_csv('Eel_DEM_clip_baseline_channel_info.csv')
terraceIDs = sorted(list(set(list(terraces.TerraceID))))
xTerraces = [] # along-channel x
zTerraces = []
yTerraces = [] # cross-channel y
XTerraces = [] # UTM x
YTerraces = [] # UTM y
XbarTerraces = [] # UTM x -- mean
YbarTerraces = [] # UTM y -- mean
zbarTerraces = [] # z mean
ATerraces = [] # Areas -- 1 point / m**2
dips = []
dip_directions = []
for terraceID in terraceIDs:
    _terrace_subset = (terraces.TerraceID.values == terraceID)
    _x = terraces['DistAlongBaseline'].values[_terrace_subset]
    _y = terraces['DistToBaseline'].values[_terrace_subset]
    _z = terraces['Elevation'].values[_terrace_subset]
    _X = terraces['X'].values[_terrace_subset]
    _Y = terraces['Y'].values[_terrace_subset]
    if len(_z) > 40000:
        # fit a plane to these points
        # form: Z = C[0]*X + C[1]*Y + C[2]
        _XY = np.vstack((_X, _Y, np.ones(len(_Y)))).transpose()
        C,_,_,_ = linalg.lstsq(_XY, _z)
        # Get dip and dip direction
        _dip_slope = (C[0]**2 + C[1]**2)**.5
        _dip = np.arctan(_dip_slope)
        # Dip orientation is +/- this (pos/neg, not approx.)
        print C[0], _dip
        _dip_orientation_pm = np.arccos(C[0] / _dip_slope)
        # Test which solution (+/-) is required
        _isplus = (np.round(_dip * np.cos(_dip_orientation_pm),6) == np.round(C[0],6))
        #print _dip * np.cos(_dip_orientation_pm), C[0]
        # + C[1] * np.sin(_dip_orientation_pm)
        _sign = 2 * _isplus - 1
        _dip_orientation = _sign * _dip_orientation_pm
        # Test which direction, now that we have orientation
        # Which goes downwards?
        _slope = C[0] * np.cos(_dip_orientation) + C[1] * np.sin(_dip_orientation)
        #_neg = C[0] * np.cos(-_dip_orientation) + C[1] * np.sin(-_dip_orientation)
        _ispos = _slope > 0#np.abs(_pos) < 0 #> np.abs(_neg)
        # Flip with np.pi if the slope is going up -- dips go downslope
        _dip_direction = _dip_orientation + np.pi * _ispos
        # * (_ispos == 0)
        # Test for dip direction -- this is like a difference, except first
        # point is the origin
        #_dip_direction_sign = np.sign( C[0] * np.cos(_dip_orientation) + C[1] * np.sin(_dip_orientation) )
        #_dip_direction = _dip_direction_sign * _dip_orientation
        #_dip_direction = (np.pi * (_dip_orientation < 0)) + _dip_orientation
        xTerraces.append(_x)
        zTerraces.append(_z)
        XTerraces.append(_X)
        YTerraces.append(_Y)
        XbarTerraces.append(np.mean(_X))
        YbarTerraces.append(np.mean(_Y))
        zbarTerraces.append(np.mean(_z))
        ATerraces.append(len(_z))
        dips.append(_dip)
        dip_directions.append(_dip_direction)

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        #ax.scatter(_X, _Y, _z, c='blue', depthshade=True)
        _x0, _y0, _z0 = np.mean(_X), np.mean(_Y), np.mean(_z)
        _u0, _v0, _w0 = np.cos(_dip_direction), np.sin(_dip_direction), -_dip_slope
        #_u0, _v0, _w0 = -np.cos(_dip_direction), -np.sin(_dip_direction), -_dip
        #_u0, _v0, _w0 = C[0]/-_dip, C[1]/-_dip, -_dip
        ax.quiver(_x0, _y0, _z0, _u0, _v0, _w0, length=50, arrow_length_ratio=.05, edgecolor='black', linewidth=4)
        surf_xy = np.meshgrid( np.linspace(np.min(_X), np.max(_X), 20), 
        np.linspace(np.min(_Y), np.max(_Y), 20))
        surf_z = C[0] * surf_xy[0] + C[1] * surf_xy[1] + C[2]
        ax.plot_wireframe(surf_xy[0], surf_xy[1], surf_z, facecolors='k', alpha=1)
        plt.title(terraceID)
        
        # 54, 17, can I use slope?
        
    """
    _x_unique = sorted(list(set(list(_x))))
    _z_unique = []
    # Filter
    if len(_x) > 200 and len(_x_unique) > 1 and len(_x_unique) < 1000:
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
    """


# Export
dips_degrees = np.array(dips) * 180./np.pi
dip_directions_degrees = np.array(dip_directions) * 180./np.pi
dip_directions_degrees[dip_directions_degrees < 0] += 360
outarray = np.vstack((XbarTerraces, YbarTerraces, dips_degrees, dip_directions_degrees)).transpose()
_column_names = ('X', 'Y', 'dip', 'dip_azimuth')
_index = np.arange(len(dips))+1
output_pd = pd.DataFrame(data = outarray, index=_index, columns=_column_names)
output_pd.to_csv('Eel_Terraces_Dip_DipDirection.csv')

# Detailed terrace figure
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(XTerraces[0], YTerraces[0], zTerraces[0], c='blue', depthshade=True)
_x0, _y0, _z0 = (np.min(XTerraces[0]) + np.max(XTerraces[0]))/2., (np.min(YTerraces[0]) + np.max(YTerraces[0]))/2., (np.min(zTerraces[0]) + np.max(zTerraces[0]))/2.
_u0, _v0, _w0 = np.cos(dip_direction), np.sin(dip_direction), -dip
ax.quiver(_x0, _y0, _z0, _u0, _v0, _w0, length=20)
surf_xy = np.meshgrid( np.linspace(np.min(XTerraces[0]), np.max(XTerraces[0]), 20), np.linspace(np.min(YTerraces[0]), np.max(YTerraces[0]), 20))
surf_z = C[0] * surf_xy[0] + C[1] * surf_xy[1] + C[2]
ax.plot_wireframe(surf_xy[0], surf_xy[1], surf_z, facecolors='k', alpha=1)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
#ax.scatter(_X, _Y, _z, c='blue', depthshade=True)
_x0, _y0, _z0 = np.mean(_X), np.mean(_Y), np.mean(_z)
_u0, _v0, _w0 = np.cos(-_dip_direction), np.sin(-_dip_direction), -_dip_slope
ax.quiver(_x0, _y0, _z0, _u0, _v0, _w0, length=50, arrow_length_ratio=.05)
surf_xy = np.meshgrid( np.linspace(np.min(_X), np.max(_X), 20), 
np.linspace(np.min(_Y), np.max(_Y), 20))
surf_z = C[0] * surf_xy[0] + C[1] * surf_xy[1] + C[2]
ax.plot_wireframe(surf_xy[0], surf_xy[1], surf_z, facecolors='k', alpha=1)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.quiver(XbarTerraces, YbarTerraces, np.cos(dip_directions), np.sin(dip_directions))


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.quiver(XbarTerraces, YbarTerraces, zbarTerraces, np.cos(dip_directions), np.sin(dip_directions), -dip, arrow_length_ratio=0, length=100)


# Terrace long profile figure
fig = plt.figure()
for i in range(len(xTerraces)):
    plt.plot(xTerraces[i], zTerraces[i], 'o')
plt.ion()
plt.show()


x_terraces = sorted(list(set(list(lp.DistAlongBaseline))))
z_terraces = []
for x_i in x_terraces:
    z_terraces.append(np.mean(lp.Elevation.values[x_terraces == x_i]))

