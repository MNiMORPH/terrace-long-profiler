import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

lp = pd.read_csv('Rio_Toro_baseline_channel_info.csv')

x_lp = sorted(list(set(list(lp.DistAlongBaseline))))
z_lp = []
for x_i in x_lp:
    z_lp.append(np.mean(lp.Elevation.values[x_lp == x_i]))
# /usr/bin/ipython:2: VisibleDeprecationWarning: boolean index did not match indexed array along dimension 0; dimension is 3790331 but corresponding boolean dimension is 6490

plt.plot(x_lp, z_lp, 'k-', linewidth=2)

terraces = pd.read_csv('Rio_Toro_terrace_info.csv')
terraceIDs = sorted(list(set(list(terraces.TerraceID))))
xTerraces = []
zTerraces = []
for terraceID in terraceIDs:
    _terrace_subset = (terraces.TerraceID.values == terraceID)
    _x = terraces['DistAlongBaseline'].values[_terrace_subset]
    _z = terraces['Elevation'].values[_terrace_subset]
    _x_unique = sorted(list(set(list(_x))))
    _z_unique = []
    for _x_unique_i in _x_unique:
        _z_unique.append(np.mean(_z[_x_unique == _x_unique_i]))
    xTerraces.append(_x_unique)
    zTerraces.append(_z_unique)

fig = plt.figure()
for i in range(len(xTerraces)):
    plt.plot(xTerraces[i], zTerraces[i], 'ko')

x_terraces = sorted(list(set(list(lp.DistAlongBaseline))))
z_terraces = []
for x_i in x_terraces:
    z_terraces.append(np.mean(lp.Elevation.values[x_terraces == x_i]))

