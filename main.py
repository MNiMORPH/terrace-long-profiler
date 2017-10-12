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
