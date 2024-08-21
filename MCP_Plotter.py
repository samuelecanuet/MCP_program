import matplotlib.pyplot as plt
from root2mpl import *
from ROOT import *

import seaborn as sns

custom_params = {
        "xtick.direction" : "out",
        "ytick.direction" : "out",
        "lines.markeredgecolor" : "k",
        "lines.markeredgewidth" : 0.5,
        "lines.linewidth" : 1,
        "lines.markersize" : 5,
        "font.family" : "serif",
        "ytick.labelsize" : 20,
        "xtick.labelsize" : 20,
        "axes.labelsize" : 25,
        "axes.titlesize" : 20,
        "legend.fontsize" : 15,
        "text.usetex" : True,
        # 'figure.subplot.left': 0.20, 
        # 'figure.subplot.bottom': 0.15, 
        # 'figure.subplot.right': 0.95, 
        # 'figure.subplot.top': 0.90
        }

sns.set_theme(style = "ticks", rc=custom_params)

map = "plasma"
file = TFile("MCP3_1800V_100V_Grid40V_WithBeamScan_0001.root")
fig, ax = plt.subplots(figsize=(5,5))
DisplayTH2D(file.Get("h_RAW"), ax = ax, title="Raw data", xlim=(-0.6, 0.6), ylim=(-0.6, 0.6), xlabel='X [au]', ylabel="Y [au]", vmin=5, vmax = 50, color=map)
ax.title.set_size(30)
ax.set_aspect('equal')
ax.set_xticks([-0.5,-0.25, 0, 0.25, 0.5])
ax.set_yticks([-0.5,-0.25, 0, 0.25, 0.5])
plt.subplots_adjust(left=0.25, right=0.95, top=0.9, bottom=0.16)
plt.savefig("RAW.png", dpi=500)
plt.show()


fig, ax = plt.subplots(figsize=(5,5))
DisplayTH2D(file.Get("h_Image_corr"), ax = ax, title="Reconstructed data", xlim=(-8, 8), ylim=(-8, 8), xlabel='X [mm]', ylabel="Y [mm]", vmin=5, vmax = 50, color=map)
ax.title.set_size(30)
ax.set_aspect('equal')
ax.set_xticks([-6,-3, 0, 3, 6])
ax.set_yticks([-6,-3, 0, 3, 6])
plt.subplots_adjust(left=0.25, right=0.95, top=0.9, bottom=0.16)
plt.savefig("REC.png", dpi=500)
plt.show()