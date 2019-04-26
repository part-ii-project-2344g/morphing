import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
import os
import sys
from _python_master import master_main
from _python_master import IND_MESHES, IND_METHOD, IND_IS_R_X, IND_IS_R_Y, IND_IS_R_Z, IND_IS_SIM, IND_IS_VOL, IND_VALUE

# ##########################################################

axes=["M","V","mode"]
axis1 = axes[0]
axis2 = axes[1]
axis3 = axes[2]
HEADER = ",".join(axes)

# ##########################################################

def write_file(file, measurements):
    for i in range(len(measurements)):
        for j in range(i+1, len(measurements)):
            if measurements[i][IND_MESHES] == measurements[j][IND_MESHES] \
            and measurements[i][IND_METHOD] == measurements[j][IND_METHOD] \
            and measurements[i][IND_IS_R_X] == measurements[j][IND_IS_R_X] \
            and measurements[i][IND_IS_R_Y] == measurements[j][IND_IS_R_Y] \
            and measurements[i][IND_IS_R_Z] == measurements[j][IND_IS_R_Z] \
            and measurements[i][IND_IS_VOL] ^ measurements[j][IND_IS_VOL] \
            and True:
                V_val = measurements[i][IND_VALUE] if measurements[i][IND_IS_VOL] else measurements[j][IND_VALUE]
                M_val = measurements[j][IND_VALUE] if measurements[i][IND_IS_VOL] else measurements[i][IND_VALUE]
                file.write(str(M_val)+","+str(V_val)+","+str(measurements[i][IND_METHOD])+"\n")

# ##########################################################

def write_plot(data_filename, plot_filename):
    df = pd.read_csv(data_filename)
    fig, ax = plt.subplots()
    
    colormap = cm.viridis
    colorlist = [colors.rgb2hex(colormap(i)) for i in np.linspace(0, 0.9, len(set(df[axis3])))]
    
    for i,c in enumerate(colorlist):
        xs,ys = zip(*[(x,y) for (x,y,l) in zip(list(df[axis1]),list(df[axis2]),list(df[axis3])) if (l==["none","poisson","direct"][i])])
        
        ax.scatter(xs, ys, label=["none","poisson","direct"][i], s=30, linewidth=0.1, c=c)

    ax.legend()
    ax.set_xlabel(axis1)
    ax.set_ylabel(axis2)
    data_range = max(df[axis1].max(), df[axis2].max())
    print("data_range = " + str(data_range))
    plt.xlim([1.0, data_range + 0.02*(data_range-1.0)])
    plt.ylim([1.0, data_range + 0.02*(data_range-1.0)])
    plt.gca().set_aspect('equal', adjustable='box')
    #plt.show()
    fig.savefig(plot_filename, dpi=200)

# ##########################################################

if __name__=="__main__":
    master_main(sys.argv[0].split(".")[0], write_file, write_plot, HEADER)
    