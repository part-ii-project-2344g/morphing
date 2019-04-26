import pandas as pd
import matplotlib.pyplot as plt
import os
import sys
from _python_master import master_main
from _python_master import IND_MESHES, IND_METHOD, IND_IS_R_X, IND_IS_R_Y, IND_IS_R_Z, IND_IS_SIM, IND_IS_VOL, IND_VALUE

# ##########################################################

axis1 = "M"
axis2 = "V"
HEADER = axis1+","+axis2

# ##########################################################

def write_file(file, measurements):
    for i in range(len(measurements)):
        for j in range(i+1, len(measurements)):
            if measurements[i][IND_MESHES] == measurements[j][IND_MESHES] \
            and measurements[i][IND_METHOD] == measurements[j][IND_METHOD] \
            and measurements[i][IND_IS_R_X] == measurements[j][IND_IS_R_X] \
            and measurements[i][IND_IS_R_Y] == measurements[j][IND_IS_R_Y] \
            and measurements[i][IND_IS_R_Z] == measurements[j][IND_IS_R_Z] \
            and measurements[i][IND_IS_SIM] == measurements[j][IND_IS_SIM] \
            and measurements[i][IND_IS_SIM] \
            and measurements[i][IND_IS_VOL] ^ measurements[j][IND_IS_VOL] \
            and True:
                V_val = measurements[i][IND_VALUE] if measurements[i][IND_IS_VOL] else measurements[j][IND_VALUE]
                M_val = measurements[j][IND_VALUE] if measurements[i][IND_IS_VOL] else measurements[i][IND_VALUE]
                file.write(str(M_val)+","+str(V_val)+"\n")

# ##########################################################

def write_plot(data_filename, plot_filename):
    df = pd.read_csv(data_filename)
    plot = df.plot(kind='scatter',x=axis1,y=axis2) # scatter plot
    #plt.axis('equal')
    data_range = max(df[axis1].max(), df[axis2].max())
    print("data_range = " + str(data_range))
    plt.xlim([1.0, data_range + 0.02*(data_range-1.0)])
    plt.ylim([1.0, data_range + 0.02*(data_range-1.0)])
    plt.gca().set_aspect('equal', adjustable='box')
    #plt.show()
    fig = plot.get_figure()
    fig.savefig(plot_filename, dpi=200)

# ##########################################################

if __name__=="__main__":
    master_main(sys.argv[0].split(".")[0], write_file, write_plot, HEADER)
    