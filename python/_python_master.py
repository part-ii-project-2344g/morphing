import sys
import pandas as pd
import statistics
import numpy as np
    
# ##########################################################

IND_MESHES, IND_METHOD, IND_IS_R_X, IND_IS_R_Y, IND_IS_R_Z, IND_IS_SIM, IND_IS_VOL, IND_VALUE = 0, 1, 2, 3, 4, 5, 6, 7
str_true = "true"
str_false = "false"


def master_main(subprogram, write_file, write_plot, data_header):
    print("Running subprogram '" + subprogram + "'...")
    
    data_filename = subprogram+".csv"
    plot_filename = subprogram+".png"
    f = open("_measurements.csv", "r")
    f2 = open(data_filename, "w")
    f2.write(data_header+"\n")
    measurements = []
    skipped_header = False
    for line in f:
        if not skipped_header:
            skipped_header = True
            continue
        meshes = line.split(",")[0]+","+line.split(",")[1]
        method = line.split(",")[2]
        is_rotx = line.split(",")[3] == str_true
        is_roty = line.split(",")[4] == str_true
        is_rotz = line.split(",")[5] == str_true
        is_simple = line.split(",")[6] == str_true
        is_v = line.split(",")[7] == str_true
        val = float(line.split(",")[8][:-1])
        measurements.append((meshes, method, is_rotx, is_roty, is_rotz, is_simple, is_v, val))
    f.close()
    print("Running write_file()...")
    write_file(f2, measurements)
    f2.flush()
    f2.close()
    print("Running write_plot()...")
    write_plot(data_filename, plot_filename)

    # Data analysis.
    df = pd.read_csv(data_filename)
    xs = list(df.iloc[:,0])
    ys = list(df.iloc[:,1])
    diffs = [x-y for (x,y) in zip(xs,ys)]
    print("average improvement x-y: " + "{0:.3f}".format(sum(diffs)/len(diffs)))
    print("median improvement x-y: " + "{0:.3f}".format(statistics.median(diffs)))
    print("percentage of improvement x>y: " + "{0:.3f}".format(100.0*len([d for d in diffs if d > 0.0])/len(diffs)) + "%")
    print("percentage of worsening y>x: " + "{0:.3f}".format(100.0*len([d for d in diffs if d < 0.0])/len(diffs)) + "%")
    A = np.vstack([xs, np.ones(len(xs))]).T
    slope, constant = np.linalg.lstsq(A, ys, rcond=None)[0]
    print("best-fit line: y = " + "{0:.3f}".format(slope) + "x + " + "{0:.3f}".format(constant) + ".")
    pearsonr = np.corrcoef(xs,ys)[1,0]
    print("Pearson correlation coefficient: " + str(pearsonr))
    print("average x: " + "{0:.3f}".format(sum(xs)/len(xs)))
    print("average y: " + "{0:.3f}".format(sum(ys)/len(ys)))