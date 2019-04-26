import os
import os.path

str_false = "false"
str_true = "true"

if __name__ == "__main__":
    f = open("_measurements.csv", "w")
    f.write("mesh1, mesh2, method, is_rotx, is_roty, is_rotz, is_simple, is_v, value\n")
    for dirpath, dirnames, filenames in os.walk("D:\\!K\\UNIVERSITY\\YEAR3\\Dissertation\\morphing\\morphs"):
        for filename in [fname for fname in filenames if (fname.endswith(".msr") or fname.endswith(".msrv"))]:
            print("path = " + os.path.join(dirpath, filename))
            msr = open(os.path.join(dirpath, filename), "r")
            line = msr.read()
            # Final measurement result: 230.157/115.579 = 1.99134.
            msr.close()
            val = line.split("=")[1][1:-2]
            mesh1 = filename.split("@")[0][1:]
            mesh2 = filename.split("@")[1]
            method = filename.split("@")[2].split("-")[0]
            is_rotx, is_roty, is_rotz = str_false, str_false, str_false
            if mesh1.endswith("_rotx80"):
                mesh1 = mesh1[:-7]
                is_rotx = str_true
            elif mesh1.endswith("_roty80"):
                mesh1 = mesh1[:-7]
                is_roty = str_true
            elif mesh1.endswith("_rotz80"):
                mesh1 = mesh1[:-7]
                is_rotz = str_true
            is_simple = str_true if (("simple" in mesh1) or ("simple" in mesh2)) else str_false
            is_v = str_true if filename.endswith(".msrv") else str_false
            f.write(mesh1 + "," + mesh2 + "," + method + "," + is_rotx + "," + is_roty + "," + is_rotz \
            + "," + is_simple + "," + is_v + "," + val + "\n")
    f.flush()
    f.close()
