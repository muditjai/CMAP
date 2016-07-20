"""
Create training data with header. Transpose the original gene rows.
"""

import numpy as np
import os
import sys
from timeit import default_timer as timer

if __name__ == '__main__':
    inp_file, out_file = ("../training.csv", "../training.transposed.txt")

    os.chdir(os.path.dirname(sys.argv[0]))
    start = timer()
    inp_matrix = np.loadtxt(inp_file, delimiter=',')
    end = timer()
    print("Read done in {}".format(end - start))

    header1 = ",".join(["Gene_{}".format(i+1) for i in range(inp_matrix.shape[0])])

    start = timer()
    np.savetxt(out_file, inp_matrix.transpose(), delimiter=",", fmt="%.2f", header=header1, comments='')
    end = timer()
    print("Write done in {}".format(end - start))
