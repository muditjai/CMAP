import numpy as np
import os
import sys
from timeit import default_timer as timer

if __name__ == '__main__':
    os.chdir(os.path.dirname(sys.argv[0]))
    start = timer()
    inp_matrix = np.loadtxt("../training.csv", delimiter=',')
    end = timer()
    print("Read done in {}".format(end - start))

    header1 = ",".join(["Gene_{}".format(i+1) for i in range(inp_matrix.shape[0])])

    start = timer()
    np.savetxt("../training.transposed.txt", inp_matrix.transpose(), delimiter=",", fmt="%.2f", header=header1)
    end = timer()
    print("Write done in {}".format(end - start))
