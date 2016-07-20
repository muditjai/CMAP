import numpy as np
import os
import sys

if __name__ == '__main__':

    os.chdir(os.path.dirname(sys.argv[0]))

    landmarks_file = "..\Test\landmarks.csv"
    landmarks = np.loadtxt(landmarks_file, delimiter=",")

    truth_file = r"..\Test\truth.csv"
    truth = np.loadtxt(truth_file, delimiter=",")

    out_test_mat = np.hstack((landmarks.transpose(), truth.transpose()))

    print("Shapes - landmarks{0}  truth{1}  merged{2}".format(landmarks.shape, truth.shape, out_test_mat.shape))

    header = ",".join(["Gene_{}".format(i+1) for i in range(out_test_mat.shape[1])])
    out_test_file = "..\Test\landmark.truth.merge.transpose.txt"
    np.savetxt(out_test_file, out_test_mat, fmt="%.2f", delimiter=",", header=header, comments='')
