import numpy as np
import os
import sys


(model_prediction_file, final_prediction_file) = (r"..\Results\Run2_fixheadertest_Gene971\97.inst.txt",
                                                  r"..\Test\prediction_gene971_all1.txt")


def get_model_prediction():
    return np.loadtxt(model_prediction_file, delimiter="\t", skiprows=1)


if __name__ == '__main__':

    os.chdir(os.path.dirname(sys.argv[0]))

    model_output = get_model_prediction()

    truth_file = np.ones([11350, 1000])  # np.loadtxt(r"..\Test\truth.csv", delimiter=",")

    prediction = np.vstack((model_output[:, 2], truth_file[1:]))

    print("Shapes - model_output{0}  truth_file{1}  prediction{2}".format(model_output.shape, truth_file.shape,
                                                                          prediction.shape))

    np.savetxt(final_prediction_file, prediction, delimiter=",", comments="", fmt="%.2f")
