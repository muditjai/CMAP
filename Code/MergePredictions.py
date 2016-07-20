"""
1. Take 10 optional inputs
2. Find the available predictions from headers
3. Construct final prediction file with above predictions and rest as 1
4. convert this to exe
"""

import numpy as np
from typing import Tuple, Dict
import argparse

# TODO Use numRows=1650 for final test, and 1000 for local test
# Arguments - --geneStartId 971 --geneEndId 12320 --numRows 1650,1000 --genFullFile False --inputFileList MergePredictionData\Gene971.txt,MergePredictionData\Gene971.txt,MergePredictionData\Gene971.txt --outputFile MergePredictionData\Gene971_972_973.txt --outputFileColFormat MergePredictionData\Gene971_972_973_col.txt
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--geneStartId", required=True, type=int, nargs='?', default=971)
    parser.add_argument("--geneEndId", required=True, type=int, nargs='?', default=12320)
    parser.add_argument("--numRows", required=True, type=int, nargs='?', default=1650)
    parser.add_argument("--genFullFile", required=True, nargs='?')
    parser.add_argument("--inputFileList", required=True, nargs='?', default="")
    parser.add_argument("--outputFile", required=True, nargs='?', default="")
    parser.add_argument("--outputFileColFormat", required=True, nargs='?', default="")
    args = parser.parse_args()

    gene_start = args.geneStartId
    gene_end = args.geneEndId
    num_rows = args.numRows
    gen_full_file = (args.genFullFile == "True")
    output_file = args.outputFile
    output_file_col_format = args.outputFileColFormat

    # Get and read all input files
    input_file_list = args.inputFileList.split(",", )
    input_file_list = list(filter(None, input_file_list))

    file_content_dict = {}      # type: Dict[str, np.ndarray]
    geneid_filename_dict = {}       # type: Dict[int, Tuple[str,int]]

    for file_name in input_file_list:
        with open(file_name) as fp:

            # Parse header lines
            header = fp.readline()
            genes_in_file = header.split(sep=",")
            gene_col_idx = 0

            # Build a list of tuples (geneid, filename). Detect duplicate geneid
            for gene_header in genes_in_file:
                geneid = int(gene_header.split("_")[1])

                if geneid in geneid_filename_dict:
                    print("**Duplicate gene found. Geneid - {0} in both \nFile1-{1}\nFile2-{2}\nExiting.".format(
                        geneid, geneid_filename_dict[geneid][0], file_name))
                    exit()
                geneid_filename_dict[geneid] = (file_name, gene_col_idx)
                gene_col_idx += 1

        # Read all file contents into dict
        file_data_arr = np.loadtxt(file_name, delimiter=",", skiprows=1)
        file_content_dict[file_name] = file_data_arr.reshape(num_rows, -1)

    # Output to console the source of genes. Output continuous sources as 1 range.
    # Merge the genes in sequence. Use valid gene predictions from files. Rest keep as 1
    output_prediction = np.asarray([])
    output_header = ""
    prev_match_file = ""

    for geneid in range(gene_start, gene_end + 1):
        pred_list = []
        match_file = ""

        if geneid in geneid_filename_dict:
            match_file = geneid_filename_dict[geneid][0]
            match_file_col_idx = geneid_filename_dict[geneid][1]
            pred_list = (file_content_dict[match_file])[:, [match_file_col_idx]]
        elif gen_full_file:     # if full file requested, generate all 1s
            pred_list = np.ones([num_rows, 1])

        # Print gene source status
        if match_file != prev_match_file:
            if match_file == "":
                print("Gene_{0} generated as all 1 - IsFinalGen-{1}".format(geneid, gen_full_file))
            else:
                print("Gene_{0} from file- {1}".format(geneid, match_file))
            prev_match_file = match_file

        if len(pred_list) != 0:
            output_header += ",Gene_" + str(geneid)
            if (len(output_prediction)) == 0:
                output_prediction = pred_list
            else:
                output_prediction = np.hstack((output_prediction, pred_list))

    print("GeneId final value Gene_{0}".format(geneid))

    output_header = output_header.strip(',')

    # Print the shape of merge. Transpose and do file out.
    print("\nOutput Header - " + output_header)
    print("Output transpose no header shape-{0}".format(output_prediction.shape))

    np.savetxt(output_file_col_format, output_prediction, fmt="%.2f", delimiter=",", header=output_header, comments="")
    np.savetxt(output_file, output_prediction.transpose(), fmt="%.2f", delimiter=",", header="", comments="")
