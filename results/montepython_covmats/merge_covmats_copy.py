"""
Skript by Justus Schwagereit and Sefa Pamuk to "Merge" MP-Covmats. Will decompose Matrix in blocks of Same and different Parameters, eg Cosmmo and nuicance.
The block with the same parameters will get added to eachother using fischer formalism while the diffrent parameters will be placed in a block diagonal form.
"""

import numpy as np
from numpy.linalg import inv
import sys

if len(sys.argv) != 4:
    print ("Usage:\n\
\t$ python3 merge_covmats.py /path/to/input_one.covmat /path/to/input_two.covmat /path/to/output.covmat\n\n\
All parameters present in covmat two but not in covmat one will be attached to covmat one and the\n\
result will be saved in the output covmat. All correlations between the new and old parameters will\n\
be set to zero (even if the second covmat contains correlations between new parameters and parameters\n\
also present in covmat one), or else you will run into positive-definiteness errors in Monte Python.\n\n\
Exiting...")
    exit()

_, matrix_1_path, matrix_2_path, output_path = sys.argv

# read matrices

with open(matrix_1_path, 'r') as f:
    header_1 = np.array(f.readline().strip('\n').strip('#').split())
    matrix_1 = np.array([[float(num) for num in line.split()] for line in f])
    header_1[-1] = str(header_1[-1]) + ","

with open(matrix_2_path, 'r') as f:
    header_2 = np.array(f.readline().strip('\n').strip('#').split())
    matrix_2 = np.array([[float(num) for num in line.split()] for line in f])
    header_2[-1] = str(header_2[-1]) + ","

# Find decomposition of input covmats into  
AohneB_indexinA    = np.array([ iA for iA, Ai in enumerate(header_1) if not Ai in header_2], dtype=int)
BohneA_indexinB    = np.array([ iB for iB, Bi in enumerate(header_2) if not Bi in header_1], dtype=int)
AschnittB_indexinA = np.array([ iA for iA, Ai in enumerate(header_1) if     Ai in header_2], dtype=int)
AschnittB_indexinB = np.array([ iB for iB, Bi in enumerate(header_2) if     Bi in header_1], dtype=int)

# Translate Matricies into fishers 
fisher_1 = inv(matrix_1)
fisher_2 = inv(matrix_2)

CinA = fisher_1[AschnittB_indexinA].T[AschnittB_indexinA]
CinB = fisher_2[AschnittB_indexinB].T[AschnittB_indexinB]
CC = CinA + CinB

AA = fisher_1[AohneB_indexinA].T[AohneB_indexinA].T
BB = fisher_2[BohneA_indexinB].T[BohneA_indexinB].T

AC = fisher_1[AschnittB_indexinA].T[AohneB_indexinA].T
BC = fisher_2[AschnittB_indexinB].T[BohneA_indexinB].T
AB = np.zeros((AA.shape[0],BB.shape[1]))

#Create Output Array and Header
result_fisher=np.block([[CC  ,AC  , BC],
                        [AC.T,AA  , AB],
                        [BC.T,AB.T, BB]])

result_matrix= inv(result_fisher)
result_matrix[np.where(result_matrix==0)]=0

result_header = np.concatenate([header_1[AschnittB_indexinA],header_1[AohneB_indexinA],header_2[BohneA_indexinB]])

with open (output_path, 'w') as f:
    # write header
    f.write('# ' + f'{str(result_header[0]): <27}')
    for i, elem in enumerate(result_header[1:-1]):
        f.write(f'{str(elem): <28}')
    f.write(str(result_header[-1].strip(',')) + "\n")

    # write body
    for (i,j), elem in np.ndenumerate(result_matrix):

        f.write((" " if elem >= 0 else "") + str("{:.20e}".format(elem)) + "\t" + ("\n" if j == len(result_matrix)-1 else ""))

with open ('MCMC_inverse.mat', 'w') as f:
    # write header
    f.write('# ' + f'{str(result_header[0]): <27}')
    for i, elem in enumerate(result_header[1:-1]):
        f.write(f'{str(elem): <28}')
    f.write(str(result_header[-1].strip(',')) + "\n")

    # write body
    for (i,j), elem in np.ndenumerate(result_fisher):

        f.write((" " if elem >= 0 else "") + str("{:.20e}".format(elem)) + "\t" + ("\n" if j == len(result_fisher)-1 else ""))

