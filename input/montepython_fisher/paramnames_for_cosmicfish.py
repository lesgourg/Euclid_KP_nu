#!/usr/bin/python

# creates a .paramnames file matching the expectation of cosmicfish for comparison plots out of MP's log file
#
# takes one argument: the directory. Example:
#
# python3 input/montepython_fisher/paramnames_for_cosmicfish.py results/montepython_fisher/photometric/pessimistic

import sys
import re

logfile = str(sys.argv[1])+'/log.param'
paramnamesfile = str(sys.argv[1])+'/fisher.paramnames'
print ('Reading in:',logfile)
print ('Writing in:',paramnamesfile)

with open(logfile) as f:
    lines = f.readlines()

o = open(paramnamesfile,"w")

for line in lines:
    if "#" in line:
        pass
    elif "['Omega_b']" in line:
        fiducial = ((line.split("=")[1][2:-2]).split(",")[0]).strip()
        o.write(r'Omegab    \Omega_{\mathrm{b}, 0}    '+fiducial+'\n')
    elif "['h']" in line:
        fiducial = ((line.split("=")[1][2:-2]).split(",")[0]).strip()
        o.write(r'h    h    '+fiducial+'\n')
    elif "['w0_fld']" in line:
        fiducial = ((line.split("=")[1][2:-2]).split(",")[0]).strip()
        o.write(r'w0    w_0    '+fiducial+'\n')
    elif "['wa_fld']" in line:
        fiducial = ((line.split("=")[1][2:-2]).split(",")[0]).strip()
        o.write(r'wa    w_a    '+fiducial+'\n')
    elif "['n_s']" in line:
        fiducial = ((line.split("=")[1][2:-2]).split(",")[0]).strip()
        o.write(r'ns    n_\mathrm{s}    '+fiducial+'\n')
    elif "['sigma8']" in line:
        fiducial = ((line.split("=")[1][2:-2]).split(",")[0]).strip()
        o.write(r'sigma8    \sigma_8    '+fiducial+'\n')
    #elif "['N_eff_camb']" in line:
    #    fiducial = ((line.split("=")[1][2:-2]).split(",")[0]).strip()
    #    o.write(r'Neff    N_\mathrm{eff}    '+fiducial+'\n')
    #elif "['m_nu_camb']" in line:
    #    fiducial = ((line.split("=")[1][2:-2]).split(",")[0]).strip()
    #    o.write(r'mnu    M_\nu    '+fiducial+'\n')
    elif "['Omega_m_camb']" in line:
        fiducial = ((line.split("=")[1][2:-2]).split(",")[0]).strip()
        o.write(r'Omegam    \Omega_{\mathrm{m}, 0}    '+fiducial+'\n')
    elif "['bias_1']" in line:
        fiducial = ((line.split("=")[1][2:-2]).split(",")[0]).strip()
        o.write(r'b1    b_1    '+fiducial+'\n')
    elif "['bias_2']" in line:
        fiducial = ((line.split("=")[1][2:-2]).split(",")[0]).strip()
        o.write(r'b2    b_2    '+fiducial+'\n')
    elif "['bias_3']" in line:
        fiducial = ((line.split("=")[1][2:-2]).split(",")[0]).strip()
        o.write(r'b3    b_3    '+fiducial+'\n')
    elif "['bias_4']" in line:
        fiducial = ((line.split("=")[1][2:-2]).split(",")[0]).strip()
        o.write(r'b4    b_4    '+fiducial+'\n')
    elif "['bias_5']" in line:
        fiducial = ((line.split("=")[1][2:-2]).split(",")[0]).strip()
        o.write(r'b5    b_5    '+fiducial+'\n')
    elif "['bias_6']" in line:
        fiducial = ((line.split("=")[1][2:-2]).split(",")[0]).strip()
        o.write(r'b6    b_6    '+fiducial+'\n')
    elif "['bias_7']" in line:
        fiducial = ((line.split("=")[1][2:-2]).split(",")[0]).strip()
        o.write(r'b7    b_7    '+fiducial+'\n')
    elif "['bias_8']" in line:
        fiducial = ((line.split("=")[1][2:-2]).split(",")[0]).strip()
        o.write(r'b8    b_8    '+fiducial+'\n')
    elif "['bias_9']" in line:
        fiducial = ((line.split("=")[1][2:-2]).split(",")[0]).strip()
        o.write(r'b9    b_9    '+fiducial+'\n')
    elif "['bias_10']" in line:
        fiducial = ((line.split("=")[1][2:-2]).split(",")[0]).strip()
        o.write(r'b10    b_10    '+fiducial+'\n')
    elif "['aIA']" in line:
        fiducial = ((line.split("=")[1][2:-2]).split(",")[0]).strip()
        o.write(r'AIA    A_\mathrm{IA}    '+fiducial+'\n')
    elif "['etaIA']" in line:
        fiducial = ((line.split("=")[1][2:-2]).split(",")[0]).strip()
        o.write(r'etaIA    \eta_mathrm{IA}    '+fiducial+'\n')
    elif "['lnbsigma8_0']" in line:
        fiducial = ((line.split("=")[1][2:-2]).split(",")[0]).strip()
        o.write(r'lnbgs8_1    lnbgs8_1    '+fiducial+'\n')
    elif "['lnbsigma8_1']" in line:
        fiducial = ((line.split("=")[1][2:-2]).split(",")[0]).strip()
        o.write(r'lnbgs8_2    lnbgs8_2    '+fiducial+'\n')
    elif "['lnbsigma8_2']" in line:
        fiducial = ((line.split("=")[1][2:-2]).split(",")[0]).strip()
        o.write(r'lnbgs8_3    lnbgs8_3    '+fiducial+'\n')
    elif "['lnbsigma8_3']" in line:
        fiducial = ((line.split("=")[1][2:-2]).split(",")[0]).strip()
        o.write(r'lnbgs8_4    lnbgs8_4    '+fiducial+'\n')
    elif "['P_shot0']" in line:
        fiducial = ((line.split("=")[1][2:-2]).split(",")[0]).strip()
        o.write(r'Ps_1    Ps_1    '+fiducial+'\n')
    elif "['P_shot1']" in line:
        fiducial = ((line.split("=")[1][2:-2]).split(",")[0]).strip()
        o.write(r'Ps_2    Ps_2    '+fiducial+'\n')
    elif "['P_shot2']" in line:
        fiducial = ((line.split("=")[1][2:-2]).split(",")[0]).strip()
        o.write(r'Ps_3    Ps_3    '+fiducial+'\n')
    elif "['P_shot3']" in line:
        fiducial = ((line.split("=")[1][2:-2]).split(",")[0]).strip()
        o.write(r'Ps_4    Ps_4    '+fiducial+'\n')
