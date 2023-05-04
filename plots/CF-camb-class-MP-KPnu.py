import os, sys, argparse

os.chdir(os.path.dirname(os.path.realpath(__file__)))
sys.path.append('../../')
from plot_master import *

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter) 
parser.add_argument('--error-only',action='store_true',dest='error_only',
                    help='    Plot error comparions plots only',
                    default=False)
args = parser.parse_args()
error_only = args.error_only

fish_files =  [
    'sefa/CosmicFish_v1.0_wCDM+mnu+Neff_internal_camb-Optimistic-3PT_GCsp_fishermatrix.txt',
    'sefa/CosmicFish_v1.0_wCDM+mnu+Neff_internal_class-Optimistic-3PT_GCsp_fishermatrix.txt'
    #,
#    'sefa/MPfisher-GCsp.mat',
#    '../../../results/cosmicfish_internal/spectroscopic/optimistic/CosmicFish_v0.9_w0wa_internal_camb-Optimistic-own_GCsp_fishermatrix.txt'
              ]
    
labels=  [
          r'${\tt CF/int/CAMB-Pcb-9p}$',
          r'${\tt CF/int/CLASS-Pcb-9p}$'
        #,
 #         r'${\tt MP-Pcb-MnuNeff-9p}$'
#          r'${\tt CF/int/CAMB}$'
]

nuis=['lnbgs8_1', 'lnbgs8_2', 'lnbgs8_3', 'lnbgs8_4', 'Ps_1', 'Ps_2', 'Ps_3', 'Ps_4']
cosmo_names=['Omegam', 'Omegab','ns', 'h','sigma8', 'w0', 'wa', 'mnu', 'Neff']
cutnames=cosmo_names+nuis

compare_errors_dict={'legend_title':'GCsp opt'}
plotter(fish_files=fish_files,labels=labels,pars=cutnames, outfile_name='all9pars',
        error_only=error_only, compare_errors_dict=compare_errors_dict, 
        outpath='./sefa/', cosmo_names=cosmo_names, yrange=[-2, 2],
        marginalise=False)
