import os, sys, argparse

os.chdir(os.path.dirname(os.path.realpath(__file__)))
sys.path.append('../../')

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter) 
parser.add_argument('--error-only',action='store_true',dest='error_only',
                    help='    Plot error comparions plots only',
                    default=False)
args = parser.parse_args()
error_only = args.error_only

from plot_master import plotter

fish_files =  [
    '../../../results/cosmicfish_external/spectroscopic/pessimistic/CosmicFish_v0.9_varying_mnu__external_class-Pessimistic-own_HP_GCsp_fishermatrix.txt',
    '../../../results/montepython_fisher/spectroscopic/pessimistic/fisher.mat'    
        ]

labels=  [r'CF_ext_class GCsp pess',
          r'MP GCsp pess']

cutnames=['Neff', 'mnu', 'Omegam', 'Omegab','ns', 'h','sigma8', 'lnbgs8_1', 'lnbgs8_2', 'lnbgs8_3', 'lnbgs8_4', 'Ps_1', 'Ps_2', 'Ps_3', 'Ps_4']

plotter(fish_files=fish_files,labels=labels,pars=cutnames,error_only=error_only)
