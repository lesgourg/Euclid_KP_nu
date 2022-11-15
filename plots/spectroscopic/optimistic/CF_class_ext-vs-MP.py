import os, sys, argparse

os.chdir(os.path.dirname(os.path.realpath(__file__)))
sys.path.append('../../')

from plot_master import plotter

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('--error-only',action='store_true',dest='error_only',
                    help='    Plot error comparions plots only',
                    default=False)
args = parser.parse_args()
error_only = args.error_only

fish_files =  [
    '../../../results/cosmicfish_external/spectroscopic/optimistic/CosmicFish_v1.0_nulcdm_external_class-Optimistic-3PT_HP_GCsp_fishermatrix.txt',
    '../../../results/montepython_fisher/spectroscopic/optimistic_HP/fisher.mat'
              ]

labels=  [r'CF_ext_class GCsp opt',
          r'MP GCsp opt']

cutnames=['Omegam', 'Omegab','ns', 'h','sigma8', 'mnu', 'Neff',  'lnbgs8_1', 'lnbgs8_2', 'lnbgs8_3', 'lnbgs8_4', 'Ps_1', 'Ps_2', 'Ps_3', 'Ps_4']

plotter(fish_files=fish_files,labels=labels,pars=cutnames,error_only=error_only)
